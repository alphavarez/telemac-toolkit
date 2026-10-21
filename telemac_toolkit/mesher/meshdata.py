"""A solver-neutral mesh: nodes, triangles and quadrilaterals, attributes, boundary strings.

The generator works in triangles. Structured zones (channel strips, regular
grids) are made of cells split in two triangles; ``Mesh.from_result`` can put
those pairs back together as quadrilaterals for programs that accept them.
Every writer reads this one structure, so supporting a new program means one
writer, not a new mesher.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..telemac_io.mesh import boundary_loops, signed_areas
from .geometry import points_in_rings

__all__ = ["BoundaryString", "Mesh", "Frame", "pairs_to_quads", "find_quads", "split_quads"]


@dataclass
class BoundaryString:
    """Consecutive boundary nodes sharing one boundary condition."""

    name: str
    codes: Tuple[int, int, int]  # TELEMAC LIHBOR, LIUBOR, LIVBOR (2 2 2 = wall)
    nodes: np.ndarray  # 0-based, in boundary order
    value: Optional[float] = None
    number: int = 0  # liquid boundary number (TELEMAC numbering), 0 for walls

    @property
    def is_wall(self) -> bool:
        return tuple(int(c) for c in self.codes) == (2, 2, 2)


@dataclass
class Frame:
    """Coordinate frame of an export: x' = R(-angle) (x - origin)."""

    origin: Tuple[float, float] = (0.0, 0.0)
    angle_deg: float = 0.0

    @property
    def is_identity(self) -> bool:
        return self.origin == (0.0, 0.0) and self.angle_deg == 0.0

    def forward(self, x, y):
        a = math.radians(self.angle_deg)
        dx, dy = np.asarray(x, float) - self.origin[0], np.asarray(y, float) - self.origin[1]
        return dx * math.cos(a) + dy * math.sin(a), -dx * math.sin(a) + dy * math.cos(a)

    def inverse(self, u, v):
        a = math.radians(self.angle_deg)
        u, v = np.asarray(u, float), np.asarray(v, float)
        return self.origin[0] + u * math.cos(a) - v * math.sin(a), self.origin[1] + u * math.sin(a) + v * math.cos(a)

    def describe(self) -> str:
        return "origin (%.3f, %.3f), rotation %.4f deg" % (self.origin[0], self.origin[1], self.angle_deg)


def _ccw_quads(x, y, q: np.ndarray) -> np.ndarray:
    """Order each quad's four nodes counter-clockwise around its centre."""
    q = np.asarray(q, dtype=np.int64).reshape(-1, 4)
    if len(q) == 0:
        return q
    cx, cy = x[q].mean(axis=1, keepdims=True), y[q].mean(axis=1, keepdims=True)
    ang = np.arctan2(y[q] - cy, x[q] - cx)
    order = np.argsort(ang, axis=1)
    return np.take_along_axis(q, order, axis=1)


def pairs_to_quads(x, y, ikle: np.ndarray) -> np.ndarray:
    """Quads from consecutive triangle pairs (2k, 2k+1), as the structured meshers write them."""
    ik = np.asarray(ikle, dtype=np.int64)
    n = len(ik) // 2
    quads = np.empty((n, 4), dtype=np.int64)
    for k in range(n):
        nodes = list(dict.fromkeys(ik[2 * k].tolist() + ik[2 * k + 1].tolist()))
        if len(nodes) != 4:
            raise ValueError("triangles %d and %d do not form a quadrilateral" % (2 * k, 2 * k + 1))
        quads[k] = nodes
    return _ccw_quads(np.asarray(x, float), np.asarray(y, float), quads)


def find_quads(ikle: np.ndarray, quads: np.ndarray):
    """Which quads still exist as two triangles of ``ikle`` (either diagonal).

    Returns (quads found, bool mask of the triangles they use).
    """
    ik = np.asarray(ikle, dtype=np.int64)
    index = {tuple(sorted(t)): i for i, t in enumerate(ik.tolist())}
    used = np.zeros(len(ik), dtype=bool)
    found = []
    for q in np.asarray(quads, dtype=np.int64).reshape(-1, 4).tolist():
        a, b, c, d = q
        for t1, t2 in (((a, b, c), (a, c, d)), ((a, b, d), (b, c, d))):
            i = index.get(tuple(sorted(t1)))
            j = index.get(tuple(sorted(t2)))
            if i is not None and j is not None and not used[i] and not used[j]:
                used[i] = used[j] = True
                found.append(q)
                break
    return np.array(found, dtype=np.int64).reshape(-1, 4), used


def split_quads(x, y, quads: np.ndarray) -> np.ndarray:
    """Two triangles per quad, cut along the shorter diagonal, counter-clockwise."""
    q = np.asarray(quads, dtype=np.int64).reshape(-1, 4)
    if len(q) == 0:
        return np.zeros((0, 3), dtype=np.int64)
    d1 = np.hypot(x[q[:, 0]] - x[q[:, 2]], y[q[:, 0]] - y[q[:, 2]])
    d2 = np.hypot(x[q[:, 1]] - x[q[:, 3]], y[q[:, 1]] - y[q[:, 3]])
    a = np.where((d1 <= d2)[:, None], q[:, [0, 1, 2]], q[:, [0, 1, 3]])
    b = np.where((d1 <= d2)[:, None], q[:, [0, 2, 3]], q[:, [1, 2, 3]])
    t = np.vstack([a, b])
    neg = signed_areas(x, y, t) < 0
    t[neg] = t[neg][:, [0, 2, 1]]
    return t


@dataclass
class Mesh:
    x: np.ndarray
    y: np.ndarray
    triangles: np.ndarray
    quads: np.ndarray = field(default_factory=lambda: np.zeros((0, 4), dtype=np.int64))
    z: Optional[np.ndarray] = None  # bottom elevation at nodes
    node_vars: Dict[str, Tuple[str, np.ndarray]] = field(default_factory=dict)  # name -> (unit, values)
    material: Optional[np.ndarray] = None  # per element: triangles first, then quads
    material_names: Dict[int, str] = field(default_factory=dict)
    boundaries: List[BoundaryString] = field(default_factory=list)
    crs_wkt: str = ""
    frame: Frame = field(default_factory=Frame)
    title: str = "TELEMAC Toolkit mesh"

    # ------------------------------------------------------------------ build
    @classmethod
    def from_result(cls, res, structured_quads: Optional[np.ndarray] = None, keep_quads: bool = False) -> "Mesh":
        tri = np.asarray(res.ikle, dtype=np.int64)
        quads = np.zeros((0, 4), dtype=np.int64)
        if keep_quads and structured_quads is not None and len(structured_quads):
            quads, used = find_quads(tri, structured_quads)
            tri = tri[~used]
        return cls(np.asarray(res.x, float), np.asarray(res.y, float), tri, _ccw_quads(res.x, res.y, quads))

    # ------------------------------------------------------------------ queries
    @property
    def npoin(self) -> int:
        return len(self.x)

    @property
    def nelem(self) -> int:
        return len(self.triangles) + len(self.quads)

    @property
    def has_quads(self) -> bool:
        return len(self.quads) > 0

    def all_triangles(self) -> np.ndarray:
        if not self.has_quads:
            return self.triangles
        return np.vstack([self.triangles, split_quads(self.x, self.y, self.quads)])

    def element_centroids(self):
        cx = np.concatenate([self.x[self.triangles].mean(axis=1), self.x[self.quads].mean(axis=1)]) if self.has_quads else self.x[self.triangles].mean(axis=1)
        cy = np.concatenate([self.y[self.triangles].mean(axis=1), self.y[self.quads].mean(axis=1)]) if self.has_quads else self.y[self.triangles].mean(axis=1)
        return cx, cy

    def bottom(self) -> np.ndarray:
        return self.z if self.z is not None else np.zeros(self.npoin)

    def boundary_order(self) -> np.ndarray:
        loops = boundary_loops(self.all_triangles(), self.x, self.y)
        return np.concatenate(loops) if loops else np.zeros(0, dtype=np.int64)

    def boundary_loops(self) -> List[np.ndarray]:
        return boundary_loops(self.all_triangles(), self.x, self.y)

    # ------------------------------------------------------------------ edits
    def assign_zones(self, zones: Sequence[Tuple[List[np.ndarray], int, Optional[float], str]],
                     default_material: int = 1, friction_name: str = "BOTTOM FRICTION",
                     default_friction: Optional[float] = None) -> List[str]:
        """Material IDs per element and a friction variable at nodes from polygons.

        ``zones`` = [(rings, material id, friction or None, name), ...]; later
        zones win where they overlap. Returns report lines.
        """
        cx, cy = self.element_centroids()
        mat = np.full(self.nelem, int(default_material), dtype=np.int64)
        fric = np.full(self.npoin, np.nan if default_friction is None else float(default_friction))
        report = []
        names = {int(default_material): "default"}
        any_friction = default_friction is not None
        for rings, mid, f, name in zones:
            inside = points_in_rings(cx, cy, rings)
            mat[inside] = int(mid)
            names[int(mid)] = name or names.get(int(mid), "zone %d" % mid)
            if f is not None:
                any_friction = True
                fric[points_in_rings(self.x, self.y, rings)] = float(f)
            report.append("Zone %s: material %d, %d elements%s" % (name or mid, mid, int(inside.sum()),
                                                                  ", friction %g" % f if f is not None else ""))
        self.material = mat
        self.material_names = names
        if any_friction:
            if np.isnan(fric).any():
                miss = int(np.isnan(fric).sum())
                fill = np.nanmedian(fric) if np.isfinite(fric).any() else 0.0
                fric = np.where(np.isnan(fric), fill, fric)
                report.append("%d nodes outside every friction zone take %g" % (miss, fill))
            self.node_vars[friction_name] = ("", fric)
        return report

    def transformed(self, frame: Frame) -> "Mesh":
        """Copy of the mesh in a local frame (translation, then rotation)."""
        u, v = frame.forward(self.x, self.y)
        m = Mesh(u, v, self.triangles.copy(), self.quads.copy(), None if self.z is None else self.z.copy(),
                 dict(self.node_vars), None if self.material is None else self.material.copy(), dict(self.material_names),
                 list(self.boundaries), self.crs_wkt, frame, self.title)
        if not frame.is_identity:
            m.crs_wkt = ""  # local coordinates: no map CRS unless the writer stores the offset
        return m

    def summary(self) -> str:
        s = "%d nodes, %d elements" % (self.npoin, self.nelem)
        if self.has_quads:
            s += " (%d triangles, %d quadrilaterals)" % (len(self.triangles), len(self.quads))
        return s
