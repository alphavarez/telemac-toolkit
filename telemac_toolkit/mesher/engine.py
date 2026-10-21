"""Unstructured triangle mesh generator with structured parts.

1. Domain rings and guide lines are resampled to the local element size;
   breaklines keep every vertex (long segments are subdivided); line
   intersections become nodes.
2. Fixed points, breakline nodes and the boundary nodes of structured parts
   are added. Where nodes of different objects overlap, the higher priority
   one is kept: structured part > fixed point > breakline > domain > guide line.
3. The inside is filled with nodes on hexagonal lattices, one per size level,
   finest first, spaced by the size field.
4. A Delaunay triangulation connects all nodes (ghost nodes outside the domain
   keep it well conditioned); triangles outside the domain, in islands or in
   structured parts are removed. Interior nodes are kept clear of every
   constraint segment so each one is a Delaunay edge; a segment still missing is
   split and the step repeated.
5. Oversized triangles get a node at their circumcentre, interior nodes are
   smoothed (Laplacian, never inverting a triangle), structured parts are
   merged in, and triangles with three boundary nodes are fixed by an edge flip.

The Delaunay step is supplied by the caller (GEOS through QGIS, or SciPy), so
this module only needs NumPy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace as dc_replace
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..telemac_io.mesh import boundary_loops, mesh_statistics, signed_areas
from .geometry import as_ring, distance_to_segments, points_in_rings, segment_intersections, turning_angles
from .sizefield import LineDensity, SizeField

__all__ = ["MeshInput", "MeshResult", "generate", "NODE_TYPES", "PointGrid"]

# node types (internal names: OUTLINE = domain boundary, HARD* = breakline, SOFT* = guide line)
INTERIOR, OUTLINE, SUBMESH, HARDPOINT, HARDLINE, HARDISLAND, SOFTLINE, SOFTISLAND = range(8)
NODE_TYPES = {
    INTERIOR: "interior", OUTLINE: "domain boundary", SUBMESH: "structured part", HARDPOINT: "fixed point",
    HARDLINE: "breakline", HARDISLAND: "closed breakline or island", SOFTLINE: "guide line", SOFTISLAND: "closed guide line",
}
import os as _os
_DEBUG = bool(_os.environ.get("TT_MESH_DEBUG"))
GHOST = 8  # helper nodes outside the domain; never part of the result
_PRIORITY = {SUBMESH: 0, HARDPOINT: 1, HARDLINE: 2, HARDISLAND: 2, OUTLINE: 3, SOFTISLAND: 4, SOFTLINE: 4, INTERIOR: 5, GHOST: 6}

Triangulator = Callable[[np.ndarray], np.ndarray]
Progress = Optional[Callable[[float, str], bool]]


@dataclass
class Submesh:
    """A ready-made mesh stitched in node for node (structured channel strip, regular grid)."""

    x: np.ndarray
    y: np.ndarray
    ikle: np.ndarray
    quads: Optional[np.ndarray] = None  # cells as quadrilaterals (node ids), for exporters that keep them


@dataclass
class MeshInput:
    """Everything the generator needs, in one CRS.

    domain: polygons to mesh, each a list of rings [outer, hole, ...]; the
    caller has already removed islands and submesh areas (GEOS difference).
    """

    domain: List[List[np.ndarray]]
    size: SizeField
    soft_lines: List[Tuple[np.ndarray, bool]] = field(default_factory=list)
    hard_lines: List[Tuple[np.ndarray, bool]] = field(default_factory=list)
    hard_points: Optional[np.ndarray] = None
    submeshes: List[Submesh] = field(default_factory=list)
    resample_outline: bool = True
    smooth_iterations: int = 8
    submesh_density: bool = True  # submesh boundaries act as density lines (their own edge length)
    corner_angle: float = 30.0  # outline/soft-line vertices turning more than this are kept
    clearance: float = 0.55  # interior nodes stay this fraction of the size away from constraints
    spacing: float = 0.80  # minimum distance between interior nodes, fraction of the size


@dataclass
class MeshResult:
    x: np.ndarray
    y: np.ndarray
    ikle: np.ndarray
    node_type: np.ndarray
    density: np.ndarray
    stats: dict
    log: List[str]
    quads: np.ndarray = field(default_factory=lambda: np.zeros((0, 4), dtype=np.int64))  # structured cells


# --------------------------------------------------------------------------
# spatial grid for neighbour queries
# --------------------------------------------------------------------------
class PointGrid:
    def __init__(self, xy: np.ndarray, cell: float):
        self.xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        self.cell = float(cell)
        if len(self.xy) == 0:
            self.x0 = self.y0 = 0.0
            self.nx = self.ny = 1
        else:
            self.x0, self.y0 = self.xy.min(axis=0) - self.cell
            span = self.xy.max(axis=0) + self.cell - (self.x0, self.y0)
            self.nx = max(1, int(span[0] // self.cell) + 2)
            self.ny = max(1, int(span[1] // self.cell) + 2)
        ix, iy = self._ij(self.xy[:, 0], self.xy[:, 1])
        key = iy * self.nx + ix
        self.order = np.argsort(key, kind="stable")
        self.keys = key[self.order]

    def _ij(self, x, y):
        ix = np.clip(((np.asarray(x) - self.x0) // self.cell).astype(np.int64), 0, self.nx - 1)
        iy = np.clip(((np.asarray(y) - self.y0) // self.cell).astype(np.int64), 0, self.ny - 1)
        return ix, iy

    def query_rect(self, xmin, ymin, xmax, ymax) -> np.ndarray:
        if len(self.xy) == 0:
            return np.zeros(0, dtype=np.int64)
        ix0, iy0 = self._ij(xmin, ymin)
        ix1, iy1 = self._ij(xmax, ymax)
        parts = []
        for iy in range(int(iy0), int(iy1) + 1):
            a = np.searchsorted(self.keys, iy * self.nx + ix0, side="left")
            b = np.searchsorted(self.keys, iy * self.nx + ix1, side="right")
            if b > a:
                parts.append(self.order[a:b])
        return np.concatenate(parts) if parts else np.zeros(0, dtype=np.int64)

    def min_dist_ratio(self, px, py, rad_p, rad_tab) -> np.ndarray:
        """For each query point, min over grid points of d / (0.5 (rad_p + rad_q)).

        Values below 1 mean a neighbour closer than the allowed spacing.
        """
        px, py, rad_p = (np.asarray(v, dtype=np.float64) for v in (px, py, rad_p))
        out = np.full(px.shape, np.inf)
        if len(self.xy) == 0 or px.size == 0:
            return out
        ix, iy = self._ij(px, py)
        starts = np.searchsorted(self.keys, np.arange(self.nx * self.ny), side="left")
        ends = np.searchsorted(self.keys, np.arange(self.nx * self.ny), side="right")
        occ = ends - starts
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                jx, jy = ix + dx, iy + dy
                valid = (jx >= 0) & (jx < self.nx) & (jy >= 0) & (jy < self.ny)
                key = np.where(valid, jy * self.nx + jx, 0)
                n = np.where(valid, occ[key], 0)
                if n.max(initial=0) == 0:
                    continue
                for k in range(int(n.max())):
                    has = n > k
                    if not has.any():
                        break
                    idx = self.order[starts[key[has]] + k]
                    d = np.hypot(px[has] - self.xy[idx, 0], py[has] - self.xy[idx, 1])
                    r = d / (0.5 * (rad_p[has] + rad_tab[idx]))
                    out[has] = np.minimum(out[has], r)
        return out


# --------------------------------------------------------------------------
# constraint lines
# --------------------------------------------------------------------------
@dataclass
class _Line:
    coords: np.ndarray
    closed: bool
    kind: int  # node type for its nodes
    resample: bool
    fixed: np.ndarray = None  # bool per vertex: never moved/removed by resampling

    def __post_init__(self):
        if self.fixed is None:
            self.fixed = np.zeros(len(self.coords), dtype=bool)
            if not self.closed:
                self.fixed[0] = self.fixed[-1] = True


def _insert_points(line: _Line, inserts: List[Tuple[int, float, float, float]]):
    """Insert (segment, t, x, y) points into a line as fixed vertices."""
    if not inserts:
        return
    c = line.coords
    n = len(c)
    nseg = n if line.closed else n - 1
    by_seg: Dict[int, list] = {}
    for seg, t, x, y in inserts:
        by_seg.setdefault(seg, []).append((t, x, y))
    new_c, new_f = [], []
    for i in range(n):
        new_c.append(c[i])
        new_f.append(line.fixed[i])
        if i < nseg and i in by_seg:
            for t, x, y in sorted(by_seg[i]):
                if t <= 1e-9:
                    new_f[-1] = True
                    continue
                if t >= 1 - 1e-9:
                    continue  # marks the next vertex; handled below
                new_c.append(np.array([x, y]))
                new_f.append(True)
    # endpoints hit exactly (t == 1) become fixed
    for seg, t, x, y in inserts:
        if t >= 1 - 1e-9:
            j = (seg + 1) % n
            for k, p in enumerate(new_c):
                if abs(p[0] - c[j][0]) < 1e-12 and abs(p[1] - c[j][1]) < 1e-12:
                    new_f[k] = True
    line.coords = np.array(new_c)
    line.fixed = np.array(new_f, dtype=bool)


def _segments(line: _Line):
    c = line.coords
    if line.closed:
        return c, np.roll(c, -1, axis=0)
    return c[:-1], c[1:]


def _resample_piece(pts: np.ndarray, size: SizeField) -> np.ndarray:
    """Interior nodes for one piece between two fixed vertices (ends excluded)."""
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    total = float(seglen.sum())
    if total <= 0:
        return np.zeros((0, 2))
    smin = float(np.min(size(pts[:, 0], pts[:, 1])))
    step = max(min(smin, total) / 6.0, total / 20000.0)
    s = np.unique(np.concatenate([np.arange(0.0, total, step), np.cumsum(np.concatenate([[0], seglen])), [total]]))
    cum = np.concatenate([[0], np.cumsum(seglen)])
    x = np.interp(s, cum, pts[:, 0])
    y = np.interp(s, cum, pts[:, 1])
    h = size(x, y)
    inv = 1.0 / np.maximum(h, 1e-12)
    T = np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(s))])
    nseg = max(1, int(round(T[-1])))
    if nseg == 1:
        return np.zeros((0, 2))
    targets = np.linspace(0, T[-1], nseg + 1)[1:-1]
    sp = np.interp(targets, T, s)
    return np.stack([np.interp(sp, cum, pts[:, 0]), np.interp(sp, cum, pts[:, 1])], axis=1)


def _resample_line(line: _Line, size: SizeField, corner_angle: float, frozen: set) -> _Line:
    c = line.coords
    if not line.resample or len(c) < 2:
        return line
    fixed = line.fixed.copy()
    ang = turning_angles(c, line.closed)
    fixed |= ang > corner_angle
    if line.closed and not fixed.any():
        fixed[0] = True
    n = len(c)
    idx_fixed = np.flatnonzero(fixed)
    out_c, out_f = [], []
    if line.closed:
        order = list(idx_fixed) + [idx_fixed[0] + n]
    else:
        order = list(idx_fixed)
    for a, b in zip(order[:-1], order[1:]):
        piece_idx = [k % n for k in range(a, b + 1)]
        pts = c[piece_idx]
        out_c.append(pts[0])
        out_f.append(True)
        key = (tuple(np.round(pts[0], 9)), tuple(np.round(pts[-1], 9)))
        if len(pts) == 2 and (key in frozen or key[::-1] in frozen):
            continue
        for p in _resample_piece(pts, size):
            out_c.append(p)
            out_f.append(False)
    if not line.closed:
        out_c.append(c[order[-1]])
        out_f.append(True)
    return _Line(np.array(out_c), line.closed, line.kind, False, np.array(out_f, dtype=bool))


# --------------------------------------------------------------------------
# interior nodes
# --------------------------------------------------------------------------
def _interior_candidates(inp: MeshInput, rings, bbox, progress: Progress, log: List[str]):
    size = inp.size
    hmin = size.minimum
    hmax = size.default
    ratio = 1.3
    levels = [hmin]
    while levels[-1] * ratio < hmax:
        levels.append(levels[-1] * ratio)
    xmin, ymin, xmax, ymax = bbox
    g = max(size.growth - 1.0, 0.0)
    out_xy, out_h = [], []
    for li, h in enumerate(levels):
        lo = h if li > 0 else 0.0
        hi = h * ratio if li < len(levels) - 1 else np.inf
        # coarse screening grid: skip cells where the size cannot fall in [lo, hi)
        cell = max(h * 12.0, max(xmax - xmin, ymax - ymin) / 300.0)
        cx = np.arange(xmin, xmax + cell, cell)
        cy = np.arange(ymin, ymax + cell, cell)
        if cx.size * cy.size == 0:
            continue
        CX, CY = np.meshgrid(cx + cell / 2, cy + cell / 2)
        sc = size(CX.ravel(), CY.ravel())
        slack = g * cell * 0.75 + 1e-9
        keep = (sc - slack < hi) & (sc + slack >= lo)
        if not keep.any():
            continue
        spacing = h * math.sqrt(ratio)  # lattice pitch in the middle of the band
        cand = []
        for k in np.flatnonzero(keep):
            x0, y0 = CX.ravel()[k] - cell / 2, CY.ravel()[k] - cell / 2
            # a global lattice (same origin everywhere) so neighbouring cells match
            ox = xmin + math.floor((x0 - xmin) / spacing) * spacing
            oy_step = spacing * math.sqrt(3) / 2
            ky = math.floor((y0 - ymin) / oy_step)
            ys = ymin + np.arange(ky, ky + int(cell / oy_step) + 2) * oy_step
            for j, yy in zip(range(ky, ky + len(ys)), ys):
                off = spacing / 2 if (j % 2) else 0.0
                xs = np.arange(ox + off, x0 + cell + spacing, spacing)
                xs = xs[(xs >= x0) & (xs < x0 + cell)]
                if yy >= y0 and yy < y0 + cell and xs.size:
                    cand.append(np.stack([xs, np.full(xs.size, yy)], axis=1))
        if not cand:
            continue
        cand = np.vstack(cand)
        hs = size(cand[:, 0], cand[:, 1])
        band = (hs >= lo) & (hs < hi)
        cand, hs = cand[band], hs[band]
        if cand.size == 0:
            continue
        inside = points_in_rings(cand[:, 0], cand[:, 1], rings)
        out_xy.append(cand[inside])
        out_h.append(hs[inside])
        if progress and progress(0.15 + 0.25 * (li + 1) / len(levels), "Placing nodes (level %d/%d)" % (li + 1, len(levels))):
            raise RuntimeError("cancelled")
    if not out_xy:
        return np.zeros((0, 2)), np.zeros(0), levels
    log.append("Interior candidates: %d on %d size levels from %.3g to %.3g" % (sum(len(a) for a in out_xy), len(levels), hmin, levels[-1]))
    return out_xy, out_h, levels


def _filter_interior(cand_levels, h_levels, cons_xy, cons_h, seg_a, seg_b, inp: MeshInput):
    """Greedy acceptance: finest level first, keep distance to accepted nodes and to constraints."""
    acc_xy = cons_xy.copy()
    acc_h = cons_h.copy()
    n_cons = len(cons_xy)
    for cand, hs in zip(cand_levels, h_levels):
        if len(cand) == 0:
            continue
        # clearance from constraint segments (keeps them Delaunay edges)
        reach = float(np.max(hs)) * inp.clearance
        grid = PointGrid(cand, max(reach, 1e-9))
        too_close = np.zeros(len(cand), dtype=bool)
        for a, b in zip(seg_a, seg_b):
            ids = grid.query_rect(min(a[0], b[0]) - reach, min(a[1], b[1]) - reach,
                                  max(a[0], b[0]) + reach, max(a[1], b[1]) + reach)
            if ids.size == 0:
                continue
            d = distance_to_segments(cand[ids, 0], cand[ids, 1], a[None, :], b[None, :])
            seglen = math.hypot(b[0] - a[0], b[1] - a[1])
            lim = np.maximum(inp.clearance * hs[ids], 0.5 * seglen * 1.02)
            too_close[ids[d < lim]] = True
        cand, hs = cand[~too_close], hs[~too_close]
        if len(cand) == 0:
            continue
        # spacing against already accepted nodes (finer levels and constraints)
        rmax = float(max(np.max(hs), np.max(acc_h) if len(acc_h) else 0.0))
        agrid = PointGrid(acc_xy, rmax * inp.spacing + 1e-9)
        ratio = agrid.min_dist_ratio(cand[:, 0], cand[:, 1], hs, acc_h)
        ok = ratio >= inp.spacing
        acc_xy = np.vstack([acc_xy, cand[ok]])
        acc_h = np.concatenate([acc_h, hs[ok]])
    return acc_xy[n_cons:], acc_h[n_cons:]


# --------------------------------------------------------------------------
# helpers on triangulations
# --------------------------------------------------------------------------
def _edge_set(ikle: np.ndarray) -> set:
    e = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    e.sort(axis=1)
    return set(map(tuple, e.tolist()))


def _orient_ccw(x, y, ikle):
    a = signed_areas(x, y, ikle)
    ikle = ikle.copy()
    neg = a < 0
    ikle[neg] = ikle[neg][:, [0, 2, 1]]
    return ikle


def _smooth(x, y, ikle, free: np.ndarray, iterations: int):
    """Laplacian smoothing of free nodes; a move that inverts a triangle is undone."""
    if iterations <= 0 or not free.any():
        return x, y
    e = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    n = x.size
    deg = np.bincount(e.ravel(), minlength=n).astype(float)
    x, y = x.copy(), y.copy()
    base_area = np.abs(signed_areas(x, y, ikle))
    for _ in range(iterations):
        sx = np.bincount(e[:, 0], weights=x[e[:, 1]], minlength=n) + np.bincount(e[:, 1], weights=x[e[:, 0]], minlength=n)
        sy = np.bincount(e[:, 0], weights=y[e[:, 1]], minlength=n) + np.bincount(e[:, 1], weights=y[e[:, 0]], minlength=n)
        nx_ = np.where(free & (deg > 0), sx / np.maximum(deg, 1), x)
        ny_ = np.where(free & (deg > 0), sy / np.maximum(deg, 1), y)
        ox, oy = x, y
        x, y = nx_, ny_
        for _guard in range(5):
            a = signed_areas(x, y, ikle)
            bad = a <= 1e-3 * base_area
            if not bad.any():
                break
            nodes = np.unique(ikle[bad].ravel())
            x[nodes], y[nodes] = ox[nodes], oy[nodes]
    return x, y


def _fix_overconstrained(x, y, ikle, on_boundary: np.ndarray, constrained_edges: set):
    """Flip the interior edge of triangles whose three nodes are on the boundary."""
    ikle = ikle.copy()
    bad = np.flatnonzero(on_boundary[ikle].all(axis=1))
    if bad.size == 0:
        return ikle, 0
    edge_map: Dict[Tuple[int, int], List[int]] = {}
    for t, tri in enumerate(ikle):
        for k in range(3):
            a, b = int(tri[k]), int(tri[(k + 1) % 3])
            edge_map.setdefault((min(a, b), max(a, b)), []).append(t)
    fixed = 0
    for t in bad:
        tri = ikle[t]
        if not on_boundary[tri].all():
            continue
        for k in range(3):
            a, b, c = int(tri[k]), int(tri[(k + 1) % 3]), int(tri[(k + 2) % 3])
            key = (min(a, b), max(a, b))
            owners = edge_map.get(key, [])
            if len(owners) != 2 or key in constrained_edges:
                continue
            u = owners[0] if owners[1] == t else owners[1]
            d = next(int(v) for v in ikle[u] if v not in (a, b))
            # quad a-d-b-c must be convex: new triangles (c, a, d) and (c, d, b) both positive
            t1 = np.array([[c, a, d]])
            t2 = np.array([[c, d, b]])
            if signed_areas(x, y, t1)[0] > 0 and signed_areas(x, y, t2)[0] > 0:
                ikle[t] = t1[0]
                ikle[u] = t2[0]
                fixed += 1
                edge_map = {}
                for tt, tri2 in enumerate(ikle):
                    for kk in range(3):
                        p, q = int(tri2[kk]), int(tri2[(kk + 1) % 3])
                        edge_map.setdefault((min(p, q), max(p, q)), []).append(tt)
                break
    return ikle, fixed


def _clear_of_segments(px, py, hh, sa, sb, clearance):
    """True where a point is at least clearance * hh away from every segment sa-sb."""
    if len(sa) == 0 or px.size == 0:
        return np.ones(px.shape, dtype=bool)
    reach = float(np.max(hh)) * clearance
    dd = np.full(px.shape, np.inf)
    grid = PointGrid(np.stack([px, py], axis=1), max(reach, 1e-9))
    for a, b in zip(sa, sb):
        ids = grid.query_rect(min(a[0], b[0]) - reach, min(a[1], b[1]) - reach,
                              max(a[0], b[0]) + reach, max(a[1], b[1]) + reach)
        if ids.size:
            dd[ids] = np.minimum(dd[ids], distance_to_segments(px[ids], py[ids], a[None, :], b[None, :]))
    return dd >= clearance * hh


def _thin(px, py, h, frac=0.5):
    """Greedy thinning: keep points at least frac * h apart (coarsest first)."""
    if px.size == 0:
        return np.zeros(0, dtype=np.int64)
    order = np.argsort(-h, kind="stable")
    grid = PointGrid(np.stack([px, py], axis=1), float(np.max(h)) * frac + 1e-9)
    taken = np.zeros(px.size, dtype=bool)
    acc = []
    for i in order:
        r = frac * h[i]
        ids = grid.query_rect(px[i] - r, py[i] - r, px[i] + r, py[i] + r)
        if np.any(taken[ids] & (np.hypot(px[ids] - px[i], py[ids] - py[i]) < r)):
            continue
        taken[i] = True
        acc.append(i)
    return np.array(acc, dtype=np.int64)


def _ghost_nodes(xy, seg_pairs, size: SizeField, clearance: float, rings) -> np.ndarray:
    """Apex of an equilateral triangle on the outer side of each boundary segment."""
    sp = np.asarray(seg_pairs, dtype=np.int64).reshape(-1, 2)
    if len(sp) == 0:
        return np.zeros((0, 2))
    a, b = xy[sp[:, 0]], xy[sp[:, 1]]
    m = 0.5 * (a + b)
    d = b - a
    L = np.hypot(d[:, 0], d[:, 1])
    ok = L > 0
    n = np.zeros_like(d)
    n[ok] = np.stack([-d[ok, 1], d[ok, 0]], axis=1) / L[ok, None]
    cand = np.vstack([m + n * (0.866 * L)[:, None], m - n * (0.866 * L)[:, None]])
    hh = np.concatenate([L, L])
    keep = np.concatenate([ok, ok]) & ~points_in_rings(cand[:, 0], cand[:, 1], rings)
    cand, hh = cand[keep], hh[keep]
    if len(cand) == 0:
        return np.zeros((0, 2))
    clear = _clear_of_segments(cand[:, 0], cand[:, 1], hh, a, b, clearance)
    cand, hh = cand[clear], hh[clear]
    idx = _thin(cand[:, 0], cand[:, 1], hh, 0.6)
    return cand[idx]


def _gap_nodes(xy, ikle, seg_list, size: SizeField, clearance: float, rings, limit: float = 0.85) -> np.ndarray:
    """New interior nodes for triangles whose circumradius exceeds ``limit`` times the local size.

    An equilateral triangle with the target edge length has a circumradius of
    0.577 h, so 0.85 h only catches triangles with edges around 1.5 h. The
    circumcentre is used when it lies inside the domain and clear of every
    constraint segment, else the centroid if that is clear. Candidates closer
    than half the local size to each other are thinned.
    """
    if len(ikle) == 0:
        return np.zeros((0, 2))
    P = xy[ikle]
    ax, ay = P[:, 0, 0], P[:, 0, 1]
    bx, by = P[:, 1, 0], P[:, 1, 1]
    cx, cy = P[:, 2, 0], P[:, 2, 1]
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ok = np.abs(d) > 1e-300
    d = np.where(ok, d, 1.0)
    a2, b2, c2 = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    R = np.hypot(ux - ax, uy - ay)
    gx, gy = P[:, :, 0].mean(axis=1), P[:, :, 1].mean(axis=1)
    h = size(gx, gy)
    big = ok & (R > limit * h)
    if not big.any():
        return np.zeros((0, 2))
    idx = np.flatnonzero(big)
    ux, uy, gx, gy, h = ux[idx], uy[idx], gx[idx], gy[idx], h[idx]
    sp = np.array(seg_list, dtype=np.int64).reshape(-1, 2)
    sa, sb = xy[sp[:, 0]], xy[sp[:, 1]]

    def clear(px, py, hh):
        if len(sa) == 0:
            return np.ones(px.shape, dtype=bool)
        reach = float(np.max(hh)) * clearance
        dd = np.full(px.shape, np.inf)
        grid = PointGrid(np.stack([px, py], axis=1), max(reach, 1e-9))
        for a, b in zip(sa, sb):
            ids = grid.query_rect(min(a[0], b[0]) - reach, min(a[1], b[1]) - reach,
                                  max(a[0], b[0]) + reach, max(a[1], b[1]) + reach)
            if ids.size:
                dd[ids] = np.minimum(dd[ids], distance_to_segments(px[ids], py[ids], a[None, :], b[None, :]))
        return dd >= clearance * hh

    use_c = np.isfinite(ux) & np.isfinite(uy)
    use_c[use_c] = points_in_rings(ux[use_c], uy[use_c], rings)
    if use_c.any():
        use_c[use_c] = clear(ux[use_c], uy[use_c], h[use_c])
    px = np.where(use_c, ux, gx)
    py = np.where(use_c, uy, gy)
    keep = use_c | clear(gx, gy, h)
    px, py, h = px[keep], py[keep], h[keep]
    if px.size == 0:
        return np.zeros((0, 2))
    # thin candidates that crowd each other (neighbouring big triangles)
    order = np.argsort(-h)
    acc = []
    grid_pts = PointGrid(np.stack([px, py], axis=1), float(np.max(h)) * 0.5 + 1e-9)
    taken = np.zeros(px.size, dtype=bool)
    for i in order:
        r = 0.5 * h[i]
        ids = grid_pts.query_rect(px[i] - r, py[i] - r, px[i] + r, py[i] + r)
        if np.any(taken[ids] & (np.hypot(px[ids] - px[i], py[ids] - py[i]) < r)):
            continue
        taken[i] = True
        acc.append(i)
    return np.stack([px[acc], py[acc]], axis=1)


def min_angles(x, y, ikle) -> np.ndarray:
    """Smallest interior angle of each triangle, degrees."""
    P = np.stack([x[ikle], y[ikle]], axis=2)
    out = np.full(len(ikle), 180.0)
    for i in range(3):
        a = P[:, (i + 1) % 3] - P[:, i]
        b = P[:, (i + 2) % 3] - P[:, i]
        c = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-300)
        out = np.minimum(out, np.degrees(np.arccos(np.clip(c, -1, 1))))
    return out


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------
def _merge_only(subs: Sequence[Submesh], log: List[str]) -> MeshResult:
    """The mesh when structured parts cover the whole domain: their nodes merged, nothing to triangulate."""
    allx = np.concatenate([s.x for s in subs])
    ally = np.concatenate([s.y for s in subs])
    ext = max(float(np.ptp(allx)), float(np.ptp(ally)), 1.0)
    snap = max(ext * 1e-9, 1e-9)
    key: Dict[Tuple[int, int], int] = {}
    xs, ys, tris, quads = [], [], [], []
    for s in subs:
        loc = np.empty(len(s.x), dtype=np.int64)
        for k, (px, py) in enumerate(zip(s.x, s.y)):
            kk = (round(float(px) / snap), round(float(py) / snap))
            if kk not in key:
                key[kk] = len(xs)
                xs.append(px)
                ys.append(py)
            loc[k] = key[kk]
        tris.append(loc[np.asarray(s.ikle, dtype=np.int64)])
        if s.quads is not None and len(s.quads):
            quads.append(loc[np.asarray(s.quads, dtype=np.int64)])
    x, y = np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64)
    ikle = _orient_ccw(x, y, np.vstack(tris))
    stats = mesh_statistics(x, y, ikle)
    amin = min_angles(x, y, ikle)
    stats["share_min_angle_below_20"] = float(np.mean(amin < 20.0)) if amin.size else 0.0
    stats["share_min_angle_below_30"] = float(np.mean(amin < 30.0)) if amin.size else 0.0
    stats["min_angle_p05"] = float(np.percentile(amin, 5)) if amin.size else 0.0
    log.append("Structured parts cover the whole domain: %d nodes, %d triangles" % (len(x), len(ikle)))
    q = np.vstack(quads) if quads else np.zeros((0, 4), dtype=np.int64)
    return MeshResult(x, y, ikle, np.full(len(x), SUBMESH, dtype=np.int64), np.zeros(len(x)), stats, log, q)


def generate(inp: MeshInput, triangulate: Triangulator, progress: Progress = None) -> MeshResult:
    log: List[str] = []
    if not inp.domain:
        if not inp.submeshes:
            raise ValueError("the domain is empty")
        return _merge_only(inp.submeshes, log)
    if inp.submeshes and inp.submesh_density:
        extra = []
        for sm in inp.submeshes:
            for loop in boundary_loops(sm.ikle, sm.x, sm.y):
                c = np.stack([sm.x[loop], sm.y[loop]], axis=1)
                e = np.hypot(*(np.roll(c, -1, axis=0) - c).T)
                extra.append(LineDensity(c, float(np.median(e)), closed=True))
        s0 = inp.size
        inp = dc_replace(inp, size=SizeField(s0.default, s0.growth, list(s0.polygons), list(s0.points), list(s0.lines) + extra))
    size = inp.size

    def step(frac, msg):
        if progress and progress(frac, msg):
            raise RuntimeError("cancelled")

    # ---- domain rings ---------------------------------------------------
    rings = [as_ring(r) for poly in inp.domain for r in poly]
    rings = [r for r in rings if len(r) >= 3]
    if not rings:
        raise ValueError("the outline is empty")
    allr = np.vstack(rings)
    bbox = (allr[:, 0].min(), allr[:, 1].min(), allr[:, 0].max(), allr[:, 1].max())
    extent = max(bbox[2] - bbox[0], bbox[3] - bbox[1])
    snap = max(extent * 1e-9, 1e-6)

    # ---- submeshes: fixed boundary nodes and frozen edges ----------------
    sub_nodes = []
    frozen = set()
    for sm in inp.submeshes:
        for loop in boundary_loops(sm.ikle, sm.x, sm.y):
            pts = np.stack([sm.x[loop], sm.y[loop]], axis=1)
            sub_nodes.append(pts)
            for p, q in zip(pts, np.roll(pts, -1, axis=0)):
                frozen.add((tuple(np.round(p, 9)), tuple(np.round(q, 9))))
    sub_xy = np.vstack(sub_nodes) if sub_nodes else np.zeros((0, 2))
    if len(sub_xy):
        # only submesh nodes on the domain rings take part in the triangulation; the
        # others (inside the submesh hole, or where a submesh sticks out of the
        # outline) would spoil the Delaunay edges of nearby constraint segments
        ring_keys = set(map(tuple, np.round(np.vstack([as_ring(r) for poly in inp.domain for r in poly]), 9).tolist()))
        on_ring = np.array([tuple(p) in ring_keys for p in np.round(sub_xy, 9).tolist()], dtype=bool)
        sub_xy = sub_xy[on_ring]

    # ---- constraint lines -------------------------------------------------
    lines: List[_Line] = []
    for poly in inp.domain:
        for k, r in enumerate(poly):
            rr = as_ring(r)
            if len(rr) >= 3:
                lines.append(_Line(rr, True, OUTLINE if k == 0 else HARDISLAND, inp.resample_outline))
    for coords, closed in inp.soft_lines:
        c = as_ring(coords) if closed else np.asarray(coords, float)[:, :2]
        if len(c) >= 2:
            lines.append(_Line(c, closed, SOFTISLAND if closed else SOFTLINE, True))
    for coords, closed in inp.hard_lines:
        c = as_ring(coords) if closed else np.asarray(coords, float)[:, :2]
        if len(c) >= 2:
            hl = _Line(c, closed, HARDISLAND if closed else HARDLINE, True)
            hl.fixed[:] = True  # every vertex stays a node; long segments are subdivided
            lines.append(hl)
    # submesh nodes lying on outline rings are fixed vertices there
    if len(sub_xy):
        key_sub = set(map(tuple, np.round(sub_xy, 9).tolist()))
        for ln in lines:
            m = np.array([tuple(p) in key_sub for p in np.round(ln.coords, 9).tolist()])
            ln.fixed |= m
    step(0.03, "Intersecting lines")

    # intersections between all constraint lines become fixed nodes
    segs = []
    for li, ln in enumerate(lines):
        a, b = _segments(ln)
        for si in range(len(a)):
            segs.append((li, si, a[si], b[si]))
    if len(lines) > 1 and segs:
        A = np.array([s[2] for s in segs])
        B = np.array([s[3] for s in segs])
        inserts: Dict[int, list] = {}
        for i, j, t, u, x, y in segment_intersections(A, B, A, B):
            li, si = segs[i][0], segs[i][1]
            lj = segs[j][0]
            if li >= lj:
                continue  # each pair once; lines do not self-intersect here
            inserts.setdefault(li, []).append((si, t, x, y))
            inserts.setdefault(lj, []).append((segs[j][1], u, x, y))
        for li, ins in inserts.items():
            _insert_points(lines[li], ins)
        if inserts:
            log.append("Line intersections inserted: %d" % (sum(len(v) for v in inserts.values()) // 2))

    step(0.08, "Resampling lines")
    lines = [_resample_line(ln, size, inp.corner_angle, frozen) for ln in lines]

    # ---- constraint nodes and segments (with priorities) -------------------
    node_xy: List[np.ndarray] = []
    node_type: List[int] = []
    seg_pairs: List[Tuple[int, int]] = []
    index: Dict[Tuple[float, float], int] = {}

    def add_node(p, kind) -> int:
        key = (round(float(p[0]) / snap) * snap, round(float(p[1]) / snap) * snap)
        if key in index:
            j = index[key]
            if _PRIORITY[kind] < _PRIORITY[node_type[j]]:
                node_type[j] = kind
            return j
        index[key] = len(node_xy)
        node_xy.append(np.asarray(p, float))
        node_type.append(kind)
        return index[key]

    for p in sub_xy:
        add_node(p, SUBMESH)
    if inp.hard_points is not None:
        for p in np.asarray(inp.hard_points, float).reshape(-1, 2):
            add_node(p, HARDPOINT)
    for ln in sorted(lines, key=lambda l: _PRIORITY[l.kind]):
        ids = [add_node(p, ln.kind) for p in ln.coords]
        pairs = list(zip(ids[:-1], ids[1:]))
        if ln.closed:
            pairs.append((ids[-1], ids[0]))
        seg_pairs += [(a, b) for a, b in pairs if a != b]
    cons_xy = np.array(node_xy)
    cons_type = np.array(node_type, dtype=np.int64)

    # overlapping constraint nodes of different objects: drop the lower priority one
    cons_h = size(cons_xy[:, 0], cons_xy[:, 1])
    grid = PointGrid(cons_xy, float(np.max(cons_h)) * 0.3 + 1e-9)
    merge = np.arange(len(cons_xy))
    seg_set = set((min(a, b), max(a, b)) for a, b in seg_pairs)
    order = np.argsort([_PRIORITY[t] for t in cons_type], kind="stable")
    alive = np.ones(len(cons_xy), dtype=bool)
    for i in order:
        if not alive[i]:
            continue
        r = 0.25 * cons_h[i]
        ids = grid.query_rect(cons_xy[i, 0] - r, cons_xy[i, 1] - r, cons_xy[i, 0] + r, cons_xy[i, 1] + r)
        for j in ids:
            if j == i or not alive[j]:
                continue
            if _PRIORITY[cons_type[j]] <= _PRIORITY[cons_type[i]] and cons_type[j] in (SUBMESH, HARDPOINT, HARDLINE, HARDISLAND):
                continue  # never drop fixed geometry
            if (min(i, j), max(i, j)) in seg_set:
                continue  # neighbours on the same line
            if math.hypot(*(cons_xy[j] - cons_xy[i])) < r:
                alive[j] = False
                merge[j] = i
    if not alive.all():
        seg_pairs = [(merge[a], merge[b]) for a, b in seg_pairs]
        seg_pairs = [(a, b) for a, b in seg_pairs if a != b]
        log.append("Overlapping nodes merged: %d" % int((~alive).sum()))
    remap = -np.ones(len(cons_xy), dtype=np.int64)
    remap[alive] = np.arange(int(alive.sum()))
    cons_xy, cons_type, cons_h = cons_xy[alive], cons_type[alive], cons_h[alive]
    seg_pairs = sorted(set((int(min(remap[a], remap[b])), int(max(remap[a], remap[b]))) for a, b in seg_pairs))
    log.append("Constraint nodes: %d, segments: %d" % (len(cons_xy), len(seg_pairs)))

    # ---- interior nodes ------------------------------------------------------
    step(0.12, "Placing interior nodes")
    cand_levels, h_levels, levels = _interior_candidates(inp, rings, bbox, progress, log)
    if isinstance(cand_levels, np.ndarray):
        cand_levels, h_levels = [cand_levels], [h_levels]
    sp = np.array(seg_pairs, dtype=np.int64).reshape(-1, 2)
    seg_a, seg_b = cons_xy[sp[:, 0]], cons_xy[sp[:, 1]]
    step(0.42, "Spacing interior nodes")
    int_xy, int_h = _filter_interior(cand_levels, h_levels, cons_xy, cons_h, seg_a, seg_b, inp)
    log.append("Interior nodes: %d" % len(int_xy))

    # ---- triangulate, keep constraints, smooth --------------------------------------
    xy = np.vstack([cons_xy, int_xy]) if len(int_xy) else cons_xy.copy()
    ntype = np.concatenate([cons_type, np.zeros(len(int_xy), dtype=np.int64)])
    dens = np.concatenate([cons_h, int_h])
    # ghost nodes just outside the boundary keep the Delaunay step well conditioned
    # where the outline is concave or nearly straight over long distances
    ghost = _ghost_nodes(xy, sp, size, inp.clearance, rings)
    if len(ghost):
        xy = np.vstack([xy, ghost])
        ntype = np.concatenate([ntype, np.full(len(ghost), GHOST, dtype=np.int64)])
        dens = np.concatenate([dens, size(ghost[:, 0], ghost[:, 1])])
    seg_list = [tuple(s) for s in sp.tolist()]

    def conform(xy, ntype, dens, seg_list, label):
        """Delaunay + removal of outside triangles; split constraint segments until all are edges."""
        missing = []
        for rnd in range(8):
            ikle = np.asarray(triangulate(xy), dtype=np.int64)
            cx = xy[ikle, 0].mean(axis=1)
            cy = xy[ikle, 1].mean(axis=1)
            ikle = ikle[points_in_rings(cx, cy, rings)]
            edges = _edge_set(ikle)
            missing = [s for s in seg_list if (min(s), max(s)) not in edges]
            if not missing:
                return xy, ntype, dens, seg_list, ikle, 0
            if _DEBUG:
                for a, b in missing[:20]:
                    raw = np.asarray(triangulate(xy), dtype=np.int64)
                    adj = [t for t in raw.tolist() if a in t and b in t]
                    print("DEBUG missing", label, rnd, np.round(xy[a], 3).tolist(), np.round(xy[b], 3).tolist(), int(ntype[a]), int(ntype[b]),
                          "adjacent", [(np.round(xy[t].mean(axis=0), 2).tolist(), [int(ntype[v]) for v in t],
                                        bool(points_in_rings(xy[t, 0].mean(keepdims=True), xy[t, 1].mean(keepdims=True), rings)[0])) for t in adj])
            log.append("%s, pass %d: %d constraint segments missing, splitting them (first at %.6g, %.6g)" % (
                label, rnd + 1, len(missing), *(0.5 * (xy[missing[0][0]] + xy[missing[0][1]]))))
            new_pts, new_types = [], []
            drop = np.zeros(len(xy), dtype=bool)
            n0 = len(xy)
            miss = set(missing)
            new_segs = [s for s in seg_list if s not in miss]
            for a, b in missing:
                m = 0.5 * (xy[a] + xy[b])
                r = 0.5 * math.hypot(*(xy[b] - xy[a]))
                d = np.hypot(xy[:, 0] - m[0], xy[:, 1] - m[1])
                drop |= (d < r * 1.05) & ((ntype == INTERIOR) | (ntype == GHOST))
                j = n0 + len(new_pts)
                new_pts.append(m)
                # the new node belongs to the lower-priority line of the two endpoints
                new_types.append(int(ntype[a] if _PRIORITY[int(ntype[a])] >= _PRIORITY[int(ntype[b])] else ntype[b]))
                new_segs += [(a, j), (j, b)]
            new_pts = np.array(new_pts)
            xy = np.vstack([xy, new_pts])
            ntype = np.concatenate([ntype, np.array(new_types, dtype=np.int64)])
            dens = np.concatenate([dens, size(new_pts[:, 0], new_pts[:, 1])])
            keep = ~np.concatenate([drop, np.zeros(len(new_pts), dtype=bool)])
            remap = -np.ones(len(xy), dtype=np.int64)
            remap[keep] = np.arange(int(keep.sum()))
            xy, ntype, dens = xy[keep], ntype[keep], dens[keep]
            seg_list = [(int(remap[a]), int(remap[b])) for a, b in new_segs]
        return xy, ntype, dens, seg_list, ikle, len(missing)

    step(0.55, "Triangulating")
    xy, ntype, dens, seg_list, ikle, n_missing = conform(xy, ntype, dens, seg_list, "Triangulation")
    # gap filling: triangles much larger than the size field get their circumcentre
    # (or centroid) as a new interior node, kept clear of the constraint segments
    added_total = 0
    for rnd in range(4):
        new = _gap_nodes(xy, ikle, seg_list, size, inp.clearance, rings)
        if len(new) == 0:
            break
        xy = np.vstack([xy, new])
        ntype = np.concatenate([ntype, np.zeros(len(new), dtype=np.int64)])
        dens = np.concatenate([dens, size(new[:, 0], new[:, 1])])
        added_total += len(new)
        xy, ntype, dens, seg_list, ikle, n_missing = conform(xy, ntype, dens, seg_list, "Gap filling %d" % (rnd + 1))
    if added_total:
        log.append("Gap filling: %d nodes added in oversized triangles" % added_total)
    # smoothing rounds: move interior nodes (topology fixed, no inversion), then re-triangulate
    rounds = max(0, int(inp.smooth_iterations) // 3)
    for rnd in range(rounds):
        step(0.65 + 0.2 * (rnd + 1) / max(rounds, 1), "Smoothing (round %d of %d)" % (rnd + 1, rounds))
        used = np.zeros(len(xy), dtype=bool)
        used[ikle.ravel()] = True
        free = (ntype == INTERIOR) & used
        ik = _orient_ccw(xy[:, 0], xy[:, 1], ikle)
        sx, sy = _smooth(xy[:, 0], xy[:, 1], ik, free, 3)
        trial = np.stack([sx, sy], axis=1)
        t_xy, t_type, t_dens, t_segs, t_ikle, t_missing = conform(trial, ntype.copy(), dens.copy(), list(seg_list), "Smoothing round %d" % (rnd + 1))
        if t_missing:
            log.append("Smoothing round %d rejected (constraints could not be kept)" % (rnd + 1))
            break
        xy, ntype, dens, seg_list, ikle = t_xy, t_type, t_dens, t_segs, t_ikle
    if n_missing:
        log.append("Warning: %d constraint segments could not be recovered" % n_missing)

    # ---- merge submeshes ---------------------------------------------------------
    x, y = xy[:, 0].copy(), xy[:, 1].copy()
    quads_all: List[np.ndarray] = []
    if inp.submeshes:
        key = {(round(float(px) / snap), round(float(py) / snap)): i for i, (px, py) in enumerate(zip(x, y))}
        add_x, add_y, add_t, add_d, tris = [], [], [], [], []
        n = len(x)
        for sm in inp.submeshes:
            loc = np.empty(len(sm.x), dtype=np.int64)
            for k, (px, py) in enumerate(zip(sm.x, sm.y)):
                kk = (round(float(px) / snap), round(float(py) / snap))
                if kk in key:
                    loc[k] = key[kk]
                    ntype[key[kk]] = SUBMESH
                else:
                    loc[k] = n + len(add_x)
                    key[kk] = loc[k]
                    add_x.append(px)
                    add_y.append(py)
                    add_t.append(SUBMESH)
            tris.append(loc[np.asarray(sm.ikle, dtype=np.int64)])
            if sm.quads is not None and len(sm.quads):
                quads_all.append(loc[np.asarray(sm.quads, dtype=np.int64)])
        if add_x:
            x = np.concatenate([x, add_x])
            y = np.concatenate([y, add_y])
            ntype = np.concatenate([ntype, np.array(add_t, dtype=np.int64)])
            dens = np.concatenate([dens, size(np.array(add_x), np.array(add_y))])
        ikle = np.vstack([ikle] + tris)

    # ---- tidy: unused nodes, orientation, smoothing, over-constrained -----------
    used = np.zeros(len(x), dtype=bool)
    used[ikle.ravel()] = True
    remap = -np.ones(len(x), dtype=np.int64)
    remap[used] = np.arange(int(used.sum()))
    x, y, ntype, dens = x[used], y[used], ntype[used], dens[used]
    ikle = remap[ikle]
    struct_quads = remap[np.vstack(quads_all)] if quads_all else np.zeros((0, 4), dtype=np.int64)
    seg_edges = set()
    for a, b in seg_list:
        if 0 <= a < len(remap) and 0 <= b < len(remap) and remap[a] >= 0 and remap[b] >= 0:
            seg_edges.add((int(min(remap[a], remap[b])), int(max(remap[a], remap[b]))))
    ikle = _orient_ccw(x, y, ikle)
    step(0.92, "Final smoothing")
    x, y = _smooth(x, y, ikle, ntype == INTERIOR, 2)
    e = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    e = np.sort(e, axis=1)
    ue, cnt = np.unique(e, axis=0, return_counts=True)
    on_b = np.zeros(len(x), dtype=bool)
    on_b[ue[cnt == 1].ravel()] = True
    ikle, nfix = _fix_overconstrained(x, y, ikle, on_b, seg_edges)
    if nfix:
        log.append("Over-constrained triangles fixed by edge flips: %d" % nfix)
    stats = mesh_statistics(x, y, ikle)
    amin = min_angles(x, y, ikle)
    stats["share_min_angle_below_20"] = float(np.mean(amin < 20.0)) if amin.size else 0.0
    stats["share_min_angle_below_30"] = float(np.mean(amin < 30.0)) if amin.size else 0.0
    stats["min_angle_p05"] = float(np.percentile(amin, 5)) if amin.size else 0.0
    log.append("Mesh: %d nodes, %d triangles, smallest angle %.1f deg (5%% of triangles below %.1f deg, %.2f%% below 20 deg)" % (
        stats["nodes"], stats["elements"], stats["angle_min"], stats["min_angle_p05"], 100 * stats["share_min_angle_below_20"]))
    step(1.0, "Done")
    return MeshResult(x, y, ikle, ntype, dens, stats, log, struct_quads)

