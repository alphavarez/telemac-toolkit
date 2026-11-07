"""Unstructured triangle mesh generator with structured parts.

1. Domain rings and guide lines are resampled to the local element size:
   corners (turning more than ``corner_angle`` at the scale of an element) stay
   nodes, and so does any other vertex the new edges would otherwise miss by
   more than the outline tolerance (5 % of the size by default), unless it is
   closer than half an element to a node already kept; breaklines keep every
   vertex (long segments are subdivided); line intersections become nodes.
   From here on the resampled rings define the inside.
2. Fixed points, breakline nodes and the boundary nodes of structured parts
   are added. Where nodes of different objects overlap, the higher priority
   one is kept: structured part > fixed point > breakline > domain > guide line
   (two nodes of the same line are never merged).
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
    # outline/line vertices turning less are kept too where the resampled edges would stray from them by
    # more than this (metres); None = outline_fraction x the local size
    outline_tolerance: Optional[float] = None
    outline_fraction: float = 0.05
    clearance: float = 0.55  # interior nodes stay this fraction of the size away from constraints
    spacing: float = 0.80  # minimum distance between interior nodes, fraction of the size
    max_nodes: Optional[int] = 2_000_000  # refuse a mesh estimated above this (None = no limit)
    min_angle: Optional[float] = None  # quality-guaranteed method: refine until every angle reaches this (degrees)


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


def _resample_arcs(pts: np.ndarray, size: SizeField) -> np.ndarray:
    """Distances along one piece (between two fixed vertices) of its new interior nodes."""
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    total = float(seglen.sum())
    if total <= 0:
        return np.zeros(0)
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
        return np.zeros(0)
    targets = np.linspace(0, T[-1], nseg + 1)[1:-1]
    return np.interp(targets, T, s)


def _resample_piece(pts: np.ndarray, size: SizeField) -> np.ndarray:
    """Interior nodes for one piece between two fixed vertices (ends excluded)."""
    sp = _resample_arcs(pts, size)
    if not sp.size:
        return np.zeros((0, 2))
    seg = np.diff(pts, axis=0)
    cum = np.concatenate([[0], np.cumsum(np.hypot(seg[:, 0], seg[:, 1]))])
    return np.stack([np.interp(sp, cum, pts[:, 0]), np.interp(sp, cum, pts[:, 1])], axis=1)


def _chord_deviation(pts: np.ndarray, cum: np.ndarray, nodes: np.ndarray, ncum: np.ndarray) -> np.ndarray:
    """Distance of each interior vertex of ``pts`` (at distances ``cum`` along it) from the chord of
    ``nodes`` (at ``ncum`` along it, ends included) that spans it."""
    v = pts[1:-1]
    j = np.clip(np.searchsorted(ncum, cum[1:-1], side="right") - 1, 0, len(nodes) - 2)
    a = nodes[j]
    d = nodes[j + 1] - a
    L2 = np.einsum("ij,ij->i", d, d)
    t = np.clip(np.einsum("ij,ij->i", v - a, d) / np.maximum(L2, 1e-300), 0.0, 1.0)
    q = v - a - t[:, None] * d
    return np.hypot(q[:, 0], q[:, 1])


def _follow_piece(pts: np.ndarray, size: SizeField, tol, stats: dict) -> List[np.ndarray]:
    """Resample one piece; while a drawn vertex strays more than ``tol`` from the new chords, keep the
    worst one as a node and resample both halves (Douglas-Peucker on the resampled line).

    A vertex closer than half an element (along the line) to a node already kept is not kept itself:
    following a detail smaller than that would need tiny elements (digitizing zigzags, noise); those
    stay off by more than ``tol`` and are counted in ``stats``. Returns the sub-pieces' node lists in
    order, each without its last node."""
    out = []
    stack = [pts]
    while stack:
        p = stack.pop()
        seg = np.diff(p, axis=0)
        cum = np.concatenate([[0], np.cumsum(np.hypot(seg[:, 0], seg[:, 1]))])
        sp = _resample_arcs(p, size)
        ncum = np.concatenate([[0.0], sp, [cum[-1]]])
        nodes = np.stack([np.interp(ncum, cum, p[:, 0]), np.interp(ncum, cum, p[:, 1])], axis=1)
        nodes[0], nodes[-1] = p[0], p[-1]
        if len(p) > 2:
            dev = _chord_deviation(p, cum, nodes, ncum)
            v = p[1:-1]
            lim = tol(v) if callable(tol) else float(tol)
            gap = 0.5 * size(v[:, 0], v[:, 1])
            room = (cum[1:-1] >= gap) & (cum[-1] - cum[1:-1] >= gap)
            over = dev - lim
            cand = np.where(room, over, -np.inf)
            k = int(np.argmax(cand))
            if cand[k] > 0:
                stack.append(p[k + 1:])  # popped second: comes after
                stack.append(p[:k + 2])
                continue
            j = int(np.argmax(dev))
            if dev[j] > stats.get("dev", 0.0):
                stats["dev"], stats["at"] = float(dev[j]), (float(v[j, 0]), float(v[j, 1]))
            stats["over"] = stats.get("over", 0) + int(np.sum(over > 0))
        out.append(nodes[:-1])
    return out


def _turning_at_scale(c: np.ndarray, closed: bool, idx: np.ndarray, half: np.ndarray) -> np.ndarray:
    """Change of direction (degrees) at vertices ``idx`` between the chords to the points ``half`` before
    and after them along the line."""
    pts = np.vstack([c, c[:1]]) if closed else c
    seg = np.diff(pts, axis=0)
    cum = np.concatenate([[0.0], np.cumsum(np.hypot(seg[:, 0], seg[:, 1]))])
    total = cum[-1]
    s0 = cum[idx]
    sa, sb = s0 - half, s0 + half
    if closed and total > 0:
        sa, sb = np.mod(sa, total), np.mod(sb, total)
    else:
        sa, sb = np.clip(sa, 0.0, total), np.clip(sb, 0.0, total)
    pa = np.stack([np.interp(sa, cum, pts[:, 0]), np.interp(sa, cum, pts[:, 1])], axis=1)
    pb = np.stack([np.interp(sb, cum, pts[:, 0]), np.interp(sb, cum, pts[:, 1])], axis=1)
    v1 = c[idx] - pa
    v2 = pb - c[idx]
    h1 = np.arctan2(v1[:, 1], v1[:, 0])
    h2 = np.arctan2(v2[:, 1], v2[:, 0])
    out = np.degrees(np.abs((h2 - h1 + np.pi) % (2 * np.pi) - np.pi))
    # an end of an open line, or a degenerate window: keep the vertex's own turning
    bad = (np.hypot(v1[:, 0], v1[:, 1]) <= 0) | (np.hypot(v2[:, 0], v2[:, 1]) <= 0)
    out[bad] = 180.0
    return out


def _resample_line(line: _Line, size: SizeField, corner_angle: float, frozen: set, tol=None,
                   stats: Optional[dict] = None) -> _Line:
    """Nodes along a line at the local size. Vertices turning more than ``corner_angle`` stay; others go,
    unless the new edges would stray more than ``tol`` (metres, or a function of (n, 2) points) from them."""
    c = line.coords
    if not line.resample or len(c) < 2:
        return line
    fixed = line.fixed.copy()
    ang = turning_angles(c, line.closed)
    corner = np.flatnonzero(ang > corner_angle)
    if corner.size:
        # a corner at the scale of the elements, not a zigzag of the digitizing: the line still turns
        # that much between half an element before the vertex and half an element after it
        h = size(c[corner, 0], c[corner, 1])
        corner = corner[_turning_at_scale(c, line.closed, corner, 0.5 * h) > corner_angle]
        fixed[corner] = True
    if line.closed and not fixed.any():
        fixed[0] = True
    n = len(c)
    idx_fixed = np.flatnonzero(fixed)
    out_c, out_f = [], []
    if line.closed:
        order = list(idx_fixed) + [idx_fixed[0] + n]
    else:
        order = list(idx_fixed)
    kept = 0
    for a, b in zip(order[:-1], order[1:]):
        piece_idx = [k % n for k in range(a, b + 1)]
        pts = c[piece_idx]
        key = (tuple(np.round(pts[0], 9)), tuple(np.round(pts[-1], 9)))
        if len(pts) == 2 and (key in frozen or key[::-1] in frozen):
            out_c.append(pts[0])
            out_f.append(True)
            continue
        parts = (_follow_piece(pts, size, tol, stats if stats is not None else {}) if tol is not None
                 else [np.vstack([pts[:1], _resample_piece(pts, size)])])
        kept += len(parts) - 1
        for part in parts:
            out_c.append(part[0])
            out_f.append(True)
            for p in part[1:]:
                out_c.append(p)
                out_f.append(False)
    if not line.closed:
        out_c.append(c[order[-1]])
        out_f.append(True)
    res = _Line(np.array(out_c), line.closed, line.kind, False, np.array(out_f, dtype=bool))
    res.kept = kept
    return res


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


def _fix_overconstrained(x, y, ikle, on_boundary: np.ndarray, constrained_edges: set, angle_floor: float = 20.0):
    """Flip the interior edge of triangles whose three nodes are on the boundary.

    Neighbours are found in a sorted table of all edges built once with NumPy; the few triangles a flip
    changes are kept aside and looked up directly, so the cost does not grow with the number of flips.
    """
    ikle = ikle.copy()
    bad = np.flatnonzero(on_boundary[ikle].all(axis=1))
    if bad.size == 0:
        return ikle, 0
    n = int(max(x.size, int(ikle.max()) + 1))
    a0 = np.concatenate([ikle[:, 0], ikle[:, 1], ikle[:, 2]]).astype(np.int64)
    b0 = np.concatenate([ikle[:, 1], ikle[:, 2], ikle[:, 0]]).astype(np.int64)
    keys = np.minimum(a0, b0) * n + np.maximum(a0, b0)
    owner = np.tile(np.arange(len(ikle), dtype=np.int64), 3)
    order = np.argsort(keys, kind="stable")
    keys, owner = keys[order], owner[order]
    changed: set = set()

    def tri_keys(t):
        tri = ikle[t]
        return {min(int(tri[k]), int(tri[(k + 1) % 3])) * n + max(int(tri[k]), int(tri[(k + 1) % 3])) for k in range(3)}

    def owners(key):
        lo = int(np.searchsorted(keys, key, side="left"))
        hi = int(np.searchsorted(keys, key, side="right"))
        out = [int(t) for t in owner[lo:hi] if int(t) not in changed]
        out += [t for t in changed if key in tri_keys(t)]
        return out

    fixed = 0
    for t in bad:
        t = int(t)
        tri = ikle[t]
        if not on_boundary[tri].all():
            continue
        for k in range(3):
            a, b, c = int(tri[k]), int(tri[(k + 1) % 3]), int(tri[(k + 2) % 3])
            key = min(a, b) * n + max(a, b)
            if (min(a, b), max(a, b)) in constrained_edges:
                continue
            own = owners(key)
            if len(own) != 2:
                continue
            u = own[0] if own[1] == t else own[1]
            d = next(int(v) for v in ikle[u] if v not in (a, b))
            # quad a-d-b-c must be convex: new triangles (c, a, d) and (c, d, b) both positive
            t1 = np.array([[c, a, d]])
            t2 = np.array([[c, d, b]])
            if signed_areas(x, y, t1)[0] > 0 and signed_areas(x, y, t2)[0] > 0:
                # and the flip must not make the worst angle of the pair worse (at sharp corners it would)
                old_min = float(min_angles(x, y, ikle[[t, u]]).min())
                new_min = float(min_angles(x, y, np.vstack([t1, t2])).min())
                if new_min < min(old_min, angle_floor):
                    continue
                ikle[t] = t1[0]
                ikle[u] = t2[0]
                changed.update((t, u))
                fixed += 1
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


def _circumcentres(xy: np.ndarray, ikle: np.ndarray):
    """Circumcentres and circumradii of triangles."""
    A, B, C = xy[ikle[:, 0]], xy[ikle[:, 1]], xy[ikle[:, 2]]
    b = B - A
    c = C - A
    d = 2.0 * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
    d = np.where(np.abs(d) < 1e-300, 1e-300, d)
    b2 = (b ** 2).sum(axis=1)
    c2 = (c ** 2).sum(axis=1)
    ux = (c[:, 1] * b2 - b[:, 1] * c2) / d
    uy = (b[:, 0] * c2 - c[:, 0] * b2) / d
    return A + np.stack([ux, uy], axis=1), np.hypot(ux, uy)


def _quality_points(xy, ntype, ikle, seg_list, size: SizeField, rings, target: float, floor: float = 0.3):
    """One round of Delaunay refinement (Ruppert): new nodes for triangles with an angle below ``target``.

    Each bad triangle, worst first, gets its circumcentre, unless that point is outside the domain,
    or within the diametral circle of a constraint segment: then the segment is split at its middle
    instead (segments of structured parts are never split; their triangle is left as it is). A
    triangle whose circumradius is below ``floor`` times the local size is not refined, so refinement
    stops near small angles of the input (sharp outline corners) instead of running forever.
    Returns (new points, their types, segments to split as {segment index: True}, number of bad triangles).
    """
    amin = min_angles(xy[:, 0], xy[:, 1], ikle)
    bad = np.flatnonzero(amin < target)
    if bad.size == 0:
        return np.zeros((0, 2)), [], set(), 0
    cc, R = _circumcentres(xy, ikle[bad])
    h = size(cc[:, 0], cc[:, 1])
    keep = R >= floor * h
    order = np.argsort(amin[bad][keep])
    cc, R, h = cc[keep][order], R[keep][order], h[keep][order]
    n_bad = int(bad.size)
    if len(cc) == 0:
        return np.zeros((0, 2)), [], set(), n_bad
    inside = points_in_rings(cc[:, 0], cc[:, 1], rings)
    sp = np.array(seg_list, dtype=np.int64).reshape(-1, 2)
    sa, sb = xy[sp[:, 0]], xy[sp[:, 1]]
    mid = 0.5 * (sa + sb)
    rad = 0.5 * np.hypot(*(sb - sa).T)
    frozen = (ntype[sp[:, 0]] == SUBMESH) & (ntype[sp[:, 1]] == SUBMESH)
    grid = PointGrid(mid, float(np.max(rad)) + 1e-9 if len(rad) else 1.0)
    split: set = set()
    accepted: List[np.ndarray] = []
    rmax = float(np.max(rad)) if len(rad) else 0.0

    def encroached(p):
        ids = grid.query_rect(p[0] - rmax, p[1] - rmax, p[0] + rmax, p[1] + rmax) if len(rad) else []
        return [j for j in ids if math.hypot(*(p - mid[j])) < rad[j] * (1 - 1e-9)]

    for i in range(len(cc)):
        p = cc[i]
        enc = encroached(p)
        if enc:
            # split the encroached segments; edges of structured parts are never split (their triangle stays)
            for j in enc:
                if not frozen[j] and rad[j] >= 0.25 * floor * h[i]:
                    split.add(int(j))
            continue
        if not inside[i]:
            continue
        # keep new nodes apart from each other (several bad triangles can share a circumcentre region)
        if accepted:
            A = np.array(accepted)
            if np.min(np.hypot(A[:, 0] - p[0], A[:, 1] - p[1])) < 0.5 * min(h[i], 2 * R[i]):
                continue
        accepted.append(p)
    return (np.array(accepted) if accepted else np.zeros((0, 2))), [INTERIOR] * len(accepted), split, n_bad


def _ring_area(r) -> float:
    r = as_ring(r)
    x, y = r[:, 0], r[:, 1]
    return 0.5 * abs(float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def _ring_length(r) -> float:
    r = as_ring(r)
    d = np.roll(r, -1, axis=0) - r
    return float(np.hypot(d[:, 0], d[:, 1]).sum())


def estimate_node_count(inp: MeshInput, rings, bbox, grid: int = 300) -> int:
    """Rough number of nodes the size field asks for (order of magnitude, before any meshing).

    Two estimates, the larger one is kept: the size field integrated on a grid over the domain
    (an equilateral triangle of edge h takes 0.866 h² per node), and the sum of what each
    refinement source asks for on its own (its area, its growth zone), which also sees sources
    smaller than a grid cell.
    """
    size = inp.size
    k = 0.8660254037844386
    x0, y0, x1, y1 = bbox
    gx = (np.arange(grid) + 0.5) / grid * (x1 - x0) + x0
    gy = (np.arange(grid) + 0.5) / grid * (y1 - y0) + y0
    X, Y = np.meshgrid(gx, gy)
    X, Y = X.ravel(), Y.ravel()
    inside = points_in_rings(X, Y, rings)
    cell = (x1 - x0) * (y1 - y0) / (grid * grid)
    est_grid = 0.0
    if inside.any():
        h = size(X[inside], Y[inside])
        est_grid = float(np.sum(cell / (k * np.maximum(h, 1e-300) ** 2)))

    area = sum(_ring_area(poly[0]) - sum(_ring_area(r) for r in poly[1:]) for poly in inp.domain if len(poly))
    d = size.default
    g = size.growth - 1.0
    est_src = area / (k * d * d)
    for pd in size.polygons:
        if pd.size <= 0:
            continue
        a = sum(_ring_area(r) for r in pd.rings[:1]) - sum(_ring_area(r) for r in pd.rings[1:])
        est_src += max(a, 0.0) / (k * pd.size ** 2)
        if g > 0 and d > pd.size:
            est_src += sum(_ring_length(r) for r in pd.rings) / (k * g) * (1.0 / pd.size - 1.0 / d)
    for pt in size.points:
        for sz in np.asarray(pt.size, dtype=float).ravel():
            if sz > 0 and g > 0 and d > sz:
                est_src += 2.0 * math.pi / (k * g * g) * (math.log(d / sz) + sz / d - 1.0)
    for ln in size.lines:
        if ln.size > 0 and g > 0 and d > ln.size:
            c = np.asarray(ln.coords, dtype=float)
            length = _ring_length(c) if ln.closed else float(np.hypot(*np.diff(c, axis=0).T).sum())
            est_src += 2.0 * length / (k * g) * (1.0 / ln.size - 1.0 / d)
    return int(min(max(est_grid, est_src), 1e15))


def generate(inp: MeshInput, triangulate: Triangulator, progress: Progress = None, interior: Optional[Callable] = None) -> MeshResult:
    """Mesh ``inp``. ``interior``: another mesher for the inside (see gmsh_engine.gmsh_interior); the domain rings,
    constraint lines, structured parts and the tidying are this module's in every case."""
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
        inp = dc_replace(inp, size=SizeField(s0.default, s0.growth, list(s0.polygons), list(s0.points), list(s0.lines) + extra,
                                             list(s0.grids)))
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

    # ---- is the mesh feasible? Refuse before spending hours on it ------------
    est = estimate_node_count(inp, rings, bbox)
    log.append("Estimated nodes: about %s" % format(int(est), ","))
    if inp.max_nodes and est > inp.max_nodes:
        raise ValueError(
            "The mesh would have about %s nodes, more than the limit of %s. The domain is %.6g x %.6g (map units), "
            "the default element size is %.6g and the smallest size asked for is %.6g. Increase the default size or the "
            "sizes of the refinement areas, points and lines (or make the refined areas smaller), then generate again."
            % (format(int(est), ","), format(int(inp.max_nodes), ","), bbox[2] - bbox[0], bbox[3] - bbox[1],
               size.default, size.minimum)
        )

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
    n_ring_lines = len(lines)  # the domain's own rings come first; resampling keeps the order
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
    if inp.outline_tolerance is not None and inp.outline_tolerance > 0:
        tol = float(inp.outline_tolerance)
    else:
        frac = float(inp.outline_fraction)
        tol = (lambda q: frac * size(q[:, 0], q[:, 1])) if frac > 0 else None
    follow: Dict[str, object] = {}
    lines = [_resample_line(ln, size, inp.corner_angle, frozen, tol, follow) for ln in lines]
    # from here on "inside" means inside the outline as resampled: its edges are chords of the drawn one,
    # and triangles between the two must neither be kept nor thrown away by the drawn rings
    rings = [as_ring(ln.coords) for ln in lines[:n_ring_lines]]
    kept = sum(getattr(ln, "kept", 0) for ln in lines)
    if tol is not None:
        within = "%g m" % tol if not callable(tol) else "%g%% of the element size" % (100 * frac)
        if kept:
            log.append("Drawn vertices kept to follow the outline and lines within %s: %d" % (within, kept))
        if follow.get("dev"):
            log.append("Largest distance between the drawn lines and the mesh edges: %.3g m (at %.6g, %.6g)"
                       % (follow["dev"], *follow["at"]))
        if follow.get("over"):
            log.append("%d drawn vertices are more than %s off the mesh edges: they are details smaller than half "
                       "an element (zigzags of the digitizing?). Smaller elements there, or a smoother line, "
                       "follow them closer." % (follow["over"], within))

    # ---- constraint nodes and segments (with priorities) -------------------
    node_xy: List[np.ndarray] = []
    node_type: List[int] = []
    seg_pairs: List[Tuple[int, int]] = []
    index: Dict[Tuple[float, float], int] = {}

    node_line: List[int] = []  # the line a node was first made for (-1: a submesh or fixed point)

    def add_node(p, kind, line=-1) -> int:
        key = (round(float(p[0]) / snap) * snap, round(float(p[1]) / snap) * snap)
        if key in index:
            j = index[key]
            if _PRIORITY[kind] < _PRIORITY[node_type[j]]:
                node_type[j] = kind
            return j
        index[key] = len(node_xy)
        node_xy.append(np.asarray(p, float))
        node_type.append(kind)
        node_line.append(line)
        return index[key]

    for p in sub_xy:
        add_node(p, SUBMESH)
    if inp.hard_points is not None:
        for p in np.asarray(inp.hard_points, float).reshape(-1, 2):
            add_node(p, HARDPOINT)
    ring_raw: List[Tuple[int, int]] = []
    for li, ln in sorted(enumerate(lines), key=lambda t: _PRIORITY[t[1].kind]):
        ids = [add_node(p, ln.kind, li) for p in ln.coords]
        pairs = list(zip(ids[:-1], ids[1:]))
        if ln.closed:
            pairs.append((ids[-1], ids[0]))
        seg_pairs += [(a, b) for a, b in pairs if a != b]
        if li < n_ring_lines:
            ring_raw += [(a, b) for a, b in pairs if a != b]
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
            if node_line[i] >= 0 and node_line[i] == node_line[j]:
                continue  # two nodes of one line: merging would cut out the ones between them
            if math.hypot(*(cons_xy[j] - cons_xy[i])) < r:
                alive[j] = False
                merge[j] = i
    if not alive.all():
        seg_pairs = [(merge[a], merge[b]) for a, b in seg_pairs]
        seg_pairs = [(a, b) for a, b in seg_pairs if a != b]
        ring_raw = [(merge[a], merge[b]) for a, b in ring_raw if merge[a] != merge[b]]
        log.append("Overlapping nodes merged: %d" % int((~alive).sum()))
    remap = -np.ones(len(cons_xy), dtype=np.int64)
    remap[alive] = np.arange(int(alive.sum()))
    cons_xy, cons_type, cons_h = cons_xy[alive], cons_type[alive], cons_h[alive]
    seg_pairs = sorted(set((int(min(remap[a], remap[b])), int(max(remap[a], remap[b]))) for a, b in seg_pairs))
    ring_set = set((int(min(remap[a], remap[b])), int(max(remap[a], remap[b]))) for a, b in ring_raw)
    ring_mask = np.array([sg in ring_set for sg in seg_pairs], dtype=bool)
    log.append("Constraint nodes: %d, segments: %d" % (len(cons_xy), len(seg_pairs)))

    ext_quads = np.zeros((0, 4), dtype=np.int64)
    if interior is not None:
        step(0.15, "Meshing the inside")

        def sub_step(frac, msg):
            step(0.15 + 0.7 * frac, msg)

        xy, ntype, dens, ikle, ext_quads, ilog = interior(cons_xy, cons_type, cons_h, seg_pairs, ring_mask, size, sub_step)
        log += list(ilog)
        seg_list = [tuple(sg) for sg in seg_pairs]
        edges = _edge_set(ikle)
        n_missing = sum(1 for a, b in seg_list if (min(a, b), max(a, b)) not in edges)
    else:
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
                    if ntype[a] == SUBMESH and ntype[b] == SUBMESH:
                        # an edge of a fixed part (structured strip, existing mesh) is never split: a node in
                        # the middle would hang on the other side; clearing its circle is enough
                        new_segs.append((a, b))
                        continue
                    j = n0 + len(new_pts)
                    new_pts.append(m)
                    # the new node belongs to the lower-priority line of the two endpoints
                    new_types.append(int(ntype[a] if _PRIORITY[int(ntype[a])] >= _PRIORITY[int(ntype[b])] else ntype[b]))
                    new_segs += [(a, j), (j, b)]
                new_pts = np.array(new_pts).reshape(-1, 2)
                xy = np.vstack([xy, new_pts])
                ntype = np.concatenate([ntype, np.array(new_types, dtype=np.int64)])
                dens = np.concatenate([dens, size(new_pts[:, 0], new_pts[:, 1])]) if len(new_pts) else dens
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
        if inp.min_angle:
            target = float(inp.min_angle)
            added = 0
            for rnd in range(15):
                step(0.85 + 0.06 * (rnd + 1) / 15, "Quality refinement (round %d)" % (rnd + 1))
                new, new_t, split, n_bad = _quality_points(xy, ntype, ikle, seg_list, size, rings, target)
                if len(new) == 0 and not split:
                    break
                n0 = len(xy)
                mids, mid_t, new_segs = [], [], []
                sp = list(seg_list)
                for j, (a, b) in enumerate(sp):
                    if j in split:
                        k = n0 + len(new) + len(mids)
                        mids.append(0.5 * (xy[a] + xy[b]))
                        mid_t.append(int(ntype[a] if _PRIORITY[int(ntype[a])] >= _PRIORITY[int(ntype[b])] else ntype[b]))
                        new_segs += [(a, k), (k, b)]
                    else:
                        new_segs.append((a, b))
                pts = [p for p in (new, np.array(mids).reshape(-1, 2)) if len(p)]
                allp = np.vstack(pts)
                xy = np.vstack([xy, allp])
                ntype = np.concatenate([ntype, np.array(list(new_t) + mid_t, dtype=np.int64)])
                dens = np.concatenate([dens, size(allp[:, 0], allp[:, 1])])
                seg_list = new_segs
                added += len(allp)
                xy, ntype, dens, seg_list, ikle, n_missing = conform(xy, ntype, dens, seg_list, "Quality round %d" % (rnd + 1))
            amin = min_angles(xy[:, 0], xy[:, 1], ikle)
            below = amin < target
            by_struct = below & ((ntype[ikle] == SUBMESH).sum(axis=1) >= 2)
            log.append("Quality refinement (target %.0f deg): %d nodes added; smallest angle of the triangles %.1f deg, "
                       "%d below the target" % (target, added, float(amin.min()) if amin.size else 0.0, int(below.sum()))
                       + (" (%d of them on an edge of a structured part, which stays fixed)" % int(by_struct.sum())
                          if by_struct.any() else ""))
    if n_missing:
        log.append("Warning: %d constraint segments could not be recovered" % n_missing)

    # helper nodes outside the domain never belong to the mesh; a triangle with one is outside it
    if len(ikle):
        ghost_tri = (ntype[ikle] == GHOST).any(axis=1)
        if ghost_tri.any():
            log.append("Warning: %d triangles outside the outline removed (at %.6g, %.6g)" % (
                int(ghost_tri.sum()), *xy[ikle[ghost_tri][0]].mean(axis=0)))
            ikle = ikle[~ghost_tri]

    # ---- merge submeshes ---------------------------------------------------------
    x, y = xy[:, 0].copy(), xy[:, 1].copy()
    quads_all: List[np.ndarray] = [ext_quads] if len(ext_quads) else []
    if inp.submeshes:
        key = {(round(float(px) / snap), round(float(py) / snap)): i for i, (px, py) in enumerate(zip(x, y))}
        add_x, add_y, add_t, tris = [], [], [], []
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
    n_free = len(ikle) - sum(len(t) for t in tris) if inp.submeshes else len(ikle)

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
    if not inp.min_angle and interior is None:  # smoothing would undo the angle guarantee (Gmsh optimises its own)
        x, y = _smooth(x, y, ikle, ntype == INTERIOR, 2)
    e = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    e = np.sort(e, axis=1)
    ue, cnt = np.unique(e, axis=0, return_counts=True)
    on_b = np.zeros(len(x), dtype=bool)
    on_b[ue[cnt == 1].ravel()] = True
    # the elements of fixed parts (structured strips, grids, the old elements around a re-meshed area) are
    # never flipped: every edge of theirs counts as constrained
    if n_free < len(ikle):
        fe = np.sort(np.concatenate([ikle[n_free:, [0, 1]], ikle[n_free:, [1, 2]], ikle[n_free:, [2, 0]]]), axis=1)
        seg_edges |= set(map(tuple, fe.tolist()))
    ikle, nfix = _fix_overconstrained(x, y, ikle, on_b, seg_edges, max(20.0, float(inp.min_angle or 0.0)))
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
