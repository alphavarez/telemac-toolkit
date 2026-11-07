"""Checks of a drawing before it is meshed, measured against the element size wanted there.

A mesher follows the drawing: two vertices 2 cm apart on a breakline give elements 2 cm wide, a domain
outline that comes back within 1 m of itself gives slivers where the elements are 5 m, a line that ends
3 cm short of the outline leaves a thin strip to fill. None of this is wrong geometry for a GIS, so QGIS
does not flag it. This module finds it, with NumPy only:

``vertex_problems``  near-duplicate vertices (closer than 1 % of the element size), spikes (the line turns
                     back on itself), corners closer together than a quarter of the element size on an
                     outline, and breakline vertices closer than a tenth of it (every breakline vertex is a
                     node);
``clean_vertices``   the same line without its near-duplicates and its small spikes (the safe fixes);
``narrow_gaps``      places where an outline comes within half an element of itself (a neck) or of another
                     ring (an island almost touching the outline or another island);
``near_ends``        line ends and vertices lying close to an outline but not on it, and the point on the
                     outline they should be snapped to.

Each problem is a ``Problem``: what, where, how far, the element size there, and whether a safe fix exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import turning_angles

__all__ = ["Problem", "vertex_problems", "clean_vertices", "narrow_gaps", "near_ends", "DUPLICATE_FRACTION",
           "GAP_FRACTION", "SNAP_FRACTION"]

#: vertices closer than this fraction of the element size are near-duplicates
DUPLICATE_FRACTION = 0.01
#: corners of an outline closer than this fraction of the element size force elements of that size
CORNER_FRACTION = 0.25
#: breakline vertices closer than this fraction of the element size force elements of that size
BREAKLINE_FRACTION = 0.1
#: parts of outlines closer than this fraction of the element size make slivers
GAP_FRACTION = 0.5
#: line ends closer than this fraction of the element size to an outline (but not on it) should be snapped
SNAP_FRACTION = 0.25
#: a vertex where the line turns by more than this (degrees) is a spike: it goes back on itself
SPIKE_DEG = 170.0

SizeFn = Callable[[np.ndarray, np.ndarray], np.ndarray]

MESSAGES = {
    "duplicate": "two vertices %.3g m apart (elements here %.3g m): a near-duplicate",
    "spike": "the line turns back on itself (%.0f deg) at a vertex %.3g m from the next (elements here %.3g m)",
    "sharp": "a very sharp corner (%.0f deg turn): the triangles there will be thin",
    "corners": "two corners %.3g m apart (elements here %.3g m): the elements shrink to that size there",
    "breakline": "breakline vertices %.3g m apart (elements here %.3g m): every breakline vertex is a node, so the "
                 "elements shrink to that size; remove vertices or make it a guide line",
    "neck": "the outline comes within %.3g m of itself (elements here %.3g m): slivers; widen the gap, close it, or "
            "set a smaller size there",
    "gap": "%s within %.3g m of %s (elements here %.3g m): slivers; widen the gap, join them, or set a smaller size "
           "there",
    "near": "%s %.3g m from the outline (elements here %.3g m): it should touch it; snap it",
}


@dataclass
class Problem:
    code: str                 # duplicate, spike, sharp, corners, breakline, neck, gap, near (QGIS adds invalid, ...)
    x: float
    y: float
    message: str
    distance: float = math.nan  # the length at stake (m)
    size: float = math.nan      # the element size there (m)
    fixable: bool = False
    index: int = -1             # vertex index in the line or ring
    target: Optional[Tuple[float, float]] = None  # where to snap (near)
    other: str = ""             # the other object (gap)

    @property
    def severity(self) -> float:
        """Smaller is worse: the length at stake over the element size."""
        if not (self.distance >= 0 and self.size > 0):
            return 0.0
        return float(self.distance / self.size)


def _seg_lengths(c: np.ndarray, closed: bool) -> np.ndarray:
    nxt = np.roll(c, -1, axis=0) if closed else c[1:]
    d = nxt - (c if closed else c[:-1])
    return np.hypot(d[:, 0], d[:, 1])


def vertex_problems(coords, closed: bool, size: SizeFn, kind: str = "outline") -> List[Problem]:
    """Problems at the vertices of one line or ring (closed: without its closing vertex).

    ``kind``: "outline" (domain, island, grid zone: corners matter), "breakline" (every vertex is a node) or
    "line" (guide lines, channel banks: only duplicates and spikes).
    """
    c = np.asarray(coords, dtype=np.float64)[:, :2]
    n = len(c)
    out: List[Problem] = []
    if n < 2:
        return out
    L = _seg_lengths(c, closed)
    a = c if closed else c[:-1]
    b = np.roll(c, -1, axis=0) if closed else c[1:]
    mid = 0.5 * (a + b)
    h = np.asarray(size(mid[:, 0], mid[:, 1]), dtype=np.float64)
    dup = L < DUPLICATE_FRACTION * h
    for k in np.flatnonzero(dup):
        out.append(Problem("duplicate", float(mid[k, 0]), float(mid[k, 1]), MESSAGES["duplicate"] % (L[k], h[k]),
                           float(L[k]), float(h[k]), True, int(k)))
    if n >= 3:
        ang = turning_angles(c, closed)
        hv = np.asarray(size(c[:, 0], c[:, 1]), dtype=np.float64)
        # edge before and after each vertex
        before = np.roll(L, 1) if closed else np.concatenate([[np.inf], L])
        after = L if closed else np.concatenate([L, [np.inf]])
        short_side = np.minimum(before, after)
        for k in np.flatnonzero(ang > SPIKE_DEG):
            if dup[k % len(L)] or dup[(k - 1) % len(L)]:
                continue  # the duplicate is the problem; removing it settles the turn
            small = short_side[k] < hv[k]
            if small:
                out.append(Problem("spike", float(c[k, 0]), float(c[k, 1]), MESSAGES["spike"] % (ang[k], short_side[k], hv[k]),
                                   float(short_side[k]), float(hv[k]), True, int(k)))
            else:
                out.append(Problem("sharp", float(c[k, 0]), float(c[k, 1]), MESSAGES["sharp"] % ang[k],
                                   float(short_side[k]), float(hv[k]), False, int(k)))
        if kind == "outline":
            # corners as the mesher sees them: turning at the vertex and at the scale of an element (the zigzags
            # of fast digitizing are not corners, the resampled outline smooths them away)
            corner = ang > 30.0
            idx = np.flatnonzero(corner)
            if idx.size:
                from .engine import _turning_at_scale

                corner[idx] = _turning_at_scale(c, closed, idx, 0.5 * hv[idx]) > 30.0
            spike = ang > SPIKE_DEG
            ia = np.arange(len(L))
            ib = (ia + 1) % n
            both = corner[ia] & corner[ib] & ~spike[ia] & ~spike[ib] & ~dup & (L < CORNER_FRACTION * h)
            for k in np.flatnonzero(both):
                out.append(Problem("corners", float(mid[k, 0]), float(mid[k, 1]), MESSAGES["corners"] % (L[k], h[k]),
                                   float(L[k]), float(h[k]), False, int(k)))
    if kind == "breakline":
        for k in np.flatnonzero(~dup & (L < BREAKLINE_FRACTION * h)):
            out.append(Problem("breakline", float(mid[k, 0]), float(mid[k, 1]), MESSAGES["breakline"] % (L[k], h[k]),
                               float(L[k]), float(h[k]), False, int(k)))
    return out


def clean_vertices(coords, closed: bool, size: SizeFn, rounds: int = 10) -> np.ndarray:
    """The line without its near-duplicate vertices and its small spikes (the end vertices of an open line stay)."""
    c = np.asarray(coords, dtype=np.float64)[:, :2].copy()
    for _ in range(rounds):
        n = len(c)
        if n < (3 if closed else 2):
            return c
        keep = np.ones(n, dtype=bool)
        changed = False
        for p in vertex_problems(c, closed, size, kind="line"):
            if not p.fixable:
                continue
            if p.code == "duplicate":
                j = (p.index + 1) % n  # drop the second vertex of the pair, unless it ends an open line
                if not closed and j == n - 1:
                    j = p.index
                if not closed and j == 0:
                    continue
            else:  # spike
                j = p.index
                if not closed and j in (0, n - 1):
                    continue
            if keep[j]:
                keep[j] = False
                changed = True
        if not changed:
            break
        c = c[keep]
        if closed and len(c) < 3:
            break
    return c


def _sample(coords, closed: bool, size: SizeFn, frac: float):
    """Points along a line at most ``frac`` x the local size apart, with their distance along it."""
    c = np.asarray(coords, dtype=np.float64)[:, :2]
    pts = np.vstack([c, c[:1]]) if closed else c
    seg = np.diff(pts, axis=0)
    L = np.hypot(seg[:, 0], seg[:, 1])
    mid = 0.5 * (pts[1:] + pts[:-1])
    h = np.asarray(size(mid[:, 0], mid[:, 1]), dtype=np.float64)
    nsub = np.maximum(1, np.ceil(L / np.maximum(frac * h, 1e-12)).astype(np.int64))
    nsub = np.minimum(nsub, 10000)
    cum = np.concatenate([[0.0], np.cumsum(L)])
    seg_id = np.repeat(np.arange(len(L)), nsub)
    t = (np.arange(int(nsub.sum())) - np.repeat(np.cumsum(nsub) - nsub, nsub)) / np.repeat(nsub, nsub)
    xy = pts[seg_id] + t[:, None] * seg[seg_id]
    s = cum[seg_id] + t * L[seg_id]
    if not closed:
        xy = np.vstack([xy, pts[-1:]])
        s = np.concatenate([s, [cum[-1]]])
    return xy, s, float(cum[-1])


def _pairs_within(xy: np.ndarray, cell: float):
    """All pairs of points (i < j) in the same or neighbouring grid cells of size ``cell`` (vectorised)."""
    if len(xy) < 2 or not cell > 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    ix = np.floor((xy[:, 0] - xy[:, 0].min()) / cell).astype(np.int64)
    iy = np.floor((xy[:, 1] - xy[:, 1].min()) / cell).astype(np.int64)
    W = int(iy.max()) + 3
    key = (ix + 1) * W + (iy + 1)
    order = np.argsort(key, kind="stable")
    ks = key[order]
    I, J = [], []
    for dx, dy in ((0, 0), (1, -1), (1, 0), (1, 1), (0, 1)):
        target = key + dx * W + dy
        lo = np.searchsorted(ks, target, side="left")
        hi = np.searchsorted(ks, target, side="right")
        cnt = hi - lo
        if not cnt.any():
            continue
        i = np.repeat(np.arange(len(xy)), cnt)
        start = np.repeat(lo, cnt)
        off = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        j = order[start + off]
        if dx == 0 and dy == 0:
            m = j > i
            i, j = i[m], j[m]
        I.append(i)
        J.append(j)
    if not I:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return np.concatenate(I), np.concatenate(J)


def narrow_gaps(rings: Sequence[np.ndarray], names: Sequence[str], size: SizeFn, frac: float = GAP_FRACTION,
                limit: int = 50) -> List[Problem]:
    """Where rings (closed, without the closing vertex) come within ``frac`` x the element size of themselves
    (a neck: far apart along the ring) or of each other. One problem per place, the narrowest first."""
    xs, ss, owner = [], [], []
    Lr = np.zeros(len(rings))
    for k, r in enumerate(rings):
        r = np.asarray(r, dtype=np.float64)[:, :2]
        if len(r) < 3:
            continue
        xy, s, L = _sample(r, True, size, 0.5 * frac)
        xs.append(xy)
        ss.append(s)
        owner.append(np.full(len(xy), k, dtype=np.int64))
        Lr[k] = L
    if not xs:
        return []
    xy = np.vstack(xs)
    s = np.concatenate(ss)
    own = np.concatenate(owner)
    h = np.asarray(size(xy[:, 0], xy[:, 1]), dtype=np.float64)
    g = frac * h
    i, j = _pairs_within(xy, float(g.max()))
    if i.size == 0:
        return []
    d = np.hypot(xy[i, 0] - xy[j, 0], xy[i, 1] - xy[j, 1])
    lim = np.minimum(g[i], g[j])
    close = d < lim
    same = own[i] == own[j]
    ds = np.abs(s[i] - s[j])
    ds = np.where(same, np.minimum(ds, Lr[own[i]] - ds), np.inf)
    # along the same ring, points closer than this along it are neighbours, not a neck
    far_along = ds > 3.0 * np.maximum(lim, d) + 1e-12
    hit = close & (~same | far_along) & (d > 0)
    i, j, d = i[hit], j[hit], d[hit]
    if i.size == 0:
        return []
    order = np.argsort(d / np.minimum(h[i], h[j]))
    out: List[Problem] = []
    taken: List[Tuple[float, float, float]] = []
    for k in order:
        a, b = int(i[k]), int(j[k])
        x, y = 0.5 * (xy[a] + xy[b])
        hh = float(min(h[a], h[b]))
        if any(math.hypot(x - tx, y - ty) < 3.0 * max(th, hh) for tx, ty, th in taken):
            continue
        taken.append((x, y, hh))
        if own[a] == own[b]:
            out.append(Problem("neck", float(x), float(y), MESSAGES["neck"] % (d[k], hh), float(d[k]), hh,
                               other=names[own[a]]))
        else:
            out.append(Problem("gap", float(x), float(y), MESSAGES["gap"] % (names[own[a]], d[k], names[own[b]], hh),
                               float(d[k]), hh, other=names[own[b]]))
        if len(out) >= limit:
            break
    return out


def _nearest_on_segments(p: np.ndarray, a: np.ndarray, b: np.ndarray, chunk: int = 2048):
    """For each point, the distance to the nearest segment a-b and the nearest point on it."""
    best = np.full(len(p), np.inf)
    q = np.zeros_like(p)
    if len(a) == 0:
        return best, q
    d = b - a
    L2 = np.maximum(np.einsum("ij,ij->i", d, d), 1e-300)
    for lo in range(0, len(p), chunk):
        P = p[lo:lo + chunk]
        t = np.clip(((P[:, None, 0] - a[None, :, 0]) * d[None, :, 0] + (P[:, None, 1] - a[None, :, 1]) * d[None, :, 1])
                    / L2[None, :], 0.0, 1.0)
        qx = a[None, :, 0] + t * d[None, :, 0]
        qy = a[None, :, 1] + t * d[None, :, 1]
        dist = np.hypot(P[:, None, 0] - qx, P[:, None, 1] - qy)
        k = np.argmin(dist, axis=1)
        r = np.arange(len(P))
        best[lo:lo + chunk] = dist[r, k]
        q[lo:lo + chunk, 0] = qx[r, k]
        q[lo:lo + chunk, 1] = qy[r, k]
    return best, q


def near_ends(points, labels: Sequence[str], outline_rings: Sequence[np.ndarray], size: SizeFn,
              frac: float = SNAP_FRACTION) -> List[Problem]:
    """Points (line ends, vertices) closer than ``frac`` x the element size to an outline but not on it."""
    p = np.asarray(points, dtype=np.float64).reshape(-1, 2)
    if len(p) == 0 or not outline_rings:
        return []
    A, B = [], []
    for r in outline_rings:
        r = np.asarray(r, dtype=np.float64)[:, :2]
        if len(r) < 2:
            continue
        A.append(r)
        B.append(np.roll(r, -1, axis=0))
    a, b = np.vstack(A), np.vstack(B)
    dist, q = _nearest_on_segments(p, a, b)
    h = np.asarray(size(p[:, 0], p[:, 1]), dtype=np.float64)
    scale = max(1.0, float(np.abs(a).max()))
    eps = max(1e-9 * scale, 1e-6 * float(np.min(h)))
    out = []
    for k in np.flatnonzero((dist > eps) & (dist < frac * h)):
        out.append(Problem("near", float(p[k, 0]), float(p[k, 1]), MESSAGES["near"] % (labels[k], dist[k], h[k]),
                           float(dist[k]), float(h[k]), True, int(k), (float(q[k, 0]), float(q[k, 1]))))
    return out
