"""Vectorised 2D geometry helpers (NumPy only)."""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np

__all__ = [
    "as_ring",
    "ring_area",
    "points_in_rings",
    "distance_to_segments",
    "segment_intersections",
    "turning_angles",
    "polyline_length",
    "snap_rings_to_nodes",
    "snap_line_to_nodes",
]


def as_ring(coords) -> np.ndarray:
    """(n, 2) array without the closing duplicate vertex."""
    a = np.asarray(coords, dtype=np.float64)[:, :2]
    if len(a) > 1 and np.allclose(a[0], a[-1]):
        a = a[:-1]
    return a


def ring_area(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def points_in_rings(px: np.ndarray, py: np.ndarray, rings: Sequence[np.ndarray]) -> np.ndarray:
    """Even-odd point-in-polygon over several rings (outer rings and holes alike)."""
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    inside = np.zeros(px.shape, dtype=bool)
    for ring in rings:
        r = as_ring(ring)
        x0, y0 = r[:, 0], r[:, 1]
        x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
        xmin, xmax = r[:, 0].min(), r[:, 0].max()
        ymin, ymax = r[:, 1].min(), r[:, 1].max()
        cand = (px >= xmin) & (px <= xmax) & (py >= ymin) & (py <= ymax)
        idx = np.flatnonzero(cand)
        if idx.size == 0:
            continue
        qx, qy = px[idx], py[idx]
        hit = np.zeros(idx.size, dtype=bool)
        for k in range(len(x0)):
            ya, yb = y0[k], y1[k]
            if ya == yb:
                continue
            cross = (ya > qy) != (yb > qy)
            if not cross.any():
                continue
            xint = (x1[k] - x0[k]) * (qy - ya) / (yb - ya) + x0[k]
            hit ^= cross & (qx < xint)
        inside[idx] ^= hit
    return inside


def distance_to_segments(px, py, a: np.ndarray, b: np.ndarray, cutoff: Optional[float] = None) -> np.ndarray:
    """Minimum distance from points to segments a[i]-b[i]; inf beyond ``cutoff``."""
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    out = np.full(px.shape, np.inf)
    if len(a) == 0 or px.size == 0:
        return out
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    big = cutoff if cutoff is not None else np.inf
    for k in range(len(a)):
        ax, ay = a[k]
        bx, by = b[k]
        if np.isfinite(big):
            sel = np.flatnonzero((px >= min(ax, bx) - big) & (px <= max(ax, bx) + big)
                                 & (py >= min(ay, by) - big) & (py <= max(ay, by) + big))
            if sel.size == 0:
                continue
            qx, qy = px[sel], py[sel]
        else:
            sel = slice(None)
            qx, qy = px, py
        dx, dy = bx - ax, by - ay
        l2 = dx * dx + dy * dy
        if l2 == 0:
            d = np.hypot(qx - ax, qy - ay)
        else:
            t = np.clip(((qx - ax) * dx + (qy - ay) * dy) / l2, 0.0, 1.0)
            d = np.hypot(qx - (ax + t * dx), qy - (ay + t * dy))
        out[sel] = np.minimum(out[sel], d)
    return out


def segment_intersections(p: np.ndarray, q: np.ndarray, r: np.ndarray, s: np.ndarray, eps: float = 1e-12):
    """Proper intersections between segments p-q (set A) and r-s (set B).

    Returns (i, j, t, u, x, y): segment i of A meets segment j of B at
    parameters t (on A) and u (on B), point (x, y). Collinear overlaps are
    ignored (handled by snapping).
    """
    p, q, r, s = (np.asarray(v, dtype=np.float64) for v in (p, q, r, s))
    out = []
    if len(p) == 0 or len(r) == 0:
        return out
    bmin = np.minimum(r, s)
    bmax = np.maximum(r, s)
    for i in range(len(p)):
        lo = np.minimum(p[i], q[i])
        hi = np.maximum(p[i], q[i])
        cand = np.flatnonzero((bmax[:, 0] >= lo[0]) & (bmin[:, 0] <= hi[0]) & (bmax[:, 1] >= lo[1]) & (bmin[:, 1] <= hi[1]))
        if cand.size == 0:
            continue
        d1 = q[i] - p[i]
        d2 = s[cand] - r[cand]
        den = d1[0] * d2[:, 1] - d1[1] * d2[:, 0]
        w = r[cand] - p[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w[:, 0] * d2[:, 1] - w[:, 1] * d2[:, 0]) / den
            u = (w[:, 0] * d1[1] - w[:, 1] * d1[0]) / den
        ok = (np.abs(den) > eps) & (t >= -1e-9) & (t <= 1 + 1e-9) & (u >= -1e-9) & (u <= 1 + 1e-9)
        for k in np.flatnonzero(ok):
            tt, uu = float(np.clip(t[k], 0, 1)), float(np.clip(u[k], 0, 1))
            x, y = p[i] + tt * d1
            out.append((i, int(cand[k]), tt, uu, float(x), float(y)))
    return out


def turning_angles(coords: np.ndarray, closed: bool) -> np.ndarray:
    """Absolute change of direction at each vertex, degrees (0 = straight)."""
    a = np.asarray(coords, dtype=np.float64)
    n = len(a)
    ang = np.zeros(n)
    if n < 3:
        return ang
    prev = np.roll(a, 1, axis=0)
    nxt = np.roll(a, -1, axis=0)
    v1 = a - prev
    v2 = nxt - a
    h1 = np.arctan2(v1[:, 1], v1[:, 0])
    h2 = np.arctan2(v2[:, 1], v2[:, 0])
    d = np.degrees(np.abs((h2 - h1 + np.pi) % (2 * np.pi) - np.pi))
    if not closed:
        d[0] = d[-1] = 0.0
    return d


def polyline_length(coords: np.ndarray, closed: bool = False) -> float:
    a = np.asarray(coords, dtype=np.float64)
    if closed:
        a = np.vstack([a, a[:1]])
    return float(np.sum(np.hypot(*np.diff(a, axis=0).T)))


def _snap_vertices(r: np.ndarray, nodes: np.ndarray, tol: float, segments=None, on_tol: float = 0.0) -> np.ndarray:
    r = np.asarray(r, dtype=np.float64)[:, :2].copy()
    if len(nodes) == 0 or len(r) == 0:
        return r
    cand = np.arange(len(r))
    if segments is not None:
        a, b = segments
        d = distance_to_segments(r[:, 0], r[:, 1], a, b, cutoff=max(on_tol, 1e-12) * 2)
        cand = np.flatnonzero(d <= on_tol)
    for i in cand:
        d = np.hypot(nodes[:, 0] - r[i, 0], nodes[:, 1] - r[i, 1])
        k = int(np.argmin(d))
        if 0 < d[k] <= tol:
            r[i] = nodes[k]
    return r


def _drop_repeats(r: np.ndarray, closed: bool) -> np.ndarray:
    keep = np.ones(len(r), dtype=bool)
    keep[1:] = np.any(r[1:] != r[:-1], axis=1)
    r = r[keep]
    if closed and len(r) > 1 and np.all(r[0] == r[-1]):
        r = r[:-1]
    return r


def snap_rings_to_nodes(rings: Sequence[np.ndarray], nodes: np.ndarray, tol: float,
                        segments=None, on_tol: float = 0.0) -> List[np.ndarray]:
    """Move ring vertices within ``tol`` of a node onto it; drop repeated vertices.

    Used after a GEOS difference between the outline and submesh areas: where
    the outline crosses a submesh boundary GEOS inserts a vertex in the middle
    of a submesh edge, which would leave a hanging node. Snapping it to the
    nearest submesh node keeps the two meshes conforming. With ``segments``
    (arrays a, b), only vertices within ``on_tol`` of those segments move.
    """
    nodes = np.asarray(nodes, dtype=np.float64).reshape(-1, 2)
    out: List[np.ndarray] = []
    for ring in rings:
        r = _drop_repeats(_snap_vertices(as_ring(ring), nodes, tol, segments, on_tol), True)
        if len(r) >= 3 and abs(ring_area(r)) > 0:
            out.append(r)
    return out


def snap_line_to_nodes(coords: np.ndarray, nodes: np.ndarray, tol: float, segments=None, on_tol: float = 0.0) -> np.ndarray:
    """Same as :func:`snap_rings_to_nodes` for an open polyline."""
    nodes = np.asarray(nodes, dtype=np.float64).reshape(-1, 2)
    return _drop_repeats(_snap_vertices(coords, nodes, tol, segments, on_tol), False)
