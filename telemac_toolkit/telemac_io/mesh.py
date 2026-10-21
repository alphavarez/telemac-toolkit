"""Mesh utilities for triangular TELEMAC meshes (pure NumPy).

* :func:`boundary_loops`  ordered boundary node loops (outer first, then islands)
* :class:`TriangleLocator` fast point-in-triangle search with barycentric weights
* :func:`sample_polyline` regular points along a line, with chainage
* :func:`mesh_statistics` edge lengths, areas and angles
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

__all__ = [
    "unique_edges",
    "boundary_edges",
    "boundary_loops",
    "TriangleLocator",
    "sample_polyline",
    "mesh_statistics",
    "signed_areas",
]


def signed_areas(x: np.ndarray, y: np.ndarray, ikle: np.ndarray) -> np.ndarray:
    """Signed area of every triangle (positive when counter-clockwise)."""
    a, b, c = ikle[:, 0], ikle[:, 1], ikle[:, 2]
    return 0.5 * ((x[b] - x[a]) * (y[c] - y[a]) - (x[c] - x[a]) * (y[b] - y[a]))


def unique_edges(ikle: np.ndarray):
    """Unique undirected edges and how many triangles share each one."""
    e = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    e.sort(axis=1)
    edges, counts = np.unique(e, axis=0, return_counts=True)
    return edges, counts


def boundary_edges(ikle: np.ndarray) -> np.ndarray:
    """Directed boundary edges, oriented as in their triangle.

    For counter-clockwise triangles the outer boundary comes out
    counter-clockwise and islands clockwise.
    """
    d = np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]])
    s = np.sort(d, axis=1)
    _, inv, counts = np.unique(s, axis=0, return_inverse=True, return_counts=True)
    return d[counts[inv.ravel()] == 1]


def boundary_loops(ikle: np.ndarray, x: Optional[np.ndarray] = None, y: Optional[np.ndarray] = None) -> list:
    """Ordered loops of boundary nodes.

    Each loop is an array of node indices (not closed: first != last). When x
    and y are given, triangles are oriented counter-clockwise first and the
    loop with the largest enclosed area (the outer boundary) comes first.
    Each loop starts at its lowest node index.
    """
    ikle = np.asarray(ikle, dtype=np.int64)
    if ikle.size == 0:
        return []
    if x is not None and y is not None:
        ar = signed_areas(x, y, ikle)
        ikle = ikle.copy()
        neg = ar < 0
        ikle[neg] = ikle[neg][:, [0, 2, 1]]
    edges = boundary_edges(ikle)
    if edges.size == 0:
        return []
    nxt = {}
    for a, b in edges:
        # pinch nodes (two boundary loops touching) keep several successors
        nxt.setdefault(int(a), []).append(int(b))
    loops = []
    remaining = {(int(a), int(b)) for a, b in edges}
    while remaining:
        start = min(a for a, _ in remaining)
        loop = [start]
        cur = start
        while True:
            succ = [b for b in nxt[cur] if (cur, b) in remaining]
            if not succ:
                break
            b = succ[0]
            remaining.discard((cur, b))
            if b == start:
                break
            loop.append(b)
            cur = b
        loops.append(np.array(loop, dtype=np.int64))
    if x is not None and y is not None and len(loops) > 1:
        def area(lp):
            xs, ys = x[lp], y[lp]
            return 0.5 * np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys)
        loops.sort(key=lambda lp: -area(lp))
    return loops


class TriangleLocator:
    """Find the triangle containing a point and its barycentric weights.

    Uses a uniform bin grid over the mesh extent; each triangle is registered
    in every bin its bounding box touches. Build cost is O(ntri), query cost is
    a handful of triangles per point.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, ikle: np.ndarray, bins_per_tri: float = 1.0):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.ikle = np.asarray(ikle, dtype=np.int64)[:, :3]
        nt = self.ikle.shape[0]
        if nt == 0:
            raise ValueError("mesh has no triangles")
        tx = self.x[self.ikle]
        ty = self.y[self.ikle]
        self.xmin, self.xmax = float(tx.min()), float(tx.max())
        self.ymin, self.ymax = float(ty.min()), float(ty.max())
        w = max(self.xmax - self.xmin, 1e-12)
        h = max(self.ymax - self.ymin, 1e-12)
        nb = max(1, int(nt * bins_per_tri))
        cell = np.sqrt(w * h / nb)
        self.nx = max(1, int(np.ceil(w / cell)))
        self.ny = max(1, int(np.ceil(h / cell)))
        self.dx = w / self.nx
        self.dy = h / self.ny

        ix0 = self._ix(tx.min(axis=1))
        ix1 = self._ix(tx.max(axis=1))
        iy0 = self._iy(ty.min(axis=1))
        iy1 = self._iy(ty.max(axis=1))
        cx = ix1 - ix0 + 1
        cy = iy1 - iy0 + 1
        cnt = cx * cy
        tri = np.repeat(np.arange(nt), cnt)
        # position k inside each triangle's block of bins
        start = np.repeat(np.cumsum(cnt) - cnt, cnt)
        k = np.arange(tri.size) - start
        bx = ix0[tri] + k % cx[tri]
        by = iy0[tri] + k // cx[tri]
        b = by * self.nx + bx
        order = np.argsort(b, kind="stable")
        self._tri = tri[order]
        self._ptr = np.zeros(self.nx * self.ny + 1, dtype=np.int64)
        np.add.at(self._ptr, b + 1, 1)
        np.cumsum(self._ptr, out=self._ptr)
        # precompute barycentric transforms
        a = self.ikle
        x0, y0 = self.x[a[:, 0]], self.y[a[:, 0]]
        self._x0, self._y0 = x0, y0
        self._e1 = np.stack([self.x[a[:, 1]] - x0, self.y[a[:, 1]] - y0], axis=1)
        self._e2 = np.stack([self.x[a[:, 2]] - x0, self.y[a[:, 2]] - y0], axis=1)
        self._det = self._e1[:, 0] * self._e2[:, 1] - self._e1[:, 1] * self._e2[:, 0]

    def _ix(self, v):
        return np.clip(((np.asarray(v) - self.xmin) / self.dx).astype(np.int64), 0, self.nx - 1)

    def _iy(self, v):
        return np.clip(((np.asarray(v) - self.ymin) / self.dy).astype(np.int64), 0, self.ny - 1)

    def locate(self, px: float, py: float, eps: float = 1e-9):
        """(triangle index, nodes[3], weights[3]) or None when outside the mesh."""
        if not (self.xmin <= px <= self.xmax and self.ymin <= py <= self.ymax):
            return None
        b = int(self._iy(py)) * self.nx + int(self._ix(px))
        cand = self._tri[self._ptr[b]:self._ptr[b + 1]]
        if cand.size == 0:
            return None
        dx = px - self._x0[cand]
        dy = py - self._y0[cand]
        det = self._det[cand]
        ok = np.abs(det) > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            l1 = (dx * self._e2[cand, 1] - dy * self._e2[cand, 0]) / det
            l2 = (self._e1[cand, 0] * dy - self._e1[cand, 1] * dx) / det
        l0 = 1.0 - l1 - l2
        inside = ok & (l0 >= -eps) & (l1 >= -eps) & (l2 >= -eps)
        hit = np.flatnonzero(inside)
        if hit.size == 0:
            return None
        j = hit[0]
        t = int(cand[j])
        w = np.clip(np.array([l0[j], l1[j], l2[j]]), 0.0, 1.0)
        w /= w.sum()
        return t, self.ikle[t].copy(), w

    def locate_many(self, px: Sequence[float], py: Sequence[float]):
        """Vectorised wrapper: (tri[n], nodes[n,3], weights[n,3]); tri = -1 outside."""
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        n = px.size
        tri = np.full(n, -1, dtype=np.int64)
        nodes = np.zeros((n, 3), dtype=np.int64)
        wts = np.zeros((n, 3))
        for i in range(n):
            r = self.locate(px[i], py[i])
            if r is not None:
                tri[i], nodes[i], wts[i] = r
        return tri, nodes, wts

    def interpolate(self, values: np.ndarray, px, py) -> np.ndarray:
        """Linear interpolation of node values at points (NaN outside)."""
        tri, nodes, w = self.locate_many(np.atleast_1d(px), np.atleast_1d(py))
        out = np.einsum("ij,ij->i", np.asarray(values)[nodes], w)
        out[tri < 0] = np.nan
        return out

    def nearest_node(self, px: float, py: float) -> int:
        d2 = (self.x - px) ** 2 + (self.y - py) ** 2
        return int(np.argmin(d2))


def sample_polyline(xy: Sequence[Sequence[float]], step: Optional[float] = None, n: Optional[int] = None):
    """Points along a polyline at a regular step (or n points).

    Returns (px, py, chainage). Every vertex of the line is kept, so profile
    breaks at the vertices are not smoothed away.
    """
    pts = np.asarray(xy, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[0] < 2:
        raise ValueError("a polyline needs at least two vertices")
    seg = np.diff(pts, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    total = cum[-1]
    if total <= 0:
        raise ValueError("polyline has zero length")
    if step is None:
        n = n or 200
        step = total / max(n - 1, 1)
    s = np.arange(0.0, total, step)
    s = np.unique(np.concatenate([s, cum, [total]]))
    px = np.interp(s, cum, pts[:, 0])
    py = np.interp(s, cum, pts[:, 1])
    return px, py, s


def mesh_statistics(x: np.ndarray, y: np.ndarray, ikle: np.ndarray) -> dict:
    """Basic quality figures for a triangle mesh."""
    ikle = np.asarray(ikle, dtype=np.int64)[:, :3]
    ar = signed_areas(x, y, ikle)
    edges, counts = unique_edges(ikle)
    el = np.hypot(x[edges[:, 0]] - x[edges[:, 1]], y[edges[:, 0]] - y[edges[:, 1]])
    # interior angles
    P = np.stack([x[ikle], y[ikle]], axis=2)  # (nt, 3, 2)
    ang = []
    for i in range(3):
        a = P[:, (i + 1) % 3] - P[:, i]
        b = P[:, (i + 2) % 3] - P[:, i]
        cosv = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-300)
        ang.append(np.degrees(np.arccos(np.clip(cosv, -1, 1))))
    ang = np.stack(ang, axis=1)
    on_boundary = np.zeros(x.size, dtype=bool)
    on_boundary[edges[counts == 1].ravel()] = True
    over = int(np.sum(on_boundary[ikle].all(axis=1)))
    return {
        "nodes": int(x.size),
        "elements": int(ikle.shape[0]),
        "area_total": float(np.abs(ar).sum()),
        "edge_min": float(el.min()) if el.size else 0.0,
        "edge_mean": float(el.mean()) if el.size else 0.0,
        "edge_max": float(el.max()) if el.size else 0.0,
        "angle_min": float(ang.min()) if ang.size else 0.0,
        "angle_max": float(ang.max()) if ang.size else 0.0,
        "clockwise_elements": int(np.sum(ar < 0)),
        "degenerate_elements": int(np.sum(np.abs(ar) < 1e-12)),
        "boundary_nodes": int(on_boundary.sum()),
        "overconstrained_elements": over,
    }
