"""Breaklines from a DEM: crests (levee and bank tops, road embankments) and toes (bank feet, channel edges).

A crest is where the terrain bends down on both sides (strongly negative curvature across it), a toe where it
bends up (strongly positive). A line is worth drawing where an element of the default size ``size`` would miss
the bend by more than the tolerance ``tol``: |curvature| > 8 tol / size² (the same estimate as the terrain
sizes, mesher.demsize). The cells beyond that threshold are thinned to one-cell-wide skeletons, traced into
polylines, simplified, and the short ones dropped.

Only NumPy is used; the DEM is any callable ``dem(x, y) -> z``.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Tuple

import numpy as np

__all__ = ["breaklines_from_dem", "thin", "trace", "simplify"]

_NB = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]  # P2..P9, clockwise from north


def thin(mask: np.ndarray, max_iter: int = 200) -> np.ndarray:
    """Zhang-Suen thinning of a boolean image to one-pixel-wide lines."""
    img = np.pad(mask.astype(np.uint8), 1)
    for _ in range(max_iter):
        changed = False
        for step in (0, 1):
            P = [np.roll(np.roll(img, -di, axis=0), -dj, axis=1) for di, dj in _NB]
            B = sum(p.astype(np.int16) for p in P)
            seq = P + [P[0]]
            A = sum(((seq[k] == 0) & (seq[k + 1] == 1)).astype(np.int16) for k in range(8))
            p2, p4, p6, p8 = P[0], P[2], P[4], P[6]
            if step == 0:
                c = (p2 * p4 * p6 == 0) & (p4 * p6 * p8 == 0)
            else:
                c = (p2 * p4 * p8 == 0) & (p2 * p6 * p8 == 0)
            rm = (img == 1) & (B >= 2) & (B <= 6) & (A == 1) & c
            rm[0, :] = rm[-1, :] = False
            rm[:, 0] = rm[:, -1] = False
            if rm.any():
                img[rm] = 0
                changed = True
        if not changed:
            break
    return img[1:-1, 1:-1].astype(bool)


def trace(skel: np.ndarray) -> List[np.ndarray]:
    """Polylines (row, col arrays) of a one-pixel skeleton, split at junctions."""
    pts = set(zip(*np.nonzero(skel)))
    if not pts:
        return []

    def nbrs(p):
        i, j = p
        return [(i + di, j + dj) for di, dj in _NB if (i + di, j + dj) in pts]

    deg = {p: len(nbrs(p)) for p in pts}
    nodes = {p for p, d in deg.items() if d != 2}
    used_edges = set()
    lines: List[np.ndarray] = []

    def walk(start, nxt):
        path = [start, nxt]
        used_edges.add(frozenset((start, nxt)))
        prev, cur = start, nxt
        while cur not in nodes:
            cand = [q for q in nbrs(cur) if q != prev and frozenset((cur, q)) not in used_edges]
            if not cand:
                break
            q = cand[0]
            used_edges.add(frozenset((cur, q)))
            path.append(q)
            prev, cur = cur, q
            if cur == start:
                break
        return np.array(path)

    for s in nodes:
        for q in nbrs(s):
            if frozenset((s, q)) not in used_edges:
                lines.append(walk(s, q))
    # closed loops without any end or junction
    rest = [p for p in pts if p not in nodes and not any(frozenset((p, q)) in used_edges for q in nbrs(p))]
    seen = set()
    for p in rest:
        if p in seen:
            continue
        nb = nbrs(p)
        if not nb:
            continue
        path = walk(p, nb[0])
        seen.update(map(tuple, path.tolist()))
        lines.append(path)
    return lines


def simplify(c: np.ndarray, tol: float) -> np.ndarray:
    """Douglas-Peucker simplification of a polyline (n, 2)."""
    if len(c) < 3:
        return c
    keep = np.zeros(len(c), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(c) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        p, q = c[a], c[b]
        d = q - p
        L = math.hypot(*d)
        seg = c[a + 1:b]
        if L == 0:
            dist = np.hypot(*(seg - p).T)
        else:
            dist = np.abs(d[0] * (seg[:, 1] - p[1]) - d[1] * (seg[:, 0] - p[0])) / L
        k = int(np.argmax(dist))
        if dist[k] > tol:
            m = a + 1 + k
            keep[m] = True
            stack += [(a, m), (m, b)]
    return c[keep]


def breaklines_from_dem(dem: Callable, bbox: Tuple[float, float, float, float], tol: float, size: float,
                        cell: float, min_length: float = None, max_cells: int = 4_000_000,
                        smooth: int = 1) -> Dict[str, List[np.ndarray]]:
    """{"crest": [polylines], "toe": [polylines]} in map coordinates.

    ``tol``: vertical tolerance (m); ``size``: the element size the lines must help (usually the default size);
    ``cell``: grid spacing (about the DEM resolution); ``min_length``: shortest line kept (default 3 x size).
    """
    x0, y0, x1, y1 = bbox
    w, h = max(x1 - x0, 1e-9), max(y1 - y0, 1e-9)
    c = max(float(cell), math.sqrt(w * h / max_cells))
    nx = int(math.ceil(w / c)) + 1
    ny = int(math.ceil(h / c)) + 1
    X, Y = np.meshgrid(x0 + np.arange(nx) * c, y0 + np.arange(ny) * c)
    z = np.asarray(dem(X.ravel(), Y.ravel()), dtype=np.float64).reshape(ny, nx)
    # signed principal curvature of largest magnitude
    zxx = np.zeros(z.shape)
    zyy = np.zeros(z.shape)
    zxy = np.zeros(z.shape)
    zxx[:, 1:-1] = (z[:, 2:] - 2 * z[:, 1:-1] + z[:, :-2]) / c ** 2
    zyy[1:-1, :] = (z[2:, :] - 2 * z[1:-1, :] + z[:-2, :]) / c ** 2
    zxy[1:-1, 1:-1] = (z[2:, 2:] - z[2:, :-2] - z[:-2, 2:] + z[:-2, :-2]) / (4 * c ** 2)
    tr = 0.5 * (zxx + zyy)
    disc = np.sqrt(np.maximum(0.25 * (zxx - zyy) ** 2 + zxy ** 2, 0.0))
    l1, l2 = tr + disc, tr - disc
    lam = np.where(np.abs(l1) >= np.abs(l2), l1, l2)
    lam = np.where(np.isfinite(lam), lam, 0.0)
    for _ in range(max(int(smooth), 0)):  # a 3 x 3 mean keeps single noisy DEM cells from making lines
        p = np.pad(lam, 1, mode="edge")
        lam = sum(p[1 + di:1 + di + ny, 1 + dj:1 + dj + nx] for di in (-1, 0, 1) for dj in (-1, 0, 1)) / 9.0
    lam[0, :] = lam[-1, :] = 0
    lam[:, 0] = lam[:, -1] = 0
    kt = 8.0 * tol / size ** 2
    min_length = 3.0 * size if min_length is None else min_length
    out: Dict[str, List[np.ndarray]] = {}
    for name, mask in (("crest", lam < -kt), ("toe", lam > kt)):
        lines = []
        for path in trace(thin(mask)):
            xy = np.stack([x0 + path[:, 1] * c, y0 + path[:, 0] * c], axis=1)
            xy = simplify(xy, c)
            if len(xy) >= 2 and float(np.hypot(*np.diff(xy, axis=0).T).sum()) >= min_length:
                lines.append(xy)
        out[name] = lines
    return out
