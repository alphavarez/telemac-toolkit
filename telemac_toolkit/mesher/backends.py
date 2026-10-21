"""Delaunay triangulation back-ends for the mesh generator."""

from __future__ import annotations

import struct

import numpy as np

__all__ = ["unique_points", "triangulate_geos", "triangulate_scipy", "default_triangulator"]


def unique_points(xy: np.ndarray):
    """Unique rows and the inverse map (duplicates confuse every Delaunay code)."""
    u, inv = np.unique(np.asarray(xy, dtype=np.float64), axis=0, return_inverse=True)
    return u, inv.ravel()


def _first_index(inv: np.ndarray, n_unique: int) -> np.ndarray:
    first = np.full(n_unique, -1, dtype=np.int64)
    first[inv[::-1]] = np.arange(len(inv) - 1, -1, -1)
    return first


def _map_back(tri_xy: np.ndarray, xy: np.ndarray) -> np.ndarray:
    """Indices of triangle vertex coordinates in xy (unique rows, exact match)."""
    order = np.lexsort((xy[:, 1], xy[:, 0]))
    sx, sy = xy[order, 0], xy[order, 1]
    flat = tri_xy.reshape(-1, 2)
    lo = np.searchsorted(sx, flat[:, 0], side="left")
    hi = np.searchsorted(sx, flat[:, 0], side="right")
    idx = np.full(len(flat), -1, dtype=np.int64)
    simple = (hi - lo) == 1
    idx[simple] = lo[simple]
    for k in np.flatnonzero(~simple):  # several points share x: match y inside the run
        if hi[k] <= lo[k]:
            continue
        run = np.arange(lo[k], hi[k])
        j = np.searchsorted(sy[run], flat[k, 1])
        if j < run.size:
            idx[k] = run[j]
    safe = np.clip(idx, 0, len(sx) - 1)
    ok = (idx >= 0) & (sx[safe] == flat[:, 0]) & (sy[safe] == flat[:, 1])
    if not ok.all():
        raise RuntimeError("triangulation returned %d unknown vertices" % int((~ok).sum()))
    return order[idx].reshape(-1, 3)


def repair_degenerate(xy: np.ndarray, ikle: np.ndarray, rel: float = 1e-10, max_rounds: int = 50) -> np.ndarray:
    """Remove zero-area triangles that some Delaunay codes emit along collinear points.

    A flat triangle (p, m, q) with m between p and q is merged with its
    neighbour (p, q, d) across the long edge into (p, m, d) and (m, q, d); a
    flat triangle on the convex hull is simply dropped. Neighbours are found
    with a sorted edge-key table, so the cost is O(n log n) whatever the
    number of flat triangles.
    """
    ik = np.asarray(ikle, dtype=np.int64).copy()
    npts = max(len(xy), 1)
    for _ in range(max_rounds):
        if len(ik) == 0:
            return ik
        P = xy[ik]
        e = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 1], P[:, 0] - P[:, 2]], axis=1)
        L2 = np.einsum("ijk,ijk->ij", e, e)
        cross = e[:, 0, 0] * (P[:, 2, 1] - P[:, 0, 1]) - e[:, 0, 1] * (P[:, 2, 0] - P[:, 0, 0])
        flat = np.abs(cross) <= rel * L2.max(axis=1)
        if not flat.any():
            return ik
        a = np.stack([ik[:, 0], ik[:, 1], ik[:, 2]], axis=1)
        b = np.stack([ik[:, 1], ik[:, 2], ik[:, 0]], axis=1)
        keys = (np.minimum(a, b) * npts + np.maximum(a, b)).ravel()
        order = np.argsort(keys, kind="stable")
        skeys = keys[order]
        gone = np.zeros(len(ik), dtype=bool)
        new = []
        for t in np.flatnonzero(flat):
            if gone[t]:
                continue
            tri = ik[t].tolist()
            k = int(np.argmax(L2[t]))  # longest edge is tri[k] -> tri[k+1]
            p, q, m = tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]
            key = min(p, q) * npts + max(p, q)
            lo = int(np.searchsorted(skeys, key, side="left"))
            hi = int(np.searchsorted(skeys, key, side="right"))
            others = [int(order[j]) // 3 for j in range(lo, hi) if int(order[j]) // 3 != t]
            if not others:
                gone[t] = True  # on the hull: just drop it
                continue
            nb = [u for u in others if not gone[u] and not flat[u]]
            if not nb:
                continue  # neighbour flat too, or rebuilt in this round: next round
            u = nb[0]
            d = [w for w in ik[u].tolist() if w not in (p, q)][0]
            gone[t] = gone[u] = True
            new += [[p, m, d], [m, q, d]]
        if not gone.any():
            # no progress (chains of flat triangles only): drop them, they cover no area
            return ik[~flat]
        ik = np.vstack([ik[~gone]] + ([np.array(new, dtype=np.int64)] if new else []))
        P = xy[ik]
        cr = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 1, 1] - P[:, 0, 1]) * (P[:, 2, 0] - P[:, 0, 0])
        ik[cr < 0] = ik[cr < 0][:, [0, 2, 1]]
    return ik


def _jittered(u: np.ndarray, rel: float = 1e-7, seed: int = 20240611) -> np.ndarray:
    """Points moved to a local origin plus a tiny deterministic jitter.

    GEOS's incremental Delaunay mishandles exactly collinear points (long
    straight outline edges resampled at even spacing): it emits flat
    triangles and leaves holes. A jitter of 1e-7 of the extent removes exact
    collinearity without changing which edges are Delaunay in practice; the
    caller maps the result back to the original indices.
    """
    lo = u.min(axis=0)
    ext = float(np.max(u.max(axis=0) - lo)) or 1.0
    rng = np.random.default_rng(seed)
    return (u - lo) + rng.uniform(-1.0, 1.0, size=u.shape) * rel * ext


def triangulate_geos(xy: np.ndarray) -> np.ndarray:
    """Delaunay triangulation with GEOS through QGIS (always available in QGIS)."""
    from qgis.core import QgsGeometry, QgsPointXY

    u, inv = unique_points(xy)
    j = _jittered(u)
    mp = QgsGeometry.fromMultiPointXY([QgsPointXY(float(a), float(b)) for a, b in j])
    tri = mp.delaunayTriangulation(0.0, False)
    wkb = bytes(tri.asWkb())
    endian = "<" if wkb[0] == 1 else ">"
    n = struct.unpack_from(endian + "I", wkb, 5)[0]
    rec = 1 + 4 + 4 + 4 + 4 * 16  # polygon, 1 ring, 4 points (closed)
    arr = np.frombuffer(wkb, dtype=np.uint8, count=n * rec, offset=9).reshape(n, rec)
    coords = arr[:, 13:13 + 48].copy().view(endian + "f8").reshape(n, 3, 2)
    ikle_u = _map_back(coords, j)
    # orientation in the true coordinates, then any triangle still flat there
    P = u[ikle_u]
    cr = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 1, 1] - P[:, 0, 1]) * (P[:, 2, 0] - P[:, 0, 0])
    ikle_u[cr < 0] = ikle_u[cr < 0][:, [0, 2, 1]]
    ikle_u = repair_degenerate(u, ikle_u)
    return _first_index(inv, len(u))[ikle_u]


def triangulate_scipy(xy: np.ndarray) -> np.ndarray:
    from scipy.spatial import Delaunay

    u, inv = unique_points(xy)
    simp = Delaunay(u).simplices.astype(np.int64)
    return _first_index(inv, len(u))[simp]


def default_triangulator():
    try:
        import qgis.core  # noqa: F401

        return triangulate_geos
    except ImportError:
        return triangulate_scipy
