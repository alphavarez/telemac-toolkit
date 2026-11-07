"""Measures of a triangular mesh, to compare meshing methods (or meshes made elsewhere) on equal terms.

``mesh_metrics`` returns one dictionary per mesh:

Element shape
    ``angle_min``, ``angle_p05`` (degrees), ``share_below_30``, ``share_below_20`` (fractions),
    ``angle_max``, ``share_above_120``; ``aspect_p95``, ``aspect_max``: circumradius over twice the
    inradius (1 for an equilateral triangle); ``edge_min``, ``edge_median``, ``edge_max`` (m).
Size changes
    ``grading_p95``, ``grading_max``: area ratio of the two elements sharing each interior edge
    (1 = equal neighbours). Sudden jumps make a model less accurate and slower to converge.
Outline
    ``outline_dev_max``, ``outline_dev_mean`` (m): distance from the drawn outline to the mesh
    boundary, sampled along the outline. Zero when the mesh follows it exactly.
Terrain (with a DEM)
    The mesh carries the bed at its nodes and TELEMAC interpolates it linearly in each element. Seven
    points per element compare that plane with the DEM: ``bed_rms``, ``bed_p95``, ``bed_max`` (m).
    With a water ``level``: ``volume_error`` (fraction) of the water volume under that level.
Time step
    ``dt_courant1`` (s): the largest time step at a Courant number of 1, from each element's smallest
    altitude and the celerity |U| + sqrt(g h), with h from ``level`` minus the bed (wet elements only)
    or a given ``depth``. Small elements in deep fast water limit it.
Size
    ``nodes``, ``elements``, ``area`` (m²), ``nodes_per_km2``.
"""

from __future__ import annotations

from typing import Callable, Dict, Optional, Sequence

import numpy as np

from ..telemac_io.mesh import boundary_loops

__all__ = ["mesh_metrics", "METRIC_INFO", "best_is_low"]

#: (label, unit, lower is better?)
METRIC_INFO = {
    "nodes": ("Nodes", "", None),
    "elements": ("Elements", "", None),
    "seconds": ("Meshing time", "s", True),
    "angle_min": ("Smallest angle", "°", False),
    "angle_p05": ("Angle, 5th percentile", "°", False),
    "share_below_30": ("Elements below 30°", "%", True),
    "share_below_20": ("Elements below 20°", "%", True),
    "angle_max": ("Largest angle", "°", True),
    "share_above_120": ("Elements above 120°", "%", True),
    "aspect_p95": ("Aspect ratio, 95th percentile", "", True),
    "aspect_max": ("Aspect ratio, largest", "", True),
    "edge_min": ("Shortest edge", "m", None),
    "edge_median": ("Median edge", "m", None),
    "edge_max": ("Longest edge", "m", None),
    "grading_p95": ("Neighbour area ratio, 95th percentile", "", True),
    "grading_max": ("Neighbour area ratio, largest", "", True),
    "outline_dev_max": ("Outline deviation, largest", "m", True),
    "outline_dev_mean": ("Outline deviation, mean", "m", True),
    "bed_rms": ("Bed error vs DEM, RMS", "m", True),
    "bed_p95": ("Bed error vs DEM, 95th percentile", "m", True),
    "bed_max": ("Bed error vs DEM, largest", "m", True),
    "volume_error": ("Water volume error under the level", "%", True),
    "dt_courant1": ("Time step at Courant 1", "s", False),
    "area": ("Area", "m²", None),
    "nodes_per_km2": ("Nodes per km²", "", None),
}


def best_is_low(key: str) -> Optional[bool]:
    return METRIC_INFO.get(key, ("", "", None))[2]


def _angles(P):
    out = []
    for i in range(3):
        a = P[:, (i + 1) % 3] - P[:, i]
        b = P[:, (i + 2) % 3] - P[:, i]
        c = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-300)
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.stack(out, axis=1)


def _dist_to_segments(px, py, ax, ay, bx, by, chunk=2048):
    out = np.empty(len(px))
    ex, ey = bx - ax, by - ay
    l2 = np.maximum(ex * ex + ey * ey, 1e-300)
    for s in range(0, len(px), chunk):
        X = px[s:s + chunk, None] - ax[None, :]
        Y = py[s:s + chunk, None] - ay[None, :]
        t = np.clip((X * ex + Y * ey) / l2, 0, 1)
        d = (X - t * ex) ** 2 + (Y - t * ey) ** 2
        out[s:s + chunk] = np.sqrt(d.min(axis=1))
    return out


# seven sample points per element (barycentric) with equal weights
_BARY = np.array([[1 / 3, 1 / 3, 1 / 3], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5],
                  [2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])


def mesh_metrics(x, y, ikle, outline: Optional[Sequence[np.ndarray]] = None, dem: Optional[Callable] = None,
                 z: Optional[np.ndarray] = None, level: Optional[float] = None, depth: Optional[float] = None,
                 velocity: float = 1.0, g: float = 9.81) -> Dict[str, float]:
    """Measures of the mesh (see the module text). ``z``: bed at the nodes (default: the DEM at the nodes)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ik = np.asarray(ikle, dtype=np.int64)
    P = np.stack([x[ik], y[ik]], axis=2)  # (n, 3, 2)
    e = [P[:, (i + 1) % 3] - P[:, i] for i in range(3)]
    L = np.stack([np.hypot(v[:, 0], v[:, 1]) for v in e], axis=1)
    area = 0.5 * np.abs(e[0][:, 0] * (-e[2][:, 1]) - e[0][:, 1] * (-e[2][:, 0]))
    ang = _angles(P)
    amin, amax = ang.min(axis=1), ang.max(axis=1)
    s = L.sum(axis=1) / 2
    inr = area / np.maximum(s, 1e-300)
    circ = L.prod(axis=1) / np.maximum(4 * area, 1e-300)
    aspect = circ / np.maximum(2 * inr, 1e-300)
    out: Dict[str, float] = {
        "nodes": float(len(x)), "elements": float(len(ik)),
        "angle_min": float(amin.min()), "angle_p05": float(np.percentile(amin, 5)),
        "share_below_30": float(np.mean(amin < 30.0)), "share_below_20": float(np.mean(amin < 20.0)),
        "angle_max": float(amax.max()), "share_above_120": float(np.mean(amax > 120.0)),
        "aspect_p95": float(np.percentile(aspect, 95)), "aspect_max": float(aspect.max()),
        "edge_min": float(L.min()), "edge_median": float(np.median(L)), "edge_max": float(L.max()),
        "area": float(area.sum()),
    }
    out["nodes_per_km2"] = out["nodes"] / max(out["area"] / 1e6, 1e-300)
    # neighbour area ratio across interior edges
    ed = np.sort(np.concatenate([ik[:, [0, 1]], ik[:, [1, 2]], ik[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(ik)), 3)
    order = np.lexsort((ed[:, 1], ed[:, 0]))
    ed, owner = ed[order], owner[order]
    same = np.all(ed[1:] == ed[:-1], axis=1)
    a1, a2 = area[owner[:-1][same]], area[owner[1:][same]]
    ratio = np.maximum(a1, a2) / np.maximum(np.minimum(a1, a2), 1e-300)
    out["grading_p95"] = float(np.percentile(ratio, 95)) if ratio.size else 1.0
    out["grading_max"] = float(ratio.max()) if ratio.size else 1.0
    # outline
    if outline:
        loops = boundary_loops(ik, x, y)
        bx = np.concatenate([x[lp] for lp in loops])
        by = np.concatenate([y[lp] for lp in loops])
        nx = np.concatenate([x[np.roll(lp, -1)] for lp in loops])
        ny = np.concatenate([y[np.roll(lp, -1)] for lp in loops])
        step = 0.25 * out["edge_median"]
        sx, sy = [], []
        for r in outline:
            r = np.asarray(r, dtype=np.float64)[:, :2]
            if np.hypot(*(r[0] - r[-1])) > 0:
                r = np.vstack([r, r[:1]])
            for a, b in zip(r[:-1], r[1:]):
                n = max(1, int(np.hypot(*(b - a)) / step))
                t = np.arange(n) / n
                sx.append(a[0] + t * (b[0] - a[0]))
                sy.append(a[1] + t * (b[1] - a[1]))
        sx, sy = np.concatenate(sx), np.concatenate(sy)
        d = _dist_to_segments(sx, sy, bx, by, nx, ny)
        out["outline_dev_max"] = float(d.max())
        out["outline_dev_mean"] = float(d.mean())
    # terrain
    zb = None
    if z is not None:
        zb = np.asarray(z, dtype=np.float64)
    elif dem is not None:
        zb = np.asarray(dem(x, y), dtype=np.float64)
    if dem is not None and zb is not None:
        sx = np.einsum("kj,nj->nk", _BARY, x[ik])
        sy = np.einsum("kj,nj->nk", _BARY, y[ik])
        zl = np.einsum("kj,nj->nk", _BARY, zb[ik])
        zd = np.asarray(dem(sx.ravel(), sy.ravel()), dtype=np.float64).reshape(sx.shape)
        ok = np.isfinite(zd) & np.isfinite(zl)
        err = np.abs(zl - zd)[ok]
        if err.size:
            out["bed_rms"] = float(np.sqrt(np.mean(err ** 2)))
            out["bed_p95"] = float(np.percentile(err, 95))
            out["bed_max"] = float(err.max())
        if level is not None:
            w = (area / 7.0)[:, None] * np.ones_like(zl)
            v_mesh = float(np.sum(np.where(ok, w * np.maximum(level - zl, 0), 0)))
            v_dem = float(np.sum(np.where(ok, w * np.maximum(level - zd, 0), 0)))
            if v_dem > 0:
                out["volume_error"] = (v_mesh - v_dem) / v_dem
    # time step at Courant 1
    alt = 2 * area / np.maximum(L.max(axis=1), 1e-300)
    h = None
    if level is not None and zb is not None:
        h = np.maximum(level - zb[ik].mean(axis=1), 0.0)
        wet = h > 0.05
    elif depth is not None:
        h = np.full(len(ik), float(depth))
        wet = np.ones(len(ik), dtype=bool)
    if h is not None and wet.any():
        c = abs(velocity) + np.sqrt(g * h[wet])
        out["dt_courant1"] = float(np.min(alt[wet] / c))
    return out
