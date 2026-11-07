"""Element sizes from the terrain: small where linear elements would misrepresent the DEM, large where it is flat.

TELEMAC carries the bed at the nodes and interpolates it linearly inside each element. Over an element of
size h the linear interpolation of a surface misses its curvature by about h² |κ| / 8, with |κ| the largest
absolute curvature of the DEM (the largest eigenvalue of its Hessian). For a vertical tolerance ``tol`` the
size is therefore

    h = sqrt(8 tol / |κ|),   clipped to [hmin, hmax].

Banks, levees, road embankments and channel edges get small elements; flat floodplain keeps the default size.
The sizes are then limited so that they do not grow faster than the growth rate from one place to the next
(a size may exceed a neighbour's by at most (growth - 1) times the distance between them).

The DEM is any callable ``dem(x, y) -> z`` (NaN outside); the result is a :class:`GridDensity` the size
field takes like the other refinement sources. Only NumPy is used.
"""

from __future__ import annotations

import math
from typing import Callable, Optional, Tuple

import numpy as np

from .sizefield import GridDensity

__all__ = ["terrain_sizes", "limit_gradation", "curvature", "apply_floor", "field_sizes"]


def curvature(z: np.ndarray, cell: float) -> np.ndarray:
    """Largest absolute eigenvalue of the Hessian of ``z`` (grid spacing ``cell``); NaN where z is missing."""
    zxx = np.full(z.shape, np.nan)
    zyy = np.full(z.shape, np.nan)
    zxy = np.full(z.shape, np.nan)
    zxx[:, 1:-1] = (z[:, 2:] - 2 * z[:, 1:-1] + z[:, :-2]) / cell ** 2
    zyy[1:-1, :] = (z[2:, :] - 2 * z[1:-1, :] + z[:-2, :]) / cell ** 2
    zxy[1:-1, 1:-1] = (z[2:, 2:] - z[2:, :-2] - z[:-2, 2:] + z[:-2, :-2]) / (4 * cell ** 2)
    # edges: copy the nearest interior value
    for a in (zxx, zyy, zxy):
        a[:, 0], a[:, -1] = a[:, 1], a[:, -2]
        a[0, :], a[-1, :] = a[1, :], a[-2, :]
    tr = 0.5 * (zxx + zyy)
    disc = np.sqrt(np.maximum(0.25 * (zxx - zyy) ** 2 + zxy ** 2, 0.0))
    return np.maximum(np.abs(tr + disc), np.abs(tr - disc))


def limit_gradation(h: np.ndarray, cell: float, growth: float, passes: int = 2) -> np.ndarray:
    """Sizes no larger than a neighbour's plus (growth - 1) times the distance (along rows and columns).

    Each pass is h[i] = min_j (h[j] + g |i - j| cell) along every row and column, computed with running
    minima, so it is linear in the grid size.
    """
    g = max(growth - 1.0, 0.0) * cell
    if g <= 0:
        return np.full(h.shape, np.nanmin(h))
    h = h.copy()
    for _ in range(passes):
        for axis in (1, 0):
            n = h.shape[axis]
            ramp = g * np.arange(n, dtype=np.float64)
            shape = [1, 1]
            shape[axis] = n
            ramp = ramp.reshape(shape)
            fwd = np.minimum.accumulate(h - ramp, axis=axis) + ramp
            rev = np.flip(np.minimum.accumulate(np.flip(h + ramp, axis=axis), axis=axis), axis=axis) - ramp
            h = np.minimum(h, np.minimum(fwd, rev))
    return h


def apply_floor(gd: GridDensity, seg_a: np.ndarray, seg_b: np.ndarray, sizes: np.ndarray, growth: float) -> GridDensity:
    """Raise the sizes near fixed edges (a structured channel's boundary) to those edges' length.

    A size s on a fixed edge holds the grid at s - (growth - 1) d at distance d, so the elements next to the
    edge have its length and shrink towards the terrain sizes no faster than the growth rate (without it,
    small terrain sizes against a long fixed edge make thin triangles).
    """
    from .geometry import distance_to_segments

    v = np.asarray(gd.values, dtype=np.float64).copy()
    ny, nx = v.shape
    g = max(growth - 1.0, 1e-6)
    X, Y = np.meshgrid(gd.x0 + np.arange(nx) * gd.dx, gd.y0 + np.arange(ny) * gd.dy)
    px, py = X.ravel(), Y.ravel()
    sizes = np.asarray(sizes, dtype=np.float64)
    floor = np.zeros(px.size)
    for s in np.unique(np.round(sizes, 6)):
        sel = np.isclose(sizes, s)
        reach = s / g
        d = distance_to_segments(px, py, seg_a[sel], seg_b[sel], cutoff=reach)
        floor = np.maximum(floor, np.where(np.isfinite(d), s - g * d, 0.0))
    v = np.maximum(v, floor.reshape(ny, nx))
    return GridDensity(gd.x0, gd.y0, gd.dx, gd.dy, v)


def terrain_sizes(dem: Callable, bbox: Tuple[float, float, float, float], tol: float, hmin: float, hmax: float,
                  growth: float = 1.2, cell: Optional[float] = None, max_cells: int = 4_000_000,
                  smooth: int = 1) -> GridDensity:
    """A :class:`GridDensity` of terrain-driven sizes over ``bbox`` (x0, y0, x1, y1).

    ``cell``: grid spacing (default: half the smallest size, coarsened to stay under ``max_cells``).
    ``smooth``: passes of a 3 x 3 mean on the curvature, which keeps single noisy DEM cells from refining.
    """
    if tol <= 0:
        raise ValueError("the terrain tolerance must be positive")
    if not 0 < hmin <= hmax:
        raise ValueError("terrain sizes need 0 < smallest size <= largest size")
    x0, y0, x1, y1 = bbox
    w, hgt = max(x1 - x0, 1e-9), max(y1 - y0, 1e-9)
    c = float(cell) if cell else 0.5 * hmin
    c = max(c, math.sqrt(w * hgt / max_cells))
    nx = int(math.ceil(w / c)) + 1
    ny = int(math.ceil(hgt / c)) + 1
    gx = x0 + np.arange(nx) * c
    gy = y0 + np.arange(ny) * c
    X, Y = np.meshgrid(gx, gy)
    z = np.asarray(dem(X.ravel(), Y.ravel()), dtype=np.float64).reshape(ny, nx)
    k = curvature(z, c)
    for _ in range(max(int(smooth), 0)):
        p = np.pad(k, 1, mode="edge")
        k = sum(p[1 + di:1 + di + ny, 1 + dj:1 + dj + nx] for di in (-1, 0, 1) for dj in (-1, 0, 1)) / 9.0
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.sqrt(8.0 * tol / k)
    h = np.where(np.isfinite(h), h, hmax)
    h = np.clip(h, hmin, hmax)
    h = limit_gradation(h, c, growth)
    return GridDensity(float(x0), float(y0), c, c, h)


def field_sizes(x, y, ikle, f: np.ndarray, tol: float, hmin: float, hmax: float) -> np.ndarray:
    """Sizes at the nodes of a mesh from the curvature of a field on it (solution-adaptive refinement).

    The same estimate as for the terrain, h = sqrt(8 tol / |kappa|), with the Hessian of the field recovered
    on the mesh (gradients per element, averaged to the nodes, twice). ``tol`` is in the field's unit: 0.02 m
    on the free surface, 0.1 m/s on the velocity magnitude.
    """
    from ..telemac_io.calculator import MeshOperators

    ops = MeshOperators(x, y, ikle)
    f = np.asarray(f, dtype=np.float64)
    fx, fy = ops.ddx(f), ops.ddy(f)
    fxx, fxy1 = ops.ddx(fx), ops.ddy(fx)
    fxy2, fyy = ops.ddx(fy), ops.ddy(fy)
    fxy = 0.5 * (fxy1 + fxy2)
    tr = 0.5 * (fxx + fyy)
    disc = np.sqrt(np.maximum(0.25 * (fxx - fyy) ** 2 + fxy ** 2, 0.0))
    k = np.maximum(np.abs(tr + disc), np.abs(tr - disc))
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.sqrt(8.0 * tol / k)
    h = np.where(np.isfinite(h), h, hmax)
    return np.clip(h, hmin, hmax)
