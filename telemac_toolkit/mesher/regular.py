"""Structured meshes: regular grids, each cell written as two consecutive triangles.

A rectangle given by an origin, a rotation angle, cell sizes and cell counts is
split into right triangles. ``pattern`` chooses the diagonals: "alternate"
(no directional bias, the default), "right" or "left".
"""

from __future__ import annotations

import math

import numpy as np

from ..telemac_io.mesh import signed_areas

__all__ = ["regular_grid_mesh", "structured_quad_mesh"]


def regular_grid_mesh(x0: float, y0: float, dx: float, dy: float, nx: int, ny: int,
                      angle_deg: float = 0.0, pattern: str = "alternate"):
    """Return (x, y, ikle) for nx by ny cells starting at (x0, y0) (lower-left corner)."""
    if nx < 1 or ny < 1 or dx <= 0 or dy <= 0:
        raise ValueError("cell counts must be >= 1 and cell sizes > 0")
    i, j = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1))
    u = (i * dx).ravel().astype(float)
    v = (j * dy).ravel().astype(float)
    a = math.radians(angle_deg)
    x = x0 + u * math.cos(a) - v * math.sin(a)
    y = y0 + u * math.sin(a) + v * math.cos(a)
    return x, y, _split_cells(nx, ny, pattern, x, y)


def structured_quad_mesh(left, right, bottom_count: int, side_count: int):
    """Structured mesh between two opposite polylines (transfinite in one direction).

    ``left`` and ``right`` are resampled to ``side_count`` nodes; each pair of
    matching nodes is joined by ``bottom_count`` nodes. This is the channel
    mesher without a thalweg, exposed for any four-sided region.
    """
    from .channel import channel_mesh

    L = np.asarray(left, float)
    length = float(np.sum(np.hypot(*np.diff(L, axis=0).T)))
    return channel_mesh(left, right, None, cross_nodes=bottom_count, interval=length / max(side_count - 1, 1))


def _split_cells(nx, ny, pattern, x, y):
    tris = []
    w = nx + 1
    for jj in range(ny):
        for ii in range(nx):
            p00 = jj * w + ii
            p10 = p00 + 1
            p01 = p00 + w
            p11 = p01 + 1
            right = pattern == "right" or (pattern == "alternate" and (ii + jj) % 2 == 0)
            if right:
                tris += [[p00, p10, p11], [p00, p11, p01]]
            else:
                tris += [[p00, p10, p01], [p10, p11, p01]]
    ikle = np.array(tris, dtype=np.int64)
    neg = signed_areas(x, y, ikle) < 0
    ikle[neg] = ikle[neg][:, [0, 2, 1]]
    return ikle
