"""Element sizes suggested from a DEM: what its pixels can show, the terrain, the rivers and a node budget.

* **Smallest size**: two pixels. A finer element only interpolates the same pixels: it adds nodes, not terrain.
* **Terrain tolerance**: about the vertical accuracy of a DEM with such pixels (LiDAR at 1 m: 0.05 m; 5 m
  DEMs: 0.2 m; 15 m DEMs: 0.5 m). Elements are refined where a larger one would miss the terrain by more
  (mesher.demsize: h = sqrt(8 tol / |curvature|)).
* **Default size** (floodplain): the size that the median terrain inside the domain asks for at that
  tolerance, at least four pixels, at most a twentieth of the domain's extent, and coarse enough that the
  default size alone stays within half the node budget.
* **Rivers** (widths of the project's channels): a spacing across of a eighth of the width (7 nodes or more
  across), a spacing along of 1.5 times that.
* **Nodes**: estimated with the terrain-driven sizes over the domain (1.155 / h² nodes per unit area for
  equilateral triangles); if over the budget the tolerance is relaxed until it fits.
* **Time step**: Courant number 1 on the smallest element for a typical depth and velocity.

Only NumPy; the DEM is any callable ``dem(x, y) -> z`` (NaN where there is no data).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

from .demsize import terrain_sizes

__all__ = ["SizeSuggestion", "suggest_sizes", "nice"]

_STEPS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)


def nice(v: float, up: bool = False) -> float:
    """A round number near ``v`` (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8 times a power of ten), below it unless ``up``."""
    if not (v > 0 and math.isfinite(v)):
        return v
    e = 10.0 ** math.floor(math.log10(v))
    m = v / e
    if up:
        s = next(s for s in _STEPS if s >= m - 1e-9)
    else:
        s = [s for s in _STEPS if s <= m + 1e-9][-1]
    return round(s * e, 10)


def tolerance_for_pixel(p: float) -> float:
    """About the vertical accuracy of a DEM with ``p`` m pixels."""
    for limit, tol in ((1.0, 0.05), (2.5, 0.1), (6.0, 0.2), (20.0, 0.5)):
        if p <= limit:
            return tol
    return 1.0


@dataclass
class SizeSuggestion:
    pixel: float
    smallest: float
    tolerance: float
    default: float
    nodes: int
    nodes_default_only: int
    time_step: float
    area: float
    channel_width: Optional[float] = None
    channel_across: Optional[float] = None
    channel_along: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    def report(self) -> List[str]:
        out = ["DEM pixels %.3g m; area %.4g km²" % (self.pixel, self.area / 1e6),
               "Smallest size: %.3g m (two pixels: finer elements add nodes, not terrain)" % self.smallest,
               "Terrain tolerance: %.3g m (about the vertical accuracy of such a DEM)" % self.tolerance,
               "Default size (floodplain): %.3g m" % self.default]
        if self.channel_across:
            out.append("Rivers about %.3g m wide: spacing across %.3g m (%d nodes across), along %.3g m"
                       % (self.channel_width, self.channel_across, int(round(self.channel_width / self.channel_across)) + 1,
                          self.channel_along))
        out.append("Nodes: about %s with the terrain refinement (%s with the default size alone)"
                   % (format(self.nodes, ","), format(self.nodes_default_only, ",")))
        out.append("Time step: about %.3g s (Courant number 1 on the smallest elements)" % self.time_step)
        return out + self.notes


def suggest_sizes(dem: Callable, bbox: Tuple[float, float, float, float], pixel: float,
                  inside: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None,
                  widths: Sequence[float] = (), node_budget: int = 300_000, depth: float = 2.0, velocity: float = 1.5,
                  growth: float = 1.2, max_cells: int = 1_000_000) -> SizeSuggestion:
    """Suggested sizes for a domain in ``bbox`` (``inside(x, y)`` tells which points belong to it; default: where
    the DEM has data) on a DEM with ``pixel`` m pixels; ``widths``: the widths of its rivers (m)."""
    if not pixel > 0:
        raise ValueError("the DEM pixel size must be positive")
    x0, y0, x1, y1 = bbox
    w, h = max(x1 - x0, 1e-9), max(y1 - y0, 1e-9)
    notes: List[str] = []
    c = max(float(pixel), math.sqrt(w * h / max_cells))
    nx, ny = int(math.ceil(w / c)) + 1, int(math.ceil(h / c)) + 1
    X, Y = np.meshgrid(x0 + np.arange(nx) * c, y0 + np.arange(ny) * c)
    z = np.asarray(dem(X.ravel(), Y.ravel()), dtype=np.float64).reshape(ny, nx)
    mask = np.isfinite(z)
    if inside is not None:
        mask &= np.asarray(inside(X.ravel(), Y.ravel()), dtype=bool).reshape(ny, nx)
    if not mask.any():
        raise ValueError("the DEM has no data inside the domain")
    area = float(mask.sum()) * c * c
    smallest = nice(2.0 * pixel, up=True)
    tol = tolerance_for_pixel(pixel)
    from .demsize import curvature

    zf = np.where(np.isfinite(z), z, np.nanmedian(z[mask]))
    k = curvature(zf, c)
    kin = k[mask & (k > 0) & np.isfinite(k)]
    typical = math.sqrt(8.0 * tol / float(np.median(kin))) if kin.size else math.sqrt(area)
    budget_h = math.sqrt(2.0 * 1.155 * area / max(node_budget, 1))
    upper = max(math.sqrt(area) / 20.0, 4.0 * pixel)
    default = min(max(typical, 4.0 * pixel, budget_h), upper)
    default = max(nice(default), smallest)
    # rivers
    width = across = along = None
    ws = [float(v) for v in widths if v and v > 0]
    if ws:
        width = float(np.median(ws))
        across = nice(width / 8.0)
        if width / across < 6.0:
            across = nice(width / 6.0)
        along = min(nice(1.5 * across), default)
        if across < pixel:
            notes.append("The rivers are narrower than 8 pixels of the DEM: it cannot show the shape of their cross-sections; "
                         "take the bed from surveyed sections or a formula (step 3).")
    else:
        notes.append("No river in the project yet: for a river W m wide use a spacing across of W / 8 (7 nodes or more "
                     "across) and 1.5 times that along.")
    # nodes with the terrain refinement, relaxing the tolerance to fit the budget

    def count(t):
        gd = terrain_sizes(lambda a, b: np.asarray(dem(a, b), dtype=np.float64), bbox, t, smallest, default, growth,
                           cell=c, max_cells=max_cells)
        hh = gd.values
        m = mask if hh.shape == mask.shape else np.isfinite(hh)
        return int(1.155 * float(np.sum(c * c / hh[m] ** 2)))

    nodes = count(tol)
    tries = 0
    while nodes > node_budget and tries < 6:
        tol = nice(tol * 1.6, up=True)
        nodes = count(tol)
        tries += 1
    if tries:
        notes.append("The terrain tolerance was relaxed to %.3g m to stay within %s nodes." % (tol, format(node_budget, ",")))
    n_default = int(1.155 * area / default ** 2)
    hmin = min(v for v in (smallest, across, default) if v)
    dt = hmin / (velocity + math.sqrt(9.81 * depth))
    return SizeSuggestion(float(pixel), float(smallest), float(tol), float(default), int(nodes), n_default, float(dt),
                          float(area), width, across, along, notes)
