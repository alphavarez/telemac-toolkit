"""Target element size field.

Each refinement source (area, point or line) fixes an element size at its
location; away from it the size grows with the growth rate g. After n elements of geometric
growth the distance covered is d = s (g^n - 1) / (g - 1), so the size reached
at distance d is exactly s + (g - 1) d. The field is the minimum over all
sources, capped by the default edge length.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np

from .geometry import as_ring, distance_to_segments, points_in_rings

__all__ = ["SizeField", "PolygonDensity", "PointDensity", "LineDensity"]


@dataclass
class PolygonDensity:
    """Constant edge length inside a polygon (rings: outer + holes)."""

    rings: List[np.ndarray]
    size: float


@dataclass
class PointDensity:
    """Edge length at points (hard points with a value, density points)."""

    xy: np.ndarray
    size: np.ndarray


@dataclass
class LineDensity:
    """Edge length along a line (soft or hard line with an edge_len value)."""

    coords: np.ndarray
    size: float
    closed: bool = False


@dataclass
class SizeField:
    default: float
    growth: float = 1.2
    polygons: List[PolygonDensity] = field(default_factory=list)
    points: List[PointDensity] = field(default_factory=list)
    lines: List[LineDensity] = field(default_factory=list)

    def __post_init__(self):
        if self.default <= 0:
            raise ValueError("default edge length must be positive")
        if self.growth < 1.0:
            raise ValueError("edge growth ratio must be >= 1")

    @property
    def minimum(self) -> float:
        vals = [self.default]
        vals += [p.size for p in self.polygons if p.size > 0]
        vals += [float(np.min(p.size[p.size > 0])) for p in self.points if np.any(p.size > 0)]
        vals += [ln.size for ln in self.lines if ln.size > 0]
        return float(min(vals))

    def _reach(self, s: float) -> float:
        """Distance beyond which a source no longer lowers the default size."""
        g = self.growth - 1.0
        if g <= 0:
            return np.inf
        return max(self.default - s, 0.0) / g

    def __call__(self, px, py) -> np.ndarray:
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        out = np.full(px.shape, self.default)
        g = self.growth - 1.0
        for poly in self.polygons:
            if poly.size <= 0:
                continue
            reach = self._reach(poly.size)
            inside = points_in_rings(px, py, poly.rings)
            out[inside] = np.minimum(out[inside], poly.size)
            if g > 0 and reach > 0:
                rest = np.flatnonzero(~inside)
                if rest.size:
                    a, b = [], []
                    for r in poly.rings:
                        rr = as_ring(r)
                        a.append(rr)
                        b.append(np.roll(rr, -1, axis=0))
                    d = distance_to_segments(px[rest], py[rest], np.vstack(a), np.vstack(b), cutoff=reach)
                    out[rest] = np.minimum(out[rest], poly.size + g * d)
        for pts in self.points:
            xy = np.asarray(pts.xy, dtype=np.float64)
            sz = np.asarray(pts.size, dtype=np.float64)
            for k in range(len(xy)):
                if sz[k] <= 0:
                    continue
                reach = self._reach(sz[k])
                sel = np.flatnonzero((np.abs(px - xy[k, 0]) <= reach) & (np.abs(py - xy[k, 1]) <= reach))
                if sel.size:
                    d = np.hypot(px[sel] - xy[k, 0], py[sel] - xy[k, 1])
                    out[sel] = np.minimum(out[sel], sz[k] + g * d)
        for ln in self.lines:
            if ln.size <= 0:
                continue
            c = np.asarray(ln.coords, dtype=np.float64)
            a, b = c[:-1], c[1:]
            if ln.closed:
                a, b = c, np.roll(c, -1, axis=0)
            d = distance_to_segments(px, py, a, b, cutoff=self._reach(ln.size))
            out = np.minimum(out, ln.size + g * d)
        return out
