"""Integrals on a TELEMAC mesh: over areas, along lines and through sections (fluxes).

Section integrals are exact for the finite-element fields: the line is cut at
every mesh edge it crosses, so between two cuts the depth and the velocity
are linear and their product (the unit discharge h u) is quadratic, which
Simpson's rule integrates exactly.

Sign of a flux: positive from the right of the line to its left, looking
along the line as drawn. Drawing a section from the left bank to the right
bank (looking downstream, as for HEC-RAS cross-sections) makes downstream
flow positive.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .mesh import TriangleLocator

__all__ = ["Section", "unique_edges", "element_mask_in_rings"]


def unique_edges(ikle: np.ndarray) -> np.ndarray:
    ik = np.asarray(ikle, dtype=np.int64)
    k = ik.shape[1]
    e = np.concatenate([ik[:, [i, (i + 1) % k]] for i in range(k)])
    return np.unique(np.sort(e, axis=1), axis=0)


@dataclass
class Section:
    """A polyline cut at the mesh edges, ready to integrate nodal fields along it."""

    s: np.ndarray        # cut positions along the line (breakpoints)
    nx: np.ndarray       # left normal of each interval
    ny: np.ndarray
    nodes: np.ndarray    # (2 * nint + 1, 3) nodes and weights of breakpoints and midpoints
    w: np.ndarray
    inside: np.ndarray   # (2 * nint + 1,)
    length: float
    wet_length_nodes: Optional[np.ndarray] = None

    @classmethod
    def build(cls, line_xy, x, y, ikle, locator: Optional[TriangleLocator] = None) -> "Section":
        P = np.asarray(line_xy, dtype=np.float64)[:, :2]
        if len(P) < 2:
            raise ValueError("a section line needs at least two vertices")
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        edges = unique_edges(ikle)
        ax, ay = x[edges[:, 0]], y[edges[:, 0]]
        bx, by = x[edges[:, 1]], y[edges[:, 1]]
        seg = np.diff(P, axis=0)
        seglen = np.hypot(seg[:, 0], seg[:, 1])
        cum = np.concatenate([[0.0], np.cumsum(seglen)])
        cuts = [cum]
        for k in range(len(seg)):
            if seglen[k] == 0:
                continue
            px, py = P[k]
            dx, dy = seg[k]
            ex, ey = bx - ax, by - ay
            den = dx * ey - dy * ex
            ok = np.abs(den) > 1e-14 * (seglen[k] * np.hypot(ex, ey) + 1e-300)
            with np.errstate(divide="ignore", invalid="ignore"):
                t = ((ax - px) * ey - (ay - py) * ex) / den
                u = ((ax - px) * dy - (ay - py) * dx) / den
            hit = ok & (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1)
            cuts.append(cum[k] + t[hit] * seglen[k])
        s = np.unique(np.concatenate(cuts))
        # drop near-duplicates (a cut through a node is found on several edges)
        tol = 1e-9 * max(cum[-1], 1.0)
        s = s[np.concatenate([[True], np.diff(s) > tol])]
        mid = 0.5 * (s[1:] + s[:-1])
        allp = np.empty(2 * len(s) - 1)
        allp[0::2] = s
        allp[1::2] = mid
        qx = np.interp(allp, cum, P[:, 0])
        qy = np.interp(allp, cum, P[:, 1])
        k = np.clip(np.searchsorted(cum, mid, side="right") - 1, 0, len(seg) - 1)
        tx, ty = seg[k, 0] / np.maximum(seglen[k], 1e-300), seg[k, 1] / np.maximum(seglen[k], 1e-300)
        loc = locator or TriangleLocator(x, y, ikle)
        tri, nodes, w = loc.locate_many(qx, qy)
        return cls(s, -ty, tx, nodes, w, tri >= 0, float(cum[-1]))

    # ---- integration -----------------------------------------------------------
    def _values(self, f) -> np.ndarray:
        v = np.einsum("ij,ij->i", np.asarray(f, dtype=np.float64)[self.nodes], self.w)
        v[~self.inside] = 0.0
        return v

    def _simpson(self, v: np.ndarray) -> float:
        a, m, b = v[0:-1:2], v[1::2], v[2::2]
        return float(np.sum((self.s[1:] - self.s[:-1]) / 6.0 * (a + 4.0 * m + b)))

    def integral(self, f) -> float:
        """Line integral of a nodal field (outside the mesh counts as 0)."""
        return self._simpson(self._values(f))

    def covered_length(self, mask_nodes=None) -> float:
        """Length of the line inside the mesh (or where a nodal 0/1 mask is 1, interpolated)."""
        if mask_nodes is None:
            v = self.inside.astype(np.float64)
        else:
            v = self._values(mask_nodes)
        return self._simpson(v)

    def discharge(self, h, u, v, weight=None) -> float:
        """Integral of h (u, v) . n with h and the velocity interpolated separately, as in the model:
        exact, since their product is quadratic between two cuts."""
        hv = self._values(h)
        if weight is not None:
            hv = hv * self._values(weight)
        return self._flux_values(hv * self._values(u), hv * self._values(v))

    def flux(self, qx, qy, weight=None) -> float:
        """Integral of (qx, qy) . n along the line; ``weight`` multiplies the flux (a concentration)."""
        vx, vy = self._values(qx), self._values(qy)
        if weight is not None:
            wv = self._values(weight)
            vx, vy = vx * wv, vy * wv
        return self._flux_values(vx, vy)

    def _flux_values(self, vx, vy) -> float:
        a = slice(0, -1, 2)
        m = slice(1, None, 2)
        b = slice(2, None, 2)
        nx, ny = self.nx, self.ny
        fa = vx[a] * nx + vy[a] * ny
        fm = vx[m] * nx + vy[m] * ny
        fb = vx[b] * nx + vy[b] * ny
        return float(np.sum((self.s[1:] - self.s[:-1]) / 6.0 * (fa + 4.0 * fm + fb)))


def element_mask_in_rings(x, y, ikle, rings: Sequence[np.ndarray]) -> np.ndarray:
    """Elements whose centroid lies inside the rings (even-odd, holes allowed)."""
    ik = np.asarray(ikle, dtype=np.int64)
    cx = np.asarray(x, dtype=np.float64)[ik].mean(axis=1)
    cy = np.asarray(y, dtype=np.float64)[ik].mean(axis=1)
    inside = np.zeros(cx.size, dtype=bool)
    for ring in rings:
        r = np.asarray(ring, dtype=np.float64)[:, :2]
        if len(r) > 1 and np.allclose(r[0], r[-1]):
            r = r[:-1]
        if len(r) < 3:
            continue
        x0, y0 = r[:, 0], r[:, 1]
        x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
        cand = np.flatnonzero((cx >= x0.min()) & (cx <= x0.max()) & (cy >= y0.min()) & (cy <= y0.max()))
        if cand.size == 0:
            continue
        qx, qy = cx[cand], cy[cand]
        hit = np.zeros(cand.size, dtype=bool)
        for k in range(len(x0)):
            if y0[k] == y1[k]:
                continue
            cross = (y0[k] > qy) != (y1[k] > qy)
            xint = (x1[k] - x0[k]) * (qy - y0[k]) / (y1[k] - y0[k]) + x0[k]
            hit ^= cross & (qx < xint)
        inside[cand] ^= hit
    return inside
