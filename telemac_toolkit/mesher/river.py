"""River reaches for structured channel strips: from a centreline or from two bank lines.

A reach is a list of cross-sections along the river: each one is a station on
the axis, a unit normal (pointing to the left bank, looking downstream) and the
distances from the axis to the left and right banks. The banks and the axis
then have matched vertices, so ``channel_mesh`` uses them as the strip's
cross-sections, normal to the flow, also in bends.

Sources of the banks:

* a centreline and a constant width;
* a centreline and a DEM: at every station the terrain profile across the
  river is read, the bed level found near the centreline, and each bank put
  where the terrain first rises ``bankfull_depth`` above the bed;
* two bank lines, drawn or taken from another layer.

With a DEM the reach is also turned so it runs downstream (from the higher
bed to the lower one).

Tight bends: cross-sections normal to a jagged or tightly curved centreline
can cross each other on the inner bank. The normals are smoothed along the
river until neighbouring sections no longer cross; sections that still cross
are reported.

Only NumPy is used; the DEM is any callable ``dem(x, y) -> z`` (NaN where
there is no data).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np

__all__ = ["Reach", "resample", "smooth_line", "reach_from_centreline", "reach_from_banks", "end_boundary_lines"]

Dem = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass
class Reach:
    axis: np.ndarray    # (n, 2) stations on the axis (midway between the banks)
    normal: np.ndarray  # (n, 2) unit normals, towards the left bank
    left_w: np.ndarray  # (n,) distance from the axis to the left bank
    right_w: np.ndarray  # (n,) distance from the axis to the right bank
    log: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    reversed: bool = False

    @property
    def left(self) -> np.ndarray:
        return self.axis + self.left_w[:, None] * self.normal

    @property
    def right(self) -> np.ndarray:
        return self.axis - self.right_w[:, None] * self.normal

    @property
    def width(self) -> np.ndarray:
        return self.left_w + self.right_w

    @property
    def s(self) -> np.ndarray:
        return _arc(self.axis)

    @property
    def length(self) -> float:
        return float(self.s[-1])

    def outline(self) -> np.ndarray:
        """Closed ring of the channel: right bank downstream, left bank back upstream."""
        ring = np.concatenate([self.right, self.left[::-1]])
        return np.concatenate([ring, ring[:1]])

    def crossings(self) -> np.ndarray:
        """Indices i where cross-section i crosses cross-section i + 1."""
        return _crossing_sections(self.left, self.right)

    def trimmed(self, keep: np.ndarray) -> "Reach":
        k = np.asarray(keep)
        return Reach(self.axis[k], self.normal[k], self.left_w[k], self.right_w[k], list(self.log), list(self.warnings),
                     self.reversed)


# ---------------------------------------------------------------------------
# polylines
# ---------------------------------------------------------------------------
def _arc(c: np.ndarray) -> np.ndarray:
    c = np.asarray(c, dtype=np.float64)
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])


def _clean(line) -> np.ndarray:
    c = np.asarray(line, dtype=np.float64)[:, :2]
    if len(c) < 2:
        raise ValueError("a river line needs at least two vertices")
    keep = np.concatenate([[True], np.hypot(*np.diff(c, axis=0).T) > 0])
    c = c[keep]
    if len(c) < 2:
        raise ValueError("a river line needs two distinct vertices")
    return c


def resample(line, step: float) -> np.ndarray:
    """Vertices evenly spaced along a polyline (the ends are kept)."""
    c = _clean(line)
    s = _arc(c)
    n = max(1, int(round(s[-1] / float(step))))
    t = np.linspace(0.0, s[-1], n + 1)
    return np.stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])], axis=1)


def smooth_line(c: np.ndarray, sigma: float, step: float) -> np.ndarray:
    """Gaussian smoothing along an evenly spaced polyline; the two end points do not move.

    The line is extended past each end by point reflection about the end point, so the
    smoothed line keeps its ends and its end directions.
    """
    c = np.asarray(c, dtype=np.float64)
    k = sigma / max(step, 1e-300)
    if k < 0.5 or len(c) < 3:
        return c.copy()
    r = int(min(math.ceil(3 * k), len(c) - 1))
    if r < 1:
        return c.copy()
    pad0 = 2 * c[0] - c[r:0:-1]
    pad1 = 2 * c[-1] - c[-2:-r - 2:-1]
    ext = np.concatenate([pad0, c, pad1])
    w = np.exp(-0.5 * (np.arange(-r, r + 1) / k) ** 2)
    w /= w.sum()
    out = np.stack([np.convolve(ext[:, 0], w, mode="valid"), np.convolve(ext[:, 1], w, mode="valid")], axis=1)
    out[0], out[-1] = c[0], c[-1]
    return out


def _normals(c: np.ndarray) -> np.ndarray:
    """Unit normals to the left of the direction of travel; the end normals are those of the end segments."""
    t = np.gradient(c, axis=0)
    t[0] = c[1] - c[0]
    t[-1] = c[-1] - c[-2]
    t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-300)[:, None]
    return np.stack([-t[:, 1], t[:, 0]], axis=1)


def _smooth_normals(nrm: np.ndarray, sigma_pts: float) -> np.ndarray:
    """Smooth the direction of the normals along the river (ends kept)."""
    if sigma_pts < 0.5 or len(nrm) < 3:
        return nrm
    r = int(min(math.ceil(3 * sigma_pts), len(nrm) - 1))
    w = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma_pts) ** 2)
    ext = np.concatenate([np.repeat(nrm[:1], r, axis=0), nrm, np.repeat(nrm[-1:], r, axis=0)])
    out = np.stack([np.convolve(ext[:, 0], w, mode="valid"), np.convolve(ext[:, 1], w, mode="valid")], axis=1)
    out /= np.maximum(np.hypot(out[:, 0], out[:, 1]), 1e-300)[:, None]
    out[0], out[-1] = nrm[0], nrm[-1]
    return out


def _orient(a, b, c):
    return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])


def _crossing_sections(L: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Indices i where segment R[i]-L[i] crosses R[i+1]-L[i+1], or where a bank runs backwards."""
    if len(L) < 2:
        return np.zeros(0, dtype=np.int64)
    p1, q1, p2, q2 = R[:-1], L[:-1], R[1:], L[1:]
    d1 = _orient(p1, q1, p2)
    d2 = _orient(p1, q1, q2)
    d3 = _orient(p2, q2, p1)
    d4 = _orient(p2, q2, q1)
    cross = (d1 * d2 < 0) & (d3 * d4 < 0)
    return np.flatnonzero(cross)


# ---------------------------------------------------------------------------
# DEM
# ---------------------------------------------------------------------------
def _profile_banks(axis, nrm, dem: Dem, depth: float, max_half: float, min_half: float, dt: float):
    """Left and right bank distances where the terrain first rises ``depth`` above the bed.

    Returns (left, right, bed, found_left, found_right); distances measured from the axis,
    the bed is the lowest terrain within ``min_half`` (at least a quarter of ``max_half``)
    of the axis.
    """
    t = np.arange(-max_half, max_half + 0.5 * dt, dt)
    X = axis[:, 0][:, None] + t[None, :] * nrm[:, 0][:, None]
    Y = axis[:, 1][:, None] + t[None, :] * nrm[:, 1][:, None]
    Z = np.asarray(dem(X.ravel(), Y.ravel()), dtype=np.float64).reshape(X.shape)
    core = np.abs(t) <= max(min_half, 0.25 * max_half)
    Zc = np.where(core[None, :], Z, np.nan)
    with np.errstate(all="ignore"):
        bad = np.all(np.isnan(Zc), axis=1)
        i0 = np.nanargmin(np.where(bad[:, None], 0.0, Zc), axis=1)
    n = len(axis)
    bed = np.where(bad, np.nan, Z[np.arange(n), i0])
    level = bed + depth
    left = np.full(n, max_half)
    right = np.full(n, max_half)
    fl = np.zeros(n, dtype=bool)
    fr = np.zeros(n, dtype=bool)
    for i in range(n):
        if bad[i]:
            continue
        z = Z[i]
        # left: from the bed towards +t
        above = np.flatnonzero((np.arange(len(t)) > i0[i]) & (z >= level[i]))
        if above.size:
            j = above[0]
            z0, z1 = z[j - 1], z[j]
            f = 0.0 if not np.isfinite(z0) or z1 == z0 else (level[i] - z0) / (z1 - z0)
            left[i] = t[j - 1] + float(np.clip(f, 0, 1)) * dt
            fl[i] = True
        below = np.flatnonzero((np.arange(len(t)) < i0[i]) & (z >= level[i]))
        if below.size:
            j = below[-1]
            z0, z1 = z[j + 1], z[j]
            f = 0.0 if not np.isfinite(z0) or z1 == z0 else (level[i] - z0) / (z1 - z0)
            right[i] = -(t[j + 1] - float(np.clip(f, 0, 1)) * dt)
            fr[i] = True
    return left, right, bed, fl, fr


def _running_median(v: np.ndarray, half: int) -> np.ndarray:
    if half < 1 or len(v) < 3:
        return v.copy()
    pad = np.concatenate([np.repeat(v[:1], half), v, np.repeat(v[-1:], half)])
    win = np.lib.stride_tricks.sliding_window_view(pad, 2 * half + 1)
    return np.median(win, axis=1)


def _gauss1d(v: np.ndarray, sigma_pts: float) -> np.ndarray:
    if sigma_pts < 0.5 or len(v) < 3:
        return v.copy()
    r = int(min(math.ceil(3 * sigma_pts), len(v) - 1))
    w = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma_pts) ** 2)
    w /= w.sum()
    pad = np.concatenate([np.repeat(v[:1], r), v, np.repeat(v[-1:], r)])
    return np.convolve(pad, w, mode="valid")


def _downstream(axis: np.ndarray, dem: Dem) -> Tuple[bool, float, float]:
    """(needs_reversing, bed at start, bed at end): mean terrain over the first and last tenth of the axis."""
    n = len(axis)
    k = max(1, n // 10)
    z = np.asarray(dem(axis[:, 0], axis[:, 1]), dtype=np.float64)
    with np.errstate(all="ignore"):
        z0 = float(np.nanmedian(z[:k]))
        z1 = float(np.nanmedian(z[-k:]))
    if not (np.isfinite(z0) and np.isfinite(z1)):
        return False, z0, z1
    return z0 < z1, z0, z1


# ---------------------------------------------------------------------------
# building a reach
# ---------------------------------------------------------------------------
def _square_ends(reach: Reach) -> Reach:
    """End cross-sections normal to the end segments of the axis, as the flat ends of a buffer
    of the axis (the domain the helper draws), so the end sections lie on the domain outline."""
    for i, j in ((0, 1), (-1, -2)):
        t = reach.axis[i] - reach.axis[j] if i == -1 else reach.axis[j] - reach.axis[i]
        t = t / max(float(np.hypot(*t)), 1e-300)
        reach.normal[i] = (-t[1], t[0])
    return reach


def _bad_sections(reach: Reach, min_gap: float) -> np.ndarray:
    """Intervals i (between sections i and i + 1) that cross, or whose cells on a bank are shorter than ``min_gap``."""
    bad = set(reach.crossings().tolist())
    if min_gap > 0 and len(reach.axis) > 1:
        for B in (reach.left, reach.right):
            bad.update(np.flatnonzero(np.hypot(*np.diff(B, axis=0).T) < min_gap).tolist())
    return np.array(sorted(bad), dtype=np.int64)


def _fix_folds(reach: Reach, sigma0: float, step: float, rounds: int = 8, max_ratio: float = 4.0) -> Reach:
    """Turn the cross-sections in tight bends until they no longer cross and no bank cell is shorter than
    1 / ``max_ratio`` of the section spacing. Only the sections around the bad ones change (blended into the
    smoothed normals over two smoothing lengths), so the others stay normal to the axis."""
    min_gap = step / max_ratio if max_ratio and max_ratio > 0 else 0.0
    bad = _bad_sections(reach, min_gap)
    if bad.size == 0:
        return reach
    nrm0 = reach.normal.copy()
    sig = max(sigma0 / step, 1.0)
    idx = np.arange(len(nrm0), dtype=np.float64)
    near = bad.astype(np.float64)
    for _ in range(rounds):
        sig *= 1.6
        sm = _smooth_normals(nrm0, sig)
        d = np.full(len(idx), np.inf)
        for b in (near, near + 1):  # distance in sections to the nearest end of a bad interval
            pos = np.clip(np.searchsorted(b, idx), 1, len(b) - 1) if len(b) > 1 else np.zeros(len(idx), dtype=np.int64)
            d = np.minimum(d, np.minimum(np.abs(idx - b[pos]), np.abs(idx - b[np.maximum(pos - 1, 0)])))
        w = np.clip((4.0 * sig - d) / (2.0 * sig), 0.0, 1.0)
        nrm = (1.0 - w)[:, None] * nrm0 + w[:, None] * sm
        reach.normal = nrm / np.maximum(np.hypot(nrm[:, 0], nrm[:, 1]), 1e-300)[:, None]
        reach.normal[0], reach.normal[-1] = nrm0[0], nrm0[-1]
        left = _bad_sections(reach, min_gap)
        if left.size == 0:
            reach.log.append("Normals smoothed over %.3g m around %d tight bend section(s), so cross-sections do not cross "
                             "and inner-bank cells stay at least 1/%g of the spacing; the others stay normal to the axis"
                             % (sig * step, bad.size, max_ratio))
            return reach
        near = np.union1d(near, left.astype(np.float64))
    s = reach.s
    crossing = reach.crossings()
    if crossing.size:
        where = ", ".join("%.0f" % s[i] for i in crossing[:8]) + (" …" if crossing.size > 8 else "")
        reach.warnings.append("%d cross-sections still cross their neighbour (at s = %s m): the bend is tighter than "
                              "half the width; edit the banks there or use a wider spacing" % (crossing.size, where))
    else:
        short = _bad_sections(reach, min_gap)
        where = ", ".join("%.0f" % s[i] for i in short[:8]) + (" …" if short.size > 8 else "")
        reach.warnings.append("%d cells on an inner bank are shorter than 1/%g of the spacing (at s = %s m): the bend is "
                              "tight for this width" % (short.size, max_ratio, where))
    return reach


def reach_from_centreline(line, step: float, width: Optional[float] = None, dem: Optional[Dem] = None,
                          bankfull_depth: Optional[float] = None, max_width: Optional[float] = None,
                          min_width: Optional[float] = None, smoothing: Optional[float] = None,
                          orient: bool = True, max_ratio: float = 4.0) -> Reach:
    """A reach from a centreline: constant ``width``, or banks found in the ``dem`` with ``bankfull_depth``.

    ``step``: spacing of cross-sections along the river. ``smoothing``: length (m) over which the
    centreline is smoothed before the normals are taken (default: the width, or half ``max_width``).
    ``max_ratio``: in bends, no bank cell shorter than ``step / max_ratio`` (the sections are turned there).
    """
    if step <= 0:
        raise ValueError("the cross-section spacing must be positive")
    use_dem = dem is not None and bankfull_depth is not None and bankfull_depth > 0
    if not use_dem and not (width and width > 0):
        raise ValueError("give a width, or a DEM with a bankfull depth")
    if use_dem and not (max_width and max_width > 0):
        raise ValueError("with banks from the DEM, give the largest width to search")
    log: List[str] = []
    c = resample(line, step)
    rev = False
    if orient and dem is not None:
        rev, z0, z1 = _downstream(c, dem)
        if rev:
            c = c[::-1].copy()
            log.append("Centreline reversed to run downstream (bed %.3f m at the drawn start, %.3f m at the end)" % (z0, z1))
        elif np.isfinite(z0) and np.isfinite(z1):
            log.append("Centreline runs downstream (bed %.3f m to %.3f m)" % (z0, z1))
    ref = float(width) if width else 0.5 * float(max_width)
    sm = float(smoothing) if smoothing is not None else ref
    ax = smooth_line(c, 0.5 * sm, step)
    nrm = _normals(ax)
    warnings: List[str] = []
    if use_dem:
        half_max = 0.5 * float(max_width)
        half_min = 0.5 * float(min_width or 0.0)
        dt = max(half_max / 150.0, 1e-6)
        lw, rw, bed, fl, fr = _profile_banks(ax, nrm, dem, float(bankfull_depth), half_max, half_min, dt)
        miss = ~(fl & fr)
        if miss.all():
            raise ValueError("no bank found in the DEM along the reach: the terrain never rises %g m above the bed "
                             "within %g m of the centreline (check the DEM, the depth and the largest width)"
                             % (bankfull_depth, half_max))
        if miss.any():
            warnings.append("%d of %d cross-sections without a bank in the DEM on one side; there the bank is "
                            "interpolated from its neighbours" % (int(miss.sum()), len(miss)))
        s = _arc(ax)
        for v, ok in ((lw, fl), (rw, fr)):
            if (~ok).any() and ok.any():
                v[~ok] = np.interp(s[~ok], s[ok], v[ok])
        lw = _gauss1d(_running_median(lw, 2), max(sm / step / 2, 0.0))
        rw = _gauss1d(_running_median(rw, 2), max(sm / step / 2, 0.0))
        tot = lw + rw
        if half_min > 0:
            grow = np.maximum(2 * half_min - tot, 0.0) / 2
            lw, rw = lw + grow, rw + grow
        # the strip axis is midway between the banks
        shift = 0.5 * (lw - rw)
        ax = ax + shift[:, None] * nrm
        half = 0.5 * (lw + rw)
        lw, rw = half.copy(), half.copy()
        nrm = _normals(smooth_line(ax, 0.5 * sm, step)) if len(ax) > 2 else nrm
        log.append("Banks from the DEM (%g m above the bed): width %.3g to %.3g m, median %.3g m"
                   % (bankfull_depth, float(np.min(lw + rw)), float(np.max(lw + rw)), float(np.median(lw + rw))))
    else:
        lw = np.full(len(ax), 0.5 * float(width))
        rw = lw.copy()
        if dem is None:
            log.append("No DEM: the reach runs in the drawn direction (draw the centreline downstream)")
    reach = Reach(ax, nrm, lw, rw, log, warnings, rev)
    reach = _fix_folds(_square_ends(reach), 0.5 * sm, step, max_ratio=max_ratio)
    reach.log.append("Reach: %.4g m along the axis, %d cross-sections every %.3g m" % (reach.length, len(ax), step))
    return reach


def _ray_hits(origin: np.ndarray, direction: np.ndarray, poly: np.ndarray, side: int, chunk: int = 512) -> np.ndarray:
    """Distance along each ray, on one ``side`` (+1 forwards, -1 backwards), to the nearest crossing
    of polyline ``poly``: positive for +1, negative for -1; NaN where the ray misses it."""
    A, B = poly[:-1], poly[1:]
    E = B - A
    out = np.full(len(origin), np.nan)
    for a in range(0, len(origin), chunk):
        Og = origin[a:a + chunk][:, None, :]
        D = direction[a:a + chunk][:, None, :]
        den = D[..., 0] * E[None, :, 1] - D[..., 1] * E[None, :, 0]
        W = A[None, :, :] - Og
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (W[..., 0] * E[None, :, 1] - W[..., 1] * E[None, :, 0]) / den
            u = (W[..., 0] * D[..., 1] - W[..., 1] * D[..., 0]) / den
        ok = (np.abs(den) > 0) & (u >= 0) & (u <= 1) & (side * t > 0)
        d = np.where(ok, np.abs(t), np.inf)
        k = np.argmin(d, axis=1)
        r = np.arange(d.shape[0])
        best = d[r, k]
        out[a:a + chunk] = np.where(np.isfinite(best), side * best, np.nan)
    return out


def reach_from_banks(left_line, right_line, step: float, dem: Optional[Dem] = None,
                     smoothing: Optional[float] = None, orient: bool = True, max_ratio: float = 4.0) -> Reach:
    """A reach from two bank lines (any direction, either labelled left or right).

    The centreline is the smoothed middle of the banks; each cross-section is the normal to it,
    cut by both banks, so the banks keep their drawn shape and the sections stay normal to the flow.
    """
    if step <= 0:
        raise ValueError("the cross-section spacing must be positive")
    L = _clean(left_line)
    R = _clean(right_line)
    if np.hypot(*(L[0] - R[0])) + np.hypot(*(L[-1] - R[-1])) > np.hypot(*(L[0] - R[-1])) + np.hypot(*(L[-1] - R[0])):
        R = R[::-1].copy()
    log: List[str] = []
    n = max(3, int(round(0.5 * (_arc(L)[-1] + _arc(R)[-1]) / step)) + 1)
    f = np.linspace(0.0, 1.0, n)

    def at(c):
        s = _arc(c)
        t = f * s[-1]
        return np.stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])], axis=1)

    mid = 0.5 * (at(L) + at(R))
    width0 = float(np.median(np.hypot(*(at(L) - at(R)).T)))
    rev = False
    if orient and dem is not None:
        rev, z0, z1 = _downstream(mid, dem)
        if rev:
            mid, L, R = mid[::-1].copy(), R[::-1].copy(), L[::-1].copy()
            log.append("Banks reversed to run downstream (bed %.3f m at the drawn start, %.3f m at the end)" % (z0, z1))
    sm = float(smoothing) if smoothing is not None else width0
    ax = smooth_line(resample(mid, step), 0.5 * sm, step)
    nrm = _normals(ax)
    # the bank called left must be on the left of the flow: nearest crossing side, by majority
    pos = _ray_hits(ax, nrm, L, +1)
    neg = -_ray_hits(ax, nrm, L, -1)
    pos_n = np.where(np.isfinite(pos), pos, np.inf)
    neg_n = np.where(np.isfinite(neg), neg, np.inf)
    if np.sum(neg_n < pos_n) > np.sum(pos_n < neg_n):
        L, R = R, L
        log.append("Bank lines swapped: the line given as left is on the right looking downstream")
    lw = _ray_hits(ax, nrm, L, +1)
    rw = -_ray_hits(ax, nrm, R, -1)
    ok = np.isfinite(lw) & np.isfinite(rw) & (lw > 0) & (rw > 0)
    if not ok.any():
        raise ValueError("the bank lines could not be matched: check that they face each other along the river")
    first, last = int(np.argmax(ok)), int(len(ok) - 1 - np.argmax(ok[::-1]))
    warnings: List[str] = []
    if first > 0 or last < len(ok) - 1:
        log.append("Cross-sections beyond the end of a bank dropped: %d at the start, %d at the end"
                   % (first, len(ok) - 1 - last))
    keep = np.arange(first, last + 1)
    inner = ~ok[keep]
    if inner.any():
        warnings.append("%d cross-sections did not meet both banks; their bank points are interpolated" % int(inner.sum()))
        s = _arc(ax[keep])
        for v in (lw, rw):
            vk = v[keep]
            good = ok[keep]
            vk[~good] = np.interp(s[~good], s[good], vk[good])
            v[keep] = vk
    ax, nrm, lw, rw = ax[keep], nrm[keep], lw[keep], rw[keep]
    # put the axis midway between the banks
    Lp, Rp = ax + lw[:, None] * nrm, ax - rw[:, None] * nrm
    ax = 0.5 * (Lp + Rp)
    half = 0.5 * (lw + rw)
    reach = Reach(ax, nrm, half.copy(), half.copy(), log, warnings, rev)
    reach = _fix_folds(_square_ends(reach), 0.5 * sm, step, max_ratio=max_ratio)
    reach.log.append("Reach from bank lines: %.4g m along the axis, %d cross-sections every %.3g m, width %.3g to %.3g m"
                     % (reach.length, len(ax), step, float(reach.width.min()), float(reach.width.max())))
    return reach


# ---------------------------------------------------------------------------
# boundary lines where the reach meets the domain outline
# ---------------------------------------------------------------------------
def _ring_param(ring: np.ndarray, p: np.ndarray) -> Tuple[float, float]:
    """(arc length along the closed ring of the point nearest to p, its distance to p)."""
    A, B = ring[:-1], ring[1:]
    E = B - A
    L2 = np.maximum(np.einsum("ij,ij->i", E, E), 1e-300)
    t = np.clip(np.einsum("ij,ij->i", p[None, :] - A, E) / L2, 0.0, 1.0)
    Q = A + t[:, None] * E
    d = np.hypot(*(Q - p[None, :]).T)
    k = int(np.argmin(d))
    s = _arc(ring)
    return float(s[k] + t[k] * math.sqrt(L2[k])), float(d[k])


def _ring_piece(ring: np.ndarray, s0: float, s1: float) -> np.ndarray:
    """Part of the closed ring from arc length s0 forward to s1 (wrapping)."""
    s = _arc(ring)
    P = s[-1]
    s0, s1 = s0 % P, s1 % P
    if s1 < s0:
        s1 += P
    ss = np.concatenate([s, s[1:] + P])
    rr = np.concatenate([ring, ring[1:]])
    inside = (ss > s0) & (ss < s1)

    def pt(v):
        return np.array([np.interp(v, ss, rr[:, 0]), np.interp(v, ss, rr[:, 1])])

    return np.vstack([pt(s0), rr[inside], pt(s1)])


def _corner_limits(ring: np.ndarray, s_a: float, s_b: float, corner_deg: float) -> Tuple[float, float]:
    """Extend the span [s_a, s_b] along the ring to the nearest corner (turning more than corner_deg) on each side."""
    r = ring[:-1]
    v_in = r - np.roll(r, 1, axis=0)
    v_out = np.roll(r, -1, axis=0) - r
    ang = np.degrees(np.abs(np.arctan2(v_in[:, 0] * v_out[:, 1] - v_in[:, 1] * v_out[:, 0],
                                       np.einsum("ij,ij->i", v_in, v_out))))
    s = _arc(ring)[:-1]
    P = _arc(ring)[-1]
    corners = s[ang > corner_deg]
    if corners.size == 0:
        return s_a, s_b
    # backwards from s_a, forwards from s_b (on the ring, wrapping)
    back = (s_a - corners) % P
    fwd = (corners - s_b) % P
    lo = s_a - float(np.min(back[back > 1e-9 * P])) if np.any(back > 1e-9 * P) else s_a
    hi = s_b + float(np.min(fwd[fwd > 1e-9 * P])) if np.any(fwd > 1e-9 * P) else s_b
    return lo, hi


def end_boundary_lines(reach: Reach, ring: np.ndarray, inflow_span: str = "channel", outflow_span: str = "valley",
                       margin: float = 0.0, corner_deg: float = 45.0):
    """The inflow (upstream end) and outflow (downstream end) lines on the domain outline ``ring``.

    ``channel``: the outline between the two banks of the end cross-section. ``valley``: the
    whole side of the domain at that end, up to the nearest corners. Each line stops ``margin``
    short of its ends, so the corner and bank nodes there stay walls.
    Returns [(name, type, coords, distance of the end section from the outline), ...].
    """
    out = []
    for k, name, typ, span in ((0, "inflow", "Q", inflow_span), (-1, "outflow", "H", outflow_span)):
        line, dist = end_line(reach, k, ring, span, margin, corner_deg)
        out.append((name, typ, line, dist))
    return out


def end_line(reach: Reach, k: int, ring: np.ndarray, span: str = "channel", margin: float = 0.0,
             corner_deg: float = 45.0) -> Tuple[np.ndarray, float]:
    """The boundary line on ``ring`` at end ``k`` (0 or -1) of the reach, and the end section's distance from it.

    ``span``: "channel" (between the banks) or "valley" (the whole side, to the nearest corners).
    """
    ring = np.asarray(ring, dtype=np.float64)[:, :2]
    if np.hypot(*(ring[0] - ring[-1])) > 0:
        ring = np.vstack([ring, ring[:1]])
    P = float(_arc(ring)[-1])
    sl, dl = _ring_param(ring, reach.left[k])
    sr, dr = _ring_param(ring, reach.right[k])
    # the short way round between the two bank points crosses the channel
    fwd = (sl - sr) % P
    if fwd <= P - fwd:
        a, b = sr, sr + fwd
    else:
        a, b = sl, sl + (P - fwd)
    if span == "valley":
        a0 = a % P
        a, b = _corner_limits(ring, a0, a0 + (b - a), corner_deg)
    m = min(margin, 0.25 * (b - a))
    return _ring_piece(ring, a + m, b - m), max(dl, dr)
