"""Channel planforms for laboratory and schematic models, and channel coordinates.

A planform is a centreline (the channel axis) with its heading at every
vertex; the banks are offset half a width on each side, so cross-sections are
normal to the axis and a structured mesh built station by station is
orthogonal, also in bends:

* ``straight``  a straight flume;
* ``bend``      straight inlet, circular bend of radius R and angle (90°, 180°,
                ...), straight outlet: the classic curved-flume experiments;
* ``meanders``  a sine-generated curve (Langbein & Leopold 1966): the heading
                follows theta(s) = theta0 sin(2 pi s / lambda), the planform
                of natural meanders, with lambda the length along the axis of
                one meander wavelength and theta0 the largest deflection.

``axis_coordinates`` gives every point its channel coordinates: s, the
distance along the axis from its start, and n, the signed distance from the
axis (positive on the left looking downstream). A bed written as a formula of
s and n follows any curved channel, drawn or generated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ["Planform", "straight", "bend", "meanders", "axis_coordinates", "min_radius"]


@dataclass
class Planform:
    axis: np.ndarray     # (n, 2) centreline vertices, evenly spaced along it
    heading: np.ndarray  # (n,) radians, counter-clockwise from x
    s: np.ndarray        # (n,) distance along the axis
    width: float
    kind: str = "straight"

    def bank(self, side: int) -> np.ndarray:
        """+1 = left bank, -1 = right bank (looking downstream)."""
        h = 0.5 * self.width * side
        return self.axis + h * np.stack([-np.sin(self.heading), np.cos(self.heading)], axis=1)

    @property
    def left(self) -> np.ndarray:
        return self.bank(+1)

    @property
    def right(self) -> np.ndarray:
        return self.bank(-1)

    @property
    def length(self) -> float:
        return float(self.s[-1])

    def outline(self) -> np.ndarray:
        """Closed ring, counter-clockwise: right bank downstream, left bank back upstream."""
        ring = np.concatenate([self.right, self.left[::-1]])
        return np.concatenate([ring, ring[:1]])

    def transformed(self, origin=(0.0, 0.0), angle_deg: float = 0.0) -> "Planform":
        """Rotate by ``angle_deg`` about (0, 0), then move (0, 0) to ``origin``."""
        a = math.radians(angle_deg)
        c, s = math.cos(a), math.sin(a)
        ax = self.axis[:, 0] * c - self.axis[:, 1] * s + origin[0]
        ay = self.axis[:, 0] * s + self.axis[:, 1] * c + origin[1]
        return Planform(np.stack([ax, ay], axis=1), self.heading + a, self.s.copy(), self.width, self.kind)


def _segments(total: float, step: float) -> int:
    return max(1, int(round(total / float(step))))


def _stations(parts, step: float) -> np.ndarray:
    """Distances along the axis: each part (length) gets a whole number of near-``step`` segments, so a
    vertex falls exactly on every junction (start and end of a bend) without tiny segments."""
    out = [np.zeros(1)]
    start = 0.0
    for length in parts:
        if length <= 0:
            continue
        n = _segments(length, step)
        out.append(start + np.linspace(0.0, length, n + 1)[1:])
        start += length
    return np.concatenate(out)


def straight(length: float, width: float, step: float) -> Planform:
    n = _segments(length, step)
    s = np.linspace(0.0, length, n + 1)
    return Planform(np.stack([s, np.zeros_like(s)], axis=1), np.zeros_like(s), s, float(width), "straight")


def bend(inlet: float, radius: float, angle_deg: float, outlet: float, width: float, step: float,
         turn: str = "left") -> Planform:
    """Straight ``inlet``, a circular arc of ``radius`` (on the axis) over ``angle_deg``, straight ``outlet``.

    Starts at (0, 0) heading along +x; ``turn`` is "left" (counter-clockwise) or "right".
    """
    if radius <= 0.5 * width:
        raise ValueError("the bend radius (%g, on the axis) must exceed half the width (%g): the inner bank "
                         "would cross itself" % (radius, 0.5 * width))
    if not 0 < angle_deg <= 360:
        raise ValueError("the bend angle must be between 0 and 360 degrees")
    sign = 1.0 if turn.lower().startswith("l") else -1.0
    phi = math.radians(angle_deg)
    arc = radius * phi
    total = inlet + arc + outlet
    s = _stations([inlet, arc, outlet], step)
    x = np.empty_like(s)
    y = np.empty_like(s)
    h = np.empty_like(s)
    a = s <= inlet
    x[a], y[a], h[a] = s[a], 0.0, 0.0
    b = (s > inlet) & (s <= inlet + arc + 1e-9 * total)
    th = (s[b] - inlet) / radius
    x[b] = inlet + radius * np.sin(th)
    y[b] = sign * radius * (1.0 - np.cos(th))
    h[b] = sign * th
    c = ~(a | b)
    ex, ey = inlet + radius * math.sin(phi), sign * radius * (1.0 - math.cos(phi))
    hd = sign * phi
    d = s[c] - inlet - arc
    x[c], y[c], h[c] = ex + d * math.cos(hd), ey + d * math.sin(hd), hd
    return Planform(np.stack([x, y], axis=1), h, s, float(width), "bend")


def min_radius(wavelength: float, theta0_deg: float) -> float:
    """Smallest radius of curvature of a sine-generated curve: lambda / (2 pi theta0)."""
    t0 = math.radians(theta0_deg)
    return math.inf if t0 == 0 else wavelength / (2.0 * math.pi * t0)


def meanders(wavelength: float, theta0_deg: float, count: float, width: float, step: float,
             inlet: float = 0.0, outlet: float = 0.0) -> Planform:
    """Sine-generated meanders: ``count`` wavelengths of axis length ``wavelength`` each.

    Optional straight reaches before and after. The mean direction is +x.
    """
    rmin = min_radius(wavelength, theta0_deg)
    if rmin <= 0.5 * width:
        raise ValueError("the tightest bend has a radius of %.3g (wavelength / (2 pi theta0)), less than half the width "
                         "(%.3g): use a longer wavelength or a smaller deflection angle" % (rmin, 0.5 * width))
    t0 = math.radians(theta0_deg)
    body = wavelength * float(count)
    total = inlet + body + outlet
    # integrate the heading finely, then sample the stations
    fine = np.linspace(0.0, body, max(2000, int(body / (step / 20.0)) + 1))
    th = t0 * np.sin(2.0 * math.pi * fine / wavelength)
    ds = np.diff(fine)
    mid = 0.5 * (th[1:] + th[:-1])
    fx = np.concatenate([[0.0], np.cumsum(ds * np.cos(mid))])
    fy = np.concatenate([[0.0], np.cumsum(ds * np.sin(mid))])
    s = _stations([inlet, body, outlet], step)
    x = np.empty_like(s)
    y = np.empty_like(s)
    h = np.empty_like(s)
    a = s <= inlet
    x[a], y[a], h[a] = s[a], 0.0, 0.0
    b = (s > inlet) & (s <= inlet + body + 1e-9 * total)
    u = s[b] - inlet
    x[b] = inlet + np.interp(u, fine, fx)
    y[b] = np.interp(u, fine, fy)
    h[b] = t0 * np.sin(2.0 * math.pi * u / wavelength)
    c = ~(a | b)
    ex, ey = inlet + fx[-1], fy[-1]
    hd = t0 * math.sin(2.0 * math.pi * body / wavelength)
    d = s[c] - inlet - body
    x[c], y[c], h[c] = ex + d * math.cos(hd), ey + d * math.sin(hd), hd
    return Planform(np.stack([x, y], axis=1), h, s, float(width), "meanders")


def axis_coordinates(axis, px, py, chunk: int = 4000):
    """(s, n): distance along the polyline ``axis`` and signed offset from it (left positive).

    Points beyond the ends get s < 0 or s > length, measured along the end segments.
    """
    A = np.asarray(axis, dtype=np.float64)[:, :2]
    if len(A) < 2:
        raise ValueError("the channel axis needs at least two vertices")
    P0, P1 = A[:-1], A[1:]
    d = P1 - P0
    L2 = np.einsum("ij,ij->i", d, d)
    keep = L2 > 0
    P0, d, L2 = P0[keep], d[keep], L2[keep]
    seg = np.sqrt(L2)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    px = np.asarray(px, dtype=np.float64).ravel()
    py = np.asarray(py, dtype=np.float64).ravel()
    s_out = np.empty(px.size)
    n_out = np.empty(px.size)
    last = len(seg) - 1
    for a in range(0, px.size, chunk):
        X = px[a:a + chunk, None] - P0[None, :, 0]
        Y = py[a:a + chunk, None] - P0[None, :, 1]
        t_raw = (X * d[None, :, 0] + Y * d[None, :, 1]) / L2[None, :]
        t = np.clip(t_raw, 0.0, 1.0)
        qx = X - t * d[None, :, 0]
        qy = Y - t * d[None, :, 1]
        dist2 = qx * qx + qy * qy
        k = np.argmin(dist2, axis=1)
        r = np.arange(k.size)
        tk = t[r, k]
        # beyond the first or last vertex: extrapolate along the end segment
        tr = t_raw[r, k]
        tk = np.where((k == 0) & (tr < 0), tr, tk)
        tk = np.where((k == last) & (tr > 1), tr, tk)
        s_out[a:a + chunk] = cum[k] + tk * seg[k]
        cross = d[k, 0] * Y[r, k] - d[k, 1] * X[r, k]
        n_out[a:a + chunk] = cross / seg[k]
    return s_out, n_out
