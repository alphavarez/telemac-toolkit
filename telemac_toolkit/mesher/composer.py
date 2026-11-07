"""Laboratory and synthetic channels composed from segments, drawn centrelines or imported files.

A *composition* is a channel axis built segment by segment, with a width, a cross-section and a bed slope that
may change along it. Text, one segment per line (or separated by ";"), "#" starts a comment:

    width 1.0                  start width (m); later in the text: a sudden width change (expansion or contraction)
    straight 5                 straight reach, length along the axis (m)
    bend left 2 90             circular bend: left or right, radius on the axis (m), angle (degrees)
    elbow right 90             sharp corner (no radius), as in an L-shaped flume
    meanders 10 60 2           sine-generated meanders: wavelength along the axis, largest deflection (deg), count
    kinoshita 10 110 2 0.0156 0.0313   Kinoshita meanders: wavelength, theta0 (deg), count, skewness Js, flatness Jf
    taper 2.0 3                width changing linearly to 2.0 m over a straight 3 m
    expand 2.0                 sudden width change to 2.0 m (same as a later "width")
    slope 0.001                bed slope from here on (m/m, positive downstream)
    drop 0.05                  sudden drop of the bed (a step down; negative = a step up)
    bed 0.30                   bed level on the axis from here on (m)
    section rect               cross-section from here on, see below

Cross-sections (n: distance from the axis, left positive; b: half the local width; OFFSET moves the main channel
to the left (+) or right (-) of the axis):

    section rect
    section trapezoid BOTTOM_WIDTH SIDE_SLOPE BANK_HEIGHT [OFFSET]
    section compound MAIN_WIDTH SIDE_SLOPE FLOODPLAIN_HEIGHT FLOODPLAIN_CROSS_SLOPE [OFFSET]
    section v DEPTH [OFFSET]
    section parabolic DEPTH [OFFSET]
    section asym BOTTOM_WIDTH LEFT_SLOPE RIGHT_SLOPE BANK_HEIGHT [OFFSET]
    section custom FORMULA     any formula of n and b (and s), e.g. 0.1*(n/b)**4

The bed is a formula of the channel coordinates s, n and b, so Generate mesh evaluates it on any mesh. Pieces of
the channel between sudden width changes and sharp corners are meshed as separate structured strips; the
corners and steps between them get triangles.

A drawn or imported centreline becomes a composition too: it is cut into pieces at sharp corners (turning more
than ``corner_deg`` at a vertex), and its width may vary vertex by vertex.

Only NumPy is used.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .planform import Planform

__all__ = ["Section", "Composition", "compose", "from_centreline", "read_xy_file", "parse_section", "SEGMENT_HELP"]

SEGMENT_HELP = __doc__.split("Cross-sections")[0].split("may change along it.")[1].strip()


# ---------------------------------------------------------------------------
# cross-sections
# ---------------------------------------------------------------------------
@dataclass
class Section:
    kind: str = "rect"
    params: Tuple = ()
    custom: str = ""

    def formula(self) -> str:
        """Bed above the axis bed level, as a formula of n (across, left positive) and b (half width)."""
        p = list(self.params)
        k = self.kind
        if k == "rect":
            return "0"
        off = 0.0
        if k == "trapezoid":
            bw, m, h = p[:3]
            off = p[3] if len(p) > 3 else 0.0
            n = _n(off)
            return "clip((abs(%s) - %g)/%g, 0, %g)" % (n, 0.5 * bw, max(m, 1e-12), h) if m > 0 else \
                "where(abs(%s) > %g, %g, 0)" % (n, 0.5 * bw, h)
        if k == "compound":
            bw, m, h, fs = p[:4]
            off = p[4] if len(p) > 4 else 0.0
            n = _n(off)
            edge = 0.5 * bw + m * h
            main = "clip((abs(%s) - %g)/%g, 0, %g)" % (n, 0.5 * bw, max(m, 1e-12), h) if m > 0 else \
                "where(abs(%s) > %g, %g, 0)" % (n, 0.5 * bw, h)
            return main + (" + %g*maximum(abs(%s) - %g, 0)" % (fs, n, edge) if fs else "")
        if k == "v":
            h = p[0]
            off = p[1] if len(p) > 1 else 0.0
            return "%g*abs(%s)/b" % (h, _n(off))
        if k == "parabolic":
            h = p[0]
            off = p[1] if len(p) > 1 else 0.0
            return "%g*(%s/b)**2" % (h, _n(off))
        if k == "asym":
            bw, ml, mr, h = p[:4]
            off = p[4] if len(p) > 4 else 0.0
            n = _n(off)
            return ("where(%s > 0, clip((%s - %g)/%g, 0, %g), clip((-%s - %g)/%g, 0, %g))"
                    % (n, n, 0.5 * bw, max(ml, 1e-12), h, n, 0.5 * bw, max(mr, 1e-12), h))
        if k == "custom":
            return "(%s)" % self.custom
        raise ValueError("unknown cross-section %r" % k)

    def describe(self) -> str:
        names = {"rect": "rectangular", "trapezoid": "trapezoidal", "compound": "compound", "v": "V-shaped",
                 "parabolic": "parabolic", "asym": "asymmetric trapezoid", "custom": "custom"}
        return names.get(self.kind, self.kind) + ("" if not self.params else " " + " ".join("%g" % v for v in self.params)) + \
            (" " + self.custom if self.custom else "")


def _n(offset: float) -> str:
    return "n" if not offset else "(n - %g)" % offset


_SECTION_ARGS = {"rect": (0, 0), "rectangular": (0, 0), "trapezoid": (3, 4), "trapezoidal": (3, 4), "compound": (4, 5),
                 "v": (1, 2), "parabolic": (1, 2), "asym": (4, 5), "asymmetric": (4, 5)}


def parse_section(words: Sequence[str]) -> Section:
    """'trapezoid 0.5 1 0.2' -> Section; 'custom <formula>' keeps the rest of the text."""
    if not words:
        raise ValueError("section: give its kind (rect, trapezoid, compound, v, parabolic, asym, custom)")
    k = words[0].lower()
    if k == "custom":
        f = " ".join(words[1:]).strip()
        if not f:
            raise ValueError("section custom: give a formula of n and b, e.g. 0.1*(n/b)**4")
        from .interp import _check_formula

        _check_formula(f)
        return Section("custom", (), f)
    if k not in _SECTION_ARGS:
        raise ValueError("unknown cross-section %r (use rect, trapezoid, compound, v, parabolic, asym or custom)" % words[0])
    lo, hi = _SECTION_ARGS[k]
    vals = [_num(w, "section " + k) for w in words[1:]]
    if not lo <= len(vals) <= hi:
        raise ValueError("section %s takes %d%s numbers (%d given)" % (k, lo, "" if lo == hi else " or %d" % hi, len(vals)))
    k = {"rectangular": "rect", "trapezoidal": "trapezoid", "asymmetric": "asym"}.get(k, k)
    return Section(k, tuple(vals))


def _num(w: str, what: str) -> float:
    try:
        return float(w)
    except ValueError:
        raise ValueError("%s: %r is not a number" % (what, w))


# ---------------------------------------------------------------------------
# the composition
# ---------------------------------------------------------------------------
@dataclass
class Composition:
    pieces: List[Planform]  # stretches of continuous width and heading, meeting at steps or sharp corners
    joints: List[dict]  # between consecutive pieces: {"kind": "step" or "elbow", "angle": degrees (elbows)}
    axis: np.ndarray  # the whole axis, from the inlet to the outlet
    s: np.ndarray  # distance along it at each vertex
    half_width: Tuple[np.ndarray, np.ndarray]  # (s, b): half the width along the axis, for the formula's b
    formula: str  # bed, a formula of s, n and b
    log: List[str] = field(default_factory=list)
    kind: str = "composed"

    @property
    def length(self) -> float:
        return float(self.s[-1])

    def transformed(self, origin=(0.0, 0.0), angle_deg: float = 0.0) -> "Composition":
        return Composition([p.transformed(origin, angle_deg) for p in self.pieces], self.joints,
                           _rotate(self.axis, origin, angle_deg), self.s, self.half_width, self.formula, self.log, self.kind)


def _rotate(xy: np.ndarray, origin, angle_deg: float) -> np.ndarray:
    a = math.radians(angle_deg)
    c, s = math.cos(a), math.sin(a)
    return np.stack([xy[:, 0] * c - xy[:, 1] * s + origin[0], xy[:, 0] * s + xy[:, 1] * c + origin[1]], axis=1)


def _lines(text: str) -> List[List[str]]:
    out = []
    for raw in re.split(r"[\n;]", text or ""):
        line = raw.split("#")[0].strip()
        if line:
            out.append(line.replace(",", " ").split() if not line.lower().startswith("section custom") else
                       ["section", "custom", line[len("section custom"):].strip()])
    return out


class _Builder:
    def __init__(self, width: float, step: float):
        self.step = float(step)
        self.width = float(width)
        self.x = self.y = self.th = self.s = 0.0
        self.pieces: List[List[tuple]] = [[]]  # (x, y, heading, s, width) per station
        self.joints: List[dict] = []
        self.add_station()

    def add_station(self):
        self.pieces[-1].append((self.x, self.y, self.th, self.s, self.width))

    def run(self, length: float, dtheta, width_end: Optional[float] = None):
        """A segment of axis ``length`` whose heading is th0 + dtheta(u); the width goes to ``width_end``."""
        if length <= 0:
            raise ValueError("segment lengths must be positive")
        n = max(1, int(round(length / self.step)))
        fine = max(n * 40, 400)
        u = np.linspace(0.0, length, fine + 1)
        th = self.th + dtheta(u)
        mid = 0.5 * (th[1:] + th[:-1])
        du = np.diff(u)
        fx = self.x + np.concatenate([[0.0], np.cumsum(du * np.cos(mid))])
        fy = self.y + np.concatenate([[0.0], np.cumsum(du * np.sin(mid))])
        w0 = self.width
        w1 = w0 if width_end is None else float(width_end)
        s0 = self.s
        for us in np.linspace(0.0, length, n + 1)[1:]:
            self.x = float(np.interp(us, u, fx))
            self.y = float(np.interp(us, u, fy))
            self.th = float(np.interp(us, u, th))
            self.s = s0 + float(us)
            self.width = w0 + (w1 - w0) * us / length
            self.add_station()
        self.x, self.y, self.th, self.s, self.width = float(fx[-1]), float(fy[-1]), float(th[-1]), s0 + length, w1

    def new_piece(self, joint: dict):
        if len(self.pieces[-1]) < 2:
            raise ValueError("a %s needs a segment before it" % joint["kind"])
        self.joints.append(joint)
        self.pieces.append([])
        self.add_station()


def _tightest(width: float, radius: float, what: str):
    if radius <= 0.5 * width:
        raise ValueError("%s: the radius on the axis (%.4g m) must exceed half the width (%.4g m), or the inner bank "
                         "crosses itself" % (what, radius, 0.5 * width))


def compose(text: str, width: float, step: float, z0: float = 0.0, slope: float = 0.0,
            section: Optional[Section] = None) -> Composition:
    """The composition described by ``text`` (see the module help), starting at (0, 0) heading along +x."""
    if step <= 0:
        raise ValueError("the element size must be positive")
    cmds = _lines(text)
    if not cmds:
        raise ValueError("no segment given: write one per line, e.g. 'straight 5', 'bend left 2 90'")
    b = _Builder(width, step)
    log: List[str] = []
    bed = [(0.0, float(z0), float(slope), section or Section())]  # (s from, z at s, slope, section)
    started = False
    for words in cmds:
        k = words[0].lower()
        args = words[1:]
        try:
            if k == "width" and not started:
                b.width = _num(args[0], "width")
                b.pieces[-1][-1] = (b.x, b.y, b.th, b.s, b.width)
            elif k in ("width", "expand", "contract", "step"):
                w = _num(args[0], k)
                if w <= 0:
                    raise ValueError("the width must be positive")
                b.width = w
                b.new_piece({"kind": "step", "s": b.s})
                log.append("Sudden width change to %g m at s = %.4g m" % (w, b.s))
            elif k == "straight":
                b.run(_num(args[0], k), lambda u: np.zeros_like(u))
                started = True
            elif k == "taper":
                w, length = _num(args[0], k), _num(args[1], k)
                if w <= 0:
                    raise ValueError("the width must be positive")
                b.run(length, lambda u: np.zeros_like(u), w)
                started = True
            elif k == "bend":
                side, r, ang = args[0].lower(), _num(args[1], k), _num(args[2], k)
                sign = 1.0 if side.startswith("l") else -1.0
                _tightest(b.width, r, "bend")
                b.run(r * math.radians(abs(ang)), lambda u: sign * u / r)
                started = True
            elif k == "elbow":
                side, ang = args[0].lower(), _num(args[1], k)
                if not 0 < abs(ang) < 180:
                    raise ValueError("an elbow turns between 0 and 180 degrees")
                sign = 1.0 if side.startswith("l") else -1.0
                b.new_piece({"kind": "elbow", "angle": sign * abs(ang), "s": b.s})
                b.th += sign * math.radians(abs(ang))
                b.pieces[-1][-1] = (b.x, b.y, b.th, b.s, b.width)
            elif k in ("meanders", "meander", "sine"):
                lam, t0, cnt = (_num(a, k) for a in args[:3])
                rmin = lam / (2 * math.pi * math.radians(t0)) if t0 else math.inf
                _tightest(b.width, rmin, "meanders (tightest bend)")
                b.run(lam * cnt, lambda u: math.radians(t0) * np.sin(2 * math.pi * u / lam))
                log.append("Meanders: tightest bend radius %.4g m" % rmin)
                started = True
            elif k == "kinoshita":
                lam, t0, cnt, js, jf = (_num(a, k) for a in args[:5])
                t = math.radians(t0)
                kk = 2 * math.pi / lam

                def dth(u, t=t, kk=kk, js=js, jf=jf):
                    return (t * np.sin(kk * u) + t ** 3 * (js * np.cos(3 * kk * u) - jf * np.sin(3 * kk * u))
                            - t ** 3 * js)  # shifted so the heading starts unchanged

                uu = np.linspace(0, lam, 2001)
                curv = np.abs(np.gradient(dth(uu), uu))
                rmin = 1.0 / max(float(curv.max()), 1e-300)
                _tightest(b.width, rmin, "Kinoshita meanders (tightest bend)")
                b.run(lam * cnt, dth)
                log.append("Kinoshita meanders: tightest bend radius %.4g m" % rmin)
                started = True
            elif k == "slope":
                z = _bed_at(bed, b.s)
                bed.append((b.s, z, _num(args[0], k), bed[-1][3]))
            elif k == "drop":
                z = _bed_at(bed, b.s) - _num(args[0], k)
                bed.append((b.s, z, bed[-1][2], bed[-1][3]))
            elif k == "bed":
                bed.append((b.s, _num(args[0], k), bed[-1][2], bed[-1][3]))
            elif k == "section":
                bed.append((b.s, _bed_at(bed, b.s), bed[-1][2], parse_section(args)))
            else:
                raise ValueError("unknown segment %r (use width, straight, bend, elbow, meanders, kinoshita, taper, expand, "
                                 "slope, drop, bed or section)" % words[0])
        except IndexError:
            raise ValueError("'%s': missing numbers" % " ".join(words))
        except ValueError as exc:
            raise ValueError("'%s': %s" % (" ".join(words), exc))
    if b.s <= 0:
        raise ValueError("the channel has no length: add a straight, bend or meanders segment")
    pieces, axis, s, hw = _finish(b)
    formula = _bed_formula(bed)
    log.insert(0, "Channel of %.4g m along the axis in %d piece(s) (%s)" % (
        b.s, len(pieces), ", ".join(j["kind"] + ("%+g°" % j["angle"] if j["kind"] == "elbow" else "") for j in b.joints) or "continuous"))
    return Composition(pieces, b.joints, axis, s, hw, formula, log)


def _bed_at(bed, s: float) -> float:
    s0, z0, sl, _sec = bed[-1]
    return z0 - sl * (s - s0)


def _bed_formula(bed) -> str:
    # collapse changes at the same station (the later wins)
    pieces = []
    for item in bed:
        if pieces and abs(item[0] - pieces[-1][0]) < 1e-12:
            pieces[-1] = item
        else:
            pieces.append(item)

    def term(item):
        s0, z0, sl, sec = item
        z = "%g" % z0 if not sl else "%g - %g*(s - %g)" % (z0, sl, s0) if s0 else "%g - %g*s" % (z0, sl)
        f = sec.formula()
        return z if f == "0" else "%s + %s" % (z, f)

    expr = term(pieces[-1])
    for k in range(len(pieces) - 2, -1, -1):
        expr = "(%s) if s < %g else (%s)" % (term(pieces[k]), pieces[k + 1][0], expr)
    return expr


def _finish(b: _Builder):
    pieces = []
    for st in b.pieces:
        a = np.array(st, dtype=np.float64)
        if len(a) < 2:
            continue
        pieces.append(Planform(a[:, :2].copy(), a[:, 2].copy(), a[:, 3].copy(), a[:, 4].copy(), "composed"))
    axis = np.vstack([pieces[0].axis] + [p.axis[1:] for p in pieces[1:]])
    s = np.concatenate([pieces[0].s] + [p.s[1:] for p in pieces[1:]])
    hs, hb = [], []
    for p in pieces:
        hs.append(p.s)
        hb.append(0.5 * np.asarray(p.width, dtype=np.float64))
    sw = np.concatenate(hs)
    bw = np.concatenate(hb)
    # a sudden change has two values at the same s: nudge the second one past it so the profile stays increasing
    for i in range(1, len(sw)):
        if sw[i] <= sw[i - 1]:
            sw[i] = sw[i - 1] + 1e-9 * max(float(sw[-1]), 1.0)
    return pieces, axis, s, (sw, bw)


# ---------------------------------------------------------------------------
# drawn or imported centrelines
# ---------------------------------------------------------------------------
def from_centreline(coords, step: float, width: float = 1.0, widths: Optional[Sequence[float]] = None,
                    corner_deg: float = 30.0, z0: float = 0.0, slope: float = 0.0,
                    section: Optional[Section] = None) -> Composition:
    """A composition along a polyline: cut at vertices turning more than ``corner_deg`` (sharp corners, meshed with
    triangles), resampled every ``step`` between them; ``widths``: one width per vertex (else ``width``)."""
    c = np.asarray(coords, dtype=np.float64)[:, :2]
    keep = np.concatenate([[True], np.hypot(*np.diff(c, axis=0).T) > 0])
    c = c[keep]
    wv = np.full(len(keep), float(width)) if widths is None else np.asarray(widths, dtype=np.float64)
    wv = wv[keep]
    if len(c) < 2:
        raise ValueError("the centreline needs two distinct vertices")
    seg = np.diff(c, axis=0)
    hd = np.arctan2(seg[:, 1], seg[:, 0])
    turn = np.degrees(np.angle(np.exp(1j * (hd[1:] - hd[:-1]))))
    cuts = [0] + [i + 1 for i in np.flatnonzero(np.abs(turn) > corner_deg)] + [len(c) - 1]
    sv = np.concatenate([[0.0], np.cumsum(np.hypot(*seg.T))])
    pieces, joints = [], []
    for a, z in zip(cuts[:-1], cuts[1:]):
        part = c[a:z + 1]
        sp = sv[a:z + 1]
        n = max(1, int(round((sp[-1] - sp[0]) / step)))
        st = np.linspace(sp[0], sp[-1], n + 1)
        x = np.interp(st, sp, part[:, 0])
        y = np.interp(st, sp, part[:, 1])
        w = np.interp(st, sp, wv[a:z + 1])
        xy = np.stack([x, y], axis=1)
        t = np.gradient(xy, axis=0) if len(xy) > 2 else np.repeat(xy[1:] - xy[:1], 2, axis=0)
        t[0], t[-1] = xy[1] - xy[0], xy[-1] - xy[-2]
        th = np.unwrap(np.arctan2(t[:, 1], t[:, 0]))
        pieces.append(Planform(xy, th, st, w, "drawn"))
        if z < len(c) - 1:
            joints.append({"kind": "elbow", "angle": float(turn[z - 1]), "s": float(sv[z])})
    # a smooth heading inside each piece: the half width must stay below the radius of curvature
    for p in pieces:
        if len(p.s) > 2:
            dth = np.abs(np.diff(p.heading)) / np.maximum(np.diff(p.s), 1e-300)
            r = 1.0 / max(float(dth.max()), 1e-300)
            if r <= 0.5 * float(np.max(p.width)):
                raise ValueError("the centreline turns too tightly at about s = %.4g m (radius %.3g m, half width %.3g m): "
                                 "smooth it, or make the corner sharper than %g° so it becomes an elbow"
                                 % (float(p.s[int(np.argmax(dth))]), r, 0.5 * float(np.max(p.width)), corner_deg))
    axis = np.vstack([pieces[0].axis] + [p.axis[1:] for p in pieces[1:]])
    s = np.concatenate([pieces[0].s] + [p.s[1:] for p in pieces[1:]])
    sw = np.concatenate([p.s for p in pieces])
    bw = np.concatenate([0.5 * np.asarray(p.width, dtype=np.float64) for p in pieces])
    for i in range(1, len(sw)):
        if sw[i] <= sw[i - 1]:
            sw[i] = sw[i - 1] + 1e-9 * max(float(sw[-1]), 1.0)
    formula = _bed_formula([(0.0, float(z0), float(slope), section or Section())])
    log = ["Centreline of %.4g m in %d piece(s)%s" % (s[-1], len(pieces), "" if not joints else ", sharp corners at s = " +
                                                         ", ".join("%.4g" % j["s"] for j in joints))]
    return Composition(pieces, joints, axis, s, (sw, bw), formula, log, "drawn")


def read_xy_file(path: str) -> List[np.ndarray]:
    """Polylines from a text file: two or three numbers per line (x, y[, width]) separated by spaces, tabs, commas
    or semicolons; a line without numbers (blank, a header, 'END') ends a polyline. Comment lines start with #."""
    parts: List[List[List[float]]] = [[]]
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("#") or line.startswith(":"):
                continue
            toks = [t for t in re.split(r"[\s,;]+", line) if t]
            vals = []
            for t in toks:
                try:
                    vals.append(float(t))
                except ValueError:
                    vals = []
                    break
            if len(vals) >= 2:
                parts[-1].append(vals[:3])
            elif parts[-1]:
                parts.append([])
    out = []
    for p in parts:
        if len(p) >= 2:
            width = max(len(v) for v in p)
            out.append(np.array([v + [np.nan] * (width - len(v)) for v in p], dtype=np.float64))
    if not out:
        raise ValueError("%s: no polyline found (two or three numbers per line: x y [width])" % path)
    return out
