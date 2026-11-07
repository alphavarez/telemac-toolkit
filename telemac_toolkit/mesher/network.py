"""River networks: centrelines that join at confluences or split at bifurcations.

A line that ends on another one (within ``tol``) splits it there, and lines drawn end to end are joined, so the
network becomes a set of reaches between junctions and free ends. Every reach is a structured strip
(mesher.river, mesher.channel) that stops short of its junctions: past the last cross-section that still
touches another channel, plus a margin of a quarter of its width (at least one section).

The channel area left at a junction is bounded by the strip ends and by the banks: going round the junction,
each bank is followed from its strip end towards the junction until it meets the bank of the next branch, so
the confluence keeps its corners (where the banks of two parallel branches never meet, they are joined across
the junction instead). That area is meshed with triangles; its banks become guide lines.

Free ends: the start of a line is an inflow, its end an outflow, so lines must run downstream (drawn that way,
or turned with a DEM before they get here).

Only NumPy is used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import points_in_rings, segment_intersections
from .river import Reach, _arc, _clean, _orient, reach_from_centreline

__all__ = ["Junction", "Network", "split_network", "river_network"]


@dataclass
class Junction:
    point: np.ndarray
    members: List[Tuple[int, int]]  # (reach, end: 0 = its start, -1 = its end)
    ring: Optional[np.ndarray] = None  # closed ring of the channel area between the strip ends
    banks: List[np.ndarray] = field(default_factory=list)  # banks across the junction, strip corner to strip corner


@dataclass
class Network:
    reaches: List[Reach]  # the strips: reaches cut short of their junctions
    full: List[Reach]  # the whole reaches, junction to junction
    junctions: List[Junction]
    ends: List[Tuple[int, int, str]]  # free ends: (reach, 0 or -1, "inflow" or "outflow")
    log: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# polylines
# ---------------------------------------------------------------------------
def _project(p: np.ndarray, c: np.ndarray) -> Tuple[float, float, np.ndarray]:
    """(distance, arc length along c, nearest point) of point p on polyline c."""
    A, B = c[:-1], c[1:]
    E = B - A
    L2 = np.maximum(np.einsum("ij,ij->i", E, E), 1e-300)
    t = np.clip(np.einsum("ij,ij->i", p[None, :] - A, E) / L2, 0.0, 1.0)
    Q = A + t[:, None] * E
    d = np.hypot(*(Q - p[None, :]).T)
    k = int(np.argmin(d))
    return float(d[k]), float(_arc(c)[k] + t[k] * math.sqrt(L2[k])), Q[k].copy()


def _cut(c: np.ndarray, s0: float, s1: float) -> np.ndarray:
    """Part of polyline c between arc lengths s0 < s1."""
    s = _arc(c)

    def pt(v):
        return np.array([np.interp(v, s, c[:, 0]), np.interp(v, s, c[:, 1])])

    inner = c[(s > s0) & (s < s1)]
    return np.vstack([pt(s0), inner, pt(s1)])


class _Sets:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


# ---------------------------------------------------------------------------
# topology
# ---------------------------------------------------------------------------
def split_network(lines: Sequence, tol: float):
    """Pieces of the network: [(coords, start junction or None, end junction or None, line index)], junction points,
    and notes. Lines are split where another line ends on them; lines meeting end to end are joined."""
    lines = [_clean(c) for c in lines]
    notes: List[str] = []
    n = len(lines)
    ends = [(i, e) for i in range(n) for e in (0, -1)]
    P = np.array([lines[i][e] for i, e in ends])
    sets = _Sets(len(ends))
    for a in range(len(ends)):
        for b in range(a + 1, len(ends)):
            if ends[a][0] != ends[b][0] and math.dist(P[a], P[b]) <= tol:
                sets.union(a, b)
    groups = {}
    for a in range(len(ends)):
        groups.setdefault(sets.find(a), []).append(a)
    points: List[np.ndarray] = []
    at_end = {}  # (line, end) -> junction
    cuts: List[Tuple[int, int]] = []  # (line, junction): the junction is inside that line
    for members in groups.values():
        centre = P[members].mean(axis=0)
        inside = {ends[a][0] for a in members}
        on = []
        for j in range(n):
            if j in inside:
                continue
            d, s, q = _project(centre, lines[j])
            if d <= tol and 0.0 < s < _arc(lines[j])[-1]:
                on.append((d, j, q))
        if len(members) == 1 and not on:
            continue  # a free end
        if len(on) > 1:
            notes.append("The junction at (%.1f, %.1f) lies on %d lines; it is put on the nearest" % (*centre, len(on)))
        k = len(points)
        if on:
            d, j, q = min(on, key=lambda v: v[0])
            points.append(q)
            cuts.append((j, k))
        else:
            points.append(centre)
        for a in members:
            at_end[ends[a]] = k
    pieces = []
    for i, c0 in enumerate(lines):
        c = c0.copy()
        js, je = at_end.get((i, 0)), at_end.get((i, -1))
        if js is not None:
            c[0] = points[js]
        if je is not None:
            c[-1] = points[je]
        L = float(_arc(c)[-1])
        bounds = [(0.0, js)] + sorted((_project(points[k], c)[1], k) for j, k in cuts if j == i) + [(L, je)]
        for (sa, ja), (sb, jb) in zip(bounds[:-1], bounds[1:]):
            if sb - sa <= 1e-9 * max(L, 1.0):
                continue
            seg = _cut(c, sa, sb)
            if ja is not None:
                seg[0] = points[ja]
            if jb is not None:
                seg[-1] = points[jb]
            pieces.append([seg, ja, jb, i])
    # junctions with two piece ends are not junctions: join the pieces
    while True:
        deg = {}
        for p, (_c, ja, jb, _i) in enumerate(pieces):
            for k in (ja, jb):
                if k is not None:
                    deg.setdefault(k, []).append(p)
        two = [(k, ps) for k, ps in deg.items() if len(ps) == 2 and ps[0] != ps[1]]
        if not two:
            break
        k, (pa, pb) = two[0]
        A, B = list(pieces[pa]), list(pieces[pb])
        flipped = False
        if A[2] != k:  # A starts at k
            if B[2] == k:
                A, B = B, A  # B ends at k: it goes first
            else:
                A = [A[0][::-1], A[2], A[1], A[3]]
                flipped = True
        if B[1] != k:  # B ends at k as well
            B = [B[0][::-1], B[2], B[1], B[3]]
            flipped = True
        if flipped:
            notes.append("Two lines drawn in opposite directions were joined at (%.1f, %.1f); check the direction"
                         % tuple(points[k]))
        merged = [np.vstack([A[0], B[0][1:]]), A[1], B[2], A[3]]
        pieces = [p for q, p in enumerate(pieces) if q not in (pa, pb)] + [merged]
    # one network
    sets = _Sets(len(pieces))
    by_j = {}
    for p, (_c, ja, jb, _i) in enumerate(pieces):
        for k in (ja, jb):
            if k is not None:
                if k in by_j:
                    sets.union(p, by_j[k])
                by_j[k] = p
    comps = len({sets.find(p) for p in range(len(pieces))})
    if comps > 1:
        raise ValueError("the centrelines do not join into one network: they form %d separate groups. End each "
                         "tributary on the river it joins (snapping helps), within %.3g m" % (comps, tol))
    return [tuple(p) for p in pieces], points, notes


# ---------------------------------------------------------------------------
# junction geometry
# ---------------------------------------------------------------------------
def _outgoing(r: Reach, e: int):
    """Axis, left bank and right bank of reach r looking away from its end e."""
    if e == 0:
        return r.axis, r.left, r.right
    return r.axis[::-1], r.right[::-1], r.left[::-1]


def _seg_hits_ring(a: np.ndarray, b: np.ndarray, ring: np.ndarray) -> bool:
    if points_in_rings(np.array([a[0], b[0]]), np.array([a[1], b[1]]), [ring]).any():
        return True
    p, q = ring[:-1], ring[1:]
    d1 = _orient(a, b, p)
    d2 = _orient(a, b, q)
    d3 = _orient(p, q, a)
    d4 = _orient(p, q, b)
    return bool(np.any((d1 * d2 < 0) & (d3 * d4 < 0)))


def _bank_through(path_a: np.ndarray, path_b: np.ndarray, centre: np.ndarray, reach_max: float) -> np.ndarray:
    """Bank line from path_a[0] to path_b[0]: each path runs from a strip corner towards the junction (and on, past
    it); they are cut where they first meet, or joined across the junction if they do not meet near it."""
    hits = segment_intersections(path_a[:-1], path_a[1:], path_b[:-1], path_b[1:])
    best = None
    la = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(path_a, axis=0).T))])
    lb = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(path_b, axis=0).T))])
    for i, j, t, u, x, y in hits:
        if math.dist((x, y), centre) > reach_max:
            continue
        cost = la[i] + t * (la[i + 1] - la[i]) + lb[j] + u * (lb[j + 1] - lb[j])
        if best is None or cost < best[0]:
            best = (cost, i, j, x, y)
    if best is not None:
        _c, i, j, x, y = best
        return np.vstack([path_a[:i + 1], [[x, y]], path_b[:j + 1][::-1]])
    return np.vstack([path_a[:-1], path_b[:-1][::-1]])  # without the extensions past the junction


def _drop_repeats(c: np.ndarray, eps: float) -> np.ndarray:
    keep = np.concatenate([[True], np.hypot(*np.diff(c, axis=0).T) > eps])
    return c[keep]


def _ring_is_simple(ring: np.ndarray) -> bool:
    a, b = ring[:-1], ring[1:]
    for i, j, *_rest in segment_intersections(a, b, a, b):
        if abs(i - j) > 1 and not (i == 0 and j == len(a) - 1) and not (j == 0 and i == len(a) - 1):
            return False
    return True


# ---------------------------------------------------------------------------
# the network
# ---------------------------------------------------------------------------
def river_network(lines: Sequence, step: float, widths: Optional[Sequence[Optional[float]]] = None,
                  dem: Optional[Callable] = None, bankfull_depth: Optional[float] = None,
                  max_width: Optional[float] = None, min_width: Optional[float] = None,
                  smoothing: Optional[float] = None, tol: Optional[float] = None, max_ratio: float = 4.0) -> Network:
    """Reaches, junctions and free ends of a river network drawn as several centrelines (running downstream).

    ``widths``: one constant width per line (None where the banks come from the DEM); ``tol``: how close a line
    end must be to another line to join it (default: a quarter of the largest width, at least the spacing).
    """
    lines = [_clean(c) for c in lines]
    widths = list(widths) if widths is not None else [None] * len(lines)
    if len(widths) != len(lines):
        raise ValueError("give one width per line")
    ref = max([w for w in widths if w] + [0.5 * float(max_width or 0.0), step])
    tol = float(tol) if tol else max(step, 0.25 * ref)
    pieces, points, notes = split_network(lines, tol)
    log: List[str] = list(notes)
    warnings: List[str] = []
    full: List[Reach] = []
    for k, (c, _ja, _jb, i) in enumerate(pieces):
        w = widths[i]
        try:
            r = reach_from_centreline(c, step, width=w, dem=dem if not w else None,
                                      bankfull_depth=None if w else bankfull_depth, max_width=max_width,
                                      min_width=min_width, smoothing=smoothing, orient=False, max_ratio=max_ratio)
        except ValueError as exc:
            raise ValueError("reach %d (%.4g m): %s" % (k + 1, float(_arc(c)[-1]), exc))
        full.append(r)
        warnings += ["Reach %d: %s" % (k + 1, m) for m in r.warnings]
    junctions = []
    for k, p in enumerate(points):
        members = [(r, 0) for r, (_c, ja, _jb, _i) in enumerate(pieces) if ja == k]
        members += [(r, -1) for r, (_c, _ja, jb, _i) in enumerate(pieces) if jb == k]
        if members:
            junctions.append(Junction(np.asarray(p, dtype=np.float64), members))
    log.append("River network: %d lines, %d reaches, %d junctions" % (len(lines), len(pieces), len(junctions)))
    # how far each strip stays from its junctions
    trim = {}
    outlines = [r.outline() for r in full]
    for J in junctions:
        wmax = max(float(np.median(full[r].width[:5] if e == 0 else full[r].width[-5:])) for r, e in J.members)
        reach_max = 6.0 * wmax
        for r, e in J.members:
            A, Pb, Qb = _outgoing(full[r], e)
            s = _arc(A)
            kmax = int(min(np.searchsorted(s, reach_max), len(A) - 1))
            last = 0
            for kk in range(kmax + 1):
                if any(_seg_hits_ring(Pb[kk], Qb[kk], outlines[o]) for o, _e in J.members if o != r):
                    last = kk
            if last >= kmax and kmax < len(A) - 1:
                warnings.append("Reach %d runs next to another branch for more than %.3g m from the junction at "
                                "(%.1f, %.1f); its strip starts there" % (r + 1, reach_max, *J.point))
            w_here = float(np.median(full[r].width[:5] if e == 0 else full[r].width[-5:]))
            spacing = float(s[-1]) / max(len(A) - 1, 1)
            margin = max(1, int(math.ceil(0.25 * w_here / max(spacing, 1e-300))))
            trim[(r, e)] = last + 1 + margin
    strips = []
    for r, reach in enumerate(full):
        a = trim.get((r, 0), 0)
        b = len(reach.axis) - trim.get((r, -1), 0)
        if b - a < 3:
            raise ValueError("reach %d (%.4g m) is too short for a structured strip between its junctions; use a "
                             "smaller cross-section spacing, or move the junctions apart" % (r + 1, reach.length))
        strips.append(reach.trimmed(np.arange(a, b)) if (a > 0 or b < len(reach.axis)) else reach)
    # the channel area of each junction
    for J in junctions:
        wmax = max(float(np.max(full[r].width)) for r, _e in J.members)
        ext = 3.0 * wmax
        order = []
        for r, e in J.members:
            A, Pb, Qb = _outgoing(full[r], e)
            t = trim[(r, e)]
            s = _arc(A)
            kd = int(min(np.searchsorted(s, max(full[r].width[0 if e == 0 else -1], step)), len(A) - 1))
            d = A[max(kd, 1)] - A[0]
            order.append((math.atan2(d[1], d[0]), r, e, Pb, Qb, t))
        order.sort(key=lambda v: v[0])
        ring, banks = [], []
        for i, (_ang, r, e, Pb, Qb, t) in enumerate(order):
            _a2, r2, e2, Pb2, Qb2, t2 = order[(i + 1) % len(order)]

            def path(B, tt):
                pa = B[tt::-1]
                u = B[0] - B[1]
                u = u / max(float(np.hypot(*u)), 1e-300)
                return np.vstack([pa, [B[0] + ext * u]])

            bank = _drop_repeats(_bank_through(path(Pb, t), path(Qb2, t2), J.point, 2.0 * wmax), 1e-9 * wmax)
            banks.append(bank)
            ring.append(Qb[t][None, :])
            ring.append(bank[:-1])
        R = np.vstack(ring)
        J.ring = np.vstack([R, R[:1]])
        J.banks = [b for b in banks if len(b) >= 2 and float(np.hypot(*np.diff(b, axis=0).T).sum()) > 1e-6 * wmax]
        if not _ring_is_simple(J.ring):
            warnings.append("The channel area at the junction at (%.1f, %.1f) crosses itself: check the banks there "
                            "(widths, or branches meeting at a very small angle)" % tuple(J.point))
    ends = []
    used = {m for J in junctions for m in J.members}
    for r in range(len(full)):
        if (r, 0) not in used:
            ends.append((r, 0, "inflow"))
        if (r, -1) not in used:
            ends.append((r, -1, "outflow"))
    for r, (st, rf) in enumerate(zip(strips, full)):
        cut = rf.length - st.length
        log.append("Reach %d: %.4g m, width %.3g to %.3g m; its strip %.4g m (%.3g m left to the junctions)"
                   % (r + 1, rf.length, float(rf.width.min()), float(rf.width.max()), st.length, cut))
    return Network(strips, full, junctions, ends, log, warnings)
