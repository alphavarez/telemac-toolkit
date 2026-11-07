"""Structured channel strips.

A structured mesh of an unbranched channel between two bank lines, optionally
split by an axis (thalweg) line: ``cross_nodes`` nodes span the channel (banks
included; odd with an axis, half on each side), and cross-sections are
``interval`` apart along the channel. Each cell is written as two triangles
(consecutive in ``ikle``); ``pairs_to_quads`` gives the cells back as quads.

Variable width (``across_spacing``): the number of nodes across then follows the
width, so the spacing across stays near ``across_spacing``. It changes by one
interval (on each side of the axis) at a time, with a hysteresis of 20 % and at
least ``min_run`` rows between changes; the row between two cross-sections with
different counts is a row of triangles (a zipper), the other rows stay cells.
``channel_strip`` returns the cells as quads and the rows for quality checks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from ..telemac_io.mesh import signed_areas
from .meshdata import _ccw_quads

__all__ = ["channel_mesh", "channel_strip", "Strip", "strip_quality"]


def _arc(c: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])


def _at(c: np.ndarray, frac: np.ndarray) -> np.ndarray:
    s = _arc(c)
    t = frac * s[-1]
    return np.stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])], axis=1)


def matched_vertices(L: np.ndarray, R: np.ndarray, T=None) -> bool:
    """True when vertex k of each bank (and of the axis) belong to the same cross-section.

    That is how generated channels (bends, meanders) are drawn: the strip then
    uses the vertices as its cross-sections, normal to the axis, instead of
    spacing each bank by its own length (which skews cells in bends, where the
    outer bank is longer than the inner one).
    """
    if len(L) != len(R) or len(L) < 3 or (T is not None and len(T) != len(L)):
        return False
    mid = T if T is not None else 0.5 * (L + R)
    across = R - L
    w = np.hypot(across[:, 0], across[:, 1])
    if np.any(w <= 0):
        return False
    if T is not None:
        t = np.einsum("ij,ij->i", mid - L, across) / w ** 2
        off = np.hypot(*(L + t[:, None] * across - mid).T)
        if np.any(off > 0.05 * w) or np.any((t < 0.02) | (t > 0.98)):
            return False
    tang = np.gradient(mid, axis=0)
    tl = np.hypot(tang[:, 0], tang[:, 1])
    if np.any(tl <= 0):
        return False
    cosang = np.abs(np.einsum("ij,ij->i", tang, across)) / (tl * w)
    return bool(np.all(cosang < 0.35))  # cross-sections within about 20 degrees of normal to the channel


@dataclass
class Strip:
    """A channel strip: the mesh, its cells as quads, and its cross-sections (bank ends, nodes across)."""

    x: np.ndarray
    y: np.ndarray
    ikle: np.ndarray    # the cells (two triangles each, consecutive) first, then the transition triangles
    quads: np.ndarray   # (k, 4) the cells
    left: np.ndarray    # (m, 2) left end of each cross-section
    right: np.ndarray   # (m, 2) right end
    counts: np.ndarray  # (m,) nodes across each cross-section

    def quality(self) -> Dict[str, object]:
        return strip_quality(self.left, self.right, self.counts)


def _counts(width: np.ndarray, halves: int, spacing: float, min_run: int) -> np.ndarray:
    """Intervals per half at each cross-section: about width / (halves x spacing), one change at a time."""
    target = np.maximum(np.asarray(width, dtype=np.float64) / (halves * spacing), 1.0)

    def best(t):  # the count whose spacing is nearest the wanted one, as a ratio
        lo = max(1, int(math.floor(t)))
        return lo if abs(math.log(t / lo)) <= abs(math.log(t / (lo + 1))) else lo + 1

    m = best(float(target[0]))
    out = np.empty(len(target), dtype=np.int64)
    run = min_run
    for j, t in enumerate(target):
        want = best(float(t))
        if want != m and run >= min_run and abs(t / m - 1.0) > 0.2:
            m += 1 if want > m else -1
            run = 0
        out[j] = m
        run += 1
    return out


def _zip(a: np.ndarray, b: np.ndarray, P: np.ndarray) -> list:
    """Triangles between two chains of nodes running the same way across (a row with a count change)."""
    na, nb = len(a) - 1, len(b) - 1
    i = j = 0
    tris = []
    while i < na or j < nb:
        if i == na:
            adv_b = True
        elif j == nb:
            adv_b = False
        else:
            fa, fb = (i + 1) / na, (j + 1) / nb
            if abs(fa - fb) > 1e-9:
                adv_b = fb < fa
            else:  # both next nodes at the same place across: the shorter diagonal
                adv_b = math.dist(P[a[i]], P[b[j + 1]]) < math.dist(P[a[i + 1]], P[b[j]])
        if adv_b:
            tris.append((a[i], b[j + 1], b[j]))
            j += 1
        else:
            tris.append((a[i], a[i + 1], b[j]))
            i += 1
    return tris


def strip_quality(left, right, counts) -> Dict[str, object]:
    """Shape of a strip's cells: along/across length ratio at both banks and skew of the cross-sections.

    ``ratio_min`` and ``ratio_max`` are the cell length along the bank over the spacing across (an inner bank
    in a tight bend gives small values); ``skew_max`` is the largest departure (degrees) of a cross-section
    from normal to a bank. ``*_at`` are the (x, y) places of the worst values.
    """
    L = np.asarray(left, dtype=np.float64)
    R = np.asarray(right, dtype=np.float64)
    n = np.maximum(np.asarray(counts, dtype=np.float64) - 1.0, 1.0)
    across = np.hypot(*(R - L).T) / n
    out: Dict[str, object] = {"ratio_min": float("nan"), "ratio_max": float("nan"), "skew_max": 0.0,
                              "ratio_min_at": None, "ratio_max_at": None, "skew_at": None}
    if len(L) < 2:
        return out
    a_mid = np.maximum(0.5 * (across[:-1] + across[1:]), 1e-300)
    rl = np.hypot(*np.diff(L, axis=0).T) / a_mid
    rr = np.hypot(*np.diff(R, axis=0).T) / a_mid
    ratio = np.concatenate([rl, rr])
    where = np.concatenate([0.5 * (L[:-1] + L[1:]), 0.5 * (R[:-1] + R[1:])])
    k0, k1 = int(np.argmin(ratio)), int(np.argmax(ratio))
    # skew: how far the corners of each cell (bank segment against cross-section) are from 90 degrees
    sec = R - L
    skew = np.zeros(len(L) - 1)
    for B in (L, R):
        e = np.diff(B, axis=0)
        le = np.maximum(np.hypot(*e.T), 1e-300)
        for s in (sec[:-1], sec[1:]):
            c = np.abs(np.einsum("ij,ij->i", e, s)) / (le * np.maximum(np.hypot(*s.T), 1e-300))
            skew = np.maximum(skew, np.degrees(np.arcsin(np.clip(c, 0.0, 1.0))))
    ks = int(np.argmax(skew))
    out.update(ratio_min=float(ratio[k0]), ratio_max=float(ratio[k1]), skew_max=float(skew[ks]),
               ratio_min_at=tuple(float(v) for v in where[k0]), ratio_max_at=tuple(float(v) for v in where[k1]),
               skew_at=tuple(float(v) for v in 0.25 * (L[ks] + R[ks] + L[ks + 1] + R[ks + 1])))
    return out


def channel_mesh(left, right, thalweg=None, cross_nodes: int = 11, interval: float = 10.0, diagonals: str = "shorter",
                 across_spacing: Optional[float] = None):
    """Return (x, y, ikle) of a channel strip.

    ``diagonals``: how each cell is cut into two triangles. "shorter" (default) takes the shorter
    diagonal, ties in one direction, so interior nodes have six triangles as in an unstructured
    mesh; "alternate" flips the diagonal from cell to cell (nodes with four and eight triangles in
    turn, which can make TELEMAC's depths oscillate node to node).
    """
    s = channel_strip(left, right, thalweg, cross_nodes, interval, diagonals, across_spacing)
    return s.x, s.y, s.ikle


def channel_strip(left, right, thalweg=None, cross_nodes: int = 11, interval: float = 10.0, diagonals: str = "shorter",
                  across_spacing: Optional[float] = None, min_run: int = 2, min_nodes: int = 3) -> Strip:
    """The strip of ``channel_mesh`` with its quads and cross-sections; ``across_spacing`` makes the number
    of nodes across follow the width (``cross_nodes`` is then ignored, and ``min_nodes`` is the fewest)."""
    L = np.asarray(left, dtype=np.float64)[:, :2]
    R = np.asarray(right, dtype=np.float64)[:, :2]
    if len(L) < 2 or len(R) < 2:
        raise ValueError("each bank needs at least two vertices")
    if cross_nodes < 2:
        raise ValueError("CrossChannelNodeCount must be at least 2")
    if interval <= 0:
        raise ValueError("AlongChannelInterval must be positive")
    # banks must run the same way
    same = np.hypot(*(L[0] - R[0])) + np.hypot(*(L[-1] - R[-1]))
    swapped = np.hypot(*(L[0] - R[-1])) + np.hypot(*(L[-1] - R[0]))
    if same > swapped:
        R = R[::-1]
    T = None
    if thalweg is not None:
        T = np.asarray(thalweg, dtype=np.float64)[:, :2]
        if np.hypot(*(T[0] - L[0])) > np.hypot(*(T[-1] - L[0])):
            T = T[::-1]
        if cross_nodes % 2 == 0 and not across_spacing:
            raise ValueError("with a thalweg the cross-channel node count must be odd")
    if matched_vertices(L, R, T):
        mid = T if T is not None else 0.5 * (L + R)
        sm = _arc(mid)
        spacing = sm[-1] / (len(mid) - 1)
        if 0.75 * interval <= spacing <= 1.25 * interval:
            st = sm  # the drawn cross-sections themselves
        else:
            st = np.linspace(0.0, sm[-1], max(2, int(round(sm[-1] / float(interval))) + 1))
        k = np.clip(np.searchsorted(sm, st, side="right") - 1, 0, len(mid) - 2)
        f = ((st - sm[k]) / np.maximum(sm[k + 1] - sm[k], 1e-300))[:, None]
        pl = L[k] + f * (L[k + 1] - L[k])
        pr = R[k] + f * (R[k + 1] - R[k])
        pt_matched = None if T is None else T[k] + f * (T[k + 1] - T[k])
        m = len(st)
    else:
        length = 0.5 * (_arc(L)[-1] + _arc(R)[-1])
        m = max(2, int(round(length / float(interval))) + 1)
        frac = np.linspace(0.0, 1.0, m)
        pl, pr = _at(L, frac), _at(R, frac)
        pt_matched = None
    pt = None
    if T is not None:
        pt = pt_matched if pt_matched is not None else _at(T, frac)
    halves = 1 if T is None else 2
    if across_spacing and across_spacing > 0:
        per = _counts(np.hypot(*(pr - pl).T), halves, float(across_spacing), max(int(min_run), 1))
        per = np.maximum(per, max(1, (int(min_nodes) - 1) // halves))
    else:
        n = int(cross_nodes)
        per = np.full(m, n - 1 if T is None else (n + 1) // 2 - 1, dtype=np.int64)
    counts = halves * per + 1
    # nodes, row by row from the left bank (through the axis) to the right bank
    rows, xs, ys = [], [], []
    start = 0
    for j in range(m):
        w = np.linspace(0.0, 1.0, int(per[j]) + 1)
        if T is None:
            p = pl[j][None, :] * (1 - w)[:, None] + pr[j][None, :] * w[:, None]
        else:
            a = pl[j][None, :] * (1 - w)[:, None] + pt[j][None, :] * w[:, None]
            b = pt[j][None, :] * (1 - w)[:, None] + pr[j][None, :] * w[:, None]
            p = np.concatenate([a, b[1:]])
        rows.append(np.arange(start, start + len(p)))
        start += len(p)
        xs.append(p[:, 0])
        ys.append(p[:, 1])
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    P = np.stack([x, y], axis=1)
    pairs, quads, zips = [], [], []
    for j in range(m - 1):
        ra, rb = rows[j], rows[j + 1]
        if len(ra) == len(rb):
            nn = len(ra)
            i = np.arange(nn - 1)
            p00, p01, p10, p11 = ra[i], ra[i + 1], rb[i], rb[i + 1]
            if diagonals == "alternate":
                main = (i + j) % 2 == 0
            else:
                d1 = np.hypot(x[p11] - x[p00], y[p11] - y[p00])
                d2 = np.hypot(x[p10] - x[p01], y[p10] - y[p01])
                main = d1 <= d2 * (1 + 1e-9)
            t_main = np.stack([np.stack([p00, p01, p11], 1), np.stack([p00, p11, p10], 1)], 1)
            t_other = np.stack([np.stack([p00, p01, p10], 1), np.stack([p01, p11, p10], 1)], 1)
            pairs.append(np.where(main[:, None, None], t_main, t_other).reshape(-1, 3))
            quads.append(np.stack([p00, p01, p11, p10], axis=1))
        else:  # a count change: a row of triangles on each side of the axis
            ka, kb = int(per[j]), int(per[j + 1])
            for h in range(halves):
                zips += _zip(ra[h * ka:(h + 1) * ka + 1], rb[h * kb:(h + 1) * kb + 1], P)
    parts = pairs + ([np.array(zips, dtype=np.int64)] if zips else [])
    ikle = np.concatenate(parts).astype(np.int64) if parts else np.zeros((0, 3), dtype=np.int64)
    neg = signed_areas(x, y, ikle) < 0
    ikle[neg] = ikle[neg][:, [0, 2, 1]]
    q = _ccw_quads(x, y, np.concatenate(quads)) if quads else np.zeros((0, 4), dtype=np.int64)
    return Strip(x, y, ikle, q, pl, pr, counts)
