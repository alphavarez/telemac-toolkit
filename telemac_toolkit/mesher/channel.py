"""Structured channel strips.

A structured mesh of an unbranched channel between two bank lines, optionally
split by an axis (thalweg) line: ``cross_nodes`` nodes span the channel (banks
included; odd with an axis, half on each side), and cross-sections are
``interval`` apart along the channel. Each cell is written as two triangles
(consecutive in ``ikle``); ``pairs_to_quads`` gives the cells back as quads.
"""

from __future__ import annotations

import numpy as np

from ..telemac_io.mesh import signed_areas

__all__ = ["channel_mesh"]


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


def channel_mesh(left, right, thalweg=None, cross_nodes: int = 11, interval: float = 10.0):
    """Return (x, y, ikle) of a channel strip."""
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
        if cross_nodes % 2 == 0:
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
    n = int(cross_nodes)
    if T is None:
        w = np.linspace(0.0, 1.0, n)
        pts = pl[:, None, :] * (1 - w)[None, :, None] + pr[:, None, :] * w[None, :, None]
    else:
        pt = pt_matched if pt_matched is not None else _at(T, frac)
        half = (n + 1) // 2
        w = np.linspace(0.0, 1.0, half)
        a = pl[:, None, :] * (1 - w)[None, :, None] + pt[:, None, :] * w[None, :, None]
        b = pt[:, None, :] * (1 - w)[None, :, None] + pr[:, None, :] * w[None, :, None]
        pts = np.concatenate([a, b[:, 1:, :]], axis=1)
    x = pts[:, :, 0].ravel()
    y = pts[:, :, 1].ravel()
    tris = []
    for j in range(m - 1):
        for i in range(n - 1):
            p00 = j * n + i
            p01 = p00 + 1
            p10 = p00 + n
            p11 = p10 + 1
            if (i + j) % 2 == 0:
                tris += [[p00, p01, p11], [p00, p11, p10]]
            else:
                tris += [[p00, p01, p10], [p01, p11, p10]]
    ikle = np.array(tris, dtype=np.int64)
    neg = signed_areas(x, y, ikle) < 0
    ikle[neg] = ikle[neg][:, [0, 2, 1]]
    return x, y, ikle
