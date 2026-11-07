"""Hydraulics for boundary and initial conditions: rating curves, initial states, boundary checks.

``rating_curve`` computes the normal-depth stage-discharge curve of a cross-section with Manning's
formula, Q = sqrt(S) * sum(K_i), K_i = A_i R_i^(2/3) / n_i, over subsections split at the bank
stations and wherever the roughness changes (the usual divided-channel method, as HEC-RAS splits
the main channel and the overbanks). Each wet stretch of a subsection counts separately. Above
the lower end of the section its ends act as vertical walls, as the closed boundaries of the
TELEMAC mesh do.

Only NumPy is used.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = ["wet_geometry", "rating_curve", "initial_state", "times_to_seconds", "boundary_report"]

G = 9.81


def _clip_segment(s0, z0, s1, z1, level):
    """Wet part of one ground segment below ``level``: (wet width, wet area, wet perimeter)."""
    d0, d1 = level - z0, level - z1
    if d0 <= 0 and d1 <= 0:
        return 0.0, 0.0, 0.0
    ds = s1 - s0
    if d0 > 0 and d1 > 0:
        return ds, 0.5 * (d0 + d1) * ds, math.hypot(ds, z1 - z0)
    # partly wet: the waterline crosses the segment
    f = d0 / (d0 - d1) if d0 > 0 else d1 / (d1 - d0)
    w = abs(ds) * f
    d = d0 if d0 > 0 else d1
    return w, 0.5 * d * w, math.hypot(w, d)


def wet_geometry(st: np.ndarray, z: np.ndarray, level: float, walls: bool = True):
    """Per ground segment: wet width, area and perimeter at ``level``; with ``walls`` the two ends of the
    section are vertical walls up to the level (their wetted height is added to the end perimeters)."""
    n = len(st) - 1
    W, A, P = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        W[i], A[i], P[i] = _clip_segment(st[i], z[i], st[i + 1], z[i + 1], level)
    if walls and n > 0:
        P[0] += max(level - z[0], 0.0)
        P[-1] += max(level - z[-1], 0.0)
    return W, A, P


def _subsections(st: np.ndarray, nseg: np.ndarray, banks: Optional[Sequence[float]]):
    """Segment index ranges [a, b) of the subsections: split at bank stations and where n changes."""
    cuts = {0, len(nseg)}
    for i in range(1, len(nseg)):
        if nseg[i] != nseg[i - 1]:
            cuts.add(i)
    for b in banks or ():
        k = int(np.searchsorted(st, b))
        if 0 < k < len(st):
            cuts.add(k)
    c = sorted(cuts)
    return [(a, b) for a, b in zip(c[:-1], c[1:]) if b > a]


def rating_curve(st, z, n, slope: float, banks: Optional[Sequence[float]] = None, count: int = 40,
                 top: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Normal-depth stage-discharge curve of a cross-section.

    ``st``: stations along the section (m, increasing); ``z``: ground levels there; ``n``: Manning's n,
    one value or one per point; ``slope``: energy slope (m/m, positive). ``banks``: stations where the
    main channel meets the floodplain (optional). Levels run from the lowest point to ``top`` (default:
    the highest ground of the section), denser near the bed. Returns (Z, Q, notes): Q starts at 0 at the
    lowest point and increases strictly, as TELEMAC needs to interpolate Z(Q).
    """
    st = np.asarray(st, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    if st.size < 2 or st.size != z.size:
        raise ValueError("the section needs at least two points with a level each")
    if np.any(~np.isfinite(z)):
        raise ValueError("the section has points without a level (no data)")
    if slope <= 0:
        raise ValueError("the slope must be positive")
    order = np.argsort(st, kind="stable")
    st, z = st[order], z[order]
    nn = np.broadcast_to(np.asarray(n, dtype=np.float64), st.shape)[order]
    if np.any(nn <= 0):
        raise ValueError("Manning's n must be positive")
    nseg = np.round(0.5 * (nn[:-1] + nn[1:]), 6)
    subs = _subsections(st, nseg, banks)
    zmin = float(z.min())
    zt = float(top if top is not None else z.max())
    if zt <= zmin:
        raise ValueError("the top level must be above the lowest point of the section")
    t = np.linspace(0.0, 1.0, int(count))[1:]
    levels = zmin + (zt - zmin) * t ** 1.5
    Z, Q = [zmin], [0.0]
    sq = math.sqrt(slope)
    for lv in levels:
        W, A, P = wet_geometry(st, z, lv)
        q = 0.0
        for a, b in subs:
            # each wet stretch of the subsection separately
            wet = A[a:b] > 0
            k = a
            while k < b:
                if not wet[k - a]:
                    k += 1
                    continue
                j = k
                while j < b and wet[j - a]:
                    j += 1
                Ai, Pi = float(A[k:j].sum()), float(P[k:j].sum())
                if Ai > 0 and Pi > 0:
                    ni = float(np.average(nseg[k:j], weights=np.maximum(W[k:j], 1e-12)))
                    q += sq * Ai * (Ai / Pi) ** (2.0 / 3.0) / ni
                k = j
        if q > Q[-1] * (1 + 1e-9) + 1e-12:
            Z.append(float(lv))
            Q.append(q)
    notes = []
    lowest_end = min(float(z[0]), float(z[-1]))
    if zt > lowest_end + 1e-9:
        notes.append("Above %.3f m the section's ends act as vertical walls (the closed boundaries of the mesh)" % lowest_end)
    notes.append("Curve from %.3f m (Q = 0) to %.3f m (Q = %.4g m3/s), %d points, %d subsection(s), slope %.4g"
                 % (Z[0], Z[-1], Q[-1], len(Q), len(subs), slope))
    return np.array(Z), np.array(Q), notes


def critical_depth_check(Z, Q, st, z):
    """Froude number at each point of a curve (for a warning when the section is supercritical)."""
    st = np.asarray(st, float)
    z = np.asarray(z, float)
    fr = []
    for lv, q in zip(Z, Q):
        W, A, _P = wet_geometry(st, z, lv)
        Atot, T = float(A.sum()), float(W.sum())
        fr.append(0.0 if Atot <= 0 or T <= 0 else (q / Atot) / math.sqrt(G * Atot / T))
    return np.array(fr)


def initial_state(z: np.ndarray, level=None, depth: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Water depth and free surface at nodes: from a water level (one value or one per node), or from a
    depth above the bed. Where the level is below the bed the node is dry (depth 0, surface = bed)."""
    z = np.asarray(z, dtype=np.float64)
    if depth is not None:
        h = np.full(z.shape, max(float(depth), 0.0))
    elif level is not None:
        lv = np.broadcast_to(np.asarray(level, dtype=np.float64), z.shape)
        h = np.maximum(lv - z, 0.0)
    else:
        raise ValueError("give a water level or a depth")
    return h, z + h


_UNITS = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0}


def times_to_seconds(values, unit: str = "s"):
    """Times as seconds from the first one. ``values``: numbers in ``unit`` (s, min, h, d), or
    datetimes (``unit`` 'datetime'). Returns (seconds, first datetime or None)."""
    if unit == "datetime":
        first = values[0]
        secs = np.array([(v - first).total_seconds() for v in values], dtype=np.float64)
        return secs, first
    f = _UNITS.get(unit)
    if f is None:
        raise ValueError("unknown time unit %r (s, min, h, d or datetime)" % unit)
    v = np.asarray(values, dtype=np.float64) * f
    return v, None


def boundary_report(x, y, z, table, liq: Optional[Dict[str, np.ndarray]] = None,
                    stage: Optional[Dict[int, Tuple[np.ndarray, np.ndarray]]] = None,
                    cas: Optional[Dict[str, str]] = None):
    """Liquid boundaries of a mesh with their checks.

    ``table``: BoundaryTable with loops and coordinates (numbered as TELEMAC numbers them); ``liq``: the
    columns of a liquid boundaries file {name: values}; ``stage``: stage-discharge curves {number: (Q, Z)};
    ``cas``: steering keywords {KEYWORD: value}. Returns (boundaries, messages): each boundary is a dict
    with number, name, codes, nodes (global, 0-based), coordinates, length, bed min and max, kind.
    """
    from ..telemac_io.cli import liquid_boundaries

    x, y, z = (np.asarray(v, dtype=np.float64) for v in (x, y, z))
    out, msgs = [], []
    lb = liquid_boundaries(table)
    for b in lb:
        nodes = np.asarray(b["nodes"], dtype=np.int64)
        h, u, v = (int(c) for c in b["codes"])
        kind = {(4, 5, 5): "Q", (5, 4, 4): "H", (5, 5, 5): "QH", (4, 6, 6): "UV", (5, 6, 6): "UVH",
                (4, 4, 4): "free", (1, 1, 1): "waves"}.get((h, u, v), "%d%d%d" % (h, u, v))
        xs, ys = x[nodes], y[nodes]
        length = float(np.hypot(np.diff(xs), np.diff(ys)).sum()) if nodes.size > 1 else 0.0
        zb = z[nodes] if nodes.size else np.zeros(0)
        name = table.comments[b["start"]] if b["start"] < len(table.comments) else ""
        d = {"number": b["number"], "name": name, "codes": (h, u, v), "kind": kind, "nodes": nodes,
             "x": xs, "y": ys, "length": length, "zmin": float(zb.min()) if zb.size else float("nan"),
             "zmax": float(zb.max()) if zb.size else float("nan")}
        out.append(d)
        label = "Boundary %d%s (%s, %d %d %d)" % (d["number"], " " + name if name else "", kind, h, u, v)
        msgs.append("%s: %d nodes, %.4g m long, bed %.3f to %.3f m" % (label, nodes.size, length, d["zmin"], d["zmax"]))
        if nodes.size < 5:
            msgs.append("  Warning: only %d nodes; 5 or more are recommended" % nodes.size)
        if liq is not None:
            need = []
            if u == 5:
                need.append("Q(%d)" % d["number"])
            if h == 5 and not (stage and d["number"] in stage):
                need.append("SL(%d)" % d["number"])
            if u == 6:
                need += ["U(%d)" % d["number"], "V(%d)" % d["number"]]
            missing = [c for c in need if c not in liq]
            if missing:
                msgs.append("  Warning: the liquid boundaries file has no column %s; TELEMAC takes the steering file values"
                            % ", ".join(missing))
            sl = liq.get("SL(%d)" % d["number"])
            if sl is not None and zb.size:
                if float(np.max(sl)) <= d["zmin"]:
                    msgs.append("  Warning: the prescribed level never rises above the lowest bed of the boundary "
                                "(%.3f m): the boundary is dry" % d["zmin"])
                elif float(np.min(sl)) < d["zmin"]:
                    msgs.append("  Note: at times the prescribed level is below the whole bed of the boundary")
            q = liq.get("Q(%d)" % d["number"])
            if q is not None and float(np.min(q)) < 0:
                msgs.append("  Note: negative discharges (water leaving through a Q boundary) at some times")
        if stage and d["number"] in stage:
            if h != 5:
                msgs.append("  Warning: a stage-discharge curve needs a prescribed elevation boundary (5 4 4), not %d %d %d"
                            % (h, u, v))
            qq, zz = stage[d["number"]]
            if np.any(np.diff(qq) <= 0):
                msgs.append("  Warning: the discharges of the stage-discharge curve must increase")
            if zb.size and float(np.min(zz)) > d["zmin"] + 0.05:
                msgs.append("  Note: the curve starts at %.3f m, above the lowest bed of the boundary (%.3f m)"
                            % (float(np.min(zz)), d["zmin"]))
    # two liquid boundaries touching each other
    lih = table.records["lihbor"]
    for a, b in (table.loops or [(0, len(table))]):
        m = b - a
        for k in range(m):
            i, j = a + k, a + (k + 1) % m
            if lih[i] != 2 and lih[j] != 2 and tuple(table.codes[i]) != tuple(table.codes[j]):
                msgs.append("Warning: two liquid boundaries with different types touch at boundary nodes %d and %d; "
                            "leave a wall node between them" % (i + 1, j + 1))
    if liq is not None:
        t = liq.get("T")
        if t is None:
            msgs.append("Warning: the liquid boundaries file has no T column")
        elif np.any(np.diff(t) <= 0):
            msgs.append("Warning: the times of the liquid boundaries file must increase")
        known = {"T"} | {"%s(%d)" % (c, d["number"]) for d in out for c in ("Q", "SL", "U", "V")}
        extra = [c for c in liq if c not in known and not c.startswith(("TR(", "C("))]
        if extra:
            msgs.append("Warning: columns for boundaries the mesh does not have: %s" % ", ".join(extra))
    if cas:
        nlb = len(out)
        for key in ("PRESCRIBED FLOWRATES", "PRESCRIBED ELEVATIONS", "STAGE-DISCHARGE CURVES"):
            val = cas.get(key)
            if val is not None:
                k = len([p for p in str(val).replace(",", ";").split(";") if p.strip() != ""])
                if k != nlb:
                    msgs.append("Warning: %s has %d values for %d liquid boundaries" % (key, k, nlb))
        sdc = cas.get("STAGE-DISCHARGE CURVES")
        if sdc is not None:
            vals = [p.strip() for p in str(sdc).replace(",", ";").split(";")]
            for d, v in zip(out, vals):
                if v not in ("", "0") and d["codes"][0] != 5:
                    msgs.append("Warning: boundary %d has a stage-discharge curve but is not a prescribed elevation"
                                % d["number"])
    if not out:
        msgs.append("No liquid boundary: every boundary node is a wall.")
    return out, msgs
