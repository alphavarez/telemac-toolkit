"""TELEMAC text files for boundary conditions: liquid boundaries file (.liq) and stage-discharge curves.

Liquid boundaries file (keyword LIQUID BOUNDARIES FILE)::

    # comment
    T       Q(1)    SL(2)
    s       m3/s    m
    0.      10.     4.4
    3600.   25.     4.4

Stage-discharge curves file (keyword STAGE-DISCHARGE CURVES FILE), as TELEMAC's READ_FIC_CURVES
reads it: per curve a line naming the columns (Q(n) and Z(n), in either order), a units line, then
the rows; a curve ends at a blank line or a line starting with #::

    # outflow, normal depth
    Q(2)    Z(2)
    m3/s    m
    0.      3.20
    12.5    3.85
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = ["write_liq", "read_liq", "write_stage_discharge", "read_stage_discharge", "read_steering"]


def _fmt(v: float) -> str:
    return "%.6g" % v if abs(v) < 1e7 else "%.10g" % v


def write_liq(path: str, times: Sequence[float], columns: Sequence[Tuple[str, str, Sequence[float]]],
              comment: str = "Liquid boundaries file written by TELEMAC Toolkit") -> str:
    """``columns``: (name such as Q(1) or SL(2), unit, values), in any order; T comes first."""
    t = np.asarray(times, dtype=np.float64)
    if t.ndim != 1 or t.size < 1:
        raise ValueError("no times")
    if np.any(np.diff(t) <= 0):
        raise ValueError("the times must increase")
    for name, _u, vals in columns:
        if len(vals) != t.size:
            raise ValueError("column %s has %d values for %d times" % (name, len(vals), t.size))
    names = ["T"] + [c[0] for c in columns]
    units = ["s"] + [c[1] for c in columns]
    w = [max(12, len(n) + 2) for n in names]
    with open(path, "w", encoding="ascii", newline="\n") as f:
        for line in comment.splitlines():
            f.write("# %s\n" % line)
        f.write("".join(n.ljust(k) for n, k in zip(names, w)).rstrip() + "\n")
        f.write("".join(u.ljust(k) for u, k in zip(units, w)).rstrip() + "\n")
        for i in range(t.size):
            row = [_fmt(t[i])] + [_fmt(float(c[2][i])) for c in columns]
            f.write("".join(v.ljust(k) for v, k in zip(row, w)).rstrip() + "\n")
    return path


def read_liq(path: str) -> Dict[str, np.ndarray]:
    """Columns of a liquid boundaries file {name: values} (names as written, e.g. T, Q(1), SL(2))."""
    rows: List[List[float]] = []
    names: Optional[List[str]] = None
    units_seen = False
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if names is None:
                names = tok
                continue
            if not units_seen:
                units_seen = True
                try:
                    [float(v) for v in tok]
                except ValueError:
                    continue  # the units line
            rows.append([float(v) for v in tok[: len(names)]])
    if names is None:
        raise ValueError("%s: no column names" % path)
    data = np.array(rows, dtype=np.float64).reshape(-1, len(names))
    return {n: data[:, i] for i, n in enumerate(names)}


def write_stage_discharge(path: str, curves: Dict[int, Tuple[Sequence[float], Sequence[float]]],
                          comments: Optional[Dict[int, str]] = None) -> str:
    """``curves``: {liquid boundary number: (Q, Z)}, Q strictly increasing."""
    comments = comments or {}
    with open(path, "w", encoding="ascii", newline="\n") as f:
        for n in sorted(curves):
            q, z = (np.asarray(v, dtype=np.float64) for v in curves[n])
            if q.size != z.size or q.size < 2:
                raise ValueError("curve %d needs at least two points" % n)
            if np.any(np.diff(q) <= 0):
                raise ValueError("the discharges of curve %d must increase" % n)
            f.write("#\n# Stage-discharge curve, liquid boundary %d%s\n#\n"
                    % (n, (": " + comments[n]) if comments.get(n) else ""))
            f.write("Q(%d)        Z(%d)\n" % (n, n))
            f.write("m3/s        m\n")
            for a, b in zip(q, z):
                f.write("%-12s%s\n" % (_fmt(a), "%.4f" % b))
    return path


def read_stage_discharge(path: str) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """{boundary number: (Q, Z)} read the way TELEMAC reads the file."""
    out: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    with open(path, "r", encoding="latin-1") as f:
        lines = f.read().splitlines()
    i = 0
    pat = re.compile(r"^\s*([QZ])\((\d+)\)")
    while i < len(lines):
        s = lines[i]
        if not s.strip() or s.startswith("#"):
            i += 1
            continue
        m = pat.match(s)
        if not m:
            raise ValueError("%s line %d: expected Q(n) or Z(n), found %r" % (path, i + 1, s))
        first, num = m.group(1), int(m.group(2))
        i += 2  # names and units
        qs, zs = [], []
        while i < len(lines) and lines[i].strip() and not lines[i].startswith("#"):
            a, b = (float(v) for v in lines[i].split()[:2])
            if first == "Q":
                qs.append(a)
                zs.append(b)
            else:
                zs.append(a)
                qs.append(b)
            i += 1
        out[num] = (np.array(qs), np.array(zs))
    return out


def read_steering(path: str) -> Dict[str, str]:
    """Keywords of a steering file {KEYWORD: value text}.

    Line based: a line with ``=`` or ``:`` outside quotes starts a keyword (the part before the first
    such separator); other lines continue the value of the previous keyword (lists over several lines).
    ``/`` outside quotes starts a comment. Quotes around a whole value are removed.
    """
    res: Dict[str, str] = {}
    key = None
    with open(path, "r", encoding="latin-1") as f:
        for raw in f:
            out, q, sep = [], None, -1
            for ch in raw.rstrip("\n"):
                if q:
                    out.append(ch)
                    if ch == q:
                        q = None
                    continue
                if ch in "'\"":
                    q = ch
                elif ch == "/":
                    break
                elif ch in "=:" and sep < 0:
                    sep = len(out)
                out.append(ch)
            line = "".join(out)
            if not line.strip():
                continue
            if sep >= 0:
                key = " ".join(line[:sep].split()).upper()
                res[key] = line[sep + 1:].strip()
            elif key is not None:
                res[key] = (res[key] + " " + line.strip()).strip()
    return {k: v.strip().strip("'\"") for k, v in res.items() if k}
