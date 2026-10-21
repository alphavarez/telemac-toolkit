"""TELEMAC boundary conditions file (.cli, CONLIM table).

One line per boundary node, in boundary numbering order:

    LIHBOR LIUBOR LIVBOR HBOR UBOR VBOR AUBOR LITBOR TBOR ATBOR BTBOR N K [# comment]

Codes: 2 = wall, 4 = free, 5 = prescribed, 6 = prescribed velocity (UV),
1 = incident waves. A node's (LIHBOR, LIUBOR, LIVBOR) triple gives its type,
e.g. 2 2 2 wall, 4 5 5 prescribed Q, 5 4 4 prescribed H.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

__all__ = ["BoundaryTable", "read_cli", "write_cli", "BOUNDARY_TYPES", "liquid_boundaries"]

#: TELEMAC boundary types (LIHBOR, LIUBOR, LIVBOR) with a label and a map colour.
BOUNDARY_TYPES = {
    (2, 2, 2): ("Closed boundary (wall)", "#8c5a2b"),
    (4, 5, 5): ("Open boundary, prescribed Q", "#1f5fd1"),
    (5, 4, 4): ("Open boundary, prescribed H", "#2e9e44"),
    (5, 5, 5): ("Open boundary, prescribed Q and H", "#17b8c4"),
    (4, 6, 6): ("Open boundary, prescribed UV", "#d62728"),
    (5, 6, 6): ("Open boundary, prescribed UV and H", "#f28e1c"),
    (1, 1, 1): ("Open boundary, incident waves", "#e3c21b"),
    (4, 4, 4): ("Open boundary, free", "#7f7f7f"),
}

_DTYPE = np.dtype([
    ("lihbor", "i4"), ("liubor", "i4"), ("livbor", "i4"),
    ("hbor", "f8"), ("ubor", "f8"), ("vbor", "f8"), ("aubor", "f8"),
    ("litbor", "i4"), ("tbor", "f8"), ("atbor", "f8"), ("btbor", "f8"),
    ("n", "i8"), ("k", "i8"),
])


@dataclass
class BoundaryTable:
    records: np.ndarray                      # structured array, _DTYPE
    comments: List[str] = field(default_factory=list)
    #: optional: record ranges [a, b) of each boundary loop and the boundary node coordinates, so the
    #: liquid boundaries are numbered exactly as TELEMAC numbers them (see liquid_boundaries)
    loops: Optional[List[Tuple[int, int]]] = None
    xy: Optional[Tuple[np.ndarray, np.ndarray]] = None

    def __len__(self):
        return self.records.size

    @property
    def nodes(self) -> np.ndarray:
        """Global node numbers, 0-based."""
        return self.records["n"] - 1

    @property
    def codes(self) -> np.ndarray:
        r = self.records
        return np.stack([r["lihbor"], r["liubor"], r["livbor"]], axis=1)

    @classmethod
    def walls(cls, boundary_nodes: np.ndarray) -> "BoundaryTable":
        """Default table: every boundary node is a wall (2 2 2)."""
        n = len(boundary_nodes)
        rec = np.zeros(n, dtype=_DTYPE)
        rec["lihbor"] = rec["liubor"] = rec["livbor"] = 2
        rec["litbor"] = 2
        rec["n"] = np.asarray(boundary_nodes) + 1
        rec["k"] = np.arange(1, n + 1)
        return cls(rec, [""] * n)


def read_cli(path: str) -> BoundaryTable:
    rows, comments = [], []
    with open(path, "r", encoding="latin-1") as f:
        for lineno, line in enumerate(f, 1):
            body, _, comment = line.partition("#")
            tok = body.split()
            if not tok:
                continue
            if len(tok) < 13:
                raise ValueError("%s line %d: expected 13 values, found %d" % (path, lineno, len(tok)))
            try:
                rows.append((
                    int(tok[0]), int(tok[1]), int(tok[2]),
                    float(tok[3]), float(tok[4]), float(tok[5]), float(tok[6]),
                    int(tok[7]), float(tok[8]), float(tok[9]), float(tok[10]),
                    int(tok[11]), int(tok[12]),
                ))
            except ValueError as exc:
                raise ValueError("%s line %d: %s" % (path, lineno, exc)) from None
            extra = " ".join(tok[13:])
            comments.append((extra + " " + comment.strip()).strip())
    return BoundaryTable(np.array(rows, dtype=_DTYPE), comments)


def write_cli(path: str, table: BoundaryTable) -> str:
    r = table.records
    with open(path, "w", encoding="latin-1", newline="\n") as f:
        for i in range(r.size):
            line = "%d %d %d %.3f %.3f %.3f %.3f %d %.3f %.3f %.3f %d %d" % (
                r["lihbor"][i], r["liubor"][i], r["livbor"][i],
                r["hbor"][i], r["ubor"][i], r["vbor"][i], r["aubor"][i],
                r["litbor"][i], r["tbor"][i], r["atbor"][i], r["btbor"][i],
                r["n"][i], r["k"][i],
            )
            c = table.comments[i] if i < len(table.comments) else ""
            f.write(line + (" # " + c if c else "") + "\n")
    return path


def _telemac_start(xs: np.ndarray, ys: np.ndarray, idx: np.ndarray) -> int:
    """FRONT2's first point among ``idx``: smallest x + y, ties (within 1e-4 of the range) by smallest y."""
    ns = xs[idx] + ys[idx]
    k0 = int(np.argmin(ns))
    eps = (ns.max() - ns.min()) * 1e-4
    best, ymin = k0, ys[idx][k0]
    for k in range(len(idx)):
        if abs(ns[k0] - ns[k]) < eps and ys[idx][k] < ymin:
            best, ymin = k, ys[idx][k]
    return int(idx[best])


def liquid_boundaries(table: BoundaryTable) -> List[dict]:
    """Liquid boundaries, numbered as TELEMAC numbers them (BIEF FRONT2).

    TELEMAC walks each boundary loop from its point with the smallest x + y
    (ties: smallest y), the loops in the order of those points, and numbers
    the liquid boundaries as it meets them; two adjacent liquid runs with
    different codes are two boundaries, and a run through the start point is
    one. A table read from a .cli alone (no coordinates) is walked from its
    first record.
    """
    lih = table.records["lihbor"].astype(np.int64)
    liu = table.records["liubor"].astype(np.int64)
    codes = table.codes
    n = lih.size
    if n == 0:
        return []
    pending = [np.arange(a, b) for a, b in (table.loops or [(0, n)]) if b > a]
    order = []
    if table.xy is not None:
        xs, ys = (np.asarray(v, dtype=np.float64) for v in table.xy)
        while pending:
            first = _telemac_start(xs, ys, np.concatenate(pending))
            k = next(i for i, idx in enumerate(pending) if first in idx)
            idx = pending.pop(k)
            order.append(np.roll(idx, -int(np.flatnonzero(idx == first)[0])))
    else:
        order = pending
    found = []  # (records of one liquid boundary), in TELEMAC order
    for idx in order:
        m = idx.size
        runs = []  # [first position, last position] along idx
        start_liquid = not (lih[idx[0]] == 2 or lih[idx[1 % m]] == 2)
        if start_liquid:
            runs.append([0, 0])
        cur = start_liquid
        for j in range(1, m):
            l1, l2, l3 = lih[idx[j - 1]], lih[idx[j]], lih[idx[(j + 1) % m]]
            if l1 == 2 and l2 != 2 and l3 != 2:
                runs.append([j, j])
                cur = True
            elif l1 != 2 and l2 != 2 and l3 == 2:
                if runs:
                    runs[-1][1] = j
                cur = False
            elif l1 != 2 and l2 != 2 and l3 != 2:
                if runs:
                    runs[-1][1] = j
                if l2 != l3 or liu[idx[j]] != liu[idx[(j + 1) % m]]:
                    if j + 1 < m:
                        runs.append([j + 1, j + 1])
                    elif start_liquid:
                        cur = False  # the run ends exactly at the start point
            elif l1 == 2 and l2 != 2 and l3 == 2:
                runs.append([j, j])  # a single liquid point (TELEMAC stops on it)
                cur = False
        if cur and runs:
            if start_liquid and len(runs) > 1:
                last = runs.pop()  # continues through the start point: part of the first boundary
                runs[0] = [last[0], runs[0][1] + m]
            else:
                runs[-1][1] = m - 1 if not start_liquid else runs[-1][1]
        for a, b in runs:
            sel = idx[np.arange(a, b + 1) % m]
            found.append(sel)
    out = []
    for number, sel in enumerate(found, 1):
        out.append({"number": number, "start": int(sel[0]), "end": int(sel[-1]), "nodes": table.nodes[sel],
                    "records": sel, "codes": tuple(codes[sel[0]])})
    return out
