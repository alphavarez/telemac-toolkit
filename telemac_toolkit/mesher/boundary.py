"""Boundary conditions from lines drawn along the mesh boundary (CONLIM).

Every boundary node starts as a wall (2 2 2). Boundary nodes
within a tolerance of a boundary line take that line's type. The table is
written in TELEMAC's boundary order (outer boundary counter-clockwise from its
lowest node, then islands), the same order the geometry writer uses for IPOBO.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..telemac_io.cli import BoundaryTable, liquid_boundaries
from ..telemac_io.mesh import boundary_loops
from .geometry import distance_to_segments

__all__ = ["BOUNDARY_CODES", "BoundaryLine", "build_boundary_table", "liq_template", "parse_type"]

#: Boundary types -> TELEMAC LIHBOR, LIUBOR, LIVBOR
BOUNDARY_CODES = {
    "wall": (2, 2, 2),
    "q": (4, 5, 5),          # prescribed discharge
    "h": (5, 4, 4),          # prescribed elevation
    "qh": (5, 5, 5),         # prescribed discharge and elevation
    "uv": (4, 6, 6),         # prescribed velocities
    "uvh": (5, 6, 6),        # prescribed velocities and elevation
    "waves": (1, 1, 1),      # incident waves
    "free": (4, 4, 4),       # free (Thompson) outflow
}
ALIASES = {
    "closed": "wall", "solid": "wall", "222": "wall",
    "discharge": "q", "flow": "q", "inflow": "q", "455": "q",
    "elevation": "h", "stage": "h", "level": "h", "outflow": "h", "544": "h",
    "555": "qh", "466": "uv", "566": "uvh", "111": "waves", "444": "free",
}


def parse_type(text) -> Tuple[int, int, int]:
    t = str(text or "").strip().lower().replace(" ", "").replace("_", "")
    t = ALIASES.get(t, t)
    if t in BOUNDARY_CODES:
        return BOUNDARY_CODES[t]
    digits = [c for c in t if c.isdigit()]
    if len(digits) == 3:
        return tuple(int(c) for c in digits)
    raise ValueError("unknown boundary type %r (use wall, Q, H, QH, UV, UVH, waves, free or a code like 455)" % text)


@dataclass
class BoundaryLine:
    coords: np.ndarray
    codes: Tuple[int, int, int]
    name: str = ""


def boundary_order(x, y, ikle) -> np.ndarray:
    """Boundary nodes in TELEMAC order (outer loop counter-clockwise first)."""
    loops = boundary_loops(ikle, x, y)
    return np.concatenate(loops) if loops else np.zeros(0, dtype=np.int64)


def attach_mesh(table: BoundaryTable, x, y, ikle) -> BoundaryTable:
    """Give a table read from a .cli the coordinates and loops of its mesh, so its liquid boundaries are
    numbered as TELEMAC numbers them (see telemac_io.cli.liquid_boundaries)."""
    nodes = table.nodes
    loops = boundary_loops(ikle, x, y)
    which = {}
    for k, lp in enumerate(loops):
        for nd in lp:
            which[int(nd)] = k
    lid = np.array([which.get(int(nd), -1) for nd in nodes])
    ranges = []
    a = 0
    for i in range(1, len(lid) + 1):
        if i == len(lid) or lid[i] != lid[a]:
            ranges.append((a, i))
            a = i
    table.loops = ranges
    table.xy = (np.asarray(x, dtype=np.float64)[nodes], np.asarray(y, dtype=np.float64)[nodes])
    return table


def build_boundary_table(x, y, ikle, lines: Sequence[BoundaryLine], tolerance: Optional[float] = None):
    """Return (BoundaryTable, report lines)."""
    loops = boundary_loops(ikle, x, y)
    order = np.concatenate(loops) if loops else np.zeros(0, dtype=np.int64)
    table = BoundaryTable.walls(order)
    sizes = np.cumsum([0] + [len(lp) for lp in loops])
    table.loops = [(int(sizes[i]), int(sizes[i + 1])) for i in range(len(loops))]
    table.xy = (np.asarray(x, dtype=np.float64)[order], np.asarray(y, dtype=np.float64)[order])
    report: List[str] = []
    if len(order) == 0:
        return table, ["The mesh has no boundary."]
    bx, by = x[order], y[order]
    if tolerance is None:
        # half of the median boundary edge length
        e = np.hypot(np.diff(np.r_[bx, bx[:1]]), np.diff(np.r_[by, by[:1]]))
        tolerance = 0.5 * float(np.median(e))
    r = table.records
    assigned = np.zeros(len(order), dtype=bool)
    for ln in lines:
        c = np.asarray(ln.coords, float)
        if len(c) < 2:
            continue
        d = distance_to_segments(bx, by, c[:-1], c[1:])
        hit = d <= tolerance
        if not hit.any():
            report.append("Line %r touches no boundary node (tolerance %.3g); check that it follows the mesh boundary." % (ln.name, tolerance))
            continue
        r["lihbor"][hit], r["liubor"][hit], r["livbor"][hit] = ln.codes
        for k in np.flatnonzero(hit):
            table.comments[k] = ln.name
        assigned |= hit
    lb = liquid_boundaries(table)
    for b in lb:
        n = len(b["nodes"])
        name = table.comments[b["start"]] or "unnamed"
        report.append("Liquid boundary %d: %s, codes %d %d %d, %d nodes (boundary nodes %d to %d)" % (
            b["number"], name, *b["codes"], n, b["start"] + 1, b["end"] + 1))
    if len(lb) > 1:
        report.append("Liquid boundaries are numbered as TELEMAC numbers them: from the boundary point with the "
                      "smallest x + y, counter-clockwise.")
        if n < 5:
            report.append("  Warning: liquid boundary %d has only %d nodes; 5 or more are recommended." % (b["number"], n))
    if not lb:
        report.append("No liquid boundary: every boundary node is a wall.")
    return table, report


def liq_template(table: BoundaryTable, times: Sequence[float] = (0.0, 3600.0), values: Optional[dict] = None) -> str:
    """A liquid boundaries file (.liq) skeleton with one column per prescribed quantity."""
    values = values or {}
    cols = []
    for b in liquid_boundaries(table):
        h, u, _v = b["codes"]
        n = b["number"]
        if u == 5:
            cols.append(("Q(%d)" % n, "m3/s"))
        if h == 5:
            cols.append(("SL(%d)" % n, "m"))
        if u == 6:
            cols.append(("U(%d)" % n, "m/s"))
            cols.append(("V(%d)" % n, "m/s"))
    lines = ["# Liquid boundaries file written by TELEMAC Toolkit; replace the values with your hydrographs",
             "T " + " ".join(c for c, _ in cols), "s " + " ".join(u for _, u in cols)]
    for t in times:
        lines.append("%g " % t + " ".join("%g" % values.get(c, 0.0) for c, _ in cols))
    return "\n".join(lines) + "\n"
