"""Blue Kenue ASCII files: line sets (.i2s, .i3s) and triangular meshes (.t3s).

Each file starts with a header of lines beginning with ``:`` or ``#`` and ends it with ``:EndHeader``.

* ``.i2s`` / ``.i3s`` (line sets): for each line a record ``n value`` (the number of vertices and the
  line's value, used by Blue Kenue for edge lengths or levels), then n vertices ``x y`` (``x y z`` in .i3s).
* ``.t3s`` (T3 mesh): the header gives ``:NodeCount`` and ``:ElementCount``; then one node per line
  ``x y z`` (further columns are extra node values), then one element per line ``n1 n2 n3`` numbered from 1.

Only NumPy is used.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

__all__ = ["read_lines", "read_t3s", "write_i2s"]


def _body(path: str) -> Tuple[Dict[str, str], List[str]]:
    header: Dict[str, str] = {}
    body: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        in_header = True
        for raw in fh:
            line = raw.strip()
            if in_header:
                if line.lower().startswith(":endheader"):
                    in_header = False
                    continue
                if line.startswith(":"):
                    parts = line[1:].split(None, 1)
                    if parts:
                        header[parts[0].lower()] = parts[1] if len(parts) > 1 else ""
                    continue
                if line.startswith("#") or not line:
                    continue
                # no header at all: the data starts here
                in_header = False
            if line and not line.startswith("#"):
                body.append(line)
    return header, body


def read_lines(path: str) -> List[Tuple[np.ndarray, float]]:
    """Lines of an .i2s or .i3s file: [(vertices (n, 2) or (n, 3), value), ...]."""
    header, body = _body(path)
    out: List[Tuple[np.ndarray, float]] = []
    k = 0
    while k < len(body):
        head = body[k].split()
        k += 1
        if len(head) < 1:
            continue
        try:
            n = int(float(head[0]))
        except ValueError:
            raise ValueError("%s: expected a vertex count, found %r" % (path, body[k - 1]))
        value = float(head[1]) if len(head) > 1 else 0.0
        if k + n > len(body):
            raise ValueError("%s: a line announces %d vertices but the file ends first" % (path, n))
        pts = np.array([[float(v) for v in body[k + i].split()[:3]] for i in range(n)], dtype=np.float64)
        k += n
        if pts.ndim == 2 and len(pts) >= 2:
            out.append((pts, value))
    return out


def read_t3s(path: str):
    """(x, y, z, ikle) of a Blue Kenue T3 mesh; z is None when the nodes have no third column."""
    header, body = _body(path)
    try:
        nn = int(header.get("nodecount", "").split()[0])
        ne = int(header.get("elementcount", "").split()[0])
    except (IndexError, ValueError):
        raise ValueError("%s: the header gives no :NodeCount and :ElementCount" % path)
    if len(body) < nn + ne:
        raise ValueError("%s: %d nodes and %d elements announced, %d data lines found" % (path, nn, ne, len(body)))
    nodes = [row.split() for row in body[:nn]]
    x = np.array([float(r[0]) for r in nodes])
    y = np.array([float(r[1]) for r in nodes])
    z = np.array([float(r[2]) for r in nodes]) if all(len(r) > 2 for r in nodes) else None
    ikle = np.array([[int(v) for v in row.split()[:3]] for row in body[nn:nn + ne]], dtype=np.int64) - 1
    if ikle.min() < 0 or ikle.max() >= nn:
        raise ValueError("%s: element node numbers outside 1..%d" % (path, nn))
    return x, y, z, ikle


def write_i2s(path: str, lines, values=None) -> None:
    """Write lines [(n, 2) arrays] as an .i2s file (for checks and exchange with Blue Kenue)."""
    values = values if values is not None else [0.0] * len(lines)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(":FileType i2s  ASCII  EnSim 1.0\n:Application TELEMAC Toolkit for QGIS\n:EndHeader\n")
        for c, v in zip(lines, values):
            c = np.asarray(c, dtype=np.float64)
            fh.write("%d %s\n" % (len(c), repr(float(v))))
            for p in c:
                fh.write("%.6f %.6f\n" % (p[0], p[1]))
