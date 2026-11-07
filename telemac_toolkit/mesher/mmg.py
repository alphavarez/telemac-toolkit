"""MMG (mmg2d) as an external program: optimise or adapt an existing triangular mesh.

MMG is LGPL-licensed and distributed as standalone programs (``mmg2d_O3``) for Windows, macOS and Linux
(https://www.mmgtools.org, or the ``mmgsuite`` package of conda-forge). The plugin does not import it: it
writes the mesh and the wanted sizes as MMG's text files (Medit ``.mesh`` and ``.sol``), runs the program and
reads the result back, so nothing has to be installed into QGIS's Python.

What is kept:

* every boundary edge is a *required edge*: the outline, the islands and their nodes do not move, so the
  boundary conditions carry over node for node;
* elements in ``keep`` (for example structured channel cells) are *required triangles*, left untouched.

Sizes: one value per node (``sizes``) for adaptation, or none with ``optimise=True`` (MMG's ``-optim``:
better shapes, the same sizes). ``hmin``, ``hmax``, ``hgrad`` and ``hausd`` are MMG's own bounds.
"""

from __future__ import annotations

import os
import shutil
# subprocess runs only MMG's mmg2d program (run_mmg2d)
import subprocess  # nosec B404
import tempfile
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..telemac_io.mesh import boundary_loops

__all__ = ["find_mmg2d", "write_medit", "read_medit", "run_mmg2d"]

_NAMES = ("mmg2d_O3", "mmg2d", "mmg2d_O3.exe", "mmg2d.exe")


def find_mmg2d(explicit: Optional[str] = None) -> Optional[str]:
    """Path of the mmg2d program: ``explicit``, the MMG2D environment variable, or the PATH."""
    for cand in (explicit, os.environ.get("MMG2D")):
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    for n in _NAMES:
        p = shutil.which(n)
        if p:
            return p
    return None


def write_medit(path: str, x, y, ikle, required_edges: Optional[np.ndarray] = None,
                required_triangles: Optional[np.ndarray] = None, sizes: Optional[np.ndarray] = None) -> None:
    """Medit .mesh (and .sol next to it when ``sizes`` are given), nodes and elements numbered from 1."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ik = np.asarray(ikle, dtype=np.int64)
    with open(path, "w") as fh:
        fh.write("MeshVersionFormatted 2\nDimension 2\nVertices\n%d\n" % len(x))
        fh.write("".join("%.17g %.17g 0\n" % (a, b) for a, b in zip(x, y)))
        fh.write("Triangles\n%d\n" % len(ik))
        fh.write("".join("%d %d %d 0\n" % (a + 1, b + 1, c + 1) for a, b, c in ik))
        if required_edges is not None and len(required_edges):
            e = np.asarray(required_edges, dtype=np.int64)
            fh.write("Edges\n%d\n" % len(e))
            fh.write("".join("%d %d 1\n" % (a + 1, b + 1) for a, b in e))
            fh.write("RequiredEdges\n%d\n" % len(e))
            fh.write("".join("%d\n" % (k + 1) for k in range(len(e))))
        if required_triangles is not None and len(required_triangles):
            t = np.flatnonzero(np.asarray(required_triangles, dtype=bool))
            fh.write("RequiredTriangles\n%d\n" % len(t))
            fh.write("".join("%d\n" % (k + 1) for k in t))
        fh.write("End\n")
    if sizes is not None:
        s = np.asarray(sizes, dtype=np.float64)
        with open(os.path.splitext(path)[0] + ".sol", "w") as fh:
            fh.write("MeshVersionFormatted 2\nDimension 2\nSolAtVertices\n%d\n1 1\n" % len(s))
            fh.write("".join("%.10g\n" % v for v in s))
            fh.write("End\n")


def read_medit(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(x, y, ikle) of a Medit .mesh file (2D, triangles)."""
    with open(path) as fh:
        tok = fh.read().split()
    k = 0
    x = y = ik = None
    while k < len(tok):
        t = tok[k]
        if t == "Vertices":
            n = int(tok[k + 1])
            a = np.array(tok[k + 2:k + 2 + 3 * n], dtype=np.float64).reshape(n, 3)
            x, y = a[:, 0], a[:, 1]
            k += 2 + 3 * n
        elif t == "Triangles":
            n = int(tok[k + 1])
            ik = np.array(tok[k + 2:k + 2 + 4 * n], dtype=np.int64).reshape(n, 4)[:, :3] - 1
            k += 2 + 4 * n
        else:
            k += 1
    if x is None or ik is None:
        raise ValueError("%s: no vertices or triangles" % path)
    return x, y, ik


def run_mmg2d(x, y, ikle, sizes: Optional[np.ndarray] = None, optimise: bool = False, keep: Optional[np.ndarray] = None,
              hmin: Optional[float] = None, hmax: Optional[float] = None, hgrad: Optional[float] = 1.3,
              hausd: Optional[float] = None, exe: Optional[str] = None, timeout: float = 3600.0,
              extra: Sequence[str] = ()) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    """Run mmg2d on the mesh; returns (x, y, ikle, log lines). The boundary nodes keep their coordinates."""
    prog = find_mmg2d(exe)
    if prog is not None and not (os.path.isfile(prog) and os.access(prog, os.X_OK)):
        raise ValueError("%s is not a program that can be run (give the path of mmg2d_O3)" % prog)
    if prog is None:
        raise ValueError("mmg2d was not found: install MMG (https://www.mmgtools.org, or conda-forge 'mmgsuite') and give "
                         "the path of mmg2d_O3, or put it on the PATH")
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ik = np.asarray(ikle, dtype=np.int64)
    ox, oy = float(x.min()), float(y.min())  # local coordinates: MMG's tolerances are relative to the extent
    loops = boundary_loops(ik, x, y)
    edges = np.vstack([np.stack([lp, np.roll(lp, -1)], axis=1) for lp in loops]) if loops else np.zeros((0, 2), dtype=np.int64)
    work = tempfile.mkdtemp(prefix="tt_mmg_")
    try:
        src = os.path.join(work, "in.mesh")
        dst = os.path.join(work, "out.mesh")
        write_medit(src, x - ox, y - oy, ik, edges, keep, None if optimise else sizes)
        cmd = [prog, "-in", src, "-out", dst, "-v", "1"]
        if optimise:
            cmd.append("-optim")
        for flag, val in (("-hmin", hmin), ("-hmax", hmax), ("-hgrad", hgrad), ("-hausd", hausd)):
            if val is not None and val > 0:
                cmd += [flag, "%.10g" % val]
        cmd += [str(a) for a in extra]
        # the user's own mmg2d (checked to be an existing program file), an argument list and no shell: nothing
        # from a drawing or a file name reaches a command line
        run = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=timeout, shell=False)  # nosec B603
        log = [line for line in (run.stdout + run.stderr).splitlines() if line.strip()]
        if run.returncode != 0 or not os.path.exists(dst):
            raise ValueError("mmg2d failed (exit %d): %s" % (run.returncode, " | ".join(log[-6:])))
        nx_, ny_, nik = read_medit(dst)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return nx_ + ox, ny_ + oy, nik, log
