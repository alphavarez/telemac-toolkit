"""Gmsh as the interior mesher (optional: only used when the ``gmsh`` Python package is installed).

The plugin's own generator still does everything around the interior: domain rings, islands and
structured parts, constraint lines intersected and resampled to the size field, fixed points, and
afterwards the stitching of structured parts, the boundary conditions and the output formats. Gmsh
receives the resulting constraint nodes and segments as fixed geometry (every segment meshed with its
two end nodes only, so no node is added on a boundary, a breakline or the edge of a structured part)
and fills the inside with one of its algorithms, following the same size field.

Algorithms (Gmsh's ``Mesh.Algorithm``):

* ``delaunay``  (5) — fastest on large plane meshes, handles complex size fields well;
* ``frontal``   (6) — Frontal-Delaunay, the highest element quality;
* ``meshadapt`` (1) — local edge swaps, splits and collapses, the most robust;
* ``quads``     (8) — Frontal-Delaunay for quadrilaterals, then recombined into quadrilaterals (for
  programs that read them; TELEMAC gets each quadrilateral as two triangles).

Gmsh is GPL-2.0-or-later, compatible with this plugin's GPL-3.0-or-later. Install it into QGIS's
Python with ``python -m pip install gmsh`` (on Windows, from the OSGeo4W shell).
"""

from __future__ import annotations

from typing import Callable, Dict, List, Tuple

import numpy as np

__all__ = ["available", "version", "ALGORITHMS", "gmsh_interior"]

ALGORITHMS = {"delaunay": 5, "frontal": 6, "meshadapt": 1, "quads": 8}


def available() -> bool:
    try:
        import gmsh  # noqa: F401

        return True
    except (ImportError, OSError):  # OSError: a system library the wheel needs is missing
        return False


def version() -> str:
    try:
        import gmsh

        return str(gmsh.__version__)
    except (ImportError, OSError, AttributeError):
        return ""


def _loops_from_segments(seg: np.ndarray) -> List[List[int]]:
    """Closed loops (node sequences) from boundary segments, each node used by two of them."""
    nb: Dict[int, List[int]] = {}
    for a, b in seg.tolist():
        nb.setdefault(a, []).append(b)
        nb.setdefault(b, []).append(a)
    seen = set()
    loops = []
    for start in nb:
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            nxt = [n for n in nb[cur] if n != prev and (n not in seen or (n == start and len(loop) > 2))]
            if not nxt:
                break
            n = nxt[0]
            if n == start:
                break
            loop.append(n)
            seen.add(n)
            prev, cur = cur, n
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def _split_quads(xy: np.ndarray, qd: np.ndarray) -> np.ndarray:
    """Two triangles per quadrilateral, cut along the diagonal that gives the larger smallest angle."""
    from .engine import min_angles

    a = np.vstack([qd[:, [0, 1, 2]], qd[:, [0, 2, 3]]])
    b = np.vstack([qd[:, [0, 1, 3]], qd[:, [1, 2, 3]]])
    n = len(qd)
    ma = min_angles(xy[:, 0], xy[:, 1], a).reshape(2, n).min(axis=0)
    mb = min_angles(xy[:, 0], xy[:, 1], b).reshape(2, n).min(axis=0)
    use_b = mb > ma
    t1 = np.where(use_b[:, None], qd[:, [0, 1, 3]], qd[:, [0, 1, 2]])
    t2 = np.where(use_b[:, None], qd[:, [1, 2, 3]], qd[:, [0, 2, 3]])
    return np.vstack([t1, t2])


def _ring_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _grid_lookup(size, origin, loc, cons_h, max_cells: int = 4_000_000):
    """A fast size callback for Gmsh: the size field sampled once on a grid, read bilinearly.

    Gmsh asks for the size at every point it considers, often millions of times; evaluating the size
    field there each time is too slow in Python. Sizes grow linearly away from their sources, so the
    grid's error is at most the growth rate times a cell diagonal. The cell is the smallest size along
    the constraints (or a fraction of the extent, for very large domains)."""
    lo = loc.min(axis=0)
    hi = loc.max(axis=0)
    span = np.maximum(hi - lo, 1e-9)
    cell = max(float(np.min(cons_h)) if len(cons_h) else span.max() / 100.0, float(np.sqrt(span[0] * span[1] / max_cells)))
    nx = int(np.ceil(span[0] / cell)) + 2
    ny = int(np.ceil(span[1] / cell)) + 2
    gx = lo[0] + (np.arange(nx) - 0.5) * cell
    gy = lo[1] + (np.arange(ny) - 0.5) * cell
    GX, GY = np.meshgrid(gx, gy)
    grid = np.asarray(size(GX.ravel() + origin[0], GY.ravel() + origin[1]), dtype=np.float64).reshape(ny, nx)
    rows = grid.tolist()
    x0, y0 = float(gx[0]), float(gy[0])
    inv = 1.0 / cell

    def lookup(dim, tag, x, y, z, lc):
        fx = (x - x0) * inv
        fy = (y - y0) * inv
        i = int(fx)
        j = int(fy)
        if i < 0:
            i, fx = 0, 0.0
        elif i > nx - 2:
            i, fx = nx - 2, float(nx - 1)
        if j < 0:
            j, fy = 0, 0.0
        elif j > ny - 2:
            j, fy = ny - 2, float(ny - 1)
        tx = fx - i
        ty = fy - j
        r0 = rows[j]
        r1 = rows[j + 1]
        return (r0[i] * (1 - tx) + r0[i + 1] * tx) * (1 - ty) + (r1[i] * (1 - tx) + r1[i + 1] * tx) * ty

    return lookup


def gmsh_interior(algorithm: str = "frontal", optimize: bool = True, threads: int = 0) -> Callable:
    """An interior mesher for ``engine.generate(..., interior=...)``.

    The returned function takes (cons_xy, cons_type, cons_h, seg_pairs, ring_mask, size, step) and returns
    (xy, ntype, dens, ikle, quads, log): the constraint nodes first, unchanged and in order, then the
    nodes Gmsh added; triangles; quadrilaterals (empty unless ``algorithm == "quads"``).
    """
    if algorithm not in ALGORITHMS:
        raise ValueError("unknown Gmsh algorithm %r (use %s)" % (algorithm, ", ".join(ALGORITHMS)))

    def run(cons_xy, cons_type, cons_h, seg_pairs, ring_mask, size, step):
        import gmsh

        from .engine import INTERIOR
        from .geometry import points_in_rings

        log: List[str] = []
        xy0 = np.asarray(cons_xy, dtype=np.float64)
        sp = np.asarray(seg_pairs, dtype=np.int64).reshape(-1, 2)
        origin = xy0.min(axis=0)
        loc = xy0 - origin  # Gmsh works in local coordinates (UTM values lose digits in its tolerances)
        # segments of the domain's rings (outline, islands, edges of structured parts) as the generator marked them
        on_ring = np.asarray(ring_mask, dtype=bool)
        mid = 0.5 * (loc[sp[:, 0]] + loc[sp[:, 1]])
        ext = float(np.max(np.ptp(loc, axis=0))) or 1.0
        loops = _loops_from_segments(sp[on_ring])
        if not loops:
            raise ValueError("Gmsh: no closed outline found among the constraint segments")
        # outer rings and holes by nesting depth
        polys = [loc[lp] for lp in loops]
        areas = np.array([abs(_ring_area(p)) for p in polys])
        depth = np.zeros(len(loops), dtype=int)
        parent = -np.ones(len(loops), dtype=int)
        for i, p in enumerate(polys):
            probe = p[0] + 1e-9 * ext  # a vertex of loop i, inside any loop that contains loop i
            best = None
            for j, q in enumerate(polys):
                if j == i or areas[j] <= areas[i]:
                    continue
                if points_in_rings(np.array([probe[0]]), np.array([probe[1]]), [q])[0]:
                    depth[i] += 1
                    if best is None or areas[j] < areas[best]:
                        best = j
            parent[i] = -1 if best is None else best

        gmsh.initialize(interruptible=False)
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.NumThreads", int(threads))
            gmsh.model.add("telemac_toolkit")
            geo = gmsh.model.geo
            n = len(loc)
            for i in range(n):
                geo.addPoint(float(loc[i, 0]), float(loc[i, 1]), 0.0, float(cons_h[i]), i + 1)
            line_of: Dict[Tuple[int, int], int] = {}
            for k, (a, b) in enumerate(sp.tolist()):
                tag = geo.addLine(a + 1, b + 1)
                line_of[(a, b)] = tag
                line_of[(b, a)] = -tag
            loop_tags = []
            for lp in loops:
                cl = [line_of[(lp[i], lp[(i + 1) % len(lp)])] for i in range(len(lp))]
                loop_tags.append(geo.addCurveLoop(cl))
            surf_of = {}
            for i in range(len(loops)):
                if depth[i] % 2 == 0:  # an outer ring: its holes are the rings directly inside it
                    holes = [loop_tags[j] for j in range(len(loops)) if parent[j] == i and depth[j] % 2 == 1]
                    surf_of[i] = geo.addPlaneSurface([loop_tags[i]] + holes)
            geo.synchronize()
            # every segment keeps exactly its two end nodes
            for tag in set(abs(t) for t in line_of.values()):
                gmsh.model.mesh.setTransfiniteCurve(tag, 2)
            # inner segments (breaklines, guide lines) and lone fixed points embedded in their surface
            outer = [i for i in surf_of]

            def surface_at(p):
                for i in outer:
                    if points_in_rings(np.array([p[0]]), np.array([p[1]]), [polys[i]])[0]:
                        inner_holes = [polys[j] for j in range(len(loops)) if parent[j] == i and depth[j] % 2 == 1]
                        if not inner_holes or not points_in_rings(np.array([p[0]]), np.array([p[1]]), inner_holes)[0]:
                            return surf_of[i]
                return None

            emb: Dict[int, List[int]] = {}
            for k in np.flatnonzero(~on_ring):
                a, b = sp[k]
                s = surface_at(mid[k])
                if s is not None:
                    emb.setdefault(s, []).append(abs(line_of[(int(a), int(b))]))
            for s, curves in emb.items():
                gmsh.model.mesh.embed(1, curves, 2, s)
            used = np.zeros(n, dtype=bool)
            used[sp.ravel()] = True
            lone = np.flatnonzero(~used)
            pts: Dict[int, List[int]] = {}
            for i in lone:
                s = surface_at(loc[i])
                if s is not None:
                    pts.setdefault(s, []).append(int(i) + 1)
            for s, p in pts.items():
                gmsh.model.mesh.embed(0, p, 2, s)
            # size: the plugin's size field, nothing else
            for opt, val in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                             ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", ALGORITHMS[algorithm]),
                             ("Mesh.Optimize", 1 if optimize else 0), ("Mesh.Smoothing", 2)):
                gmsh.option.setNumber(opt, val)
            gmsh.model.mesh.setSizeCallback(_grid_lookup(size, origin, loc, cons_h))
            if algorithm == "quads":
                gmsh.option.setNumber("Mesh.RecombineAll", 1)
                gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
            if step:
                step(0.2, "Gmsh: meshing (%s)" % algorithm)
            gmsh.model.mesh.generate(2)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            gxy = coords.reshape(-1, 3)[:, :2]
            tri_n = []
            quads = []
            types, etags, enodes = gmsh.model.mesh.getElements(2)
            for et, nodes in zip(types, enodes):
                if et == 2:
                    tri_n.append(np.asarray(nodes, dtype=np.int64).reshape(-1, 3))
                elif et == 3:
                    quads.append(np.asarray(nodes, dtype=np.int64).reshape(-1, 4))
                else:
                    raise ValueError("Gmsh returned elements of type %d (only triangles and quadrilaterals are used)" % et)
        finally:
            gmsh.finalize()
        # Gmsh node tags -> indices: constraint nodes first (matched by position), then the new ones
        tag_index = {}
        key = {(round(float(x) / (1e-9 * ext)), round(float(y) / (1e-9 * ext))): i for i, (x, y) in enumerate(loc)}
        new_xy: List[np.ndarray] = []
        for t, p in zip(tags.tolist(), gxy):
            k = (round(float(p[0]) / (1e-9 * ext)), round(float(p[1]) / (1e-9 * ext)))
            if k in key:
                tag_index[int(t)] = key[k]
            else:
                tag_index[int(t)] = n + len(new_xy)
                new_xy.append(p)
        lookup = np.zeros(int(max(tags)) + 1, dtype=np.int64)
        for t, i in tag_index.items():
            lookup[t] = i
        ikle = lookup[np.vstack(tri_n)] if tri_n else np.zeros((0, 3), dtype=np.int64)
        qd = lookup[np.vstack(quads)] if quads else np.zeros((0, 4), dtype=np.int64)
        if len(qd):  # TELEMAC reads triangles: each quadrilateral also as two triangles, cut along its better diagonal
            ikle = np.vstack([ikle, _split_quads(np.vstack([loc, np.array(new_xy).reshape(-1, 2)]), qd)])
        add = np.array(new_xy).reshape(-1, 2) + origin
        xy = np.vstack([xy0, add])
        ntype = np.concatenate([np.asarray(cons_type, dtype=np.int64), np.full(len(add), INTERIOR, dtype=np.int64)])
        dens = np.concatenate([np.asarray(cons_h, dtype=np.float64), size(add[:, 0], add[:, 1]) if len(add) else np.zeros(0)])
        log.append("Gmsh %s, algorithm %s: %d nodes added inside, %d triangles%s"
                   % (version(), algorithm, len(add), len(ikle) - 2 * len(qd),
                      (", %d quadrilaterals" % len(qd)) if len(qd) else ""))
        return xy, ntype, dens, ikle, qd, log

    return run
