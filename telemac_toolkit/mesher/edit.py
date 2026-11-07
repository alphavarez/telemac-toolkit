"""Editing existing meshes: carry values and boundary conditions over to an edited mesh, and
re-mesh a region inside an existing mesh.

``transfer`` maps the nodes of an edited mesh to the original one: nodes that did not move keep
their values exactly; new or moved nodes get values interpolated linearly in the original
triangles (nearest original node outside it). ``rebuild_boundary`` gives the edited mesh a
boundary table: boundary nodes that existed keep their record (types, values, name); new ones
take the type of the original boundary edge they lie on when both its ends share it, otherwise
they are walls.

``remesh_region`` replaces the elements whose centre lies in a region by a new mesh of that
region. Where the region borders old elements its outline keeps its nodes (those edges are fixed),
so the new elements join the old ones node for node. Where the region reaches the boundary of the
mesh, nodes can be added on the boundary; ``rebuild_boundary`` gives them its conditions.

Only NumPy is used; the triangulation is supplied by the caller, as for the generator.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..telemac_io.cli import BoundaryTable
from ..telemac_io.mesh import TriangleLocator, boundary_loops, signed_areas
from .geometry import points_in_rings

__all__ = ["match_nodes", "transfer", "rebuild_boundary", "region_elements", "remesh_region"]


def _keys(x, y, tol):
    return np.round(np.asarray(x) / tol).astype(np.int64), np.round(np.asarray(y) / tol).astype(np.int64)


def match_nodes(old_x, old_y, new_x, new_y, tol: float) -> np.ndarray:
    """For each new node, the index of the old node at the same place (within ``tol``), or -1."""
    kx, ky = _keys(old_x, old_y, tol)
    table: Dict[Tuple[int, int], int] = {}
    for i, k in enumerate(zip(kx.tolist(), ky.tolist())):
        table.setdefault(k, i)
    nx, ny = _keys(new_x, new_y, tol)
    out = np.full(len(nx), -1, dtype=np.int64)
    for i, (a, b) in enumerate(zip(nx.tolist(), ny.tolist())):
        j = table.get((a, b))
        if j is None:  # a key next door (rounding at a cell edge)
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    j = table.get((a + da, b + db))
                    if j is not None:
                        break
                if j is not None:
                    break
        if j is not None and np.hypot(old_x[j] - new_x[i], old_y[j] - new_y[i]) <= tol:
            out[i] = j
    return out


def transfer(old_x, old_y, old_ikle, values: np.ndarray, new_x, new_y, match: Optional[np.ndarray] = None, tol: float = 1e-6):
    """Values (shape (k, n_old)) at the new nodes: copied where a node matches, interpolated elsewhere.

    Returns (new values (k, n_new), number copied, number interpolated, number from the nearest node).
    """
    values = np.atleast_2d(np.asarray(values, dtype=np.float64))
    if match is None:
        match = match_nodes(old_x, old_y, new_x, new_y, tol)
    out = np.zeros((values.shape[0], len(new_x)))
    same = match >= 0
    out[:, same] = values[:, match[same]]
    rest = np.flatnonzero(~same)
    n_int = n_near = 0
    if rest.size:
        loc = TriangleLocator(np.asarray(old_x), np.asarray(old_y), np.asarray(old_ikle))
        tri, nodes, w = loc.locate_many(np.asarray(new_x)[rest], np.asarray(new_y)[rest])
        inside = tri >= 0
        for k in range(values.shape[0]):
            out[k, rest[inside]] = np.einsum("ij,ij->i", values[k][nodes[inside]], w[inside])
        n_int = int(inside.sum())
        for i in rest[~inside]:
            j = loc.nearest_node(float(new_x[i]), float(new_y[i]))
            out[:, i] = values[:, j]
            n_near += 1
    return out, int(same.sum()), n_int, n_near


def rebuild_boundary(old: BoundaryTable, old_x, old_y, new_x, new_y, new_ikle, match: np.ndarray):
    """Boundary table of the edited mesh, in TELEMAC order. Returns (table, report lines)."""
    old_x, old_y, new_x, new_y = (np.asarray(v, dtype=np.float64) for v in (old_x, old_y, new_x, new_y))
    loops = boundary_loops(new_ikle, new_x, new_y)
    order = np.concatenate(loops) if loops else np.zeros(0, dtype=np.int64)
    table = BoundaryTable.walls(order)
    sizes = np.cumsum([0] + [len(lp) for lp in loops])
    table.loops = [(int(sizes[i]), int(sizes[i + 1])) for i in range(len(loops))]
    table.xy = (new_x[order], new_y[order])
    rec_of_old = {int(n): k for k, n in enumerate(old.nodes.tolist())}
    # old boundary edges (consecutive records of a loop)
    ea, eb = [], []
    for a, b in (old.loops or [(0, len(old))]):
        idx = np.arange(a, b)
        ea.append(idx)
        eb.append(np.roll(idx, -1))
    ea = np.concatenate(ea) if ea else np.zeros(0, dtype=np.int64)
    eb = np.concatenate(eb) if eb else np.zeros(0, dtype=np.int64)
    on = old.nodes
    A = np.stack([old_x[on[ea]], old_y[on[ea]]], axis=1) if ea.size else np.zeros((0, 2))
    B = np.stack([old_x[on[eb]], old_y[on[eb]]], axis=1) if eb.size else np.zeros((0, 2))
    kept = new_rec = walls = 0
    fields = [f for f in table.records.dtype.names if f not in ("n", "k")]
    for k, node in enumerate(order.tolist()):
        j = int(match[node]) if node < len(match) else -1
        if j >= 0 and j in rec_of_old:
            r = rec_of_old[j]
            for f in fields:
                table.records[f][k] = old.records[f][r]
            table.comments[k] = old.comments[r] if r < len(old.comments) else ""
            kept += 1
            continue
        new_rec += 1
        if not ea.size:
            continue
        # the nearest old boundary edge
        ax, ay = new_x[node], new_y[node]
        dd = np.hypot(*(np.stack([B[:, 0] - A[:, 0], B[:, 1] - A[:, 1]], 1).T))
        t = np.clip(((ax - A[:, 0]) * (B[:, 0] - A[:, 0]) + (ay - A[:, 1]) * (B[:, 1] - A[:, 1])) / np.maximum(dd ** 2, 1e-300), 0, 1)
        dist = np.hypot(A[:, 0] + t * (B[:, 0] - A[:, 0]) - ax, A[:, 1] + t * (B[:, 1] - A[:, 1]) - ay)
        best = int(np.argmin(dist))
        ra, rb = int(ea[best]), int(eb[best])
        ca, cb = tuple(old.codes[ra]), tuple(old.codes[rb])
        if ca == cb and ca[0] != 2:
            for f in fields:
                table.records[f][k] = old.records[f][ra]
            table.comments[k] = old.comments[ra] if ra < len(old.comments) else ""
        else:
            walls += 1
    report = ["Boundary: %d nodes (%d before); %d kept their conditions, %d new (%d of them walls)"
              % (len(order), len(old), kept, new_rec, walls)]
    return table, report


def region_elements(x, y, ikle, rings: Sequence[np.ndarray]) -> np.ndarray:
    """Elements whose centre lies inside the region (rings: outer rings and holes alike, even-odd)."""
    cx = np.asarray(x)[ikle].mean(axis=1)
    cy = np.asarray(y)[ikle].mean(axis=1)
    return points_in_rings(cx, cy, rings)


def remesh_region(x, y, ikle, inside: np.ndarray, size: Optional[float], triangulate, growth: float = 1.2,
                  min_angle: Optional[float] = None, smoothing: int = 8):
    """Replace the elements flagged ``inside`` by a new mesh of the region they cover.

    ``size``: element size inside (None = the median edge length of the region's outline). The outline
    keeps its nodes, and the size grows from its edges at ``growth``. Returns (x, y, ikle, old_index,
    report): ``old_index`` gives for each node of the result the node of the input it is (-1 = new).
    """
    from .engine import MeshInput, Submesh, generate
    from .sizefield import SizeField

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ikle = np.asarray(ikle, dtype=np.int64)
    inside = np.asarray(inside, dtype=bool)
    if not inside.any():
        raise ValueError("no element has its centre inside the region")
    sub = ikle[inside]
    loops = boundary_loops(sub, x, y)
    if not loops:
        raise ValueError("the region has no outline")
    rings = [np.stack([x[lp], y[lp]], axis=1) for lp in loops]
    # outer rings and holes: a ring inside another ring of the region is a hole of it
    areas = [abs(float(np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1]))) / 2 for r in rings]
    order = np.argsort(areas)[::-1]
    polys: List[List[np.ndarray]] = []
    for i in order:
        r = rings[i]
        host = None
        for p in polys:
            if points_in_rings(r[:1, 0], r[:1, 1], [p[0]])[0]:
                host = p
                break
        if host is None:
            polys.append([r])
        else:
            host.append(r)
    edges = np.concatenate([np.hypot(*np.diff(np.vstack([r, r[:1]]), axis=0).T) for r in rings])
    h_out = float(np.median(edges))
    target = float(size) if size and size > 0 else h_out
    # the old elements around the region (sharing a node with its outline) are a fixed part: their
    # edges on the outline are never split, so the new elements join the old ones node for node
    ring_nodes = np.concatenate(loops)
    on_ring = np.zeros(len(x), dtype=bool)
    on_ring[ring_nodes] = True
    collar = (~inside) & on_ring[ikle].any(axis=1)
    cidx = np.unique(ikle[collar])
    cmap = -np.ones(len(x), dtype=np.int64)
    cmap[cidx] = np.arange(cidx.size)
    fixed = Submesh(x[cidx], y[cidx], cmap[ikle[collar]])
    sf = SizeField(target, growth)
    res = generate(MeshInput(polys, sf, resample_outline=False, smooth_iterations=smoothing, min_angle=min_angle,
                             max_nodes=None, submeshes=[fixed]), triangulate)
    if any("could not be recovered" in m for m in res.log):
        raise ValueError("the new elements could not be joined to the old ones; try a size closer to the outline's edges")
    # join: nodes of the result that are old nodes. Candidates are the collar's nodes and every node of the
    # region's outline: where the region reaches the mesh boundary its outline nodes touch no old element
    # outside it, and without them the new elements would get duplicates of those nodes
    tol = max(1e-9 * max(np.ptp(x), np.ptp(y)), 1e-9)
    cand = np.union1d(cidx, ring_nodes)
    m = match_nodes(x[cand], y[cand], res.x, res.y, tol)
    old_of_new = np.where(m >= 0, cand[np.maximum(m, 0)], -1)
    keep_old = np.ones(len(x), dtype=bool)
    interior_old = np.setdiff1d(np.unique(sub), ring_nodes)
    keep_old[interior_old] = False
    old_ids = np.flatnonzero(keep_old)
    new_index = -np.ones(len(x), dtype=np.int64)
    new_index[old_ids] = np.arange(old_ids.size)
    fresh = np.flatnonzero(old_of_new < 0)
    fresh_index = old_ids.size + np.arange(fresh.size)
    map_res = np.empty(len(res.x), dtype=np.int64)
    map_res[old_of_new >= 0] = new_index[old_of_new[old_of_new >= 0]]
    map_res[fresh] = fresh_index
    out_x = np.concatenate([x[old_ids], res.x[fresh]])
    out_y = np.concatenate([y[old_ids], res.y[fresh]])
    kept_el = new_index[ikle[~inside & ~collar]]
    new_el = map_res[res.ikle]
    out_ikle = np.vstack([kept_el, new_el])
    neg = signed_areas(out_x, out_y, out_ikle) < 0
    out_ikle[neg] = out_ikle[neg][:, [0, 2, 1]]
    old_index = np.concatenate([old_ids, -np.ones(fresh.size, dtype=np.int64)])
    # nodes no element uses any more (an outline node the new mesh did not keep): dropped
    used = np.zeros(len(out_x), dtype=bool)
    used[out_ikle.ravel()] = True
    if not used.all():
        remap = -np.ones(len(out_x), dtype=np.int64)
        remap[used] = np.arange(int(used.sum()))
        out_x, out_y, old_index = out_x[used], out_y[used], old_index[used]
        out_ikle = remap[out_ikle]
    # the result must cover the same area once: every boundary node has as many boundary edges in as out
    oriented = out_ikle
    e = np.concatenate([oriented[:, [0, 1]], oriented[:, [1, 2]], oriented[:, [2, 0]]])
    key = np.minimum(e[:, 0], e[:, 1]) * len(out_x) + np.maximum(e[:, 0], e[:, 1])
    _u, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    be = e[cnt[inv] == 1]
    unbalanced = np.flatnonzero(np.bincount(be[:, 0], minlength=len(out_x)) != np.bincount(be[:, 1], minlength=len(out_x)))
    over = np.flatnonzero(cnt > 2)
    if unbalanced.size or over.size:
        k = int(unbalanced[0]) if unbalanced.size else int(e[np.flatnonzero(cnt[inv] > 2)[0], 0])
        raise ValueError("re-meshing this area did not join the old elements cleanly near (%.3f, %.3f); nothing was "
                         "changed. Try a size closer to the outline's edges, or a slightly different polygon."
                         % (out_x[k], out_y[k]))
    report = ["Region: %d elements and %d interior nodes replaced by %d elements and %d nodes (size %.4g inside, %.4g on the outline)"
              % (int(inside.sum()), interior_old.size, len(res.ikle) - int(collar.sum()), fresh.size, target, h_out)]
    report += [m for m in res.log if m.startswith(("Mesh:", "Quality"))]
    return out_x, out_y, out_ikle, old_index, report
