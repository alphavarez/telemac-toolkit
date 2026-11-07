"""Editing existing meshes: value transfer, boundary rebuild, re-meshing a region."""

import numpy as np
import pytest

from telemac_toolkit.mesher.backends import triangulate_scipy
from telemac_toolkit.mesher.boundary import BoundaryLine, attach_mesh, build_boundary_table, parse_type
from telemac_toolkit.mesher.edit import match_nodes, rebuild_boundary, region_elements, remesh_region, transfer
from telemac_toolkit.mesher.engine import MeshInput, generate, min_angles
from telemac_toolkit.mesher.sizefield import SizeField
from telemac_toolkit.telemac_io.cli import liquid_boundaries
from telemac_toolkit.telemac_io.mesh import boundary_loops


@pytest.fixture(scope="module")
def square():
    ring = np.array([[0, 0], [200, 0], [200, 100], [0, 100], [0, 0]], float)
    r = generate(MeshInput([[ring]], SizeField(10.0, 1.2)), triangulate_scipy)
    lines = [BoundaryLine(np.array([[0, 20], [0, 80]]), parse_type("Q"), "in"),
             BoundaryLine(np.array([[200, 20], [200, 80]]), parse_type("H"), "out")]
    table, _ = build_boundary_table(r.x, r.y, r.ikle, lines)
    return r.x, r.y, r.ikle, table


def _edges(ikle):
    e = np.sort(np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return u, c


def test_remesh_region_keeps_the_boundary_and_conforms(square):
    x, y, ik, table = square
    region = [np.array([[60, 30], [140, 30], [140, 70], [60, 70]], float)]
    inside = region_elements(x, y, ik, region)
    nx, ny, nik, old_index, rep = remesh_region(x, y, ik, inside, 3.0, triangulate_scipy)
    assert len(nik) > len(ik)  # finer inside
    _u, c = _edges(nik)
    assert c.max() == 2  # no edge shared by three elements
    # the outer boundary is the same set of points (no hanging nodes: boundary edges only on the outline)
    lp_old = boundary_loops(ik, x, y)
    lp_new = boundary_loops(nik, nx, ny)
    assert len(lp_new) == 1 and len(lp_new[0]) == len(lp_old[0])
    assert np.sum(c == 1) == len(lp_old[0])
    assert min_angles(nx, ny, nik).min() > 15
    # old nodes kept where they are
    kept = old_index >= 0
    assert np.allclose(nx[kept], x[old_index[kept]]) and np.allclose(ny[kept], y[old_index[kept]])
    assert any("replaced" in m for m in rep)
    with pytest.raises(ValueError, match="no element"):
        remesh_region(x, y, ik, np.zeros(len(ik), bool), 3.0, triangulate_scipy)


def test_transfer_and_boundary_after_remeshing(square):
    x, y, ik, table = square
    z = 5.0 + 0.01 * x - 0.02 * y  # linear: interpolation is exact
    region = [np.array([[-1, 30], [40, 30], [40, 70], [-1, 70]], float)]  # touches the inflow boundary
    inside = region_elements(x, y, ik, region)
    nx, ny, nik, old_index, _rep = remesh_region(x, y, ik, inside, 4.0, triangulate_scipy)
    m = match_nodes(x, y, nx, ny, 1e-6)
    kept = old_index >= 0
    assert np.array_equal(m[kept], old_index[kept])  # kept nodes map back (new ones may fall on removed ones)
    vals, copied, interp, near = transfer(x, y, ik, z[None, :], nx, ny, m)
    assert copied >= int(kept.sum()) and near == 0 and interp == len(nx) - copied
    assert np.allclose(vals[0], 5.0 + 0.01 * nx - 0.02 * ny)
    attach_mesh(table, x, y, ik)
    t2, rep = rebuild_boundary(table, x, y, nx, ny, nik, m)
    before = {table.comments[b["start"]]: len(b["nodes"]) for b in liquid_boundaries(table)}
    after = {t2.comments[b["start"]]: len(b["nodes"]) for b in liquid_boundaries(t2)}
    # where the region touches the mesh boundary, nodes can be added there: they take its conditions
    assert after["out"] == before["out"] and after["in"] >= before["in"]
    yin = sorted(ny[b["nodes"]].tolist() for b in liquid_boundaries(t2) if t2.comments[b["start"]] == "in")[0]
    assert min(yin) == pytest.approx(20) and max(yin) == pytest.approx(80)  # the inflow still spans the same stretch
    assert "kept their conditions" in rep[0]


def test_rebuild_boundary_for_new_boundary_nodes(square):
    x, y, ik, table = square
    attach_mesh(table, x, y, ik)
    # an "edited" mesh: the same outline with an extra node in the middle of every boundary edge
    lp = boundary_loops(ik, x, y)[0]
    mids = np.stack([(x[lp] + x[np.roll(lp, -1)]) / 2, (y[lp] + y[np.roll(lp, -1)]) / 2], axis=1)
    pts = np.vstack([np.stack([x, y], 1), mids])
    from scipy.spatial import Delaunay

    tri = Delaunay(pts).simplices
    ex, ey = pts[:, 0], pts[:, 1]
    m = match_nodes(x, y, ex, ey, 1e-6)
    t2, _rep = rebuild_boundary(table, x, y, ex, ey, tri, m)
    before = {table.comments[b["start"]]: len(b["nodes"]) for b in liquid_boundaries(table)}
    after = {t2.comments[b["start"]]: len(b["nodes"]) for b in liquid_boundaries(t2)}
    assert set(after) == set(before)
    for k in before:
        assert after[k] == 2 * before[k] - 1  # new nodes between liquid nodes are liquid


@pytest.mark.parametrize("region", [
    [[150, -10], [210, -10], [210, 110], [150, 110]],   # across the whole right end, on the boundary
    [[170, 20], [210, 20], [210, 60], [170, 60]],       # touching the boundary on one side only
    [[-10, -10], [40, -10], [40, 40], [-10, 40]],       # a corner of the mesh
])
def test_remesh_region_on_the_mesh_boundary(square, region):
    """A region reaching the mesh boundary: no duplicate or unused nodes, one boundary loop, conditions kept."""
    x, y, ik, table = square
    inside = region_elements(x, y, ik, [np.array(region, float)])
    nx, ny, nik, old_index, _rep = remesh_region(x, y, ik, inside, 4.0, triangulate_scipy)
    used = np.zeros(len(nx), bool)
    used[nik.ravel()] = True
    assert used.all()
    keys = np.round(np.stack([nx, ny], 1), 6)
    assert len(np.unique(keys, axis=0)) == len(nx)  # no two nodes at one place
    _u, c = _edges(nik)
    assert c.max() == 2
    loops = boundary_loops(nik, nx, ny)
    assert len(loops) == 1
    # boundary conditions carried over: same liquid boundaries
    attach_mesh(table, x, y, ik)
    match = match_nodes(x, y, nx, ny, 1e-6)
    t2, _r = rebuild_boundary(table, x, y, nx, ny, nik, match)
    names = sorted(t2.comments[b["start"]] for b in liquid_boundaries(t2))
    assert names == ["in", "out"]


def test_remesh_two_areas_far_from_the_origin():
    """Map coordinates (125 km, 188 km): the outline of a re-meshed area kept its last node, and the new
    elements join the old ones (the ring test once took two nodes 1 m apart for one, and the patch overlapped
    its surroundings)."""
    from telemac_toolkit.mesher import PolygonDensity
    from telemac_toolkit.telemac_io.mesh import signed_areas

    X0, Y0 = 125000.0, 188500.0
    dom = np.array([[0, 0], [80, 0], [80, 60], [0, 60]], float) + [X0, Y0]
    a = np.array([[10, 10], [35, 12], [30, 30], [12, 25]], float) + [X0, Y0]
    b = np.array([[45, 30], [70, 35], [65, 55], [50, 50]], float) + [X0, Y0]
    sf = SizeField(4.0, 1.2, polygons=[PolygonDensity([a], 0.4), PolygonDensity([b], 0.8)])
    res = generate(MeshInput(domain=[[dom]], size=sf), triangulate_scipy)
    x, y, ik = res.x.astype(np.float32).astype(float), res.y.astype(np.float32).astype(float), res.ikle  # as read from a SELAFIN
    area0 = np.abs(signed_areas(x, y, ik)).sum()
    for poly in (a, b):
        x, y, ik, _old, _rep = remesh_region(x, y, ik, region_elements(x, y, ik, [poly]), 1.5, triangulate_scipy)
        loops = boundary_loops(ik, x, y)
        assert len(loops) == 1 and len(loops[0]) == len(boundary_loops(res.ikle, res.x, res.y)[0])
        assert abs(np.abs(signed_areas(x, y, ik)).sum() - area0) < 1e-6 * area0


def test_as_ring_keeps_close_vertices_far_from_the_origin():
    from telemac_toolkit.mesher.geometry import as_ring, same_point

    r = np.array([[500000.0, 2000000.0], [500100.0, 2000000.0], [500100.0, 2000100.0], [500000.5, 2000000.5]])
    assert len(as_ring(r)) == 4  # the last vertex is 0.7 m from the first: a vertex, not the closing one
    assert len(as_ring(np.vstack([r, r[:1]]))) == 4
    assert same_point(r[0], r[0] + 1e-9) and not same_point(r[0], r[0] + 0.01)
