"""Mesher core: T3 generator, channel mesher, regular grid, interpolation, boundary table."""

import numpy as np
import pytest

from telemac_toolkit.mesher import (
    MeshInput, PolygonDensity, SizeField, Submesh, channel_mesh, generate,
)
from telemac_toolkit.mesher.boundary import BoundaryLine, build_boundary_table, liq_template, parse_type
from telemac_toolkit.mesher.engine import min_angles
from telemac_toolkit.mesher.geometry import snap_rings_to_nodes
from telemac_toolkit.mesher.interp import RasterSource, bilinear, fill_nearest, interpolate_nodes
from telemac_toolkit.mesher.regular import regular_grid_mesh
from telemac_toolkit.telemac_io.mesh import boundary_loops, mesh_statistics, signed_areas

scipy = pytest.importorskip("scipy")
from telemac_toolkit.mesher.backends import triangulate_scipy  # noqa: E402

OUTER = np.array([[0, 0], [1000, 0], [1000, 600], [0, 600]], float)


def _edges(ikle):
    e = np.sort(np.concatenate([ikle[:, [0, 1]], ikle[:, [1, 2]], ikle[:, [2, 0]]]), axis=1)
    return set(map(tuple, e.tolist()))


def _node(x, y, p):
    d = np.hypot(x - p[0], y - p[1])
    k = int(np.argmin(d))
    assert d[k] < 1e-6, "node %s missing" % (p,)
    return k


def test_t3_quality_and_constraints():
    hole = np.array([[400, 250], [500, 250], [500, 330], [400, 330]], float)
    dens = PolygonDensity([np.array([[600, 50], [950, 50], [950, 250], [600, 250]], float)], 8.0)
    hard = np.array([[100, 100], [300, 400]], float)
    inp = MeshInput(domain=[[OUTER, hole]], size=SizeField(30.0, 1.2, polygons=[dens]),
                    soft_lines=[(np.array([[50, 500], [900, 450]], float), False)],
                    hard_lines=[(hard, False)], hard_points=np.array([[700.0, 400.0]]))
    res = generate(inp, triangulate_scipy)
    st = res.stats
    assert st["clockwise_elements"] == 0 and st["degenerate_elements"] == 0
    assert st["angle_min"] > 25.0
    assert st["angle_max"] < 125.0
    assert st["overconstrained_elements"] == 0
    assert abs(st["area_total"] - (600000.0 - 8000.0)) < 1e-6
    # two boundary loops: outline and the hole
    assert len(boundary_loops(res.ikle, res.x, res.y)) == 2
    # the hard line end points and the hard point are nodes; the hard line is a chain of edges
    a = _node(res.x, res.y, hard[0])
    b = _node(res.x, res.y, hard[1])
    _node(res.x, res.y, (700.0, 400.0))
    on = np.flatnonzero((np.abs((res.x - 100) * 300 - (res.y - 100) * 200) < 1e-6 * 360)
                        & (res.x >= 100 - 1e-9) & (res.x <= 300 + 1e-9))
    order = on[np.argsort(res.x[on])]
    assert order[0] == a and order[-1] == b
    E = _edges(res.ikle)
    assert all(tuple(sorted((int(p), int(q)))) in E for p, q in zip(order[:-1], order[1:]))
    # the density polygon is refined
    fine = (res.x > 650) & (res.x < 900) & (res.y > 80) & (res.y < 220)
    coarse = (res.x < 300) & (res.y > 450)
    assert fine.sum() > 8 * coarse.sum() / 3


def test_submesh_merge_and_grading():
    L = np.array([[100, 280], [400, 300], [700, 270], [900, 290]], float)
    cx, cy, cik = channel_mesh(L, L + [0, 60], L + [0, 30], cross_nodes=7, interval=12.0)
    gx, gy, gik = regular_grid_mesh(600, 420, 10, 10, 12, 8, angle_deg=15)
    subs = [Submesh(cx, cy, cik), Submesh(gx, gy, gik)]
    holes = []
    for s in subs:
        loop = boundary_loops(s.ikle, s.x, s.y)[0]
        holes.append(np.stack([s.x[loop], s.y[loop]], axis=1))
    outer = np.array([[0, 0], [1000, 0], [1000, 650], [0, 650]], float)
    res = generate(MeshInput(domain=[[outer] + holes], size=SizeField(30.0, 1.15), submeshes=subs), triangulate_scipy)
    assert len(boundary_loops(res.ikle, res.x, res.y)) == 1  # submeshes are stitched in
    assert abs(res.stats["area_total"] - 650000.0) < 1e-6
    assert res.stats["angle_min"] > 25.0
    tri = set(tuple(sorted(t)) for t in res.ikle.tolist())
    for s in subs:
        m = np.array([_node(res.x, res.y, p) for p in zip(s.x, s.y)])
        assert all(tuple(sorted(t)) in tri for t in m[s.ikle].tolist())


def test_channel_mesh():
    L = np.array([[0, 0], [100, 0]], float)
    R = np.array([[100, 20], [0, 20]], float)  # drawn the other way: reversed automatically
    x, y, ik = channel_mesh(L, R, np.array([[0, 10], [100, 10]], float), cross_nodes=5, interval=10)
    assert x.size == 11 * 5 and ik.shape == (2 * 10 * 4, 3)
    assert np.all(signed_areas(x, y, ik) > 0)
    assert abs(mesh_statistics(x, y, ik)["area_total"] - 2000.0) < 1e-9
    with pytest.raises(ValueError):
        channel_mesh(L, R, np.array([[0, 10], [100, 10]], float), cross_nodes=4)


def test_regular_grid():
    x, y, ik = regular_grid_mesh(10, 20, 2.0, 3.0, 4, 5, angle_deg=30)
    assert x.size == 30 and ik.shape == (40, 3)
    assert np.all(signed_areas(x, y, ik) > 0)
    assert abs(mesh_statistics(x, y, ik)["area_total"] - 8 * 15) < 1e-9
    assert abs(x[0] - 10) < 1e-12 and abs(y[0] - 20) < 1e-12
    for pat in ("left", "right"):
        assert np.all(signed_areas(*regular_grid_mesh(0, 0, 1, 1, 3, 3, pattern=pat)) > 0)


def test_bilinear_reproduces_a_plane():
    # raster 50 x 40 pixels of 2 m, top-left at (1000, 2000)
    nr, nc = 40, 50
    xc = 1000 + (np.arange(nc) + 0.5) * 2.0
    yc = 2000 - (np.arange(nr) + 0.5) * 2.0
    XC, YC = np.meshgrid(xc, yc)
    Z = 3.0 + 0.01 * XC - 0.02 * YC
    rng = np.random.default_rng(1)
    px = rng.uniform(1001, 1099, 200)
    py = rng.uniform(1921, 1999, 200)
    z = bilinear(Z, 1000, 2000, 2.0, -2.0, px, py)
    assert np.allclose(z, 3.0 + 0.01 * px - 0.02 * py)
    assert np.isnan(bilinear(Z, 1000, 2000, 2.0, -2.0, np.array([900.0]), np.array([1950.0])))[0]
    # nodata pixels are skipped, priority between sources, nearest fill
    Z2 = Z.copy()
    Z2[:, :25] = -9999
    s1 = RasterSource(Z2, 1000, 2000, 2.0, -2.0, nodata=-9999, name="partial")
    s2 = RasterSource(Z * 0 + 7.0, 1000, 2000, 2.0, -2.0, name="flat")
    v, src = interpolate_nodes(np.array([1010.0, 1090.0, 5000.0]), np.array([1950.0, 1950.0, 0.0]), [s1, s2])
    assert src.tolist() == [1, 0, -1] and v[0] == 7.0
    filled = fill_nearest(np.array([1010.0, 1090.0, 5000.0]), np.array([1950.0, 1950.0, 0.0]), v)
    assert filled[2] == v[1]


def test_boundary_table_from_lines():
    x, y, ik = regular_grid_mesh(0, 0, 10, 10, 20, 6)
    up = BoundaryLine(np.array([[0, -1], [0, 61]], float), parse_type("Q"), "upstream")
    down = BoundaryLine(np.array([[200, -1], [200, 61]], float), parse_type("544"), "downstream")
    table, report = build_boundary_table(x, y, ik, [up, down])
    codes = table.codes
    nodes = table.nodes
    assert len(table) == 2 * (20 + 6)
    assert np.all(codes[x[nodes] == 0] == (4, 5, 5))
    assert np.all(codes[x[nodes] == 200] == (5, 4, 4))
    assert np.all(codes[(x[nodes] > 0) & (x[nodes] < 200)] == (2, 2, 2))
    assert any("Liquid boundary 1" in r for r in report) and any("Liquid boundary 2" in r for r in report)
    txt = liq_template(table)
    head = txt.splitlines()[1].split()
    assert head[0] == "T" and set(head[1:]) == {"Q(%d)" % n for n in (1, 2) if "Q(%d)" % n in head} | {
        "SL(%d)" % n for n in (1, 2) if "SL(%d)" % n in head}
    assert len(head) == 3
    with pytest.raises(ValueError):
        parse_type("banana")


def test_snap_rings_to_submesh_nodes():
    # a ring vertex created by GEOS where the outline crosses a submesh edge
    ring = np.array([[0, 0], [10, 0], [10, 4.0], [10, 5.2], [0, 10]], float)
    nodes = np.array([[10, 5.0], [12, 5.0], [10, 3.0]], float)
    snapped = snap_rings_to_nodes([ring], nodes, tol=0.5)[0]
    assert [10.0, 5.0] in snapped.tolist()
    assert [10.0, 5.2] not in snapped.tolist()
    assert len(snapped) == 5
