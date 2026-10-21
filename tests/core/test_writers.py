"""Neutral mesh, structured quads, local frames and the pure-Python writers."""

import os

import numpy as np
import pytest

from telemac_toolkit.mesher.boundary import BoundaryLine, build_boundary_table, parse_type
from telemac_toolkit.mesher.channel import channel_mesh
from telemac_toolkit.mesher.interp import FormulaSource, PlaneSource, PrismaticChannelSource, channel_formula
from telemac_toolkit.mesher.meshdata import Frame, Mesh, find_quads, pairs_to_quads, split_quads
from telemac_toolkit.mesher.regular import regular_grid_mesh
from telemac_toolkit.mesher.writers import FORMATS, boundary_strings, write_mesh
from telemac_toolkit.telemac_io import SelafinFile, read_cli
from telemac_toolkit.telemac_io.mesh import signed_areas


def _flume(nx=20, ny=4, dx=0.5, dy=0.25, x0=0.0, y0=0.0):
    x, y, ik = regular_grid_mesh(x0, y0, dx, dy, nx, ny)
    return x, y, ik


def _table(x, y, ik, length):
    lines = [BoundaryLine(np.array([[0, 0.3], [0, 0.7]]), parse_type("Q"), "inlet"),
             BoundaryLine(np.array([[length, 0.3], [length, 0.7]]), parse_type("H"), "outlet")]
    return build_boundary_table(x, y, ik, lines)[0]


def test_quads_from_pairs_and_back():
    x, y, ik = _flume()
    q = pairs_to_quads(x, y, ik)
    assert q.shape == (80, 4)
    # counter-clockwise, area of each quad = dx * dy
    qa = 0.5 * np.abs(np.sum(x[q] * np.roll(y[q], -1, axis=1) - np.roll(x[q], -1, axis=1) * y[q], axis=1))
    assert np.allclose(qa, 0.125)
    found, used = find_quads(ik, q)
    assert len(found) == 80 and used.all()
    t = split_quads(x, y, q)
    assert len(t) == 160 and np.all(signed_areas(x, y, t) > 0)
    # a flipped diagonal is still found
    ik2 = ik.copy()
    a, b, c, d = q[0]
    ik2[0], ik2[1] = [a, b, d], [b, c, d]
    assert find_quads(ik2, q)[1].all()


def test_frame_round_trip():
    f = Frame((512000.0, 2155000.0), 30.0)
    x = np.array([512000.0, 512010.0, 511990.0])
    y = np.array([2155000.0, 2155005.0, 2155020.0])
    u, v = f.forward(x, y)
    assert abs(u[0]) < 1e-9 and abs(v[0]) < 1e-9
    assert abs(np.hypot(u[1], v[1]) - np.hypot(10, 5)) < 1e-9
    xx, yy = f.inverse(u, v)
    assert np.allclose(xx, x) and np.allclose(yy, y)


def test_bottom_sources():
    x = np.linspace(0, 10, 11)
    y = np.full(11, 0.5)
    assert np.allclose(PlaneSource.from_slope(1.0, 0.01).sample(x, y), 1.0 - 0.01 * x)
    src = PrismaticChannelSource(1.0, 0.01, 1.0, 1.5, 0.3, 0, 0.5)
    assert np.allclose(src.sample(x, y), 1.0 - 0.01 * x)  # on the axis: bed level
    assert np.allclose(src.sample(np.zeros(1), np.array([0.5 + 0.5 + 1.5 * 0.3 + 2])), 1.3)  # bank top
    f = channel_formula(1.0, 0.01, 1.0, 1.5, 0.3, 0, 0.5, 0)
    xs, ys = np.random.default_rng(3).uniform(-2, 12, (2, 300))
    assert np.allclose(FormulaSource(f).sample(xs, ys), src.sample(xs, ys))
    for bad in ("__import__('os')", "x.real", "open('f')", "lambda: 1"):
        with pytest.raises((ValueError, SyntaxError)):
            FormulaSource(bad)


def test_boundary_strings_follow_the_boundary():
    x, y, ik = _flume()
    mesh = Mesh(x, y, ik)
    t = _table(x, y, ik, 10.0)
    s = boundary_strings(mesh, t, {"inlet": 0.05})
    assert sum(len(b.nodes) for b in s) == len(t)
    open_ = [b for b in s if not b.is_wall]
    assert sorted(b.name for b in open_) == ["inlet", "outlet"]
    inlet = [b for b in open_ if b.name == "inlet"][0]
    assert inlet.value == 0.05 and inlet.codes == (4, 5, 5) and inlet.number in (1, 2)
    assert np.all(x[inlet.nodes] == 0) and len(inlet.nodes) == 3  # corners stay walls


def test_writers_produce_consistent_files(tmp_path):
    x, y, ik = _flume()
    q = pairs_to_quads(x, y, ik)
    mesh = Mesh(x, y, ik[:80], q[40:], z=1.0 - 0.01 * x)  # half triangles, half quads
    mesh.assign_zones([([np.array([[5, -1], [11, -1], [11, 2], [5, 2]], float)], 3, 0.02, "gravel")], default_friction=0.03)
    assert set(np.unique(mesh.material)) == {1, 3}
    t = _table(x, y, mesh.all_triangles(), 10.0)
    for key in ("2dm", "gmsh", "gid", "mike21", "adcirc"):
        p = str(tmp_path / ("m" + ("_gid" if key == "gid" else "") + FORMATS[key].ext))
        files = write_mesh(key, p, mesh, t, {"inlet": 0.05})
        assert os.path.exists(files[0]) and os.path.getsize(files[0]) > 0
    txt = open(str(tmp_path / "m.2dm")).read().splitlines()
    assert txt[0] == "MESH2D"
    assert sum(l.startswith("E3T") for l in txt) == 80 and sum(l.startswith("E4Q") for l in txt) == 40
    assert sum(l.startswith("ND") for l in txt) == len(x)
    ns = [l for l in txt if l.startswith("NS")]
    assert len(ns) == 2 and all(int(l.split()[-1]) < 0 for l in ns)
    assert 'MAT 3 "gravel"' in txt
    msh = open(str(tmp_path / "m.msh")).read()
    assert "$MeshFormat\n2.2 0 8" in msh and '"inlet"' in msh and '"gravel"' in msh
    f14 = open(str(tmp_path / "m.14")).read().splitlines()
    ne, nn = map(int, f14[1].split())
    assert ne == 160 and nn == len(x)  # quads split for ADCIRC
    assert float(f14[2].split()[3]) == pytest.approx(-(1.0 - 0.01 * x[0]))
    assert any("Number of open boundaries" in l and l.startswith("2 ") for l in f14)
    mk = open(str(tmp_path / "m.mesh")).read().splitlines()
    assert mk[0].split()[:3] == ["100079", "1000", str(len(x))]
    codes = [int(l.split()[4]) for l in mk[1:1 + len(x)]]
    assert set(codes) == {0, 1, 2, 3}


def test_telemac_writer_and_local_frame(tmp_path):
    # a 10 m flume drawn in UTM, written in a local frame with the offset stored
    x, y, ik = _flume(x0=512000.0, y0=2155000.0)
    mesh = Mesh(x, y, ik, z=np.full(len(x), 5.0))
    mesh.node_vars["BOTTOM FRICTION"] = ("", np.full(len(x), 0.03))
    f = Frame((512000.0, 2155000.0), 0.0)
    m = mesh.transformed(f)
    assert m.x.min() == 0.0 and m.y.min() == 0.0
    lines = [BoundaryLine(np.array([[0, 0.05], [0, 0.95]]), parse_type("Q"), "inlet"),
             BoundaryLine(np.array([[10, 0.05], [10, 0.95]]), parse_type("H"), "outlet")]
    t = build_boundary_table(m.x, m.y, m.all_triangles(), lines)[0]
    p = str(tmp_path / "flume.slf")
    files = write_mesh("telemac", p, m, t, {}, cli_path=str(tmp_path / "flume.cli"), liq_path=str(tmp_path / "flume.liq"),
                       store_offset=True)
    assert len(files) == 3
    with SelafinFile(p) as s:
        h = s.header
        assert h.origin == (512000.0, 2155000.0)
        assert np.allclose(h.x, x) and np.allclose(h.x_stored.min(), 0.0)
        assert [n.strip() for n in h.var_names] == ["BOTTOM", "BOTTOM FRICTION"]
    assert len(read_cli(str(tmp_path / "flume.cli"))) == h.ipobo.max()
    # not stored: truly local coordinates
    p2 = str(tmp_path / "flume_local.slf")
    write_mesh("telemac", p2, m, None, store_offset=False)
    with SelafinFile(p2) as s:
        assert s.header.origin == (0.0, 0.0) and abs(s.header.x.min()) < 1e-9


def test_gmsh_and_2dm_read_by_meshio(tmp_path):
    meshio = pytest.importorskip("meshio")
    x, y, ik = _flume()
    q = pairs_to_quads(x, y, ik)
    mesh = Mesh(x, y, ik[:80], q[40:], z=1.0 - 0.01 * x)
    t = _table(x, y, mesh.all_triangles(), 10.0)
    p = str(tmp_path / "m.msh")
    write_mesh("gmsh", p, mesh, t, {})
    m = meshio.read(p)
    cells = {c.type: len(c.data) for c in m.cells}
    assert cells["triangle"] == 80 and cells["quad"] == 40 and cells["line"] > 0
    assert np.allclose(m.points[:, 2], 1.0 - 0.01 * x)
    assert {"inlet", "outlet", "wall"} <= set(m.field_data)
    # quads must be counter-clockwise for every solver
    qq = [c.data for c in m.cells if c.type == "quad"][0]
    P = m.points[qq][:, :, :2]
    area = 0.5 * np.sum(P[:, :, 0] * np.roll(P[:, :, 1], -1, axis=1) - np.roll(P[:, :, 0], -1, axis=1) * P[:, :, 1], axis=1)
    assert np.all(area > 0)
