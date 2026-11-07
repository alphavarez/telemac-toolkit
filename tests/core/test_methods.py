"""Mesh measures, Blue Kenue files, and the Gmsh interior mesher (skipped without the gmsh package)."""

import numpy as np
import pytest

from telemac_toolkit.mesher.metrics import METRIC_INFO, mesh_metrics
from telemac_toolkit.telemac_io.bluekenue import read_lines, read_t3s, write_i2s


def _grid(nx=10, ny=5):
    xs, ys = np.meshgrid(np.arange(nx + 1.0), np.arange(ny + 1.0))
    ik = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            ik += [[a, a + 1, a + nx + 2], [a, a + nx + 2, a + nx + 1]]
    return xs.ravel(), ys.ravel(), np.array(ik)


def test_metrics_on_a_known_grid():
    x, y, ik = _grid()
    plane = lambda px, py: 2.0 + 0.1 * px - 0.05 * py  # noqa: E731
    m = mesh_metrics(x, y, ik, outline=[np.array([[0, 0], [10, 0], [10, 5], [0, 5.0]])], dem=plane, level=4.0, velocity=1.0)
    assert m["angle_min"] == pytest.approx(45.0) and m["angle_max"] == pytest.approx(90.0)
    assert m["aspect_max"] == pytest.approx((np.sqrt(2) / 2) / (2 - np.sqrt(2)))  # circumradius over twice the inradius
    assert m["grading_max"] == pytest.approx(1.0)
    assert m["outline_dev_max"] == pytest.approx(0.0, abs=1e-12)
    assert m["bed_max"] == pytest.approx(0.0, abs=1e-12)  # a plane is linear in every element
    assert m["volume_error"] == pytest.approx(0.0, abs=1e-12)
    # time step: smallest altitude 1/sqrt(2), deepest element centre
    zc = plane(x[ik].mean(axis=1), y[ik].mean(axis=1))
    dt = np.min((1 / np.sqrt(2)) / (1.0 + np.sqrt(9.81 * (4.0 - zc))))
    assert m["dt_courant1"] == pytest.approx(dt)
    assert m["area"] == pytest.approx(50.0) and set(m) <= set(METRIC_INFO)


def test_metrics_see_a_curved_bed_and_a_short_outline():
    x, y, ik = _grid()
    bowl = lambda px, py: 0.1 * (px - 5) ** 2  # noqa: E731
    m = mesh_metrics(x, y, ik, outline=[np.array([[0, 0], [12, 0], [12, 5], [0, 5.0]])], dem=bowl, level=2.0)
    # linear interpolation of x² over unit elements overestimates by at most 0.1 * 1/4
    assert 0 < m["bed_max"] <= 0.025 + 1e-9 and m["volume_error"] < 0
    assert m["outline_dev_max"] == pytest.approx(2.0, abs=1e-6)  # the mesh stops 2 m short of the outline


def test_blue_kenue_lines_and_t3s(tmp_path):
    lines = [np.array([[0, 0], [10, 0], [10, 5]], float), np.array([[1, 1], [2, 1], [2, 2], [1, 1]], float)]
    p = tmp_path / "lines.i2s"
    write_i2s(str(p), lines, [5.0, 2.5])
    back = read_lines(str(p))
    assert len(back) == 2 and back[0][1] == 5.0 and np.allclose(back[1][0], lines[1])
    i3s = tmp_path / "z.i3s"
    i3s.write_text(":FileType i3s  ASCII  EnSim 1.0\n:EndHeader\n2 0\n0 0 1.5\n3 4 2.5\n")
    (c, v), = read_lines(str(i3s))
    assert c.shape == (2, 3) and c[1, 2] == 2.5
    t3s = tmp_path / "m.t3s"
    t3s.write_text(":FileType t3s  ASCII  EnSim 1.0\n:NodeCount 4\n:ElementCount 2\n:EndHeader\n"
                   "0 0 1\n1 0 2\n1 1 3\n0 1 4\n1 2 3\n1 3 4\n")
    x, y, z, ik = read_t3s(str(t3s))
    assert ik.tolist() == [[0, 1, 2], [0, 2, 3]] and z.tolist() == [1, 2, 3, 4]
    bad = tmp_path / "bad.t3s"
    bad.write_text(":NodeCount 4\n:ElementCount 2\n:EndHeader\n0 0\n")
    with pytest.raises(ValueError, match="announced"):
        read_t3s(str(bad))


def test_gmsh_interior_keeps_constraints_and_quality():
    gmsh_engine = pytest.importorskip("telemac_toolkit.mesher.gmsh_engine")
    if not gmsh_engine.available():
        pytest.skip("gmsh is not installed")
    from telemac_toolkit.mesher.backends import default_triangulator
    from telemac_toolkit.mesher.channel import channel_mesh
    from telemac_toolkit.mesher.engine import MeshInput, Submesh, generate, min_angles
    from telemac_toolkit.mesher.sizefield import PolygonDensity, SizeField
    from telemac_toolkit.telemac_io.mesh import boundary_loops

    # a 200 m x 100 m box with an island, a breakline and a structured strip along its bottom
    island = np.array([[120, 50], [140, 50], [140, 70], [120, 70]], float)
    x, y, ik = channel_mesh(np.array([[0, 20.0], [200, 20]]), np.array([[0, 0.0], [200, 0]]), None, cross_nodes=5, interval=5.0)
    dom = [[np.array([[0, 20], [200, 20], [200, 100], [0, 100]], float), island]]
    inp = MeshInput(dom, SizeField(8.0, 1.2, [PolygonDensity([np.array([[20, 40], [60, 40], [60, 80], [20, 80.0]])], 3.0)]),
                    hard_lines=[(np.array([[80, 30], [80, 90.0]]), False)], submeshes=[Submesh(x, y, ik)])
    for alg in ("delaunay", "frontal", "meshadapt", "quads"):
        res = generate(inp, default_triangulator(), None, interior=gmsh_engine.gmsh_interior(alg))
        e = np.sort(np.concatenate([res.ikle[:, [0, 1]], res.ikle[:, [1, 2]], res.ikle[:, [2, 0]]]), axis=1)
        u, c = np.unique(e, axis=0, return_counts=True)
        assert c.max() == 2, alg
        assert len(boundary_loops(res.ikle, res.x, res.y)) == 2, alg  # outline and the island
        keys = set(map(tuple, np.round(np.stack([res.x, res.y], 1), 6).tolist()))
        assert all((80.0, yy) in keys for yy in (30.0, 90.0)), alg  # breakline ends are nodes
        assert {tuple(np.round(p, 6)) for p in np.stack([x, y], 1).tolist()} <= keys, alg  # every strip node kept
        if alg != "quads":
            assert min_angles(res.x, res.y, res.ikle).min() > 25.0, alg
        else:
            assert len(res.quads) > 0
