import os

import numpy as np
import pytest

from telemac_io import SelafinFile, TriangleLocator, boundary_loops, mesh_statistics, sample_polyline
from telemac_io.cli import BoundaryTable, liquid_boundaries, read_cli, write_cli
from telemac_io.variables import canonical, display_name, lookup, vector_pairs

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
R2D = os.path.join(DATA, "r2d_tidal_flats.slf")


def test_locator_reproduces_linear_field():
    with SelafinFile(R2D) as s:
        x, y, ikle = s.x, s.y, s.ikle
    loc = TriangleLocator(x, y, ikle)
    f = 2.0 * x - 3.0 * y + 7.0
    rng = np.random.default_rng(1)
    # random points inside random triangles
    t = rng.integers(0, ikle.shape[0], 300)
    w = rng.dirichlet([1, 1, 1], 300)
    px = np.einsum("ij,ij->i", x[ikle[t]], w)
    py = np.einsum("ij,ij->i", y[ikle[t]], w)
    got = loc.interpolate(f, px, py)
    np.testing.assert_allclose(got, 2 * px - 3 * py + 7, rtol=1e-10, atol=1e-6)
    assert np.isnan(loc.interpolate(f, [x.min() - 1000], [y.min() - 1000]))[0]


def test_boundary_loops_square_with_hole():
    # 4x4 grid of nodes, hole in the middle cell
    xs, ys = np.meshgrid(np.arange(4.0), np.arange(4.0))
    x, y = xs.ravel(), ys.ravel()
    tris = []
    for j in range(3):
        for i in range(3):
            if (i, j) == (1, 1):
                continue
            a = j * 4 + i
            tris += [[a, a + 1, a + 5], [a, a + 5, a + 4]]
    loops = boundary_loops(np.array(tris), x, y)
    assert [len(lp) for lp in loops] == [12, 4]
    st = mesh_statistics(x, y, np.array(tris))
    assert st["elements"] == 16 and st["clockwise_elements"] == 0


def test_sample_polyline_keeps_vertices():
    px, py, s = sample_polyline([[0, 0], [10, 0], [10, 5]], step=3)
    assert s[-1] == pytest.approx(15)
    assert any(np.isclose(px, 10) & np.isclose(py, 0))
    assert np.all(np.diff(s) > 0)


def test_variable_catalogue():
    assert canonical("HAUTEUR D'EAU") == canonical("WATER DEPTH") == "water_depth"
    assert canonical("vitesse u") == "velocity_u"
    assert lookup("FOND").symbol == "B"
    assert display_name("WATER DEPTH", "M") == "Water depth (m)"
    assert vector_pairs(["VITESSE U", "VITESSE V", "FOND"]) == [(0, 1, "Velocity")]
    assert vector_pairs(["FOO X", "BAR", "FOO Y"]) == [(0, 2, "Foo")]
    assert canonical("MY TRACER") == "my_tracer"


def test_cli_round_trip_and_liquid_boundaries(tmp_path):
    nodes = np.arange(10)
    t = BoundaryTable.walls(nodes)
    r = t.records
    r["lihbor"][2:5], r["liubor"][2:5], r["livbor"][2:5] = 4, 5, 5   # Q inflow
    r["lihbor"][7:9], r["liubor"][7:9], r["livbor"][7:9] = 5, 4, 4   # H outflow
    t.comments[2] = "upstream"
    p = write_cli(str(tmp_path / "b.cli"), t)
    back = read_cli(p)
    np.testing.assert_array_equal(back.records, t.records)
    assert back.comments[2] == "upstream"
    lb = liquid_boundaries(back)
    assert [(b["number"], b["start"], b["end"], b["codes"]) for b in lb] == [(1, 2, 4, (4, 5, 5)), (2, 7, 8, (5, 4, 4))]
