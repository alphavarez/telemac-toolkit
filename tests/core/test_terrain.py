"""Terrain-driven sizes, breaklines from a DEM, solution-adaptive sizes and MMG (skipped without mmg2d)."""

import math

import numpy as np
import pytest

from telemac_toolkit.mesher.breaklines import breaklines_from_dem, simplify, thin, trace
from telemac_toolkit.mesher.demsize import apply_floor, curvature, field_sizes, limit_gradation, terrain_sizes
from telemac_toolkit.mesher.sizefield import GridDensity, SizeField


def levee(x, y):
    """Flat ground at 0 with a levee along y = 0: 2 m high, crest 4 m wide, sides 1:2 (8 m each)."""
    d = np.abs(np.asarray(y, dtype=float))
    return np.clip(2.0 - np.maximum(d - 2.0, 0.0) * 0.25, 0.0, 2.0)


def test_curvature_of_a_paraboloid():
    c = 0.5
    xs = np.arange(-10, 10 + c, c)
    X, Y = np.meshgrid(xs, xs)
    k = curvature(0.3 * X ** 2 + 0.05 * Y ** 2, c)
    assert np.allclose(k[5:-5, 5:-5], 0.6)  # largest second derivative: 2 * 0.3


def test_gradation_limit_holds_along_rows_and_columns():
    h = np.full((40, 60), 50.0)
    h[20, 30] = 1.0
    out = limit_gradation(h, 1.0, 1.2)
    assert out[20, 30] == 1.0
    assert out[20, 35] == pytest.approx(1.0 + 0.2 * 5)
    assert out[10, 30] == pytest.approx(1.0 + 0.2 * 10)
    dif = np.abs(np.diff(out, axis=1))
    assert dif.max() <= 0.2 + 1e-12


def test_terrain_sizes_refine_the_levee_sides_only():
    gd = terrain_sizes(levee, (-100, -200, 100, 200), tol=0.05, hmin=1.0, hmax=30.0, growth=1.2, cell=0.5)
    sf = SizeField(30.0, 1.2, grids=[gd])
    # the breaks of slope at |y| = 2 and 10 are refined; away from them the size grows at the growth rate
    # (about 1 m + 0.2 x 45 m at y = 55) up to the default size on flat ground far away
    assert sf(np.array([0.0]), np.array([55.0]))[0] == pytest.approx(1.0 + 0.2 * 45.0, rel=0.1)
    assert sf(np.array([0.0]), np.array([190.0]))[0] == pytest.approx(30.0)
    near = sf(np.array([0.0, 0.0]), np.array([2.0, 10.0]))
    assert np.all(near < 3.0)
    # sizes grow away from the levee no faster than the growth rate
    ys = np.arange(10.0, 50.0, 1.0)
    s = sf(np.zeros_like(ys), ys)
    assert np.all(np.diff(s) <= 0.2 * 1.0 + 1e-6)
    assert isinstance(gd, GridDensity) and gd.minimum >= 1.0


def test_floor_holds_sizes_next_to_fixed_edges():
    gd = GridDensity(0.0, 0.0, 1.0, 1.0, np.full((50, 50), 2.0))
    a = np.array([[10.0, 25.0]])
    b = np.array([[40.0, 25.0]])
    out = apply_floor(gd, a, b, np.array([8.0]), 1.2)
    assert out(np.array([25.0]), np.array([25.0]))[0] == pytest.approx(8.0)
    assert out(np.array([25.0]), np.array([35.0]))[0] == pytest.approx(8.0 - 0.2 * 10)
    assert out(np.array([25.0]), np.array([49.0]))[0] == pytest.approx(3.2)  # floor 8 - 0.2 * 24 = 3.2 > 2


def test_breaklines_find_the_levee_crest_and_toes():
    res = breaklines_from_dem(levee, (-60, -40, 60, 40), tol=0.05, size=20.0, cell=0.5, min_length=30.0)
    crest_y = [np.abs(c[:, 1]).mean() for c in res["crest"]]
    toe_y = [np.abs(c[:, 1]).mean() for c in res["toe"]]
    assert crest_y and toe_y
    assert all(abs(v - 2.0) < 1.0 for v in crest_y)  # the crest edges at |y| = 2 bend down
    assert all(abs(v - 10.0) < 1.0 for v in toe_y)  # the toes at |y| = 10 bend up
    assert max(float(np.ptp(c[:, 0])) for c in res["crest"]) > 80.0  # long, continuous lines


def test_thin_trace_simplify():
    img = np.zeros((20, 40), dtype=bool)
    img[8:12, 2:38] = True  # a 4-pixel-wide bar
    sk = thin(img)
    assert sk.sum() < img.sum() / 2 and sk[:, 10].sum() == 1
    lines = trace(sk)
    assert len(lines) == 1
    xy = simplify(lines[0][:, ::-1].astype(float), 0.5)
    assert len(xy) == 2


def test_field_sizes_follow_the_curvature_of_a_result():
    xs, ys = np.meshgrid(np.linspace(0, 100, 41), np.linspace(0, 50, 21))
    x, y = xs.ravel(), ys.ravel()
    ik = []
    for j in range(20):
        for i in range(40):
            a = j * 41 + i
            ik += [[a, a + 1, a + 42], [a, a + 42, a + 41]]
    ik = np.array(ik)
    f = np.tanh((x - 50.0) / 5.0)  # a sharp front at x = 50
    h = field_sizes(x, y, ik, f, tol=0.01, hmin=0.5, hmax=20.0)
    assert h[np.abs(x - 50) < 8].mean() < h[np.abs(x - 50) > 30].mean() / 3
    assert h.min() >= 0.5 and h.max() <= 20.0


def _mmg():
    from telemac_toolkit.mesher.mmg import find_mmg2d

    return find_mmg2d()


@pytest.mark.skipif(_mmg() is None, reason="mmg2d not found (set MMG2D or put it on the PATH)")
def test_mmg_keeps_the_boundary_and_kept_elements():
    from telemac_toolkit.mesher.backends import default_triangulator
    from telemac_toolkit.mesher.engine import MeshInput, generate, min_angles
    from telemac_toolkit.mesher.mmg import run_mmg2d
    from telemac_toolkit.telemac_io.mesh import boundary_loops

    ring = np.array([[0, 0], [200, 0], [200, 100], [0, 100], [0, 0]], float)
    r = generate(MeshInput([[ring]], SizeField(10.0, 1.2)), default_triangulator())
    keep = r.x[r.ikle].mean(axis=1) < 40  # leave the left end alone
    sizes = np.where(np.hypot(r.x - 120, r.y - 50) < 30, 2.5, 10.0)
    x, y, ik, log = run_mmg2d(r.x, r.y, r.ikle, sizes=sizes, keep=keep, hgrad=1.3)
    b0 = {tuple(np.round(p, 6)) for p in np.stack([r.x, r.y], 1)[np.concatenate(boundary_loops(r.ikle, r.x, r.y))].tolist()}
    b1 = {tuple(np.round(p, 6)) for p in np.stack([x, y], 1)[np.concatenate(boundary_loops(ik, x, y))].tolist()}
    assert b0 == b1  # the boundary nodes did not move
    kept = {tuple(sorted(map(tuple, np.round(np.stack([r.x[t], r.y[t]], 1), 6).tolist()))) for t in r.ikle[keep]}
    now = {tuple(sorted(map(tuple, np.round(np.stack([x[t], y[t]], 1), 6).tolist()))) for t in ik}
    assert kept <= now  # the kept elements are still there, unchanged
    assert len(x) > len(r.x)  # refined around (120, 50)
    near = np.hypot(x - 120, y - 50) < 20
    assert near.sum() > 4 * np.sum(np.hypot(r.x - 120, r.y - 50) < 20)
    assert min_angles(x, y, ik).min() > 15.0
    x2, y2, ik2, _ = run_mmg2d(r.x, r.y, r.ikle, optimise=True)
    assert abs(len(x2) - len(r.x)) < 0.2 * len(r.x)
    assert math.isfinite(float(x2.sum()))
