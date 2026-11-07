"""River reach helper: banks from a DEM or bank lines, tight bends, boundary lines on the outline."""

import numpy as np
import pytest

from telemac_toolkit.mesher.channel import channel_mesh, matched_vertices
from telemac_toolkit.mesher.geometry import drop_spikes
from telemac_toolkit.mesher.river import end_boundary_lines, reach_from_banks, reach_from_centreline, resample, smooth_line


def terrain(x, y):
    """Meandering river 3 m deep; 49 m wide measured along y (less across the flow in bends)."""
    yc = 70.0 * np.sin(2 * np.pi * x / 900.0)
    d = np.abs(y - yc)
    floor = 10.0 - 0.002 * x + 0.006 * np.maximum(d - 25.0, 0.0) + 2.5e-5 * np.maximum(d - 25.0, 0.0) ** 2
    bank = np.clip((d - 18.0) / 8.0, 0.0, 1.0)
    return floor - 3.0 * (1.0 - bank)


def _true_width(x):
    """Bank-to-bank width normal to the flow for a bankfull depth of 2.5 m: 2 * 24.67 m * cos(heading)."""
    slope = 70.0 * 2 * np.pi / 900.0 * np.cos(2 * np.pi * x / 900.0)
    return 2 * (18.0 + 8.0 * 2.5 / 3.0) / np.hypot(1.0, slope)


def test_resample_and_smooth_keep_the_ends():
    c = resample([[0, 0], [10, 0], [10, 10]], 1.0)
    assert len(c) == 21 and np.allclose(c[0], [0, 0]) and np.allclose(c[-1], [10, 10])
    s = smooth_line(c, 2.0, 1.0)
    assert np.allclose(s[0], c[0]) and np.allclose(s[-1], c[-1])
    assert np.hypot(*(s[10] - [10, 0])) > 0.5  # the corner is rounded


def test_centreline_and_dem_find_the_banks_and_the_direction():
    rng = np.random.default_rng(1)
    xs = np.linspace(2000, 0, 60)  # drawn upstream
    cl = np.stack([xs, 70 * np.sin(2 * np.pi * xs / 900) + rng.normal(0, 4, xs.size)], axis=1)
    r = reach_from_centreline(cl, 8.0, dem=terrain, bankfull_depth=2.5, max_width=120.0)
    assert r.reversed and r.axis[0, 0] < r.axis[-1, 0]
    yc = 70 * np.sin(2 * np.pi * r.axis[:, 0] / 900)
    assert np.sqrt(np.mean((r.axis[:, 1] - yc) ** 2)) < 1.0  # the axis moved onto the river despite 4 m of noise
    inner = slice(5, -5)
    err = np.abs(r.width[inner] - _true_width(r.axis[inner, 0]))
    assert np.median(err) < 1.5 and err.max() < 4.0
    assert r.crossings().size == 0 and not r.warnings
    assert matched_vertices(r.left, r.right, r.axis)


def test_constant_width_keeps_the_drawn_direction():
    cl = [[0, 0], [300, 50], [600, 0]]
    r = reach_from_centreline(cl, 10.0, width=40.0)
    assert not r.reversed and np.allclose(r.width, 40.0)
    assert any("drawn direction" in m for m in r.log)
    with pytest.raises(ValueError, match="width"):
        reach_from_centreline(cl, 10.0)
    with pytest.raises(ValueError, match="largest width"):
        reach_from_centreline(cl, 10.0, dem=terrain, bankfull_depth=2.0)


def test_bank_lines_any_order_and_direction():
    x = np.linspace(0, 2000, 400)
    yc = 70 * np.sin(2 * np.pi * x / 900)
    sl = np.gradient(yc, x)
    nx_, ny_ = -sl / np.hypot(1, sl), 1 / np.hypot(1, sl)
    left = np.stack([x + 24.6 * nx_, yc + 24.6 * ny_], 1)
    right = np.stack([x - 24.6 * nx_, yc - 24.6 * ny_], 1)
    r = reach_from_banks(right[::-1], left, 8.0, dem=terrain)  # swapped and backwards
    assert r.axis[0, 0] < r.axis[-1, 0]
    assert abs(np.median(r.width) - 49.2) < 0.2
    # left really on the left looking downstream
    t = r.axis[1] - r.axis[0]
    v = r.left[0] - r.axis[0]
    assert t[0] * v[1] - t[1] * v[0] > 0
    assert r.crossings().size == 0


def test_tight_bend_normals_are_turned_until_sections_do_not_cross():
    # a 90 degree corner drawn as a polyline: raw normals at the corner cross on the inner bank
    cl = [[0, 0], [200, 0], [200, 200]]
    r = reach_from_centreline(cl, 5.0, width=60.0, smoothing=1.0)
    assert r.crossings().size == 0
    assert any("smoothed" in m for m in r.log) or not r.warnings


def test_boundary_lines_channel_and_valley():
    cl = [[0, 0], [1000, 0]]
    r = reach_from_centreline(cl, 10.0, width=40.0)
    ring = np.array([[0, -300], [1000, -300], [1000, 300], [0, 300], [0, -300]], float)
    lines = end_boundary_lines(r, ring, "channel", "valley", margin=6.0)
    (n1, t1, c1, d1), (n2, t2, c2, d2) = lines
    assert (n1, t1, n2, t2) == ("inflow", "Q", "outflow", "H")
    assert np.allclose(c1[:, 0], 0) and abs(np.ptp(c1[:, 1]) - 28.0) < 1e-6  # 40 m channel less 2 x 6 m
    assert np.allclose(c2[:, 0], 1000) and abs(np.ptp(c2[:, 1]) - 588.0) < 1e-6  # whole end, corners kept as walls
    assert d1 < 1e-9 and d2 < 1e-9


def test_channel_mesh_uses_the_reach_sections():
    r = reach_from_centreline([[0, 0], [400, 100], [800, 0]], 10.0, width=30.0)
    x, y, ik = channel_mesh(r.left, r.right, r.axis, cross_nodes=7, interval=10.0)
    assert len(x) == 7 * len(r.axis)


def test_drop_spikes():
    # A -> B -> A' folding back along the same line, as after snapping to a channel end
    ring = np.array([[0, 0], [10, 0], [10, 5], [10, 3], [10, 10], [0, 10]], float)
    out = drop_spikes(ring)
    assert len(out) == 5 and not any(np.allclose(p, [10, 5]) for p in out)
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    assert len(drop_spikes(square)) == 4
