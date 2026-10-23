"""Calculator, spatial operators, integrals and channel planforms (NumPy only)."""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from telemac_toolkit.mesher.regular import regular_grid_mesh  # noqa: E402
from telemac_toolkit.telemac_io import SelafinFile, write_selafin  # noqa: E402
from telemac_toolkit.telemac_io.calculator import (  # noqa: E402
    Calculator, CalculatorError, MeshOperators, parse_outputs, run_calculator, variable_names)
from telemac_toolkit.telemac_io.integrals import Section, element_mask_in_rings  # noqa: E402


def irregular(n=24, size=10.0, jitter=0.25, seed=1):
    """A regular grid with its interior nodes moved at random: an unstructured-like test mesh."""
    x, y, ik = regular_grid_mesh(0.0, 0.0, size / n, size / n, n, n, 0.0, "alternate")
    rng = np.random.default_rng(seed)
    inner = (x > 1e-9) & (x < size - 1e-9) & (y > 1e-9) & (y < size - 1e-9)
    h = size / n
    x = x.copy()
    y = y.copy()
    x[inner] += rng.uniform(-jitter, jitter, inner.sum()) * h
    y[inner] += rng.uniform(-jitter, jitter, inner.sum()) * h
    return x, y, ik


def test_gradient_exact_for_linear_fields():
    x, y, ik = irregular()
    ops = MeshOperators(x, y, ik)
    f = 2.0 * x - 3.0 * y + 7.0
    assert np.allclose(ops.ddx(f), 2.0, atol=1e-10)
    assert np.allclose(ops.ddy(f), -3.0, atol=1e-10)
    assert np.allclose(ops.laplacian(f)[(x > 0.5) & (x < 9.5) & (y > 0.5) & (y < 9.5)], 0.0, atol=1e-8)


def test_rotation_curl_divergence_strain():
    x, y, ik = irregular()
    ops = MeshOperators(x, y, ik)
    om = 0.3
    u, v = -om * (y - 5.0), om * (x - 5.0)  # solid-body rotation
    assert np.allclose(ops.curl(u, v), 2 * om, atol=1e-10)
    assert np.allclose(ops.div(u, v), 0.0, atol=1e-10)
    assert np.allclose(ops.strain(u, v), 0.0, atol=1e-10)
    assert np.all(ops.okubo_weiss(u, v) < 0)  # rotation dominates
    u2, v2 = 0.5 * x, -0.5 * y  # pure strain
    assert np.allclose(ops.div(u2, v2), 0.0, atol=1e-10)
    assert np.allclose(ops.curl(u2, v2), 0.0, atol=1e-10)
    assert np.allclose(ops.strain(u2, v2), math.sqrt(2 * 0.25 + 2 * 0.25), atol=1e-10)


def test_gradient_converges():
    errs = []
    for n in (16, 32, 64):
        x, y, ik = irregular(n=n, size=math.pi)
        ops = MeshOperators(x, y, ik)
        f = np.sin(x) * np.cos(y)
        inner = (x > 0.3) & (x < math.pi - 0.3) & (y > 0.3) & (y < math.pi - 0.3)
        errs.append(np.abs(ops.ddx(f) - np.cos(x) * np.cos(y))[inner].mean())
    assert errs[1] < 0.6 * errs[0] and errs[2] < 0.6 * errs[1]  # first order on irregular meshes


def test_laplacian_quadratic_regular_grid():
    x, y, ik = regular_grid_mesh(0.0, 0.0, 0.1, 0.1, 40, 40, 0.0, "alternate")
    ops = MeshOperators(x, y, ik)
    lap = ops.laplacian(x ** 2 + y ** 2)
    inner = (x > 0.15) & (x < 3.85) & (y > 0.15) & (y < 3.85)  # div(grad) needs two node layers
    assert np.allclose(lap[inner], 4.0, atol=1e-6)


def test_area_integrals():
    x, y, ik = irregular()
    ops = MeshOperators(x, y, ik)
    assert ops.integral(np.ones_like(x)) == pytest.approx(100.0, rel=1e-12)
    assert ops.integral(x) == pytest.approx(500.0, rel=1e-12)  # exact for linear fields
    ring = np.array([[0, 0], [5, 0], [5, 10], [0, 10], [0, 0]], float)
    mask = element_mask_in_rings(x, y, ik, [ring])
    assert ops.integral(np.ones_like(x), mask) == pytest.approx(50.0, rel=0.05)


def test_section_flux_exact():
    x, y, ik = irregular()
    h = 1.0 + 0.1 * x  # linear depth
    u = 0.5 + 0.02 * y  # linear velocity: h u is quadratic along a line
    v = np.zeros_like(x)
    line = np.array([[4.3, 10.0], [4.3, 0.0]])  # left bank (y = 10) to right bank, flow along +x
    sec = Section.build(line, x, y, ik)
    q = sec.discharge(h, u, v)
    # integral over y in [0, 10] of (1 + 0.43) (0.5 + 0.02 y)
    exact = 1.43 * (0.5 * 10 + 0.01 * 100)
    assert q == pytest.approx(exact, rel=1e-9)
    assert Section.build(line[::-1], x, y, ik).discharge(h, u, v) == pytest.approx(-exact, rel=1e-9)
    # the product of the nodal values, interpolated (FLOWRATE ALONG X written by TELEMAC), is close but not exact
    assert sec.flux(h * u, h * v) == pytest.approx(exact, rel=1e-3)
    assert sec.length == pytest.approx(10.0)
    assert sec.integral(np.ones_like(x)) == pytest.approx(10.0, rel=1e-12)


def _result(tmp_path):
    x, y, ik = irregular(n=12)
    times = np.array([0.0, 10.0, 20.0, 30.0])
    H = np.stack([np.full_like(x, 0.1 * k + 0.02 * x) for k in range(4)])
    U = np.stack([-0.2 * (y - 5) for _ in range(4)])
    V = np.stack([0.2 * (x - 5) for _ in range(4)])
    B = np.stack([-0.01 * x for _ in range(4)])
    S = H + B
    p = str(tmp_path / "r.slf")
    write_selafin(p, x, y, ik, {"VELOCITY U": U, "VELOCITY V": V, "WATER DEPTH": H, "FREE SURFACE": S, "BOTTOM": B},
                  times=times, units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M", "FREE SURFACE": "M", "BOTTOM": "M"})
    return p, x, y, times, H


def test_calculator_end_to_end(tmp_path):
    p, x, y, times, H = _result(tmp_path)
    out = str(tmp_path / "c.slf")
    text = """# comment line
VORTICITY [1/S] = curl(U, V)
HMAX [M] = tmax(H)
TPEAK [S] = tpeak(H)
RISE [M/S] = ddt(S)
VOLUME RATE = tint(H)
SLOPE X = ddx(S)
SLOPE Y = ddy(S)
FR = mag(U, V) / sqrt(g * maximum(H, 1e-6))
WETVORT = where(H > 0.15, VORTICITY, 0)
"""
    rep = run_calculator(p, out, text)
    assert rep.frames == [0, 1, 2, 3]
    with SelafinFile(out) as s:
        names = [n.strip() for n in s.var_names]
        assert names[:5] == ["VELOCITY U", "VELOCITY V", "WATER DEPTH", "FREE SURFACE", "BOTTOM"]
        assert "VORTICITY" in names and "SLOPE ALONG X" in names and "SLOPE ALONG Y" in names  # a vector for QGIS
        assert np.allclose(s.read("VORTICITY", 2), 0.4, atol=1e-6)
        assert np.allclose(s.read("HMAX", 0), H[-1], atol=1e-6)
        assert np.allclose(s.read("TPEAK", 0), 30.0)
        assert np.allclose(s.read("RISE", 1), 0.01, atol=1e-6)  # 0.1 m per 10 s
        assert np.allclose(s.read("VOLUME RATE", 0), np.trapezoid(H, times, axis=0) if hasattr(np, "trapezoid")
                           else np.trapz(H, times, axis=0), atol=1e-4)
        assert np.allclose(s.read("SLOPE ALONG X", 0), 0.02 - 0.01, atol=1e-6)
        wv = s.read("WETVORT", 0)
        assert np.all((wv == 0) | (np.abs(wv - 0.4) < 1e-5))
    from telemac_toolkit.telemac_io.variables import vector_pairs

    with SelafinFile(out) as s:
        assert [lab for _i, _j, lab in vector_pairs(s.var_names)] == ["Velocity", "Slope"]
    assert rep.renamed == {"SLOPE X": "SLOPE ALONG X", "SLOPE Y": "SLOPE ALONG Y"}


def test_vector_file_names():
    from telemac_toolkit.telemac_io.calculator import vector_file_names

    assert vector_file_names(["GRAD X", "GRAD Y", "H2"]) == {"GRAD X": "GRAD ALONG X", "GRAD Y": "GRAD ALONG Y"}
    # too long for "ALONG": U/V instead (16 characters at most)
    assert vector_file_names(["WATER SLOPE X", "WATER SLOPE Y"]) == {"WATER SLOPE X": "WATER SLOPE U",
                                                                      "WATER SLOPE Y": "WATER SLOPE V"}
    # already a vector for QGIS, or no partner: unchanged
    assert vector_file_names(["Q ALONG X", "Q ALONG Y", "A U", "A V", "LONE X"]) == {}
    # a name already in the file is never reused
    assert vector_file_names(["G X", "G Y"], taken=["G ALONG X"]) == {"G X": "G U", "G Y": "G V"}


def test_calculator_static_output_and_errors(tmp_path):
    p, *_ = _result(tmp_path)
    out = str(tmp_path / "m.slf")
    rep = run_calculator(p, out, "HMAX [M] = tmax(H)", keep=False)
    assert rep.frames == [3]
    with SelafinFile(out) as s:
        assert s.nframes == 1 and [n.strip() for n in s.var_names] == ["HMAX"]
    with pytest.raises(CalculatorError):
        parse_outputs("A VERY LONG OUTPUT NAME = H")
    with SelafinFile(p) as slf:
        with pytest.raises(CalculatorError):
            Calculator(slf, parse_outputs("X1 = __import__('os')"))
        with pytest.raises(CalculatorError):
            Calculator(slf, parse_outputs("X1 = H.real"))
        with pytest.raises(CalculatorError):
            Calculator(slf, parse_outputs("X1 = QQQ + 1"))
    with pytest.raises(CalculatorError):
        run_calculator(p, out, "BOTTOM = H")  # clashes with a kept variable


def test_variable_names():
    vi = variable_names(["VELOCITY U      ", "HAUTEUR D'EAU   ", "MY TRACER       "])
    assert vi["U"] == 0 and vi["VELOCITY_U"] == 0 and vi["velocity_u"] == 0
    assert vi["H"] == 1 and vi["water_depth"] == 1
    assert vi["MY_TRACER"] == 2


def test_planforms():
    from telemac_toolkit.mesher import planform as P

    b = P.bend(3.0, 2.0, 180.0, 3.0, 1.0, 0.05)
    assert b.length == pytest.approx(6.0 + 2.0 * math.pi, rel=1e-12)
    # the arc: every axis point of the bend is 2 m from the centre (3, 2)
    arc = (b.s > 3.0 + 1e-9) & (b.s < 3.0 + 2 * math.pi - 1e-9)
    assert np.allclose(np.hypot(b.axis[arc, 0] - 3.0, b.axis[arc, 1] - 2.0), 2.0)
    assert np.allclose(b.axis[-1], [0.0, 4.0], atol=1e-9)  # outlet straight comes back along -x
    assert np.allclose(np.hypot(*(b.left - b.axis).T), 0.5)
    s, n = P.axis_coordinates(b.axis, b.left[:, 0], b.left[:, 1])
    # the axis is a polyline: a point on the inner bank is ambiguous by about n x (turn per segment) / 2
    assert np.allclose(n, 0.5, atol=1e-4) and np.allclose(s, b.s, atol=0.01)
    s, n = P.axis_coordinates(b.axis, b.right[:, 0], b.right[:, 1])
    assert np.allclose(n, -0.5, atol=1e-4) and np.allclose(s, b.s, atol=0.01)
    with pytest.raises(ValueError):
        P.bend(1.0, 0.4, 90.0, 1.0, 1.0, 0.05)
    m = P.meanders(10.0, 60.0, 2.0, 1.0, 0.05)
    assert m.length == pytest.approx(20.0)
    assert P.min_radius(10.0, 60.0) == pytest.approx(10.0 / (2 * math.pi * math.radians(60.0)))
    assert abs(m.axis[-1, 1]) < 1e-3  # whole wavelengths end on the mean line
    with pytest.raises(ValueError):
        P.meanders(2.0, 100.0, 1.0, 1.0, 0.05)


def test_matched_channel_strip_is_orthogonal_in_a_bend():
    from telemac_toolkit.mesher import planform as P
    from telemac_toolkit.mesher.channel import channel_mesh, matched_vertices

    b = P.bend(1.0, 2.0, 90.0, 1.0, 1.0, 0.1)
    assert matched_vertices(b.left, b.right, b.axis)
    x, y, ik = channel_mesh(b.left, b.right, b.axis, cross_nodes=11, interval=0.1)
    n = 11
    m = len(x) // n
    pts = np.stack([x, y], axis=1).reshape(m, n, 2)
    # each cross-section passes through the arc centre (3 -> centre (1, 2)) inside the bend: radial lines
    arc = (b.s > 1.0 + 1e-9) & (b.s < 1.0 + math.pi - 1e-9)
    for j in np.flatnonzero(arc):
        d = pts[j] - np.array([1.0, 2.0])
        ang = np.arctan2(d[:, 1], d[:, 0])
        assert np.ptp(ang) < 1e-9


def test_liquid_boundaries_numbered_like_telemac():
    """TELEMAC (BIEF FRONT2) starts from the boundary point with the smallest x + y: in this U-bend that is
    the corner of the inlet's right bank, and walking counter-clockwise it meets the outlet first."""
    from telemac_toolkit.mesher import planform as P
    from telemac_toolkit.mesher.boundary import BoundaryLine, attach_mesh, build_boundary_table
    from telemac_toolkit.mesher.channel import channel_mesh
    from telemac_toolkit.telemac_io.cli import liquid_boundaries, read_cli, write_cli

    b = P.bend(3.0, 2.0, 180.0, 3.0, 1.0, 0.05)
    x, y, ik = channel_mesh(b.left, b.right, b.axis, cross_nodes=21, interval=0.05)
    inlet = BoundaryLine([b.right[0] + 0.03 * (b.left[0] - b.right[0]), b.left[0] - 0.03 * (b.left[0] - b.right[0])], (4, 5, 5), "inlet")
    outlet = BoundaryLine([b.right[-1] + 0.03 * (b.left[-1] - b.right[-1]), b.left[-1] - 0.03 * (b.left[-1] - b.right[-1])], (5, 4, 4), "outlet")
    table, _rep = build_boundary_table(x, y, ik, [inlet, outlet])
    lb = liquid_boundaries(table)
    assert [table.comments[q["start"]] for q in lb] == ["outlet", "inlet"]
    # the same table read back from its .cli needs the mesh to be numbered right
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "t.cli")
        write_cli(p, table)
        back = attach_mesh(read_cli(p), x, y, ik)
    assert [back.comments[q["start"]].strip() for q in liquid_boundaries(back)] == ["outlet", "inlet"]
    # a straight flume from (0, 0) to (5, 1): TELEMAC starts at the corner (0, 0), goes along y = 0 and meets
    # the outlet first, as in the flume run (tests/data/demo_flume: SL(1) Q(2))
    s = P.straight(5.0, 1.0, 0.1)
    x, y, ik = channel_mesh(s.left, s.right, s.axis, cross_nodes=11, interval=0.1)
    table, _rep = build_boundary_table(x, y + 0.5, ik, [BoundaryLine([(0.0, 0.05), (0.0, 0.95)], (4, 5, 5), "inlet"),
                                                         BoundaryLine([(5.0, 0.05), (5.0, 0.95)], (5, 4, 4), "outlet")])
    assert [table.comments[q["start"]] for q in liquid_boundaries(table)] == ["outlet", "inlet"]
    # a flume running towards -x: the walk starts at the outlet's corner (-5, 0) and meets the inlet first
    table, _rep = build_boundary_table(-x, y + 0.5, ik, [BoundaryLine([(0.0, 0.05), (0.0, 0.95)], (4, 5, 5), "inlet"),
                                                          BoundaryLine([(-5.0, 0.05), (-5.0, 0.95)], (5, 4, 4), "outlet")])
    assert [table.comments[q["start"]] for q in liquid_boundaries(table)] == ["inlet", "outlet"]
