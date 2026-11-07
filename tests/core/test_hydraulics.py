"""Rating curves, initial states, boundary checks and the boundary/run text files."""

import math
import os

import numpy as np
import pytest

from telemac_toolkit.mesher.hydraulics import boundary_report, initial_state, rating_curve, times_to_seconds
from telemac_toolkit.telemac_io.liq import read_liq, read_stage_discharge, read_steering, write_liq, write_stage_discharge
from telemac_toolkit.telemac_io.steering import courant_time_step, steering_template


def test_rating_curve_matches_manning_in_a_rectangle():
    # a 20 m wide rectangular channel, banks 5 m high: Q = A R^(2/3) sqrt(S) / n with A = b h, P = b + 2 h
    st = np.array([0.0, 0.0, 20.0, 20.0])
    z = np.array([5.0, 0.0, 0.0, 5.0])
    st = np.array([0.0, 1e-9, 20.0 - 1e-9, 20.0])  # stations must increase: vertical walls as near-vertical segments
    Z, Q, notes = rating_curve(st, z, 0.03, 0.001, count=30)
    assert Z[0] == 0.0 and Q[0] == 0.0 and np.all(np.diff(Q) > 0)
    for h, q in zip(Z[1:], Q[1:]):
        A, P = 20 * h, 20 + 2 * h
        exact = A * (A / P) ** (2 / 3) * math.sqrt(0.001) / 0.03
        assert abs(q - exact) / exact < 1e-6
    assert any("Curve from" in n for n in notes)


def test_rating_curve_divides_at_roughness_changes():
    # main channel 20 m wide, 2 m deep, with a 100 m floodplain on one side
    st = np.array([0.0, 10.0, 10.0001, 30.0, 30.0001, 130.0])
    z = np.array([4.0, 2.0, 0.0, 0.0, 2.0, 2.0])
    banks = [10.0, 30.0001]  # bank stations at the top of the channel walls
    rough = np.array([0.08, 0.08, 0.03, 0.03, 0.03, 0.08])
    Zr, Qr, _ = rating_curve(st, z, rough, 0.001, banks=banks, count=60, top=3.0)
    Zs, Qs, _ = rating_curve(st, z, 0.03, 0.001, banks=banks, count=60, top=3.0)
    # below the floodplain only the main channel flows: exact Manning for a 20 m rectangle
    for lv in Zs[(Zs > 0) & (Zs < 1.9)]:
        A, P = 20 * lv, 20 + 2 * lv
        exact = A * (A / P) ** (2 / 3) * math.sqrt(0.001) / 0.03
        assert abs(float(np.interp(lv, Zs, Qs)) - exact) / exact < 1e-3
    # above it, a rough floodplain carries less than a smooth one
    assert float(np.interp(2.5, Zr, Qr)) < float(np.interp(2.5, Zs, Qs))
    with pytest.raises(ValueError):
        rating_curve(st, z, 0.03, -0.001)


def test_initial_state_and_times():
    z = np.array([0.0, 1.0, 2.0, 3.0])
    h, s = initial_state(z, level=1.5)
    assert np.allclose(h, [1.5, 0.5, 0, 0]) and np.allclose(s, [1.5, 1.5, 2, 3])
    h, s = initial_state(z, depth=0.2)
    assert np.allclose(h, 0.2) and np.allclose(s, z + 0.2)
    secs, first = times_to_seconds([0, 1.5, 3], "h")
    assert np.allclose(secs, [0, 5400, 10800]) and first is None
    import datetime as dt

    t0 = dt.datetime(2026, 9, 14, 6)
    secs, first = times_to_seconds([t0, t0 + dt.timedelta(minutes=90)], "datetime")
    assert np.allclose(secs, [0, 5400]) and first == t0


def test_text_files_round_trip(tmp_path):
    p = str(tmp_path / "a.liq")
    write_liq(p, [0, 3600], [("Q(1)", "m3/s", [10, 20]), ("SL(2)", "m", [4.4, 4.5])], comment="test")
    d = read_liq(p)
    assert list(d) == ["T", "Q(1)", "SL(2)"] and np.allclose(d["Q(1)"], [10, 20])
    with pytest.raises(ValueError, match="increase"):
        write_liq(p, [0, 0], [("Q(1)", "m3/s", [1, 2])])
    s = str(tmp_path / "r.txt")
    write_stage_discharge(s, {2: ([0, 5, 10], [1.0, 1.5, 1.8]), 1: ([0, 1], [0.0, 0.2])})
    c = read_stage_discharge(s)
    assert set(c) == {1, 2} and np.allclose(c[2][1], [1.0, 1.5, 1.8])
    # TELEMAC's own example: columns in either order
    open(s, "w").write("#\nQ(1)     Z(1)\nm3/s      m\n61.       0.\n62.       0.1\n#\nZ(2)     Q(2)\nm      m3/s\n10.       1.\n20.       2.\n")
    c = read_stage_discharge(s)
    assert np.allclose(c[1][0], [61, 62]) and np.allclose(c[2][0], [1, 2]) and np.allclose(c[2][1], [10, 20])


def test_steering_template_and_reader(tmp_path):
    cas = str(tmp_path / "run" / "t2d.cas")
    os.makedirs(os.path.dirname(cas))
    text = steering_template(cas, str(tmp_path / "geo.slf"), str(tmp_path / "geo.cli"), 2,
                             liquid_file=str(tmp_path / "run" / "q.liq"), stage_file=str(tmp_path / "run" / "r.txt"),
                             stage_boundaries=[1], prescribed_flowrates=[0, 60], prescribed_elevations=[0, 0],
                             initial_file=str(tmp_path / "run" / "init.slf"), kinds=["H", "Q"])
    open(cas, "w").write(text)
    k = read_steering(cas)
    assert k["GEOMETRY FILE"] == "../geo.slf" and k["LIQUID BOUNDARIES FILE"] == "q.liq"
    assert k["STAGE-DISCHARGE CURVES"] == "1;0" and k["VELOCITY PROFILES"] == "1;4"
    assert k["PRESCRIBED FLOWRATES"] == "0;60" and "COMPUTATION CONTINUED" not in k
    assert k["MAXIMUM NUMBER OF ITERATIONS FOR ADVECTION SCHEMES"] == "200"
    # a 10 m equilateral grid, 1 m deep: dt = height / sqrt(g h)
    x = np.array([0.0, 10.0, 5.0])
    y = np.array([0.0, 0.0, 5 * math.sqrt(3)])
    dt = courant_time_step(x, y, np.array([[0, 1, 2]]), depth=1.0, velocity=0.0)
    assert abs(dt - 5 * math.sqrt(3) / math.sqrt(9.81)) < 1e-9


def test_boundary_report():
    from telemac_toolkit.mesher.boundary import BoundaryLine, build_boundary_table, parse_type
    from telemac_toolkit.telemac_io.mesh import mesh_statistics  # noqa: F401

    # a 100 x 20 strip of squares, inflow on the left, outflow on the right
    nx, ny = 21, 5
    X, Y = np.meshgrid(np.linspace(0, 100, nx), np.linspace(0, 20, ny))
    x, y = X.ravel(), Y.ravel()
    tris = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            tris += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    ik = np.array(tris)
    z = 1.0 - 0.01 * x
    lines = [BoundaryLine(np.array([[0, 3.5], [0, 16.5]]), parse_type("Q"), "in"),
             BoundaryLine(np.array([[100, 3.5], [100, 16.5]]), parse_type("H"), "out")]
    table, _rep = build_boundary_table(x, y, ik, lines)
    b, msgs = boundary_report(x, y, z, table)
    assert sorted(bb["kind"] for bb in b) == ["H", "Q"]
    assert all(len(bb["nodes"]) == 3 for bb in b)
    assert sum("only 3 nodes" in m for m in msgs) == 2
    nums = {bb["name"]: bb["number"] for bb in b}
    liq = {"T": np.array([0.0, 10.0]), "Q(%d)" % nums["in"]: np.array([1.0, 1.0]), "SL(%d)" % nums["out"]: np.array([-1.0, -1.0])}
    _b, msgs = boundary_report(x, y, z, table, liq=liq)
    assert any("boundary is dry" in m for m in msgs)
    _b, msgs = boundary_report(x, y, z, table, cas={"PRESCRIBED FLOWRATES": "1.;0.;0."})
    assert any("3 values for 2" in m for m in msgs)


def test_fitting_time_step_ends_on_the_duration():
    """TELEMAC runs round(DURATION / DT) steps and allows 1 ms past the last time of the liquid boundaries file."""
    from telemac_toolkit.telemac_io.steering import fitting_time_step

    for dt_max, duration in ((0.6099, 10800.0), (0.55, 10800.0), (0.24, 10800.0), (1.7, 3600.0), (0.333, 1000.0),
                             (2.0, 10800.5), (0.05, 7.3), (13.0, 86400.0)):
        dt = fitting_time_step(dt_max, duration)
        n = int(duration / dt + 0.5)
        assert dt <= dt_max + 1e-12, (dt_max, duration, dt)
        assert abs(n * dt - duration) < 1e-3, (dt_max, duration, dt)
        assert abs(dt * 1e4 - round(dt * 1e4)) < 1e-6  # at most four decimals
    assert fitting_time_step(0.6099, 10800.0) == 0.6


def test_steering_declares_a_double_precision_geometry(tmp_path):
    from telemac_toolkit.telemac_io.selafin import write_selafin

    x = np.array([0.0, 1.0, 0.0, 1.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    ik = np.array([[0, 1, 2], [1, 3, 2]])
    single, double = str(tmp_path / "s.slf"), str(tmp_path / "d.slf")
    write_selafin(single, x, y, ik, {"BOTTOM": np.zeros(4)})
    write_selafin(double, x, y, ik, {"BOTTOM": np.zeros(4)}, double=True)
    cas = str(tmp_path / "t2d.cas")
    t = steering_template(cas, double, str(tmp_path / "d.cli"), 0, initial_file=double)
    assert "GEOMETRY FILE FORMAT = 'SERAFIND'" in t and "RESULTS FILE FORMAT = 'SERAFIND'" in t
    assert "PREVIOUS COMPUTATION FILE FORMAT = 'SERAFIND'" in t
    assert "FILE FORMAT" not in steering_template(cas, single, str(tmp_path / "s.cli"), 0)
