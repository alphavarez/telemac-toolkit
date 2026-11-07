"""Boundary conditions and run files in a real (headless) QGIS: rating curve, hydrograph table, check,
initial state, steering file. With TELEMAC_CMD set (a command that runs telemac2d.py), the case is run too.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_bc_qgis.py
"""

import os
import shlex
import subprocess

import numpy as np
import pytest
from qgis.core import QgsGeometry, QgsPointXY, QgsVectorLayer

from test_mesher_qgis import CRS, X0, Y0, _add
from test_river_qgis import _line_layer, dem, qgis  # noqa: F401  (fixtures)


@pytest.fixture(scope="module")
def case(dem, tmp_path_factory):  # noqa: F811
    import processing

    d = tmp_path_factory.mktemp("bc")
    xs = np.linspace(0, 2000, 60)
    cl = np.stack([xs, 70 * np.sin(2 * np.pi * xs / 900)], axis=1)
    path = str(d / "river.gpkg")
    processing.run("telemactoolkit:riverreach", {
        "NEW_PROJECT": path, "CENTRELINE": _line_layer([cl]), "WIDTH_MODE": 1, "DEM": dem, "BANKFULL_DEPTH": 2.5,
        "MAX_WIDTH": 120.0, "SPACING": 8.0, "FLOODPLAIN": 150.0, "FLOODPLAIN_SIZE": 30.0})
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "OUTPUT": str(d / "geo.slf"), "LOAD": False})
    return d, path, g["OUTPUT"], g["OUTPUT_CLI"]


def _numbers(cli, mesh):
    from telemac_toolkit.tt_processing.bc_algorithms import load_mesh
    from telemac_toolkit.telemac_io.cli import liquid_boundaries

    *_r, table, _c = load_mesh(mesh, cli)
    return {table.comments[b["start"]].split()[-1]: b["number"] for b in liquid_boundaries(table)}


def test_chain_to_steering_file(case):
    import processing
    from telemac_toolkit.telemac_io.liq import read_liq, read_stage_discharge, read_steering

    d, path, mesh, cli = case
    num = _numbers(cli, mesh)
    assert set(num) == {"inflow", "outflow"}

    # rating curve on the outflow (the only prescribed elevation)
    r = processing.run("telemactoolkit:stagedischarge", {"MESH": mesh, "SLOPE": 0.002, "MANNING": 0.035,
                                                         "OUTPUT": str(d / "rating.txt"), "TABLE": "memory:"})
    curves = read_stage_discharge(r["OUTPUT"])
    assert list(curves) == [num["outflow"]]
    q, z = curves[num["outflow"]]
    assert q[0] == 0 and np.all(np.diff(q) > 0) and abs(z[0] - 3.0) < 0.1  # starts at the channel bed
    assert r["TABLE"].featureCount() == len(q)
    z60 = float(np.interp(60.0, q, z))
    assert 3.6 < z60 < 4.8

    # hydrograph from a CSV file (hours), the outflow needs no SL column because of the curve
    csv = d / "hydrograph.csv"
    csv.write_text("time_h,flow\n0,60\n1.5,60\n3,60\n")
    tab = QgsVectorLayer("file://%s?delimiter=,&detectTypes=yes" % csv, "hyd", "delimitedtext")
    assert tab.isValid()
    lq = processing.run("telemactoolkit:liquidboundariesfile", {
        "MESH": mesh, "TABLE": tab, "TIME_FIELD": "time_h", "TIME_UNIT": 2, "VALUE_FIELDS": ["flow"], "STAGE": r["OUTPUT"],
        "OUTPUT": str(d / "flow.liq")})
    liq = read_liq(lq["OUTPUT"])
    assert list(liq) == ["T", "Q(%d)" % num["inflow"]]
    assert np.allclose(liq["T"], [0, 5400, 10800]) and np.allclose(liq["Q(%d)" % num["inflow"]], 60)

    # dates and times, columns by name
    csv2 = d / "dated.csv"
    csv2.write_text("when,q\n2026-09-14 06:00:00,10\n2026-09-14 07:30:00,40\n2026-09-14 09:00:00,25\n")
    tab2 = QgsVectorLayer("file://%s?delimiter=,&detectTypes=yes" % csv2, "dated", "delimitedtext")
    lq2 = processing.run("telemactoolkit:liquidboundariesfile", {
        "MESH": mesh, "TABLE": tab2, "TIME_FIELD": "when", "TIME_UNIT": 4, "COLUMNS": "Q(%d)=q" % num["inflow"],
        "STAGE": r["OUTPUT"], "OUTPUT": str(d / "dated.liq")})
    assert np.allclose(read_liq(lq2["OUTPUT"])["T"], [0, 5400, 10800])

    # initial state along the axis, from the curve level downstream
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(path)
    axis = next(f for f in lay["channels"].getFeatures() if f["role"] == "axis").geometry()
    ax_layer = QgsVectorLayer("LineString?crs=%s" % CRS, "axis", "memory")
    _add(ax_layer, axis)
    L = axis.length()
    ini = processing.run("telemactoolkit:initialconditions", {
        "MESH": mesh, "MODE": 1, "LEVEL": z60 + 0.0019 * L, "LEVEL_DOWN": z60, "AXIS": ax_layer, "MESH_CRS": CRS,
        "OUTPUT": str(d / "init.slf")})
    from telemac_toolkit.telemac_io import SelafinFile

    with SelafinFile(ini["OUTPUT"]) as s:
        h = s.read("WATER DEPTH", 0)
        assert 0.1 * len(h) < (h > 0).sum() < 0.6 * len(h)  # the channel is wet, the floodplain dry

    cas = processing.run("telemactoolkit:steeringfile", {
        "MESH": mesh, "LIQ": lq["OUTPUT"], "STAGE": r["OUTPUT"], "INITIAL": ini["OUTPUT"], "DURATION": 10800.0,
        "DEPTH_GUESS": 3.0, "VELOCITY_GUESS": 1.5, "OUTPUT_PERIOD": 1800.0, "FRICTION": 0.035, "TITLE": "river test",
        "OUTPUT": str(d / "t2d.cas")})
    k = read_steering(cas["OUTPUT"])
    order = sorted(num, key=num.get)
    assert k["GEOMETRY FILE"] == "geo.slf" and k["BOUNDARY CONDITIONS FILE"] == "geo.cli"
    assert k["STAGE-DISCHARGE CURVES"] == ";".join("1" if b == "outflow" else "0" for b in order)
    assert k["VELOCITY PROFILES"] == ";".join("4" if b == "inflow" else "1" for b in order)
    assert k["PREVIOUS COMPUTATION FILE"] == "init.slf"
    assert 0.1 < float(k["TIME STEP"]) < 5

    chk = processing.run("telemactoolkit:checkboundaries", {"MESH": mesh, "STEERING": cas["OUTPUT"], "MESH_CRS": CRS,
                                                            "OUTPUT": "memory:"})
    assert chk["WARNINGS"] == 0
    feats = {f["number"]: f for f in chk["OUTPUT"].getFeatures()}
    assert feats[num["inflow"]]["kind"] == "Q" and feats[num["outflow"]]["kind"] == "H"

    cmd = os.environ.get("TELEMAC_CMD")
    if cmd:
        run = subprocess.run(shlex.split(cmd) + ["telemac2d.py", "t2d.cas"], cwd=str(d), capture_output=True, text=True,
                             timeout=3000)
        log = run.stdout + run.stderr
        assert "CORRECT END OF RUN" in log, log[-3000:]
        flux = [line for line in log.splitlines() if "FLUX BOUNDARY" in line][-2:]
        vals = sorted(abs(float(f.split(":")[1].split()[0])) for f in flux)
        assert abs(vals[0] - 60) < 1.0 and abs(vals[1] - 60) < 1.0, flux
        with SelafinFile(str(d / "t2d_res.slf")) as s:
            fs = s.read("FREE SURFACE", len(s.times) - 1)
            hh = s.read("WATER DEPTH", len(s.times) - 1)
        from telemac_toolkit.tt_processing.bc_algorithms import load_mesh
        from telemac_toolkit.telemac_io.cli import liquid_boundaries

        *_r, table, _c = load_mesh(mesh, cli)
        outn = next(b for b in liquid_boundaries(table) if b["number"] == num["outflow"])["nodes"]
        wet = outn[hh[outn] > 0.05]
        assert abs(float(np.mean(fs[wet])) - z60) < 0.05  # the level follows the rating curve


def test_check_reports_problems(case):
    import processing
    from telemac_toolkit.telemac_io.liq import write_liq

    d, _path, mesh, cli = case
    num = _numbers(cli, mesh)
    bad = str(d / "bad.liq")
    # the outflow level below the whole bed, and a column for a boundary that does not exist
    write_liq(bad, [0, 3600], [("Q(%d)" % num["inflow"], "m3/s", [10, 20]), ("SL(%d)" % num["outflow"], "m", [-5, -5]),
                               ("Q(7)", "m3/s", [1, 1])])
    chk = processing.run("telemactoolkit:checkboundaries", {"MESH": mesh, "LIQ": bad, "OUTPUT": "memory:"})
    assert chk["WARNINGS"] >= 2


def test_bad_inputs(case):
    import processing
    from qgis.core import QgsProcessingException

    d, _path, mesh, _cli = case
    with pytest.raises(QgsProcessingException, match="no liquid boundary 9"):
        processing.run("telemactoolkit:stagedischarge", {"MESH": mesh, "BOUNDARY": 9, "OUTPUT": str(d / "x.txt")})
    pts = QgsVectorLayer("Point?crs=%s&field=t:double&field=q:double" % CRS, "t", "memory")
    _add(pts, QgsGeometry.fromPointXY(QgsPointXY(X0, Y0)), t=0.0, q=1.0)
    with pytest.raises(QgsProcessingException, match="at least two rows"):
        processing.run("telemactoolkit:liquidboundariesfile", {"MESH": mesh, "TABLE": pts, "TIME_FIELD": "t", "TIME_UNIT": 0,
                                                               "COLUMNS": "Q(1)=q", "OUTPUT": str(d / "y.liq")})


def test_panel_run_page(case, monkeypatch):
    import processing
    from qgis.core import QgsCoordinateReferenceSystem
    from qgis.gui import QgsMapCanvas
    from qgis.PyQt.QtWidgets import QMainWindow
    from telemac_toolkit.gui.mesh_panel import MeshPanel

    d, path, mesh, _cli = case

    class Iface:
        def __init__(self):
            self.win = QMainWindow()
            self.canvas = QgsMapCanvas(self.win)
            self.canvas.setDestinationCrs(QgsCoordinateReferenceSystem(CRS))

        def mapCanvas(self):
            return self.canvas

    calls = []
    monkeypatch.setattr(processing, "execAlgorithmDialog", lambda alg, params=None: calls.append((alg, dict(params or {}))))
    panel = MeshPanel(Iface())
    panel.set_project(path)
    panel.last_output = mesh
    assert set(panel.run_buttons) == {"stagedischarge", "liquidboundariesfile", "initialconditions", "steeringfile", "checkboundaries"}
    open(os.path.splitext(mesh)[0] + ".liq", "w").write("T Q(1)\ns m3/s\n0 1\n10 1\n")
    panel.run_buttons["steeringfile"].click()
    panel.run_buttons["stagedischarge"].click()
    alg, params = calls[0]
    assert alg == "telemactoolkit:steeringfile" and params["MESH"] == mesh and params["LIQ"].endswith("geo.liq")
    assert params["OUTPUT"].endswith("t2d_geo.cas")
    assert calls[1][1]["OUTPUT"].endswith("geo_rating.txt")
