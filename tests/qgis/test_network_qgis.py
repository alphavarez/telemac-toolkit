"""River networks, nodes across that follow the width, and strip transition checks in a real (headless) QGIS.
With TELEMAC_CMD set, the confluence also runs in TELEMAC-2D (two inflows, one outflow).

    QT_QPA_PLATFORM=offscreen TELEMAC_CMD="..." python -m pytest tests/qgis/test_network_qgis.py
"""

import math
import os

import numpy as np
import pytest
from qgis.core import QgsCoordinateReferenceSystem, QgsFeature, QgsGeometry, QgsPointXY, QgsRasterLayer, QgsVectorLayer

from test_mesher_qgis import CRS, X0, Y0, _add, _app

ANG = math.radians(60)
T0 = np.array([1000 - 700 * math.cos(ANG), 700 * math.sin(ANG)])
J = np.array([1000.0, 0.0])


def _dseg(x, y, a, b):
    e = b - a
    t = np.clip(((x - a[0]) * e[0] + (y - a[1]) * e[1]) / float(e @ e), 0, 1)
    return np.hypot(x - a[0] - t * e[0], y - a[1] - t * e[1])


def confluence(x, y):
    """A main river along y = 0 (40 m wide, 3 m deep) and a tributary (25 m wide, 2 m deep) joining it at
    x = 1000 m from the north-west at 60 degrees; the floodplain slopes 2 per mil downstream (+x)."""
    fp = 12.0 - 0.002 * x + 0.001 * np.abs(y)
    bm = np.clip((np.abs(y) - 20.0) / 6.0, 0, 1)
    bt = np.clip((_dseg(x, y, T0, J) - 12.5) / 5.0, 0, 1)
    return fp - np.maximum(3.0 * (1 - bm), 2.0 * (1 - bt))


@pytest.fixture(scope="module")
def qgis():
    return _app()


@pytest.fixture(scope="module")
def conf_dem(qgis, tmp_path_factory):
    from osgeo import gdal, osr

    res = 2.0
    nx, ny = 1100, 550  # -100..2100 x -300..800
    xs = -100 + (np.arange(nx) + 0.5) * res
    ys = 800 - (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    path = str(tmp_path_factory.mktemp("conf") / "confluence.tif")
    ds = gdal.GetDriverByName("GTiff").Create(path, nx, ny, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((X0 - 100, res, 0, Y0 + 800, 0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32614)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(confluence(X, Y).astype(np.float32))
    ds = None
    lyr = QgsRasterLayer(path, "confluence dem")
    assert lyr.isValid()
    return lyr


def _lines(items, field="width"):
    lyr = QgsVectorLayer("LineString?crs=%s&field=%s:double" % (CRS, field), "lines", "memory")
    for c, v in items:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(X0 + a, Y0 + b) for a, b in c]))
        f.setAttribute(field, v)
        lyr.dataProvider().addFeatures([f])
    return lyr


def _liquid(cli_path, mesh):
    """{name: (number, nodes, codes)}, numbered as TELEMAC numbers them (with the mesh coordinates)."""
    from telemac_toolkit.telemac_io.cli import liquid_boundaries
    from telemac_toolkit.tt_processing.bc_algorithms import load_mesh

    *_r, t, _c = load_mesh(mesh, cli_path)
    return {t.comments[b["start"]]: (b["number"], len(b["nodes"]), b["codes"]) for b in liquid_boundaries(t)}


@pytest.fixture(scope="module")
def network(conf_dem, tmp_path_factory):
    import processing

    d = tmp_path_factory.mktemp("network")
    # the main river in one line; the tributary stops 5 m short of it, the helper joins them
    cl = _lines([([(0, 0), (2000, 0)], 40.0), ([tuple(T0), (1003.0, 4.0)], 25.0)])
    path = str(d / "net.gpkg")
    processing.run("telemactoolkit:riverreach", {"NEW_PROJECT": path, "SOURCE": 0, "CENTRELINE": cl, "WIDTH_FIELD": "width",
                                                 "DEM": conf_dem, "SPACING": 8.0, "FLOODPLAIN": 150.0, "H": 7.0})
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [conf_dem], "OUTPUT": str(d / "net.slf"),
                                                       "LOAD": False})
    return d, path, g["OUTPUT"], g["OUTPUT_CLI"]


def test_confluence_project_and_mesh(network):
    from telemac_toolkit.mesher.engine import min_angles
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile
    from telemac_toolkit.telemac_io.mesh import boundary_loops

    d, path, mesh, cli = network
    lay = M.open_mesh_project(path)
    assert sorted({f["channel"] for f in lay["channels"].getFeatures()}) == ["river_1", "river_2", "river_3"]
    assert all(f["across_spacing"] == 8.0 for f in lay["channels"].getFeatures())
    guides = [f for f in lay["constraint_lines"].getFeatures()]
    assert len(guides) == 3 and all(f["mode"] == "guide" for f in guides)
    areas = [(f["name"], f["size"]) for f in lay["refinement_areas"].getFeatures()]
    assert len(areas) == 1 and 4.0 < areas[0][1] <= 8.0
    bnd = sorted((f["name"], f["type"], f["value"]) for f in lay["boundaries"].getFeatures())
    assert [b[1] for b in bnd] == ["Q", "H", "Q"] and bnd[1][2] == 7.0
    with SelafinFile(mesh) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        a = min_angles(h.x, h.y, h.ikle)
        assert a.min() > 18.0
        assert len(boundary_loops(h.ikle, h.x, h.y)) == 1
        assert np.max(np.abs(z - confluence(h.x - X0, h.y - Y0))) < 0.5
        # the junction: no node deeper than the main bed, nodes on both banks of the tributary mouth
        near = np.hypot(h.x - X0 - 1000, h.y - Y0) < 40
        assert near.sum() > 40
    liq = _liquid(cli, mesh)
    assert {k: v[2] for k, v in liq.items()} == {"river_1 inflow": (4, 5, 5), "river_2 outflow": (5, 4, 4),
                                                 "river_3 inflow": (4, 5, 5)}
    assert liq["river_1 inflow"][1] >= 3 and liq["river_3 inflow"][1] >= 1


def test_confluence_runs_in_telemac(network):
    """Main river 40 m³/s and tributary 20 m³/s in, the rating curve at the outflow: 60 m³/s leave."""
    import shlex
    import subprocess

    import processing
    from telemac_toolkit.telemac_io.liq import read_stage_discharge

    d, path, mesh, cli = network
    num = {k: v[0] for k, v in _liquid(cli, mesh).items()}
    r = processing.run("telemactoolkit:stagedischarge", {"MESH": mesh, "SLOPE": 0.002, "MANNING": 0.035,
                                                         "OUTPUT": str(d / "rating.txt"), "TABLE": "memory:"})
    q, z = read_stage_discharge(r["OUTPUT"])[num["river_2 outflow"]]
    flows = {"river_1 inflow": "main", "river_3 inflow": "trib"}
    csv = d / "hyd.csv"
    csv.write_text("t,main,trib\n0,40,20\n10800,40,20\n")
    tab = QgsVectorLayer("file://%s?delimiter=,&detectTypes=yes" % csv, "hyd", "delimitedtext")
    fields = [flows[k] for k in sorted(flows, key=lambda k: num[k])]
    lq = processing.run("telemactoolkit:liquidboundariesfile", {"MESH": mesh, "TABLE": tab, "TIME_FIELD": "t", "TIME_UNIT": 0,
                                                                "VALUE_FIELDS": fields, "STAGE": r["OUTPUT"],
                                                                "OUTPUT": str(d / "net.liq")})
    from telemac_toolkit.qgis_bridge import meshing as M

    # the main river's axis alone leaves the tributary's inflow dry: the steering tool says so
    axis = QgsVectorLayer("LineString?crs=%s" % CRS, "axis", "memory")
    _add(axis, QgsGeometry.fromPolylineXY([QgsPointXY(X0, Y0), QgsPointXY(X0 + 2000, Y0)]))
    z60 = float(np.interp(60.0, q, z))
    dry = processing.run("telemactoolkit:initialconditions", {
        "MESH": mesh, "MODE": 1, "LEVEL": z60 + 0.002 * 2000, "LEVEL_DOWN": z60, "AXIS": axis, "MESH_CRS": CRS,
        "OUTPUT": str(d / "init_dry.slf")})
    from qgis.core import QgsProcessingContext, QgsProcessingFeedback

    class _Fb(QgsProcessingFeedback):
        def __init__(self):
            super().__init__()
            self.warn = []

        def pushWarning(self, m):
            self.warn.append(m)

    fb = _Fb()
    processing.run("telemactoolkit:steeringfile", {
        "MESH": mesh, "LIQ": lq["OUTPUT"], "STAGE": r["OUTPUT"], "INITIAL": dry["OUTPUT"], "DURATION": 10800.0,
        "OUTPUT": str(d / "dry.cas")}, context=QgsProcessingContext(), feedback=fb)
    assert any("river_3 inflow" in w and "dry" in w for w in fb.warn), fb.warn
    # every branch wet: a depth above the bed of the three channel axes
    ini = processing.run("telemactoolkit:initialconditions", {
        "MESH": mesh, "MODE": 3, "DEPTH": 1.2, "AXIS": M.open_mesh_project(path)["channels"], "MESH_CRS": CRS,
        "OUTPUT": str(d / "init.slf")})
    cas = processing.run("telemactoolkit:steeringfile", {
        "MESH": mesh, "LIQ": lq["OUTPUT"], "STAGE": r["OUTPUT"], "INITIAL": ini["OUTPUT"], "DURATION": 10800.0,
        "DEPTH_GUESS": 3.0, "VELOCITY_GUESS": 1.5, "OUTPUT_PERIOD": 3600.0, "FRICTION": 0.035, "TITLE": "confluence",
        "OUTPUT": str(d / "net.cas")})
    cmd = os.environ.get("TELEMAC_CMD")
    if not cmd:
        pytest.skip("TELEMAC_CMD not set")
    run = subprocess.run(shlex.split(cmd) + ["telemac2d.py", os.path.basename(cas["OUTPUT"])], cwd=str(d),
                         capture_output=True, text=True, timeout=3000)
    log = run.stdout + run.stderr
    assert "CORRECT END OF RUN" in log, log[-3000:]
    flux = [line for line in log.splitlines() if "FLUX BOUNDARY" in line][-3:]
    v = sorted(abs(float(f.split(":")[1].split()[0])) for f in flux)
    assert v[0] == pytest.approx(20.0, rel=0.02) and v[1] == pytest.approx(40.0, rel=0.02)
    assert v[2] == pytest.approx(60.0, rel=0.05)


def test_widening_river_follows_the_width(qgis, tmp_path):
    import processing
    from telemac_toolkit.mesher.engine import min_angles
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile

    x = np.linspace(0, 1500, 151)
    w = 30 + 120 * x / 1500
    banks = _lines([([(a, b) for a, b in zip(x, w / 2)], 0.0), ([(a, b) for a, b in zip(x, -w / 2)], 0.0)])
    path = str(tmp_path / "wide.gpkg")
    processing.run("telemactoolkit:riverreach", {"NEW_PROJECT": path, "SOURCE": 1, "BANKS": banks, "SPACING": 10.0,
                                                 "FLOODPLAIN": 100.0, "NAME": "main"})
    lay = M.open_mesh_project(path)
    b = M.build_mesh_input(lay, QgsCoordinateReferenceSystem(CRS), 40.0)
    line = next(i for i in b.info if i.startswith("Channel main"))
    assert "7 to 1" in line, line  # 7 across at the narrow end (the fewest), 30 m wide; more where it reaches 150 m
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "OUTPUT": str(tmp_path / "wide.slf"), "LOAD": False})
    with SelafinFile(g["OUTPUT"]) as s:
        h = s.header
        assert min_angles(h.x, h.y, h.ikle).min() > 18.0
        xs, ys = h.x - X0, h.y - Y0
        up = (np.abs(xs - 100) < 5) & (np.abs(ys) < 0.5 * (30 + 120 * 100 / 1500) - 0.1)
        down = (np.abs(xs - 1400) < 5) & (np.abs(ys) < 0.5 * (30 + 120 * 1400 / 1500) - 0.1)
        assert down.sum() > 2 * up.sum()  # more nodes across where the river is wider


def test_strip_transition_checks(qgis, tmp_path):
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "checks.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), settings={"default_size": "30"}, add_to_project=False)

    def pts(c):
        return [QgsPointXY(X0 + a, Y0 + b) for a, b in c]

    _add(lay["domain"], QgsGeometry.fromPolygonXY([pts([(0, -100), (1000, -100), (1000, 100), (0, 100), (0, -100)])]))
    for role, c in (("left", [(0, 20), (1000, 20)]), ("right", [(0, -20), (1000, -20)]), ("axis", [(0, 0), (1000, 0)])):
        _add(lay["channels"], QgsGeometry.fromPolylineXY(pts(c)), channel="c", role=role, nodes_across=9, spacing=10.0)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(pts([(0, -15), (0, 15)])), name="in", type="Q")
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(pts([(1000, 2), (1000, 15)])), name="out", type="H")
    _add(lay["refinement_areas"], QgsGeometry.fromPolygonXY([pts([(400, 20), (500, 20), (500, 60), (400, 60), (400, 20)])]),
         name="fine", size=2.0)
    b = M.build_mesh_input(M.open_mesh_project(path), QgsCoordinateReferenceSystem(CRS), 30.0)
    text = " | ".join(b.warnings)
    assert "covers only part" in text and "'c'" in text  # the outflow line misses the right half of the end
    assert "more than 2.5 times the element size" in text  # 10 m strip edges next to 2 m elements
    assert "first cross-section" not in text  # the inflow covers its end
    lay = M.open_mesh_project(path)
    out = next(f for f in lay["boundaries"].getFeatures() if f["name"] == "out")
    lay["boundaries"].startEditing()
    lay["boundaries"].deleteFeature(out.id())
    lay["boundaries"].commitChanges()
    b = M.build_mesh_input(M.open_mesh_project(path), QgsCoordinateReferenceSystem(CRS), 30.0)
    assert any("last cross-section lies on the domain outline but no boundary line" in w for w in b.warnings)
