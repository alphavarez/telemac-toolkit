"""Editing existing meshes in a real (headless) QGIS: re-mesh inside polygons, the QGIS mesh-editing round trip,
boundary lines from a .cli, snapping preset and drawing checks.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_edit_qgis.py
"""

import os

import numpy as np
import pytest
from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsGeometry, QgsMeshLayer, QgsPoint, QgsPointXY, QgsProject, QgsVectorLayer

from test_bc_qgis import case  # noqa: F401  (fixture)
from test_mesher_qgis import CRS, X0, Y0, _add, _app, _pts
from test_river_qgis import dem, qgis  # noqa: F401  (fixtures)


def plane(x, y):
    return 5.0 + 0.002 * (x - X0) - 0.004 * (y - Y0)


@pytest.fixture(scope="module")
def mesh(tmp_path_factory):
    _app()
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M

    d = tmp_path_factory.mktemp("edit")
    path = str(d / "proj.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), settings={"default_size": "20",
                                                                                   "bottom_formula": "5.0 + 0.002*x - 0.004*y"})
    _add(lay["domain"], QgsGeometry.fromPolygonXY([_pts([(0, 0), (600, 0), (600, 300), (0, 300), (0, 0)])]), name="box")
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(0, 40), (0, 260)])), name="inflow", type="Q", value=10.0)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(600, 40), (600, 260)])), name="outflow", type="H", value=6.0)
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "FORMULA": "5.0 + 0.002*(x-%r) - 0.004*(y-%r)" % (X0, Y0),
                                                       "OUTPUT": str(d / "geo.slf"), "LOAD": False})
    return d, path, g["OUTPUT"], lay


def _liquid(cli, slf):
    from telemac_toolkit.tt_processing.bc_algorithms import load_mesh
    from telemac_toolkit.telemac_io.cli import liquid_boundaries

    *_r, table, _c = load_mesh(slf, cli)
    return {table.comments[b["start"]]: len(b["nodes"]) for b in liquid_boundaries(table)}


def test_remesh_inside_polygons(mesh):
    import processing
    from telemac_toolkit.telemac_io import SelafinFile

    d, _path, slf, _lay = mesh
    polys = QgsVectorLayer("Polygon?crs=%s&field=size:double" % CRS, "areas", "memory")
    _add(polys, QgsGeometry.fromPolygonXY([_pts([(200, 100), (400, 100), (400, 200), (200, 200), (200, 100)])]), size=4.0)
    r = processing.run("telemactoolkit:remeshregion", {"INPUT": slf, "POLYGONS": polys, "SIZE_FIELD": "size", "MESH_CRS": CRS,
                                                       "OUTPUT": str(d / "refined.slf")})
    with SelafinFile(slf) as a, SelafinFile(r["OUTPUT"]) as b:
        assert b.header.npoin > a.header.npoin + 300
        assert b.header.var_names == a.header.var_names
        z = b.read("BOTTOM", 0)
        assert np.allclose(z, plane(b.header.x, b.header.y), atol=1e-4)  # linear bottom carried exactly
    assert _liquid(r["OUTPUT_CLI"], r["OUTPUT"]) == _liquid(os.path.splitext(slf)[0] + ".cli", slf)


def test_qgis_editing_round_trip(mesh):
    import processing
    import shutil
    from telemac_toolkit.telemac_io import SelafinFile

    d, _path, slf, _lay = mesh
    copy = str(d / "editing.slf")
    shutil.copyfile(slf, copy)
    lyr = QgsMeshLayer(copy, "m", "mdal")
    ct = QgsCoordinateTransform()
    assert lyr.startFrameEditing(ct)
    n0 = lyr.meshVertexCount()
    added = lyr.meshEditor().addPointsAsVertices([QgsPoint(X0 + 301.3, Y0 + 151.7, 0.0), QgsPoint(X0 + 123.4, Y0 + 77.7, 0.0)], 1e-6)
    assert added == 2
    assert lyr.commitFrameEditing(ct, False)
    del lyr
    with SelafinFile(copy) as s:
        assert s.header.nvar == 0  # QGIS keeps only the geometry
    r = processing.run("telemactoolkit:rebuildedited", {"EDITED": copy, "ORIGINAL": slf, "OUTPUT": str(d / "rebuilt.slf")})
    with SelafinFile(r["OUTPUT"]) as s:
        h = s.header
        assert h.npoin == n0 + 2 and "BOTTOM" in [v.strip() for v in h.var_names]
        assert np.allclose(s.read("BOTTOM", 0), plane(h.x, h.y), atol=1e-4)
    assert _liquid(r["OUTPUT_CLI"], r["OUTPUT"]) == _liquid(os.path.splitext(slf)[0] + ".cli", slf)


def test_boundary_lines_from_cli(mesh):
    import processing

    _d, _path, slf, _lay = mesh
    r = processing.run("telemactoolkit:clilines", {"MESH": slf, "MESH_CRS": CRS, "OUTPUT": "memory:"})
    feats = {f["name"]: f for f in r["OUTPUT"].getFeatures()}
    assert set(feats) == {"inflow", "outflow"}
    assert feats["inflow"]["type"] == "Q" and feats["outflow"]["type"] == "H"
    xs = [p.x() - X0 for p in feats["inflow"].geometry().vertices()]
    assert max(abs(v) for v in xs) < 1e-6  # on the left side of the box


def test_snapping_and_drawing_checks(mesh):
    from qgis.gui import QgsMapCanvas
    from qgis.PyQt.QtWidgets import QMainWindow
    from telemac_toolkit.gui.mesh_panel import MeshPanel
    from telemac_toolkit.qgis_bridge import meshing as M

    d, path, slf, lay = mesh

    class Iface:
        def __init__(self):
            self.win = QMainWindow()
            self.canvas = QgsMapCanvas(self.win)
            self.canvas.setDestinationCrs(QgsCoordinateReferenceSystem(CRS))
            self.active = None

        def mapCanvas(self):
            return self.canvas

        def setActiveLayer(self, lyr):
            self.active = lyr

    panel = MeshPanel(Iface())
    panel.set_project(path)
    assert panel.snapping_preset() >= 5
    cfg = QgsProject.instance().snappingConfig()
    assert cfg.enabled() and cfg.intersectionSnapping() and QgsProject.instance().topologicalEditing()
    # editing from the panel
    panel.last_output = slf
    lyr = panel.edit_mesh_start()
    assert lyr is not None and lyr.isEditable()
    lyr.meshEditor().addPointsAsVertices([QgsPoint(X0 + 455.5, Y0 + 222.2, 0.0)], 1e-6)
    res = panel.edit_mesh_finish()
    assert res and res["OUTPUT"].endswith("_edited.slf") and os.path.exists(res["OUTPUT_CLI"])
    # a boundary line beside the outline is reported
    lay2 = M.open_mesh_project(path)
    _add(lay2["boundaries"], QgsGeometry.fromPolylineXY(_pts([(30, 280), (200, 280)])), name="off", type="H", value=1.0)
    b = M.build_mesh_input(lay2, None, 20.0)
    assert any("'off' lies up to" in w for w in b.warnings)
    assert not any("'inflow' lies" in w for w in b.warnings)


def _run_chain(d, path, mesh, cli, label):
    """Rating curve, constant hydrograph, initial state from the curve level, steering file; then TELEMAC-2D when
    TELEMAC_CMD is set. Returns (fluxes in and out, log) or None without TELEMAC."""
    import shlex
    import subprocess

    import processing
    from telemac_toolkit.telemac_io.liq import read_stage_discharge

    from test_bc_qgis import _numbers

    num = _numbers(cli, mesh)
    assert set(num) == {"inflow", "outflow"}, num
    r = processing.run("telemactoolkit:stagedischarge", {"MESH": mesh, "SLOPE": 0.002, "MANNING": 0.035,
                                                         "OUTPUT": str(d / ("%s_rating.txt" % label)), "TABLE": "memory:"})
    q, z = read_stage_discharge(r["OUTPUT"])[num["outflow"]]
    csv = d / ("%s_hyd.csv" % label)
    csv.write_text("t,flow\n0,60\n10800,60\n")
    tab = QgsVectorLayer("file://%s?delimiter=,&detectTypes=yes" % csv, "hyd", "delimitedtext")
    lq = processing.run("telemactoolkit:liquidboundariesfile", {"MESH": mesh, "TABLE": tab, "TIME_FIELD": "t", "TIME_UNIT": 0,
                                                                "VALUE_FIELDS": ["flow"], "STAGE": r["OUTPUT"],
                                                                "OUTPUT": str(d / ("%s.liq" % label))})
    # initial level along the channel axis, from the rating-curve level downstream (as in test_bc_qgis)
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(path)
    axis = next(f for f in lay["channels"].getFeatures() if f["role"] == "axis").geometry()
    ax_layer = QgsVectorLayer("LineString?crs=%s" % CRS, "axis", "memory")
    _add(ax_layer, axis)
    z60 = float(np.interp(60.0, q, z))
    ini = processing.run("telemactoolkit:initialconditions", {
        "MESH": mesh, "MODE": 1, "LEVEL": z60 + 0.0019 * axis.length(), "LEVEL_DOWN": z60, "AXIS": ax_layer, "MESH_CRS": CRS,
        "OUTPUT": str(d / ("%s_init.slf" % label))})
    cas = processing.run("telemactoolkit:steeringfile", {
        "MESH": mesh, "LIQ": lq["OUTPUT"], "STAGE": r["OUTPUT"], "INITIAL": ini["OUTPUT"], "DURATION": 10800.0,
        "DEPTH_GUESS": 3.0, "VELOCITY_GUESS": 1.5, "OUTPUT_PERIOD": 3600.0, "FRICTION": 0.035, "TITLE": label,
        "OUTPUT": str(d / ("%s.cas" % label))})
    cmd = os.environ.get("TELEMAC_CMD")
    if not cmd:
        return None
    run = subprocess.run(shlex.split(cmd) + ["telemac2d.py", os.path.basename(cas["OUTPUT"])], cwd=str(d),
                         capture_output=True, text=True, timeout=3000)
    log = run.stdout + run.stderr
    assert "CORRECT END OF RUN" in log, log[-3000:]
    flux = [line for line in log.splitlines() if "FLUX BOUNDARY" in line][-2:]
    return sorted(abs(float(f.split(":")[1].split()[0])) for f in flux), log


def test_edited_river_meshes_run_in_telemac(case):  # noqa: F811
    """A river mesh re-meshed at its downstream end (new boundary nodes on the outflow) and edited in QGIS
    (vertices added in the channel) keeps a valid .cli; with TELEMAC_CMD both run, and the discharge passes."""
    import shutil

    import processing
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, mesh, cli = case
    # re-mesh a patch across the outflow end: floodplain and channel, finer than before
    poly = QgsVectorLayer("Polygon?crs=%s" % CRS, "area", "memory")
    _add(poly, QgsGeometry.fromPolygonXY([_pts([(1700, -220), (2100, -220), (2100, 220), (1700, 220), (1700, -220)])]))
    rm = processing.run("telemactoolkit:remeshregion", {"INPUT": mesh, "POLYGONS": poly, "SIZE": 12.0, "METHOD": 1,
                                                        "MESH_CRS": CRS, "OUTPUT": str(d / "remeshed.slf")})
    with SelafinFile(mesh) as a, SelafinFile(rm["OUTPUT"]) as b:
        assert b.header.npoin > a.header.npoin
        assert b.header.var_names == a.header.var_names
    # edit in QGIS: two vertices in the channel
    copy = str(d / "editing.slf")
    shutil.copyfile(mesh, copy)
    lyr = QgsMeshLayer(copy, "m", "mdal")
    ct = QgsCoordinateTransform()
    assert lyr.startFrameEditing(ct)
    pts = [QgsPoint(X0 + x, Y0 + 70 * np.sin(2 * np.pi * x / 900) + 2.0, 0.0) for x in (600.0, 1200.0)]
    assert lyr.meshEditor().addPointsAsVertices(pts, 1e-6) == 2
    assert lyr.commitFrameEditing(ct, False)
    del lyr
    ed = processing.run("telemactoolkit:rebuildedited", {"EDITED": copy, "ORIGINAL": mesh, "OUTPUT": str(d / "edited.slf")})
    for label, res in (("remeshed", rm), ("edited", ed)):
        out = _run_chain(d, path, res["OUTPUT"], res["OUTPUT_CLI"], label)
        if out is not None:
            flux, _log = out
            assert abs(flux[0] - 60) < 1.5 and abs(flux[1] - 60) < 1.5, (label, flux)
