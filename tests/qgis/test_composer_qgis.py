"""Composed, drawn and imported laboratory channels in a real (headless) QGIS, meshed, and with TELEMAC_CMD a
sudden-expansion flume run in TELEMAC-2D.

    QT_QPA_PLATFORM=offscreen TELEMAC_CMD="..." python -m pytest tests/qgis/test_composer_qgis.py
"""

import math
import os

import numpy as np
import pytest
from qgis.core import QgsGeometry, QgsPointXY, QgsVectorLayer

from test_mesher_qgis import _add
from test_river_qgis import qgis  # noqa: F401  (fixture)


def _mesh(path, out):
    import processing
    from telemac_toolkit.mesher.engine import min_angles
    from telemac_toolkit.telemac_io import SelafinFile
    from telemac_toolkit.telemac_io.mesh import boundary_loops

    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "OUTPUT": out, "LOAD": False})
    with SelafinFile(g["OUTPUT"]) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        assert len(boundary_loops(h.ikle, h.x, h.y)) >= 1
        return g, h.x.copy(), h.y.copy(), h.ikle.copy(), z, float(min_angles(h.x, h.y, h.ikle).min())


def _template(path, **kw):
    import processing

    params = dict(NEW_PROJECT=path, CELL=0.05, BED_Z0=0.3, Q=0.05, H=0.42)
    params.update(kw)
    return processing.run("telemactoolkit:channeltemplate", params)


def test_composed_s_bend_with_taper(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    p = str(tmp_path / "s.gpkg")
    _template(p, PLANFORM=3, WIDTH=1.0, SEGMENTS="straight 3\nbend left 2 90\nstraight 1\nbend right 2 90\ntaper 1.5 2\nstraight 4",
              SECTION=3, BANK_HEIGHT=0.1)  # V-shaped
    lay = M.open_mesh_project(p)
    assert sorted({f["channel"] for f in lay["channels"].getFeatures()}) == ["flume", "flume_axis"]
    st = M.read_settings(p)
    assert st["bottom_axis"] == "flume_axis" and st["bottom_halfwidth"] and " b" in " " + st["bottom_formula"].replace("/", " ")
    g, x, y, ik, z, a = _mesh(p, str(tmp_path / "s.slf"))
    assert a > 30.0
    assert z.min() == pytest.approx(0.3, abs=1e-3) and z.max() == pytest.approx(0.4, abs=1e-3)  # V: 0.1 m higher at the walls
    names = {f["name"] for f in lay["boundaries"].getFeatures()}
    assert names == {"inlet", "outlet"}


def test_expansion_and_elbow_pieces(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    p = str(tmp_path / "e.gpkg")
    _template(p, PLANFORM=3, SEGMENTS="width 0.5\nstraight 2\nexpand 1.5\nstraight 4", SLOPE=0.001)
    lay = M.open_mesh_project(p)
    assert sorted({f["channel"] for f in lay["channels"].getFeatures()}) == ["flume_1", "flume_2", "flume_axis"]
    assert next(lay["domain"].getFeatures()).geometry().area() == pytest.approx(0.5 * 2 + 1.5 * 4, rel=1e-9)
    *_r, z, a = _mesh(p, str(tmp_path / "e.slf"))
    assert a > 25.0 and z.max() == pytest.approx(0.3, abs=1e-3) and z.min() == pytest.approx(0.294, abs=1e-3)
    p2 = str(tmp_path / "l.gpkg")
    _template(p2, PLANFORM=3, SEGMENTS="straight 3\nelbow right 90\nstraight 3", WIDTH=1.0)
    lay = M.open_mesh_project(p2)
    assert next(lay["domain"].getFeatures()).geometry().area() == pytest.approx(6.0, rel=1e-6)  # an L of two 3.5 x 1 legs
    *_r, a = _mesh(p2, str(tmp_path / "l.slf"))
    assert a > 25.0


def test_old_planforms_take_the_new_sections(qgis, tmp_path):  # noqa: F811
    p = str(tmp_path / "b.gpkg")
    _template(p, PLANFORM=1, SECTION=2, WIDTH=1.2, BOTTOM_WIDTH=0.4, SIDE_SLOPE=1.0, BANK_HEIGHT=0.1, FLOODPLAIN_SLOPE=0.05)
    *_r, z, a = _mesh(p, str(tmp_path / "b.slf"))
    # compound: main channel 0.4 m at 0.3, floodplains 0.1 m higher from 0.3 m off the axis, rising 0.05 m/m to the
    # walls 0.6 m off it
    assert z.min() == pytest.approx(0.3, abs=1e-3) and z.max() == pytest.approx(0.4 + 0.05 * 0.3, abs=2e-3)
    p = str(tmp_path / "st.gpkg")
    _template(p, PLANFORM=0, SECTION=4, WIDTH=1.0, LENGTH=5.0, BANK_HEIGHT=0.08)  # parabolic, straight
    *_r, z, a = _mesh(p, str(tmp_path / "st.slf"))
    assert z.max() == pytest.approx(0.38, abs=1e-3) and a > 30.0


def test_drawn_and_imported_shapes(qgis, tmp_path):  # noqa: F811
    from osgeo import ogr
    from telemac_toolkit.qgis_bridge import meshing as M

    # a drawn centreline with a sharp corner, in local coordinates
    lyr = QgsVectorLayer("LineString?field=w:double", "cl", "memory")
    lyr.setCrs(M.local_crs())
    _add(lyr, QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(4, 0), QgsPointXY(4, 3)]), w=0.8)
    p = str(tmp_path / "d.gpkg")
    _template(p, PLANFORM=4, CENTRELINE=lyr, WIDTH_FIELD="w")
    lay = M.open_mesh_project(p)
    assert sorted({f["channel"] for f in lay["channels"].getFeatures()}) == ["flume_1", "flume_2", "flume_axis"]
    assert next(lay["domain"].getFeatures()).geometry().area() == pytest.approx(0.8 * (4.4 + 3.4) - 0.8 * 0.8, rel=1e-6)
    *_r, a = _mesh(p, str(tmp_path / "d.slf"))
    assert a > 25.0
    # drawn on a UTM map: the new project takes the line's CRS, so the flume keeps its size
    utm = QgsVectorLayer("LineString?crs=EPSG:32614", "utm", "memory")
    _add(utm, QgsGeometry.fromPolylineXY([QgsPointXY(500000, 2000000), QgsPointXY(500004, 2000000), QgsPointXY(500004, 2000003)]))
    p = str(tmp_path / "u.gpkg")
    _template(p, PLANFORM=4, CENTRELINE=utm, WIDTH=0.8)
    lay = M.open_mesh_project(p)
    assert lay["domain"].crs().authid() == "EPSG:32614"
    assert next(lay["domain"].getFeatures()).geometry().area() == pytest.approx(5.6, rel=1e-6)
    # an imported CSV centreline with widths
    csv = tmp_path / "cl.csv"
    csv.write_text("x;y;w\n" + "\n".join("%g;%g;%g" % (x, 0.3 * math.sin(x / 2), 0.6 + 0.05 * x) for x in np.linspace(0, 8, 33)))
    p = str(tmp_path / "c.gpkg")
    _template(p, PLANFORM=5, IMPORT_FILE=str(csv), IMPORT_AS=0)
    lay = M.open_mesh_project(p)
    ax = next(f for f in lay["channels"].getFeatures() if f["channel"] == "flume" and f["role"] == "left")
    assert ax["across_spacing"] == pytest.approx(0.05)
    *_r, a = _mesh(p, str(tmp_path / "c.slf"))
    assert a > 25.0
    # an imported outline from a DXF: a basin, an island and a wall
    dxf = str(tmp_path / "basin.dxf")
    ds = ogr.GetDriverByName("DXF").CreateDataSource(dxf)
    dl = ds.CreateLayer("entities")
    for ring in ([(0, 0), (6, 0), (6, 3), (0, 3), (0, 0)], [(2, 1), (3, 1), (3, 2), (2, 2), (2, 1)], [(4.5, 0.5), (4.5, 2.5)]):
        f = ogr.Feature(dl.GetLayerDefn())
        g = ogr.Geometry(ogr.wkbLineString)
        for x, y in ring:
            g.AddPoint_2D(x, y)
        f.SetGeometry(g)
        dl.CreateFeature(f)
    ds = None
    p = str(tmp_path / "x.gpkg")
    r = _template(p, PLANFORM=5, IMPORT_FILE=dxf, IMPORT_AS=1, CELL=0.1)
    lay = M.open_mesh_project(p)
    assert lay["domain"].featureCount() == 1 and lay["islands"].featureCount() == 1
    assert lay["constraint_lines"].featureCount() == 1 and r["BED_FORMULA"] == "0.3"
    g, x, y, ik, z, a = _mesh(p, str(tmp_path / "x.slf"))
    assert a > 25.0 and not np.any((x > 2.01) & (x < 2.99) & (y > 1.01) & (y < 1.99))  # nothing inside the island


def test_sudden_expansion_runs_in_telemac(qgis, tmp_path):  # noqa: F811
    """A 0.5 m flume opening into a 1.5 m one: 0.05 m³/s in, the level held at the outlet; it passes."""
    import shlex
    import subprocess

    import processing
    from telemac_toolkit.telemac_io import read_cli
    from telemac_toolkit.telemac_io.cli import liquid_boundaries
    from telemac_toolkit.tt_processing.bc_algorithms import load_mesh

    p = str(tmp_path / "exp.gpkg")
    _template(p, PLANFORM=3, SEGMENTS="width 0.5\nstraight 2\nexpand 1.5\nstraight 4", SLOPE=0.0, Q=0.05, H=0.42)
    g, *_r = _mesh(p, str(tmp_path / "exp.slf"))
    mesh, cli = g["OUTPUT"], g["OUTPUT_CLI"]
    *_m, table, _c = load_mesh(mesh, cli)
    num = {table.comments[b["start"]]: b["number"] for b in liquid_boundaries(table)}
    assert set(num) == {"inlet", "outlet"}
    assert read_cli(cli).records["litbor"].max() == 5  # inflow: prescribed tracer
    csv = tmp_path / "bc.csv"
    csv.write_text("t,q,h\n0,0.05,0.42\n600,0.05,0.42\n")
    tab = QgsVectorLayer("file://%s?delimiter=,&detectTypes=yes" % csv, "bc", "delimitedtext")
    fields = ["q", "h"] if num["inlet"] < num["outlet"] else ["h", "q"]
    lq = processing.run("telemactoolkit:liquidboundariesfile", {"MESH": mesh, "TABLE": tab, "TIME_FIELD": "t", "TIME_UNIT": 0,
                                                                "VALUE_FIELDS": fields, "OUTPUT": str(tmp_path / "exp.liq")})
    ini = processing.run("telemactoolkit:initialconditions", {"MESH": mesh, "MODE": 0, "LEVEL": 0.42,
                                                              "OUTPUT": str(tmp_path / "init.slf")})
    cas = processing.run("telemactoolkit:steeringfile", {
        "MESH": mesh, "LIQ": lq["OUTPUT"], "INITIAL": ini["OUTPUT"], "DURATION": 600.0, "DEPTH_GUESS": 0.12,
        "VELOCITY_GUESS": 0.5, "OUTPUT_PERIOD": 60.0, "FRICTION": 0.012, "TITLE": "expansion", "OUTPUT": str(tmp_path / "exp.cas")})
    cmd = os.environ.get("TELEMAC_CMD")
    if not cmd:
        pytest.skip("TELEMAC_CMD not set")
    run = subprocess.run(shlex.split(cmd) + ["telemac2d.py", os.path.basename(cas["OUTPUT"])], cwd=str(tmp_path),
                         capture_output=True, text=True, timeout=3000)
    log = run.stdout + run.stderr
    assert "CORRECT END OF RUN" in log, log[-3000:]
    flux = [line for line in log.splitlines() if "FLUX BOUNDARY" in line][-2:]
    v = sorted(abs(float(f.split(":")[1].split()[0])) for f in flux)
    assert v[0] == pytest.approx(0.05, rel=0.03) and v[1] == pytest.approx(0.05, rel=0.03)


def test_double_precision_mesh_loads_and_runs(qgis, tmp_path):  # noqa: F811
    """A SERAFIND mesh: QGIS reads it, the steering file declares it, TELEMAC runs on it."""
    import shlex
    import subprocess

    import processing
    from qgis.core import QgsMesh, QgsMeshLayer
    from telemac_toolkit.telemac_io import SelafinFile
    from telemac_toolkit.telemac_io.liq import read_steering

    p = str(tmp_path / "dbl.gpkg")
    _template(p, PLANFORM=0, LENGTH=6.0, WIDTH=0.5, SLOPE=0.001, Q=0.02, H=0.35, BED_Z0=0.3)
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": p, "OUTPUT": str(tmp_path / "dbl.slf"), "DOUBLE": True,
                                                       "LOAD": False})
    with SelafinFile(g["OUTPUT"]) as s:
        assert s.header.is_double
        n = s.header.npoin
    ml = QgsMeshLayer(g["OUTPUT"], "dbl", "mdal")
    qm = QgsMesh()
    assert ml.isValid()
    ml.dataProvider().populateMesh(qm)
    assert qm.vertexCount() == n
    ini = processing.run("telemactoolkit:initialconditions", {"MESH": g["OUTPUT"], "MODE": 0, "LEVEL": 0.35,
                                                              "OUTPUT": str(tmp_path / "init.slf")})
    with SelafinFile(ini["OUTPUT"]) as s:
        assert s.header.is_double
    cas = processing.run("telemactoolkit:steeringfile", {
        "MESH": g["OUTPUT"], "INITIAL": ini["OUTPUT"], "DURATION": 60.0, "DEPTH_GUESS": 0.05, "VELOCITY_GUESS": 0.5,
        "OUTPUT_PERIOD": 30.0, "FRICTION": 0.012, "TITLE": "double", "OUTPUT": str(tmp_path / "dbl.cas")})
    k = read_steering(cas["OUTPUT"])
    assert k["GEOMETRY FILE FORMAT"] == "SERAFIND" and k["PREVIOUS COMPUTATION FILE FORMAT"] == "SERAFIND"
    cmd = os.environ.get("TELEMAC_CMD")
    if not cmd:
        pytest.skip("TELEMAC_CMD not set")
    run = subprocess.run(shlex.split(cmd) + ["telemac2d.py", os.path.basename(cas["OUTPUT"])], cwd=str(tmp_path),
                         capture_output=True, text=True, timeout=3000)
    log = run.stdout + run.stderr
    assert "CORRECT END OF RUN" in log, log[-3000:]
    res = [f for f in os.listdir(str(tmp_path)) if f.endswith("_res.slf")]
    with SelafinFile(os.path.join(str(tmp_path), res[0])) as s:
        assert s.header.is_double and s.header.npoin == n
