"""Mesher in a real (headless) QGIS: mesh project, GEOS backend, Processing algorithms, workflow panel.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis
"""

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
_plugins = os.path.join(sys.prefix, "share", "qgis", "python", "plugins")
if os.path.isdir(_plugins) and _plugins not in sys.path:
    sys.path.append(_plugins)

from qgis.core import (  # noqa: E402
    QgsApplication,
    QgsCoordinateReferenceSystem,
    QgsFeature,
    QgsGeometry,
    QgsMesh,
    QgsMeshLayer,
    QgsPointXY,
    QgsProject,
    QgsRasterLayer,
    QgsVectorFileWriter,
    QgsVectorLayer,
)

X0, Y0 = 512000.0, 2155000.0
_PROVIDER = None
_APP = None
CRS = "EPSG:32614"


def _app():
    global _APP
    app = QgsApplication.instance()
    if app is None:
        # keep a reference for the whole session: if the Python wrapper is collected,
        # QGIS is torn down under the feet of the next test module
        app = _APP = QgsApplication([], True)
        app.initQgis()
    from processing.core.Processing import Processing

    Processing.initialize()
    reg = QgsApplication.processingRegistry()
    if reg.providerById("telemactoolkit") is None:
        from telemac_toolkit.tt_processing.provider import TelemacProvider

        global _PROVIDER
        _PROVIDER = TelemacProvider()  # keep a Python reference, or its algorithms are garbage collected
        reg.addProvider(_PROVIDER)
    return app


@pytest.fixture(scope="module")
def qgis():
    return _app()


def _add(layer, geom, **attrs):
    f = QgsFeature(layer.fields())
    f.setGeometry(geom)
    for k, v in attrs.items():
        f.setAttribute(k, v)
    layer.startEditing()
    assert layer.addFeature(f)
    assert layer.commitChanges(), layer.commitErrors()


def _pts(pts):
    return [QgsPointXY(X0 + a, Y0 + b) for a, b in pts]


def _plane(x, y):
    return 5.0 + 0.002 * (x - X0) - 0.004 * (y - Y0)


def _dem(path, res=4.0):
    from osgeo import gdal, osr

    nx, ny = 300, 200  # covers -100..1100 x -100..700
    xs = X0 - 100 + (np.arange(nx) + 0.5) * res
    ys = Y0 + 700 - (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    ds = gdal.GetDriverByName("GTiff").Create(path, nx, ny, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((X0 - 100, res, 0, Y0 + 700, 0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32614)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(_plane(X, Y).astype(np.float32))
    ds = None
    return path


def _read(path):
    lyr = QgsMeshLayer(path, "m", "mdal")
    qm = QgsMesh()
    if lyr.isValid():
        lyr.dataProvider().populateMesh(qm)
    return lyr, qm


@pytest.fixture(scope="module")
def project(qgis, tmp_path_factory):
    from telemac_toolkit.qgis_bridge import meshing as M

    d = tmp_path_factory.mktemp("meshproj")
    path = str(d / "proj.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), settings={"default_size": "30"})
    _add(lay["domain"], QgsGeometry.fromPolygonXY([_pts([(0, 0), (1000, 0), (1000, 600), (0, 600), (0, 0)])]), name="valley")
    _add(lay["islands"], QgsGeometry.fromPolygonXY([_pts([(150, 450), (250, 450), (250, 520), (150, 520), (150, 450)])]), name="isle")
    # a channel crossing the whole domain (its ends stick out), with an axis line
    for role, dy in (("left", 0), ("right", 60), ("axis", 30)):
        _add(lay["channels"], QgsGeometry.fromPolylineXY(_pts([(-20, 260 + dy), (300, 280 + dy), (700, 250 + dy), (1020, 270 + dy)])),
             channel="river", role=role, nodes_across=7, spacing=12.0)
    _add(lay["grid_zones"], QgsGeometry.fromPolygonXY([_pts([(600, 420), (720, 452), (700, 527), (580, 495), (600, 420)])]),
         name="bridge", dx=8.0)
    _add(lay["refinement_areas"], QgsGeometry.fromPolygonXY([_pts([(50, 50), (400, 50), (400, 200), (50, 200), (50, 50)])]), size=10.0)
    _add(lay["constraint_lines"], QgsGeometry.fromPolylineXY(_pts([(800, 50), (950, 200)])), name="dike", mode="breakline")
    _add(lay["constraint_lines"], QgsGeometry.fromPolylineXY(_pts([(450, 100), (450, 500)])), name="road", mode="guide", size=15.0)
    _add(lay["fixed_points"], QgsGeometry.fromPointXY(QgsPointXY(X0 + 900, Y0 + 500)), name="gauge")
    _add(lay["roughness_zones"], QgsGeometry.fromPolygonXY([_pts([(-50, 240), (1050, 240), (1050, 360), (-50, 360), (-50, 240)])]),
         name="river bed", material=2, friction=0.03)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(0, 20), (0, 240)])), name="inflow", type="Q", value=50.0)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(1000, 350), (1000, 580)])), name="outflow", type="544", value=4.0)
    return d, path, lay


def test_project_template(project):
    from telemac_toolkit.qgis_bridge import meshing as M

    d, path, lay = project
    assert set(lay) == set(M.LAYER_KEYS)
    for name, geom, fields, _t, _old in M.PROJECT_LAYERS:
        assert lay[name].isValid()
        assert [f.name() for f in lay[name].fields()][-len(fields):] == [f[0] for f in fields]
    group = [g for g in QgsProject.instance().layerTreeRoot().children() if g.name().startswith("Mesh project")][0]
    assert len(group.children()) == len(M.PROJECT_LAYERS)
    assert group.customProperty("telemac_toolkit/mesh_project") == path
    st = M.read_settings(path)
    assert st["default_size"] == "30" and st["format"] == "telemac"
    M.write_settings(path, {"format": "2dm", "new_key": "x"})
    st = M.read_settings(path)
    assert st["format"] == "2dm" and st["new_key"] == "x"
    M.write_settings(path, {"format": "telemac"})


def test_old_project_layers_still_read(qgis, tmp_path):
    """A 0.2.0 project (outline, density, hard_lines, soft_lines, ...) meshes the same way."""
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "old.gpkg")
    ctx = QgsProject.instance().transformContext()
    old = [("outline", "Polygon", "name:string"), ("density", "Polygon", "edge_len:double"),
           ("hard_lines", "LineString", "name:string"), ("soft_lines", "LineString", "edge_len:double")]
    for i, (name, geom, fld) in enumerate(old):
        mem = QgsVectorLayer("%s?crs=%s&field=%s" % (geom, CRS, fld), name, "memory")
        f = QgsFeature(mem.fields())
        if name == "outline":
            f.setGeometry(QgsGeometry.fromPolygonXY([_pts([(0, 0), (300, 0), (300, 200), (0, 200), (0, 0)])]))
        elif name == "density":
            f.setGeometry(QgsGeometry.fromPolygonXY([_pts([(10, 10), (100, 10), (100, 100), (10, 100), (10, 10)])]))
            f.setAttribute("edge_len", 5.0)
        elif name == "hard_lines":
            f.setGeometry(QgsGeometry.fromPolylineXY(_pts([(150, 20), (150, 180)])))
        else:
            f.setGeometry(QgsGeometry.fromPolylineXY(_pts([(200, 20), (280, 180)])))
        mem.dataProvider().addFeatures([f])
        opts = QgsVectorFileWriter.SaveVectorOptions()
        opts.driverName = "GPKG"
        opts.layerName = name
        opts.actionOnExistingFile = (QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile if i == 0
                                     else QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer)
        assert QgsVectorFileWriter.writeAsVectorFormatV3(mem, path, ctx, opts)[0] == QgsVectorFileWriter.WriterError.NoError
    lay = M.open_mesh_project(path)
    assert {"domain", "refinement_areas", "hard_lines", "soft_lines"} <= set(lay)
    b = M.build_mesh_input(lay, None, 20.0)
    assert any(i.startswith("Breaklines: 1") for i in b.info) and any(i.startswith("Guide lines: 1") for i in b.info)
    assert b.inp.size.minimum == 5.0
    res = M.run_mesher(b)
    assert res.stats["angle_min"] > 25.0
    x = res.x[np.abs(res.x - (X0 + 150)) < 1e-6]
    assert len(x) >= 2  # the breakline is a chain of nodes


def test_chain_every_format(project):
    """One drawing, every program: TELEMAC first, the others from the same mesh."""
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile, read_cli

    d, path, lay = project
    dem = QgsRasterLayer(_dem(str(d / "plane.tif")), "plane")
    src = M.bottom_sources_from([dem], crs=QgsCoordinateReferenceSystem(CRS))
    counts = {}
    for fmt in ("telemac", "2dm", "gmsh", "gid", "ugrid", "mike21", "adcirc", "hecras"):
        r = M.generate_and_export(lay, str(d / "out"), M.ExportOptions(fmt=fmt, quads=True), None, 30.0, 1.2, 8, src,
                                  default_friction=0.05)
        assert os.path.exists(r.path), fmt
        counts[fmt] = r
        if r.path.endswith((".slf", ".2dm", ".nc", ".mesh")):
            lyr, qm = _read(r.path)
            assert lyr.isValid(), fmt
            assert qm.vertexCount() == r.mesh.npoin, fmt
            sizes = sorted(set(len(qm.face(i)) for i in range(qm.faceCount())))
            assert sizes == ([3] if fmt == "telemac" else [3, 4]), fmt
            z = np.array([qm.vertex(i).z() for i in range(qm.vertexCount())]) if fmt != "telemac" else None
            if z is not None:
                xs = np.array([qm.vertex(i).x() for i in range(qm.vertexCount())])
                ys = np.array([qm.vertex(i).y() for i in range(qm.vertexCount())])
                assert np.allclose(z, _plane(xs, ys), atol=1e-3), fmt
    r = counts["telemac"]
    assert r.mesh.has_quads  # structured cells found again (TELEMAC gets them split)
    assert any("Zone river bed: material 2" in l for l in r.report)
    with SelafinFile(r.path) as s:
        names = [n.strip() for n in s.header.var_names]
        assert names == ["BOTTOM", "BOTTOM FRICTION"]
        fr = s.read("BOTTOM FRICTION", 0)
        assert set(np.round(np.unique(fr), 3)) == {0.03, 0.05}
        assert np.allclose(s.read("BOTTOM", 0), _plane(s.header.x, s.header.y), atol=2e-4)
        assert s.header.ipobo.max() == len(read_cli(str(d / "out.cli")))
    assert any(l.startswith("PRESCRIBED FLOWRATES") for l in r.cas)
    txt = open(str(d / "out.2dm")).read()
    ends = [l for l in txt.splitlines() if l.startswith("NS ") and int(l.split()[-1]) < 0]
    assert len(ends) == 2 and 'MAT 2 "river bed"' in txt  # one nodestring per open boundary
    assert os.path.exists(str(d / "out_nodestrings.csv"))
    assert os.path.exists(str(d / "out_gid.msh"))
    hec = str(d / "out.gpkg")
    for name, n in (("perimeter", 1), ("breaklines", 4), ("refinement_regions", 2), ("bc_lines", 2)):
        v = QgsVectorLayer("%s|layername=%s" % (hec, name), name, "ogr")
        assert v.isValid() and v.featureCount() == n, (name, v.featureCount())


def test_local_frames(project):
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, lay = project
    src = M.bottom_sources_from(formula="5.0")
    # TELEMAC, shift stored in the header: QGIS shows it in place, TELEMAC computes near (0, 0)
    r = M.generate_and_export(lay, str(d / "loc1"), M.ExportOptions(frame="lower_left", store_offset=True), None, 30.0,
                              bottom_sources=src)
    with SelafinFile(r.path) as s:
        assert s.header.origin == (X0 - 20.0, Y0)  # the channel sticks out 20 m to the west
        assert abs(s.header.x_stored.min()) < 1e-6 and abs(s.header.y_stored.min()) < 1e-6
        assert abs(s.header.x.min() - (X0 - 20.0)) < 1e-6
    assert os.path.exists(str(d / "loc1.prj"))
    # a chosen point at (0, 0) and a rotation: a truly local file plus a note to go back
    r = M.generate_and_export(lay, str(d / "loc2"), M.ExportOptions(fmt="2dm", frame="point", frame_point=(X0, Y0), angle=90.0),
                              None, 30.0, bottom_sources=src)
    lyr, qm = _read(r.path)
    xs = np.array([qm.vertex(i).x() for i in range(qm.vertexCount())])
    ys = np.array([qm.vertex(i).y() for i in range(qm.vertexCount())])
    # rotating the axes by 90 degrees: x_local = y, y_local = -x
    assert abs(ys.min() + 1020.0) < 1e-6 and abs(xs.min()) < 1e-6
    assert os.path.exists(str(d / "loc2_frame.txt")) and not os.path.exists(str(d / "loc2.prj"))


def test_flume_template_in_local_coordinates(qgis, tmp_path):
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile, read_cli
    from telemac_toolkit.telemac_io.cli import liquid_boundaries

    path = str(tmp_path / "flume.gpkg")
    res = processing.run("telemactoolkit:channeltemplate", {
        "NEW_PROJECT": path, "LENGTH": 12.0, "WIDTH": 0.6, "CELL": 0.05, "STRUCTURED": True, "BED_Z0": 0.3,
        "SLOPE": 0.002, "SECTION": 1, "BOTTOM_WIDTH": 0.3, "SIDE_SLOPE": 1.0, "BANK_HEIGHT": 0.1, "Q": 0.02, "H": 0.35})
    assert "clip(" in res["BED_FORMULA"]
    lay = M.open_mesh_project(path)
    assert M.is_local_crs(lay["domain"].crs())
    assert lay["domain"].featureCount() == 1 and lay["boundaries"].featureCount() == 2 and lay["grid_zones"].featureCount() == 1
    out = str(tmp_path / "flume.slf")
    r = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "OUTPUT": out, "OUTPUT_LIQ": str(tmp_path / "flume.liq"),
                                                       "LOAD": False})
    assert r["OUTPUT"] == out and r["OUTPUT_CLI"].endswith("flume.cli")
    with SelafinFile(out) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        assert h.npoin == 241 * 13  # 240 x 12 cells of 5 cm
        assert abs(h.x.min()) < 1e-9 and abs(h.y.min()) < 1e-9 and abs(h.x.max() - 12.0) < 1e-9
        axis = np.abs(h.y - 0.3) < 1e-9
        assert np.allclose(z[axis], 0.3 - 0.002 * h.x[axis])
        bank = np.abs(h.y) < 1e-9
        assert np.allclose(z[bank], 0.3 - 0.002 * h.x[bank] + 0.1)
    t = read_cli(r["OUTPUT_CLI"])
    names = sorted(t.comments[b["start"]] for b in liquid_boundaries(t))
    assert names == ["inlet", "outlet"]
    liq = open(str(tmp_path / "flume.liq")).read()
    assert "0.02" in liq and "0.35" in liq


def test_export_cut_quality_algorithms(project):
    import processing
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, lay = project
    r = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "FORMAT": 1, "QUADS": True,
                                                       "OUTPUT": str(d / "alg.2dm"), "LOAD": True, "FORMULA": "1.5"})
    assert r["OUTPUT"].endswith("alg.2dm")
    loaded = [l for l in QgsProject.instance().mapLayers().values() if isinstance(l, QgsMeshLayer)
              and l.source().strip('"').endswith("alg.2dm")]
    assert loaded and loaded[0].crs().authid() == CRS
    # convert the 2DM (with quads) to TELEMAC, with the boundary lines
    e = processing.run("telemactoolkit:exportmesh", {"INPUT": r["OUTPUT"], "BOUNDARIES": lay["boundaries"], "FORMAT": 0,
                                                     "OUTPUT": str(d / "conv.slf"), "LOAD": False})
    with SelafinFile(e["OUTPUT"]) as s:
        assert s.header.ndp == 3 and np.allclose(s.read("BOTTOM", 0), 1.5)
    assert os.path.exists(str(d / "conv.cli"))
    # cut the river out
    c = processing.run("telemactoolkit:cutmesh", {"INPUT": e["OUTPUT"], "POLYGONS": lay["roughness_zones"], "FORMAT": 0,
                                                  "OUTPUT": str(d / "river.slf"), "LOAD": False})
    with SelafinFile(c["OUTPUT"]) as s, SelafinFile(e["OUTPUT"]) as full:
        assert 0 < s.header.nelem < full.header.nelem
        assert s.header.y.min() >= Y0 + 230.0
    q = processing.run("telemactoolkit:meshquality", {"INPUT": r["OUTPUT"], "OUTPUT": "memory:"})["OUTPUT"]
    feats = list(q.getFeatures())
    assert len(feats) == q.featureCount() and {f["nodes"] for f in feats} == {3, 4}
    assert sum(f["flagged"] for f in feats) == 0


def test_regular_and_channel_algorithms(project):
    import processing
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, lay = project
    gr = processing.run("telemactoolkit:regulargrid", {"ORIGIN": "0,0", "DX": 0.1, "DY": 0.1, "NX": 30, "NY": 5, "CRS": CRS,
                                                        "FORMULA": "1 - 0.01*x", "OUTPUT": str(d / "gr.slf"), "LOAD": False})
    with SelafinFile(gr["OUTPUT"]) as s:
        assert s.header.npoin == 31 * 6 and s.header.nelem == 300
        assert np.allclose(s.read("BOTTOM", 0), 1 - 0.01 * s.header.x)
    q = processing.run("telemactoolkit:regulargrid", {"ORIGIN": "0,0", "DX": 0.1, "DY": 0.1, "NX": 30, "NY": 5, "CRS": CRS,
                                                       "FORMAT": 1, "QUADS": True, "OUTPUT": str(d / "gr.2dm"), "LOAD": False})
    assert open(q["OUTPUT"]).read().count("E4Q") == 150
    ch = processing.run("telemactoolkit:channelmesh", {"CHANNELS": lay["channels"], "OUTPUT": str(d / "ch.slf"), "LOAD": False})
    with SelafinFile(ch["OUTPUT"]) as s:
        assert s.header.npoin == 7 * 88 and s.header.nelem == 2 * 6 * 87


def test_create_project_algorithm(qgis, tmp_path):
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M

    out = str(tmp_path / "new.gpkg")
    n0 = len(QgsProject.instance().mapLayers())
    processing.run("telemactoolkit:createmeshproject", {"LOCAL": True, "DEFAULT_SIZE": 0.1, "OUTPUT": out})
    assert os.path.exists(out)
    assert len(QgsProject.instance().mapLayers()) == n0 + len(M.PROJECT_LAYERS)
    assert M.read_settings(out)["default_size"] == "0.1"
    assert M.is_local_crs(M.open_mesh_project(out)["domain"].crs())


def test_workflow_panel(project):
    from qgis.gui import QgsMapCanvas
    from qgis.PyQt.QtWidgets import QMainWindow
    from telemac_toolkit.gui.mesh_panel import MeshPanel
    from telemac_toolkit.qgis_bridge import meshing as M

    class Iface:
        def __init__(self):
            self.win = QMainWindow()
            self.canvas = QgsMapCanvas(self.win)
            self.canvas.setDestinationCrs(QgsCoordinateReferenceSystem(CRS))

        def mapCanvas(self):
            return self.canvas

    d, path, lay = project
    iface = Iface()
    panel = MeshPanel(iface)
    assert panel.autodetect() is not None  # a mesh project on the map is picked up
    panel.set_project(path)
    assert panel.table.rowCount() == len(M.LAYER_KEYS)
    assert panel.size.value() == 30.0
    panel.formula.setText("2.0")
    panel.fmt.setCurrentIndex(list(M.FORMATS).index("telemac"))
    panel.output.setFilePath(str(d / "panel.slf"))
    assert panel.check_drawing()
    res = panel.generate(wait=True)
    assert res is not None and os.path.exists(str(d / "panel.slf")) and os.path.exists(str(d / "panel.cli"))
    assert "Done" in panel.log.toPlainText()
    st = M.read_settings(path)
    assert st["bottom_formula"] == "2.0" and st["mesh_file"].endswith("panel.slf")
    pt = panel.find(False, 1)
    assert pt is not None and "Node 1" in panel.found.text()
    assert panel.find(True, 10 ** 8) is None
    q = panel.quality()
    assert q is not None and q.featureCount() > 0
    # switching the format renames the output file
    panel.fmt.setCurrentIndex(list(M.FORMATS).index("gid"))
    assert panel.output.filePath().endswith("panel_gid.msh")
    panel.cleanup()
