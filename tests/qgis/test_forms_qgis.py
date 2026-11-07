"""Project layer forms and hover tips, boundary tracer and custom codes in the .cli, CRS without an EPSG code,
sizes suggested from a DEM, and the panel's automatic snapping, in a real (headless) QGIS.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_forms_qgis.py
"""

import sqlite3

import numpy as np
import pytest
from qgis.core import (QgsCoordinateReferenceSystem, QgsExpression, QgsExpressionContext, QgsExpressionContextUtils,
                       QgsGeometry, QgsPointXY, QgsProject, QgsVectorLayer)

from test_mesher_qgis import CRS, X0, Y0, _add
from test_river_qgis import dem, qgis  # noqa: F401  (fixtures)

INEGI_WKT = ('PROJCS["UTM_Zone_14_Northern_Hemisphere",GEOGCS["GCS_ITRF_1992",DATUM["D_ITRF_1992",SPHEROID["GRS_1980",'
             '6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
             'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],'
             'PARAMETER["Central_Meridian",-99.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],'
             'UNIT["Meter",1.0]]')


def _pts(c):
    return [QgsPointXY(X0 + a, Y0 + b) for a, b in c]


def _tip(layer, feature) -> str:
    ctx = QgsExpressionContext(QgsExpressionContextUtils.globalProjectLayerScopes(layer))
    ctx.setFeature(feature)
    return QgsExpression.replaceExpressionText(layer.mapTipTemplate(), ctx)


@pytest.fixture()
def project(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "forms.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), add_to_project=False, settings={"default_size": "50"})
    _add(lay["domain"], QgsGeometry.fromPolygonXY([_pts([(0, 0), (1000, 0), (1000, 500), (0, 500), (0, 0)])]), name="valley", size=10.0)
    _add(lay["refinement_areas"], QgsGeometry.fromPolygonXY([_pts([(100, 100), (200, 100), (200, 200), (100, 200), (100, 100)])]),
         name="bridge", size=2.0)
    _add(lay["constraint_lines"], QgsGeometry.fromPolylineXY(_pts([(300, 50), (300, 450)])), name="levee", mode="guide", size=5.0)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(0, 50), (0, 450)])), name="in", type="Q", value=25.0)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(1000, 50), (1000, 200)])), name="out", type="H", value=7.5)
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(1000, 260), (1000, 450)])), name="spill", type="custom",
         codes="5 5 5 4")
    _add(lay["boundaries"], QgsGeometry.fromPolylineXY(_pts([(200, 500), (800, 500)])), name="side", type="Q", value=1.0,
         tracer="free")
    return path


def test_hover_tips_and_forms(project):
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(project)
    dom = next(lay["domain"].getFeatures())
    tip = _tip(lay["domain"], dom)
    assert "valley" in tip and "Area" in tip and "m²" in tip and "50.00 ha" in tip and "10.00 m" in tip
    assert "about 5,775" in tip  # 1.155 x 500 000 m² / 10²
    ref = _tip(lay["refinement_areas"], next(lay["refinement_areas"].getFeatures()))
    assert "bridge" in ref and "Element size" in ref and "about 2,888" in ref  # 1.155 x 10 000 / 2²
    line = _tip(lay["constraint_lines"], next(lay["constraint_lines"].getFeatures()))
    assert "guide line" in line and "Length" in line and "400" in line
    tips = {f["name"]: _tip(lay["boundaries"], f) for f in lay["boundaries"].getFeatures()}
    assert "Open boundary with prescribed Q (4 5 5)" in tips["in"] and "25" in tips["in"] and "m³/s" in tips["in"]
    assert "Open boundary with prescribed H (5 4 4)" in tips["out"] and "from the boundary code" in tips["out"]
    assert "Custom" in tips["spill"] and "5 5 5 4" in tips["spill"]
    assert "Open boundary with free tracer" in tips["side"]
    # drop-down lists as Blue Kenue's CONLIM editor
    b = lay["boundaries"]
    w = b.editorWidgetSetup(b.fields().indexOf("type"))
    assert w.type() == "ValueMap"
    labels = [list(d.keys())[0] for d in w.config()["map"]]
    assert labels[:3] == ["Closed boundary (wall)", "Open boundary with prescribed Q", "Open boundary with prescribed H"]
    assert b.editorWidgetSetup(b.fields().indexOf("tracer")).type() == "ValueMap"
    assert lay["constraint_lines"].editorWidgetSetup(lay["constraint_lines"].fields().indexOf("mode")).type() == "ValueMap"
    assert b.attributeDisplayName(b.fields().indexOf("value")).startswith("Value (m³/s for Q")


def test_tracer_and_custom_codes_in_the_cli(project, tmp_path):
    import processing
    from telemac_toolkit.telemac_io import read_cli

    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": project, "OUTPUT": str(tmp_path / "m.slf"), "LOAD": False})
    t = read_cli(g["OUTPUT_CLI"])
    r = t.records
    by = {}
    for k, c in enumerate(t.comments):
        if c:
            by.setdefault(c.split()[-1], set()).add((int(r["lihbor"][k]), int(r["liubor"][k]), int(r["livbor"][k]), int(r["litbor"][k])))
    assert by["in"] == {(4, 5, 5, 5)}  # an inflow: prescribed tracer
    assert by["out"] == {(5, 4, 4, 4)}  # an outflow: free tracer
    assert by["spill"] == {(5, 5, 5, 4)}  # custom codes, tracer given
    assert by["side"] == {(4, 5, 5, 4)}  # tracer set to free
    walls = [k for k, c in enumerate(t.comments) if not c]
    assert walls and all(int(r["litbor"][k]) == 2 for k in walls)


def test_old_project_gets_the_new_fields(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "old.gpkg")
    M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), add_to_project=False)
    # a project made before 0.7: no tracer, codes or across_spacing columns
    con = sqlite3.connect(path)
    for table, col in (("boundaries", "tracer"), ("boundaries", "codes"), ("channels", "across_spacing")):
        con.execute("ALTER TABLE %s DROP COLUMN %s" % (table, col))
    con.commit()
    con.close()
    lay = M.open_mesh_project(path)
    for key, col in (("boundaries", "tracer"), ("boundaries", "codes"), ("channels", "across_spacing")):
        assert lay[key].fields().indexOf(col) >= 0, (key, col)
    assert QgsVectorLayer("%s|layername=boundaries" % path, "b", "ogr").fields().indexOf("codes") >= 0  # written to the file


def test_crs_without_a_code(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    custom = QgsCoordinateReferenceSystem.fromWkt(INEGI_WKT)
    assert custom.isValid() and custom.authid() == ""  # listed as "name []" by QGIS
    cand = M.crs_with_code(custom)
    assert cand is not None and cand.authid() == "EPSG:4487"  # Mexico ITRF92 / UTM zone 14N
    assert M.crs_with_code(QgsCoordinateReferenceSystem(CRS)) is None  # it has a code already
    path = str(tmp_path / "custom.gpkg")
    lay = M.create_mesh_project(path, custom, add_to_project=False)
    _add(lay["domain"], QgsGeometry.fromPolygonXY([_pts([(0, 0), (100, 0), (100, 50), (0, 0)])]), name="d")
    assert M.open_mesh_project(path)["domain"].crs().authid() == ""
    assert M.assign_project_crs(path, cand) >= 10
    lay = M.open_mesh_project(path)
    assert all(v.crs().authid() == "EPSG:4487" for v in lay.values())
    g = next(lay["domain"].getFeatures()).geometry()
    assert g.boundingBox().xMinimum() == pytest.approx(X0) and g.boundingBox().yMaximum() == pytest.approx(Y0 + 50)


def test_sizes_suggested_from_the_dem(dem, tmp_path):  # noqa: F811
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M

    # without a project: the DEM's extent
    r = processing.run("telemactoolkit:suggestsizes", {"DEM": dem})
    assert r["SMALLEST"] == pytest.approx(4.0) and r["TOLERANCE"] == pytest.approx(0.1)  # 2 m pixels
    assert 8.0 <= r["DEFAULT_SIZE"] <= 200.0 and r["NODES"] > 1000 and 0.1 < r["TIME_STEP"] < 5
    # the river project: its channel width gives the spacing across and along
    xs = np.linspace(0, 2000, 60)
    cl = QgsVectorLayer("LineString?crs=%s" % CRS, "cl", "memory")
    _add(cl, QgsGeometry.fromPolylineXY(_pts(np.stack([xs, 70 * np.sin(2 * np.pi * xs / 900)], 1))))
    path = str(tmp_path / "river.gpkg")
    processing.run("telemactoolkit:riverreach", {"NEW_PROJECT": path, "CENTRELINE": cl, "WIDTH_MODE": 1, "DEM": dem,
                                                 "BANKFULL_DEPTH": 2.5, "MAX_WIDTH": 120.0, "SPACING": 8.0, "FLOODPLAIN": 150.0})
    r = processing.run("telemactoolkit:suggestsizes", {"DEM": dem, "PROJECT": path, "APPLY": True})
    assert 4.0 <= r["CHANNEL_ACROSS"] <= 8.0 and r["CHANNEL_ALONG"] >= r["CHANNEL_ACROSS"]  # banks about 44-52 m apart
    st = M.read_settings(path)
    assert float(st["default_size"]) == pytest.approx(r["DEFAULT_SIZE"]) and float(st["terrain_tolerance"]) == pytest.approx(0.1)


def test_panel_sets_snapping_when_a_project_opens(project, dem):  # noqa: F811
    from qgis.gui import QgsMapCanvas
    from qgis.PyQt.QtWidgets import QAction, QDockWidget, QMainWindow
    from telemac_toolkit.gui.mesh_panel import MeshPanel

    class Iface:
        def __init__(self):
            self.win = QMainWindow()
            self.canvas = QgsMapCanvas(self.win)
            self.tracing = QAction("tracing", self.win)
            self.tracing.setObjectName("mActionEnableTracing")
            self.tracing.setCheckable(True)
            self.win.addAction(self.tracing)
            self.tips = QAction("tips", self.win)
            self.tips.setCheckable(True)
            self.cad = QDockWidget("cad", self.win)
            self.cad.setObjectName("AdvancedDigitizingTools")
            self.cad.hide()

        def mainWindow(self):
            return self.win

        def mapCanvas(self):
            return self.canvas

        def actionMapTips(self):
            return self.tips

    QgsProject.instance().snappingConfig()
    iface = Iface()
    panel = MeshPanel(iface)
    panel.auto_snap.setChecked(True)
    panel.set_project(project)
    cfg = QgsProject.instance().snappingConfig()
    assert cfg.enabled() and cfg.intersectionSnapping() and QgsProject.instance().topologicalEditing()
    assert iface.tracing.isChecked() and iface.tips.isChecked() and not iface.cad.isHidden()
    # a DEM on the map: the suggestion fills the sizes when applied
    QgsProject.instance().addMapLayer(dem, False)
    sug = panel.suggest_sizes(dem.id(), quiet=True)
    assert sug is not None and sug.smallest == pytest.approx(4.0)
    assert any("Smallest size" in line for line in panel.log.toPlainText().splitlines())


def test_overwriting_a_mesh_loaded_in_qgis(project, tmp_path):
    """Windows cannot overwrite a file QGIS has open: Generate mesh removes its layer first (the new mesh replaces it)."""
    import processing
    from qgis.core import QgsMeshLayer
    from telemac_toolkit.qgis_bridge import meshing as M

    out = str(tmp_path / "loaded.slf")
    processing.run("telemactoolkit:generatemesh", {"PROJECT": project, "OUTPUT": out, "LOAD": False})
    lyr = QgsMeshLayer(out, "loaded", "mdal")
    assert lyr.isValid()
    QgsProject.instance().addMapLayer(lyr)
    assert [x.name() for x in M.layers_using([out])] == ["loaded"]
    assert M.free_name(out).endswith("loaded_2.slf")
    processing.run("telemactoolkit:generatemesh", {"PROJECT": project, "OUTPUT": out, "LOAD": False})
    assert not M.layers_using([out])  # removed, then written again


def test_unwritable_file_says_why(tmp_path):
    from telemac_toolkit.telemac_io.selafin import SelafinWriter

    bad = str(tmp_path / "no_such_folder" / "m.slf")
    with pytest.raises(OSError, match="cannot be written"):
        SelafinWriter(bad, ["BOTTOM"], ["M"], np.zeros(3), np.array([0.0, 1.0, 0.0]), np.array([[0, 1, 2]]))


def test_local_crs_has_a_name_and_a_code(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    c = M.local_crs()
    assert c.authid().startswith("USER:") and M.is_local_crs(c)
    path = str(tmp_path / "lab.gpkg")
    M.create_mesh_project(path, c, add_to_project=False)
    assert M.open_mesh_project(path)["domain"].crs().authid() == c.authid()


def test_refinement_area_without_size_says_it_is_ignored(project):
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(project)
    _add(lay["refinement_areas"], QgsGeometry.fromPolygonXY([_pts([(500, 100), (600, 100), (600, 200), (500, 100)])]), name="nosize")
    lay = M.open_mesh_project(project)
    f = next(x for x in lay["refinement_areas"].getFeatures() if x["name"] == "nosize")
    assert "not set: this area is ignored" in _tip(lay["refinement_areas"], f)


def test_hover_tips_readable_on_a_dark_theme(project):
    """QGIS leaves the tip text black on a dark frame: the tip carries its own colours, from the theme."""
    from qgis.PyQt.QtGui import QColor, QPalette
    from qgis.PyQt.QtWidgets import QApplication
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(project)
    old = QApplication.palette()
    dark = QPalette(old)
    dark.setColor(QPalette.ColorRole.Base, QColor("#202124"))
    dark.setColor(QPalette.ColorRole.Text, QColor("#eeeeee"))
    try:
        QApplication.setPalette(dark)
        M.refresh_map_tips(lay)
        t = lay["refinement_areas"].mapTipTemplate()
        assert "background-color:#202124" in t and "color:#eeeeee" in t
    finally:
        QApplication.setPalette(old)
    M.refresh_map_tips(lay)
    bg, fg = M.tip_colors()
    assert "color:%s" % fg in lay["domain"].mapTipTemplate()
    tip = _tip(lay["domain"], next(lay["domain"].getFeatures()))
    assert "Area" in tip and "</div>" in tip
