"""Checking the drawing of a mesh project: problems on the map, safe fixes, the Processing tool and the panel.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_drawcheck_qgis.py
"""

import pytest
from qgis.core import QgsCoordinateReferenceSystem, QgsGeometry, QgsPointXY, QgsProject

from test_mesher_qgis import CRS, X0, Y0, _add
from test_river_qgis import qgis  # noqa: F401  (fixture)


def _pts(c):
    return [QgsPointXY(X0 + a, Y0 + b) for a, b in c]


@pytest.fixture()
def project(qgis, tmp_path):  # noqa: F811
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "check.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), add_to_project=False, settings={"default_size": "30"})
    # a near-duplicate vertex (2 cm) on the bottom edge and a spike on the right edge
    ring = [(0, 0), (500, 0), (500.02, 0), (1000, 0), (1000, 300), (1004, 300.2), (1000, 300.4), (1000, 600), (0, 600), (0, 0)]
    _add(lay["domain"], QgsGeometry.fromPolygonXY([_pts(ring)]), name="valley")
    # an island 2 m from the top edge (elements of 30 m there): slivers
    _add(lay["islands"], QgsGeometry.fromPolygonXY([_pts([(200, 498), (300, 498), (300, 598), (200, 598), (200, 498)])]),
         name="pier")
    # a breakline ending 3 m short of the left edge
    _add(lay["constraint_lines"], QgsGeometry.fromPolylineXY(_pts([(3, 200), (400, 250)])), name="levee", mode="breakline")
    # a bow-tie roughness zone (not a valid polygon)
    _add(lay["roughness_zones"], QgsGeometry.fromPolygonXY([_pts([(600, 100), (700, 200), (700, 100), (600, 200), (600, 100)])]),
         name="bowtie", material=2, friction=0.03)
    return path


def _check(path):
    from telemac_toolkit.qgis_bridge import drawing as D
    from telemac_toolkit.qgis_bridge import meshing as M

    lay = M.open_mesh_project(path)
    b = M.build_mesh_input(lay, None, 30.0, max_nodes=None)
    return lay, b, D.check_drawing(lay, b.inp.size, b.crs)


def test_problems_found_and_the_safe_ones_fixed(project):
    from telemac_toolkit.qgis_bridge import drawing as D

    lay, b, chk = _check(project)
    codes = {i.problem.code for i in chk.issues}
    assert {"duplicate", "spike", "gap", "near", "invalid"} <= codes, [i.text for i in chk.issues]
    gap = [i for i in chk.issues if i.problem.code == "gap"][0]
    assert "pier" in gap.problem.message and abs(gap.problem.distance - 2.0) < 1.0
    near = [i for i in chk.issues if i.problem.code == "near"][0]
    assert "levee" in near.label and abs(near.problem.distance - 3.0) < 1e-6
    assert "Drawing check:" in chk.summary()[0]
    pts = D.issues_layer(chk)
    assert pts.featureCount() == len(chk.issues) and pts.crs().authid() == CRS
    report = D.fix_drawing(lay, chk, b.inp.size, commit=True)
    assert any("snapped" in r for r in report) and any("removed" in r for r in report), report
    lay2, b2, chk2 = _check(project)
    left = {i.problem.code for i in chk2.issues}
    assert left == {"gap"}, [i.text for i in chk2.issues]  # the island is the modeller's call
    # the breakline now ends on the outline, and the domain lost only its two bad vertices
    line = next(lay2["constraint_lines"].getFeatures()).geometry().asPolyline()
    assert abs(line[0].x() - X0) < 1e-6
    assert len(next(lay2["domain"].getFeatures()).geometry().asPolygon()[0]) == 10 - 2
    assert next(lay2["roughness_zones"].getFeatures()).geometry().isGeosValid()


def test_processing_tool_and_panel(project):
    import processing
    from qgis.gui import QgsMapCanvas
    from qgis.PyQt.QtWidgets import QMainWindow
    from telemac_toolkit.gui.mesh_panel import MeshPanel
    from telemac_toolkit.qgis_bridge import drawing as D

    r = processing.run("telemactoolkit:checkdrawing", {"PROJECT": project, "OUTPUT": "memory:"})
    assert r["PROBLEMS"] >= 5 and r["FIXABLE"] >= 4 and r["OUTPUT"].featureCount() == r["PROBLEMS"]

    class Iface:
        def __init__(self):
            self.win = QMainWindow()
            self.canvas = QgsMapCanvas(self.win)
            self.canvas.setDestinationCrs(QgsCoordinateReferenceSystem(CRS))

        def mapCanvas(self):
            return self.canvas

    panel = MeshPanel(Iface())
    panel.set_project(project, add_to_map=True)
    assert panel.check_drawing()
    checks = [lyr for lyr in QgsProject.instance().mapLayers().values()
              if lyr.customProperty("telemac_toolkit/drawing_check") == "1"]
    assert len(checks) == 1 and checks[0].featureCount() >= 5 and panel.btn_fix.isEnabled()
    panel.fix_drawing()
    log = panel.log.toPlainText()
    assert "snapped onto the outline" in log and "edit buffers" in log
    checks = [lyr for lyr in QgsProject.instance().mapLayers().values()
              if lyr.customProperty("telemac_toolkit/drawing_check") == "1"]
    assert len(checks) == 1 and {f["problem"] for f in checks[0].getFeatures()} == {"gap"}
    assert not panel.btn_fix.isEnabled()
    # mesh type: the choice drives the method and the quadrilateral option
    panel.fmt.setCurrentIndex(panel.fmt.findData("2dm"))
    panel.mesh_type.setCurrentIndex(panel.mesh_type.findData("rectangles"))
    assert panel.quads.isChecked() and panel.method.currentData() != "gmsh_quads"
    assert "Nothing structured" in panel.mesh_type_note.text()
    panel.mesh_type.setCurrentIndex(panel.mesh_type.findData("triangles"))
    assert not panel.quads.isChecked()
    panel.quads.setChecked(True)
    assert panel.mesh_type.currentData() == "rectangles"
    panel.fmt.setCurrentIndex(panel.fmt.findData("telemac"))
    assert "triangles only" in panel.mesh_type_note.text()
    assert panel.settings()["mesh_type"] == "rectangles"
    for lyr in list(QgsProject.instance().mapLayers().values()):
        if lyr.customProperty("telemac_toolkit/drawing_check") == "1" or D.CHECK_LAYER_NAME in lyr.name():
            QgsProject.instance().removeMapLayer(lyr.id())
    panel.cleanup()
