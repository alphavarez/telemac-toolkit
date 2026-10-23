"""0.4.0 in a real (headless) QGIS: arrows, CRS choice and visibility, calculator, integrals, node editing,
curved channel templates.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis
"""

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_plugins = os.path.join(sys.prefix, "share", "qgis", "python", "plugins")
if os.path.isdir(_plugins) and _plugins not in sys.path:
    sys.path.append(_plugins)

from qgis.core import (  # noqa: E402
    QgsCoordinateReferenceSystem,
    QgsFeature,
    QgsGeometry,
    QgsMapRendererSequentialJob,
    QgsMapSettings,
    QgsMeshLayer,
    QgsPointXY,
    QgsProject,
    QgsRectangle,
    QgsVectorLayer,
)
from qgis.PyQt.QtCore import QSize  # noqa: E402
from qgis.PyQt.QtGui import QColor  # noqa: E402

from test_mesher_qgis import _app  # noqa: E402

from telemac_toolkit.mesher.regular import regular_grid_mesh  # noqa: E402
from telemac_toolkit.telemac_io import SelafinFile, write_selafin  # noqa: E402


@pytest.fixture(scope="module")
def qgis():
    return _app()


def _uniform_result(path, x0=0.0, y0=0.0, n=20, d=1.0, speed=1.0):
    x, y, ik = regular_grid_mesh(x0, y0, d, d, n, n, 0.0, "alternate")
    write_selafin(path, x, y, ik, {"VELOCITY U": np.full_like(x, speed), "VELOCITY V": np.zeros_like(x),
                                   "WATER DEPTH": np.full_like(x, 0.5)},
                  times=[0.0], units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M"})
    return path


def _render(layer, extent, size=600):
    ms = QgsMapSettings()
    ms.setLayers([layer])
    ms.setOutputSize(QSize(size, size))
    ms.setExtent(extent)
    ms.setDestinationCrs(layer.crs())
    ms.setOutputDpi(96)
    ms.setBackgroundColor(QColor("white"))
    job = QgsMapRendererSequentialJob(ms)
    job.start()
    job.waitForFinished()
    img = job.renderedImage()
    w, h = img.width(), img.height()
    ptr = img.constBits()
    try:
        ptr.setsize(img.sizeInBytes())
    except AttributeError:
        pass
    arr = np.frombuffer(bytes(ptr), dtype=np.uint8).reshape(h, img.bytesPerLine())[:, : 4 * w].reshape(h, w, 4)
    return arr[:, :, :3].mean(axis=2) < 100  # dark pixels


def test_arrow_style_fits_the_grid(qgis, tmp_path):
    from telemac_toolkit.qgis_bridge import styles
    from telemac_toolkit.qgis_bridge.results import load_result
    from qgis.core import QgsMeshRendererVectorArrowSettings

    p = _uniform_result(str(tmp_path / "u.slf"))
    rs = load_result(p, crs=QgsCoordinateReferenceSystem("EPSG:32614"), add_to_project=False)
    st = styles.ArrowStyle(placement="grid", spacing=60, length="speed", longest=0.0, min_value=0.0)
    styles.apply_arrow_style(rs, st)
    g = next(iter(rs.group_of_pair.values()))
    vs = rs.layer.rendererSettings().vectorSettings(g)
    assert vs.isOnUserDefinedGrid() and vs.userGridCellWidth() == 60
    assert vs.arrowSettings().shaftLengthMethod() == QgsMeshRendererVectorArrowSettings.ArrowScalingMethod.Scaled
    # the fastest arrow (1 m/s) is 90 % of the spacing: 54 px
    assert vs.arrowSettings().scaleFactor() == pytest.approx(54.0, rel=1e-6)
    s = rs.layer.rendererSettings()
    s.setActiveScalarDatasetGroup(-1)
    s.setActiveVectorDatasetGroup(g)
    rs.layer.setRendererSettings(s)
    dark = _render(rs.layer, QgsRectangle(0, 0, 20, 20))
    # one row of arrows: the longest dark run along a row stays below the spacing (no overlap)
    runs = []
    for row in dark:
        cols = np.flatnonzero(row)
        if cols.size:
            splits = np.split(cols, np.flatnonzero(np.diff(cols) > 1) + 1)
            runs += [c[-1] - c[0] + 1 for c in splits]
    assert max(runs) < 60 and max(runs) > 40
    # fixed scale and nodes
    st2 = styles.ArrowStyle(placement="nodes", length="scale", scale=5.0)
    styles.apply_arrow_style(rs, st2)
    vs = rs.layer.rendererSettings().vectorSettings(g)
    assert not vs.isOnUserDefinedGrid()
    assert vs.arrowSettings().scaleFactor() == pytest.approx(5.0 * styles.PX_PER_MM)
    assert st2.length_of(0.5, 1.0) == pytest.approx(2.5)


def _arrow_direction(dark, tail):
    """Unit direction (map axes: x east, y north) of the arrow drawn from pixel ``tail``.

    QGIS starts each grid arrow at its grid point: the arrow is the dark component that
    touches the tail, and its farthest pixel from the tail is the tip.
    """
    r0, c0 = tail
    near = [(r, c) for r in range(r0 - 3, r0 + 4) for c in range(c0 - 3, c0 + 4) if dark[r, c]]
    assert near, "no arrow at the grid point"
    seen, todo = set(near), list(near)
    while todo:
        r, c = todo.pop()
        for rr in (r - 1, r, r + 1):
            for cc in (c - 1, c, c + 1):
                if (rr, cc) not in seen and 0 <= rr < dark.shape[0] and 0 <= cc < dark.shape[1] and dark[rr, cc]:
                    seen.add((rr, cc))
                    todo.append((rr, cc))
    r, c = max(seen, key=lambda p: (p[0] - r0) ** 2 + (p[1] - c0) ** 2)
    dx, dy = c - c0, -(r - r0)  # image rows grow downwards
    n = math.hypot(dx, dy)
    assert n > 20, "arrow too short"
    return dx / n, dy / n


@pytest.mark.parametrize("lang", ["en", "fr"])
def test_arrows_follow_the_vector_shown(qgis, tmp_path, lang):
    from telemac_toolkit.qgis_bridge import styles
    from telemac_toolkit.qgis_bridge.results import load_result

    x, y, ik = regular_grid_mesh(0.0, 0.0, 1.0, 1.0, 10, 10, 0.0, "alternate")
    one = np.ones_like(x)
    # velocity to the east, unit discharge to the north, wind to the south-west
    if lang == "en":
        names = ("VELOCITY U", "VELOCITY V", "WATER DEPTH", "FLOWRATE ALONG X", "FLOWRATE ALONG Y",
                 "WIND ALONG X", "WIND ALONG Y")
    else:
        names = ("VITESSE U", "VITESSE V", "HAUTEUR D'EAU", "DEBIT SUIVANT X", "DEBIT SUIVANT Y", "VENT X", "VENT Y")
    vals = (1.0, 0.0, 0.3, 0.0, 0.3, -3.0, -3.0)
    p = str(tmp_path / ("v_%s.slf" % lang))
    write_selafin(p, x, y, ik, {n: v * one for n, v in zip(names, vals)})
    rs = load_result(p, crs=QgsCoordinateReferenceSystem("EPSG:32614"), add_to_project=False)
    vel, q = rs.group_of_pair[(0, 1)], rs.group_of_pair[(3, 4)]
    if lang == "en":
        wind = rs.group_of_pair[(5, 6)]
    else:  # QGIS does not pair "VENT X/Y": no arrows for it, the plugin must not invent a group
        assert (5, 6) not in rs.group_of_pair
    assert len(set(rs.group_of_pair.values())) == len(rs.group_of_pair)
    # opening: depth coloured, velocity arrows
    assert styles.arrow_group(rs) == vel
    # a scalar keeps the arrows; a vector moves them
    styles.set_active_variable(rs, 2)
    assert styles.arrow_group(rs) == vel
    styles.set_active_variable(rs, (3, 4))
    assert styles.arrow_group(rs) == q
    styles.set_active_variable(rs, 2)
    assert styles.arrow_group(rs) == q
    styles.set_arrows_visible(rs, False)
    assert styles.arrow_group(rs) is None
    styles.set_active_variable(rs, (0, 1))  # arrows off stay off
    assert styles.arrow_group(rs) is None
    styles.set_arrows_visible(rs, True)
    assert styles.arrow_group(rs) == q  # the vector shown last
    # each vector is drawn in its own direction
    st = styles.ArrowStyle(placement="grid", spacing=100, length="fixed", longest=15.0, min_value=0.0, width=1.0)
    styles.apply_arrow_style(rs, st)
    expected = [(vel, (1.0, 0.0)), (q, (0.0, 1.0))] + ([(wind, (-0.7071, -0.7071))] if lang == "en" else [])
    for g, (ex, ey) in expected:
        styles.set_arrow_group(rs, g)
        s = rs.layer.rendererSettings()
        s.setActiveScalarDatasetGroup(-1)
        rs.layer.setRendererSettings(s)
        dx, dy = _arrow_direction(_render(rs.layer, QgsRectangle(1, 1, 9, 9), size=300), tail=(100, 100))
        assert dx * ex + dy * ey > 0.9, (g, dx, dy)


def test_crs_suggestions(qgis, tmp_path):
    from telemac_toolkit.qgis_bridge.crs import Extent, is_local_crs, plausible, suggest_crs

    prj = QgsProject.instance()
    prj.clear()
    utm = QgsCoordinateReferenceSystem("EPSG:32614")
    lab = Extent(0.0, 0.0, 20.0, 1.0)
    field = Extent(512000.0, 2155000.0, 514000.0, 2156000.0)
    assert plausible(utm, field) and not plausible(utm, lab)
    assert not plausible(QgsCoordinateReferenceSystem("EPSG:4326"), field)
    p = str(tmp_path / "a.slf")
    # empty project (its CRS is only a default): a laboratory file gets the local grid as suggestion, asked
    prj.setCrs(QgsCoordinateReferenceSystem("EPSG:4326"))
    c = suggest_crs(p, lab, prj)
    assert c.ask and is_local_crs(c.crs)
    # an empty project with a projected CRS: a file that fits it takes it
    prj.setCrs(utm)
    assert suggest_crs(p, field, prj).ask is False
    # a map project in UTM: a UTM file takes it silently, a laboratory one is asked about
    v = QgsVectorLayer("Point?crs=EPSG:32614", "v", "memory")
    prj.addMapLayer(v)
    c = suggest_crs(p, field, prj)
    assert not c.ask and c.crs == utm
    c = suggest_crs(p, lab, prj)
    assert c.ask and is_local_crs(c.crs)
    # a .prj next to the file wins
    with open(str(tmp_path / "a.prj"), "w") as f:
        f.write(utm.toWkt())
    c = suggest_crs(p, lab, prj)
    assert not c.ask and c.crs == utm
    prj.clear()


def test_show_layer_zooms_and_switches(qgis, tmp_path):
    from qgis.gui import QgsMapCanvas

    from telemac_toolkit.qgis_bridge.crs import is_local_crs, local_crs, show_layer

    prj = QgsProject.instance()
    prj.clear()
    canvas = QgsMapCanvas()
    canvas.resize(600, 400)
    # a mesh project in local coordinates is on screen
    lab = QgsVectorLayer("Polygon?crs=EPSG:4326", "lab", "memory")
    lab.setCrs(local_crs())
    prj.addMapLayer(lab)
    prj.setCrs(local_crs())
    canvas.setDestinationCrs(local_crs())
    canvas.setExtent(QgsRectangle(-1, -1, 21, 2))
    # a result in UTM: the project switches and the map zooms to it
    p = _uniform_result(str(tmp_path / "f.slf"), 512000.0, 2155000.0)
    lyr = QgsMeshLayer(p, "f", "mdal")
    lyr.setCrs(QgsCoordinateReferenceSystem("EPSG:32614"))
    prj.addMapLayer(lyr)
    done = show_layer(lyr, canvas, prj)
    assert "CRS" in done and "zoomed" in done
    assert prj.crs().authid() == "EPSG:32614"
    assert canvas.extent().intersects(QgsRectangle(512000, 2155000, 512020, 2155020))
    # same CRS, out of view: zoom only
    canvas.setExtent(QgsRectangle(600000, 2300000, 600100, 2300100))
    done = show_layer(lyr, canvas, prj)
    assert done == "zoomed to the layer"
    # in view: nothing to do
    assert show_layer(lyr, canvas, prj) == ""
    assert not is_local_crs(prj.crs())
    prj.clear()


def test_calculator_and_integrals_algorithms(qgis, tmp_path):
    import processing

    x, y, ik = regular_grid_mesh(0.0, 0.0, 0.25, 0.25, 40, 8, 0.0, "alternate")  # 10 m x 2 m
    H = np.stack([0.2 + 0.0 * x, 0.3 + 0.0 * x])
    U = np.stack([0.5 + 0.0 * x, 0.4 + 0.0 * x])
    V = np.zeros_like(U)
    p = str(tmp_path / "r.slf")
    write_selafin(p, x, y, ik, {"VELOCITY U": U, "VELOCITY V": V, "WATER DEPTH": H}, times=[0.0, 60.0],
                  units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M"})
    out = str(tmp_path / "c.slf")
    r = processing.run("telemactoolkit:calculator", {"INPUT": p, "EXPRESSIONS": "Q [M2/S] = H * U\nHMAX [M] = tmax(H)",
                                                     "KEEP": True, "OUTPUT": out, "LOAD": False})
    with SelafinFile(r["OUTPUT"]) as s:
        assert np.allclose(s.read("Q", 1), 0.12)
        assert np.allclose(s.read("HMAX", 0), 0.3)
    # qgis_process loads the plugin with no QGIS window (iface None): LOAD must not try to open the result on a map
    import qgis.utils

    class Headless:
        iface = None
        opened = []

        def open_path(self, path):
            self.opened.append(path)

    old = qgis.utils.plugins.get("telemac_toolkit")
    qgis.utils.plugins["telemac_toolkit"] = Headless()
    try:
        processing.run("telemactoolkit:calculator", {"INPUT": p, "EXPRESSIONS": "HMAX [M] = tmax(H)",
                                                     "KEEP": False, "OUTPUT": str(tmp_path / "c2.slf"), "LOAD": True})
        assert Headless.opened == []
    finally:
        if old is None:
            qgis.utils.plugins.pop("telemac_toolkit", None)
        else:
            qgis.utils.plugins["telemac_toolkit"] = old
    # discharge through a section drawn from the left bank (y = 2) to the right bank (y = 0)
    lines = QgsVectorLayer("LineString?crs=EPSG:32614&field=name:string", "sec", "memory")
    f = QgsFeature(lines.fields())
    f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(4.1, 2.0), QgsPointXY(4.1, 0.0)]))
    f.setAttribute("name", "S1")
    lines.dataProvider().addFeatures([f])
    r = processing.run("telemactoolkit:integrals", {"INPUT": p, "MODE": 0, "LINES": lines, "OUTPUT": "memory:"})
    feats = sorted(r["OUTPUT"].getFeatures(), key=lambda q: q["frame"])
    assert [q["discharge"] for q in feats] == pytest.approx([0.2, 0.24], rel=1e-6)  # single-precision file
    assert feats[0]["length"] == pytest.approx(2.0)
    # volume of water over the whole mesh
    r = processing.run("telemactoolkit:integrals", {"INPUT": p, "MODE": 3, "EXPRESSION": "H", "OUTPUT": "memory:"})
    vols = sorted((q["frame"], q["integral"]) for q in r["OUTPUT"].getFeatures())
    assert vols[0][1] == pytest.approx(4.0) and vols[1][1] == pytest.approx(6.0)


def test_edit_node_values(qgis, tmp_path):
    import processing

    x, y, ik = regular_grid_mesh(0.0, 0.0, 1.0, 1.0, 20, 20, 0.0, "alternate")
    p = str(tmp_path / "geo.slf")
    write_selafin(p, x, y, ik, {"BOTTOM": 10.0 + 0.0 * x}, units={"BOTTOM": "M"}, origin=(0, 0))
    with open(str(tmp_path / "geo.cli"), "w") as f:
        f.write("2 2 2 0 0 0 0 2 0 0 0 1 1\n")
    with SelafinFile(p) as s:
        ipobo = s.header.ipobo.copy()
    area = QgsVectorLayer("Polygon?crs=EPSG:32614", "a", "memory")
    f = QgsFeature()
    f.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(5, 5), QgsPointXY(15, 5), QgsPointXY(15, 15), QgsPointXY(5, 15), QgsPointXY(5, 5)]]))
    area.dataProvider().addFeatures([f])
    out = str(tmp_path / "geo_edit.slf")
    processing.run("telemactoolkit:interpolatebottom", {
        "INPUT": p, "MESH_CRS": "EPSG:32614", "VARIABLE": "BOTTOM", "CONSTANT": 8.0, "OPERATION": 3, "AREA": area,
        "BLEND": 2.0, "OUTPUT": out, "LOAD": False})
    with SelafinFile(out) as s:
        z = s.read("BOTTOM", 0)
        assert np.array_equal(s.header.ipobo, ipobo)
    centre = (x == 10) & (y == 10)
    edge = (x == 6) & (y == 10)  # 1 m inside the edge: half way
    outside = (x == 2) & (y == 2)
    assert z[centre][0] == pytest.approx(8.0) and z[edge][0] == pytest.approx(9.0) and z[outside][0] == 10.0
    assert os.path.exists(str(tmp_path / "geo_edit.cli"))
    # a new friction variable on the whole mesh
    out2 = str(tmp_path / "geo_fr.slf")
    processing.run("telemactoolkit:interpolatebottom", {"INPUT": out, "MESH_CRS": "EPSG:32614", "VARIABLE": "BOTTOM FRICTION",
                                                        "CONSTANT": 0.03, "OUTPUT": out2, "LOAD": False})
    with SelafinFile(out2) as s:
        assert [n.strip() for n in s.var_names] == ["BOTTOM", "BOTTOM FRICTION"]
        assert np.allclose(s.read("BOTTOM FRICTION", 0), 0.03)


def test_bend_template_and_bed_along_the_axis(qgis, tmp_path):
    import processing

    from telemac_toolkit.qgis_bridge import meshing as M

    g = str(tmp_path / "bend.gpkg")
    r = processing.run("telemactoolkit:channeltemplate", {
        "NEW_PROJECT": g, "PLANFORM": 1, "WIDTH": 1.0, "INLET": 2.0, "OUTLET": 2.0, "RADIUS": 2.0, "BEND_ANGLE": 90.0,
        "TURN": 1, "CELL": 0.1, "STRUCTURED": True, "BED_Z0": 1.0, "SLOPE": 0.01, "SECTION": 0})
    assert "s" in r["BED_FORMULA"]
    st = M.read_settings(g)
    assert st["bottom_axis"] == "flume"
    slf = str(tmp_path / "bend.slf")
    processing.run("telemactoolkit:generatemesh", {"PROJECT": g, "OUTPUT": slf, "LOAD": False})
    with SelafinFile(slf) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        from telemac_toolkit.telemac_io import mesh_statistics

        stt = mesh_statistics(h.x, h.y, h.ikle)
    L = 4.0 + 2.0 * math.pi / 2
    assert z.max() == pytest.approx(1.0, abs=2e-3) and z.min() == pytest.approx(1.0 - 0.01 * L, abs=2e-3)
    assert stt["angle_min"] > 35  # cross-sections normal to the axis: no skewed cells in the bend
    # turning right: the outlet ends below the start line
    assert h.y.min() < -2.0
    layers = M.open_mesh_project(g)
    roles = sorted(f["role"] for f in layers["channels"].getFeatures())
    assert roles == ["axis", "left", "right"]


def test_legend_categories(qgis, tmp_path):
    from telemac_toolkit.qgis_bridge import meshing as M

    g = str(tmp_path / "p.gpkg")
    layers = M.create_mesh_project(g, M.local_crs(), add_to_project=False)
    assert len(layers["channels"].renderer().categories()) == 3
    assert len(layers["constraint_lines"].renderer().categories()) == 2
    assert M.auto_default_size({"domain": layers["domain"]}) == 0.0


def test_calculator_xy_outputs_get_arrows(qgis, tmp_path):
    """'NAME X' / 'NAME Y' from the calculator open as one vector with arrows (QGIS pairs ALONG X/Y only)."""
    from telemac_toolkit.qgis_bridge import styles
    from telemac_toolkit.qgis_bridge.results import load_result
    from telemac_toolkit.telemac_io.calculator import run_calculator

    p = _uniform_result(str(tmp_path / "u.slf"), n=10)
    out = str(tmp_path / "slope.slf")
    run_calculator(p, out, "SLOPE X = 0.0 * H + 0.5\nSLOPE Y = 0.0 * H - 0.5")
    rs = load_result(out, crs=QgsCoordinateReferenceSystem("EPSG:32614"), add_to_project=False)
    labels = {lab: (i, j) for i, j, lab in rs.pairs}
    assert set(labels) == {"Velocity", "Slope"}
    g = rs.group_of_pair[labels["Slope"]]
    assert g != rs.group_of_pair[labels["Velocity"]]
    styles.set_active_variable(rs, labels["Slope"])
    assert styles.arrow_group(rs) == g
