"""Integration tests in a real (headless) QGIS: run with QGIS's Python.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis

Covers loading, styles, the clock, live cursor, probes, profiles, 3D files,
Processing algorithms and plugin unload, on QGIS 3.44 (Qt5) and 4.x (Qt6).
"""

import os
import shutil
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, "tests", "data")
OUT = os.environ.get("TT_SCREENSHOTS", "")

from qgis.core import (  # noqa: E402
    Qgis,
    QgsApplication,
    QgsCoordinateReferenceSystem,
    QgsDateTimeRange,
    QgsFeature,
    QgsGeometry,
    QgsMapRendererParallelJob,
    QgsPointXY,
    QgsProject,
    QgsTemporalNavigationObject,
    QgsVectorLayer,
)
from qgis.gui import QgsMapCanvas, QgsMapMouseEvent, QgsMessageBar  # noqa: E402
from qgis.PyQt.QtCore import QEvent, QPoint, QPointF, Qt  # noqa: E402
from qgis.PyQt.QtGui import QMouseEvent  # noqa: E402
from qgis.PyQt.QtWidgets import QMainWindow, QMenu  # noqa: E402

_app = None


def app():
    global _app
    if _app is None:
        _app = QgsApplication.instance()  # another test module may have started QGIS already
        if _app is None:
            _app = QgsApplication([], True)
            _app.initQgis()
        from processing.core.Processing import Processing

        Processing.initialize()
    return _app


class FakeIface:
    """The parts of QgisInterface the plugin uses."""

    def __init__(self):
        self.win = QMainWindow()
        self.win.resize(1400, 900)
        self.canvas = QgsMapCanvas(self.win)
        self.canvas.resize(900, 700)
        self.nav = QgsTemporalNavigationObject(self.canvas)
        self.canvas.setTemporalController(self.nav)
        self.win.setCentralWidget(self.canvas)
        self.bar = QgsMessageBar(self.win)
        self.menu = QMenu()
        self.handlers = []

    def mainWindow(self):
        return self.win

    def mapCanvas(self):
        return self.canvas

    def messageBar(self):
        return self.bar

    def addToolBar(self, name):
        return self.win.addToolBar(name)

    def addPluginToMenu(self, name, action):
        self.menu.addAction(action)

    def removePluginMenu(self, name, action):
        self.menu.removeAction(action)

    def addDockWidget(self, area, dock):
        self.win.addDockWidget(area, dock)

    def removeDockWidget(self, dock):
        self.win.removeDockWidget(dock)

    def registerCustomDropHandler(self, h):
        self.handlers.append(h)

    def unregisterCustomDropHandler(self, h):
        self.handlers.remove(h)


@pytest.fixture(scope="module")
def env(tmp_path_factory):
    app()
    tmp = tmp_path_factory.mktemp("tt")
    for f in os.listdir(DATA):
        if f.endswith(".slf"):
            shutil.copy(os.path.join(DATA, f), tmp / f)
    QgsProject.instance().clear()
    QgsProject.instance().setCrs(QgsCoordinateReferenceSystem("EPSG:32614"))
    iface = FakeIface()
    from telemac_toolkit import classFactory

    plugin = classFactory(iface)
    plugin.interactive = False  # no CRS question in tests: the best guess is used
    plugin.initGui()
    yield plugin, iface, tmp
    plugin.unload()
    QgsProject.instance().clear()


def mouse(canvas, map_xy, kind=QEvent.Type.MouseMove, button=Qt.MouseButton.NoButton):
    px = canvas.getCoordinateTransform().transform(QgsPointXY(*map_xy))
    qe = QMouseEvent(kind, QPointF(px.x(), px.y()), QPointF(px.x(), px.y()), button, button, Qt.KeyboardModifier.NoModifier)
    return QgsMapMouseEvent(canvas, qe)


def render(iface, rs, path):
    s = iface.canvas.mapSettings()
    s.setOutputSize(iface.canvas.size())
    s.setLayers([rs.layer])
    s.setExtent(rs.layer.extent().buffered(rs.layer.extent().width() * 0.03))
    s.setIsTemporal(True)
    s.setTemporalRange(iface.canvas.temporalRange())
    job = QgsMapRendererParallelJob(s)
    job.start()
    job.waitForFinished()
    img = job.renderedImage()
    if path:
        img.save(path)
    return img


def test_open_2d_result(env):
    plugin, iface, tmp = env
    rs = plugin.open_path(str(tmp / "r2d_tidal_flats.slf"))
    assert rs is not None and rs.layer.isValid()
    assert rs.nframes == 17
    root = rs.layer.datasetGroupTreeRootItem()
    names = sorted(root.child(i).name() for i in range(root.childCount()))
    assert names == ["Bottom (m)", "Free surface (m)", "Velocity (m/s)", "Water depth (m)"]
    s = rs.layer.rendererSettings()
    assert s.activeScalarDatasetGroup() == rs.group_of_var[rs.var_index("water_depth")]
    assert s.activeVectorDatasetGroup() == rs.group_of_pair[(0, 1)]
    d = plugin.dock
    assert d.rs is rs and d.slider.maximum() == 16
    assert d.map_var.count() == 4  # depth, surface, bottom + velocity (U, V are one vector group)
    assert d.pr_var.itemData(0) == d.WATER_PROFILE


def test_clock_drives_layer_dataset(env):
    plugin, iface, tmp = env
    rs = plugin.dock.rs
    for target in (5, 16, 0, 9):
        plugin.dock.set_frame(target)
        assert plugin.timesync.frame == target
        assert plugin.dock.slider.value() == target
        rng = iface.canvas.temporalRange()
        idx = rs.layer.activeScalarDatasetAtTime(rng)
        assert idx.dataset() == target, (target, idx.dataset(), rng.begin().toString())
    # moving the QGIS controller moves the plugin
    iface.nav.setCurrentFrameNumber(12)
    assert plugin.timesync.frame == 12
    assert "12 / 17" not in plugin.dock.time_label.text() and "13 / 17" in plugin.dock.time_label.text()


def test_live_cursor_values(env):
    plugin, iface, tmp = env
    rs = plugin.dock.rs
    iface.canvas.setDestinationCrs(rs.layer.crs())
    iface.canvas.setExtent(rs.layer.extent())
    plugin.toggle_live_cursor()
    assert iface.canvas.mapTool() is plugin.cursor_tool and plugin.a_cursor.isChecked()
    node = 200
    x, y = float(rs.slf.x[node]), float(rs.slf.y[node])
    plugin.cursor_tool.canvasMoveEvent(mouse(iface.canvas, (x, y)))
    t = plugin.dock.table
    vals = {t.item(r, 0).text(): float(t.item(r, 1).text()) for r in range(t.rowCount())}
    fr = plugin.timesync.frame
    # mouse position snaps to a pixel, so compare with interpolation at the snapped point
    cx, cy = plugin.dock.cursor_pt
    r = rs.locate(cx, cy)
    expect = rs.slf.read("WATER DEPTH", fr)[r[1]] @ r[2]
    assert vals["Water depth"] == pytest.approx(expect, rel=1e-3, abs=1e-3)
    assert "Velocity magnitude" in vals and "Velocity direction (bearing)" in vals
    plugin.toggle_live_cursor()
    assert iface.canvas.mapTool() is None or iface.canvas.mapTool() is not plugin.cursor_tool


def test_probe_series(env):
    plugin, iface, tmp = env
    rs = plugin.dock.rs
    plugin.toggle_probe()
    node = 321
    x, y = float(rs.slf.x[node]), float(rs.slf.y[node])
    ev = mouse(iface.canvas, (x, y), QEvent.Type.MouseButtonRelease, Qt.MouseButton.LeftButton)
    plugin.probe_tool.canvasReleaseEvent(ev)
    d = plugin.dock
    assert len(d.probes) == 1 and len(plugin.probe_tool.markers) == 1
    k = d.ts_var.findData(rs.var_index("free_surface"))
    d.ts_var.setCurrentIndex(k)
    ser = d.ts_plot.series[0]
    px, py = d.probes[0]
    r = rs.locate(px, py)
    np.testing.assert_allclose(ser.y, rs.slf.series_weighted("FREE SURFACE", r[1], r[2]))
    assert ser.x.size == 17
    csv = d.ts_plot.to_csv()
    assert csv.count("\n") == 18
    # peak button jumps to the probe's maximum
    d.map_var.setCurrentIndex(d.map_var.findData(rs.var_index("free_surface")))
    d.jump_to_peak()
    assert plugin.timesync.frame == int(np.argmax(ser.y))


def test_profile(env):
    plugin, iface, tmp = env
    rs = plugin.dock.rs
    e = rs.layer.extent()
    line = [(e.xMinimum() + 0.05 * e.width(), e.center().y()), (e.xMaximum() - 0.05 * e.width(), e.center().y())]
    plugin.dock.set_profile_line(line)
    p = plugin.dock.pr_plot
    labels = [s.label for s in p.series]
    assert labels == ["Bed", "Free surface"]
    assert p.band is not None
    bed = p.series[0].y
    assert np.isfinite(bed).mean() > 0.95
    before = p.series[1].y.copy()
    plugin.dock.set_frame((plugin.timesync.frame + 4) % rs.nframes)
    after = plugin.dock.pr_plot.series[1].y
    assert not np.allclose(np.nan_to_num(before), np.nan_to_num(after))  # profile follows the clock


def test_screenshots(env):
    plugin, iface, tmp = env
    if not OUT:
        pytest.skip("set TT_SCREENSHOTS to a folder to save images")
    os.makedirs(OUT, exist_ok=True)
    rs = plugin.dock.rs
    plugin.dock.set_frame(4)
    img = render(iface, rs, os.path.join(OUT, "map_%s.png" % Qgis.version().split("-")[0]))
    assert not img.isNull()
    plugin.dock.resize(430, 820)
    plugin.dock.show()
    for k in (0, 1, 2):
        plugin.dock.tabs.setCurrentIndex(k)
        plugin.dock.widget().grab().save(os.path.join(OUT, "dock_tab%d_%s.png" % (k, Qgis.version().split("-")[0])))


def test_open_3d_result(env):
    plugin, iface, tmp = env
    rs = plugin.open_path(str(tmp / "r3d_bump.slf"), plane="average")
    assert rs is not None and rs.layer.isValid()
    assert rs.derived_from_3d == "average"
    assert os.path.exists(str(tmp / "r3d_bump_2d_average.slf"))
    assert rs.nframes == 3
    assert plugin.dock.rs is rs
    assert plugin.dock.result_combo.count() == 2


def test_history_file_is_rejected_cleanly(env):
    plugin, iface, tmp = env
    n = len(plugin.results)
    assert plugin.open_path(str(tmp / "r1d_tomsail.slf")) is None
    assert len(plugin.results) == n


def test_layer_removal_detaches(env):
    plugin, iface, tmp = env
    rs = plugin.dock.rs
    QgsProject.instance().removeMapLayer(rs.layer_id)
    assert rs.layer_id not in plugin.results
    assert plugin.dock.result_combo.count() == 1


def test_processing_algorithms(env):
    plugin, iface, tmp = env
    import processing

    r2d = str(tmp / "r2d_tidal_flats.slf")
    info = processing.run("telemactoolkit:selafininfo", {"INPUT": r2d, "OUTPUT": str(tmp / "info.html")})
    assert "Water depth" in open(info["OUTPUT"], encoding="utf-8").read()

    out = processing.run("telemactoolkit:extract2d", {"INPUT": str(tmp / "r3d_tidal_flats.slf"), "MODE": 3, "PLANE": 20,
                                                      "OUTPUT": str(tmp / "top.slf")})
    from telemac_toolkit.telemac_io import SelafinFile

    with SelafinFile(out["OUTPUT"]) as s, SelafinFile(str(tmp / "r3d_tidal_flats.slf")) as s3:
        np.testing.assert_allclose(s.read("MUD", -1), s3.plane("MUD", -1, 20), atol=1e-6)

    sub = processing.run("telemactoolkit:subset", {"INPUT": r2d, "FRAMES": "0:17:4", "VARIABLES": "water depth, bottom",
                                                   "OUTPUT": str(tmp / "sub.slf")})
    with SelafinFile(sub["OUTPUT"]) as s:
        assert s.nframes == 5 and s.var_names == ["WATER DEPTH", "BOTTOM"]

    with SelafinFile(r2d) as s:
        pts = QgsVectorLayer("Point?crs=EPSG:32614", "pts", "memory")
        feats = []
        for node in (10, 50):
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(float(s.x[node]), float(s.y[node]))))
            feats.append(f)
        f = QgsFeature()
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(-1e7, -1e7)))  # outside
        feats.append(f)
        pts.dataProvider().addFeatures(feats)
        ts = processing.run("telemactoolkit:probetimeseries", {"INPUT": r2d, "POINTS": pts, "VARIABLES": "WATER DEPTH",
                                                               "OUTPUT": "memory:"})["OUTPUT"]
        assert ts.featureCount() == 2 * 17
        got = [f["water_depth"] for f in ts.getFeatures() if f["point_id"] == feats[0].id() or True][:17]
        np.testing.assert_allclose(got, s.series("WATER DEPTH", 10), atol=1e-5)

    b = processing.run("telemactoolkit:meshboundary", {"INPUT": r2d, "CRS": "EPSG:32614", "OUTPUT": "memory:"})["OUTPUT"]
    assert b.featureCount() == 1
    f = next(b.getFeatures())
    assert f["kind"] == "outer" and f["nodes"] == 264


def test_valley_flood_result(env):
    """A real TELEMAC-2D v9.1 run (tests/data/demo_valley) with wetting and drying and a Froude variable."""
    plugin, iface, tmp = env
    src = os.path.join(DATA, "demo_valley", "r2d_valley.slf")
    rs = plugin.open_path(src)  # Froude styling used to crash QGIS: regression test
    assert rs is not None and rs.nframes == 31
    assert rs.has_calendar_date and rs.header.start_datetime.isoformat() == "2026-09-14T06:00:00"
    d = plugin.dock
    labels = [d.ts_var.itemText(k) for k in range(d.ts_var.count())]
    assert len(labels) == len(set(labels)), labels  # no duplicate velocity magnitude
    # origin offset applied by both QGIS and the plugin
    e = rs.layer.extent()
    assert 511990 < e.xMinimum() < 512010 and 2154690 < e.yMinimum() < 2154710
    # time of the flood peak (area-weighted mean depth)
    assert rs.peak_frame(rs.var_index("water_depth")) == 13
    # cursor rows have unique labels
    rows = rs.values_at(512000 + 1050, 2155000 + 120, 13)
    names = [r[0] for r in rows]
    assert len(names) == len(set(names))
    # dry threshold is applied to the depth style (values below are clipped)
    g = rs.group_of_var[rs.var_index("water_depth")]
    sh = rs.layer.rendererSettings().scalarSettings(g).colorRampShader()
    assert sh.clip() and sh.minimumValue() == pytest.approx(0.01)
    plugin.dock.dry.setValue(0.05)
    plugin.dock._on_dry()
    sh = rs.layer.rendererSettings().scalarSettings(g).colorRampShader()
    assert sh.minimumValue() == pytest.approx(0.05)


def test_running_simulation_refresh(env):
    """A result still being written: open with 10 frames, append the rest, press Refresh."""
    plugin, iface, tmp = env
    from telemac_toolkit.telemac_io import SelafinFile

    full = os.path.join(DATA, "demo_valley", "r2d_valley.slf")
    with SelafinFile(full) as s:
        head, fs = s.header.header_size, s.header.frame_size
    data = open(full, "rb").read()
    part = str(tmp / "running.slf")
    with open(part, "wb") as f:
        f.write(data[: head + 10 * fs + 123])  # 10 frames + part of the 11th
    rs = plugin.open_path(part)
    assert rs.nframes == 10 and rs.warnings()
    assert rs.layer.dataProvider().datasetCount(0) in (10, 11)
    with open(part, "ab") as f:
        f.write(data[head + 10 * fs + 123:])
    plugin.dock.refresh()
    assert rs.nframes == 31
    assert plugin.dock.slider.maximum() == 30
    assert plugin.timesync.frame == 30  # jumps to the newest frame
    assert rs.layer.dataProvider().datasetCount(0) == 31


def test_sidecar_prj_sets_crs(env):
    plugin, iface, tmp = env
    src = str(tmp / "with_prj.slf")
    shutil.copy(os.path.join(DATA, "demo_valley", "r2d_valley.slf"), src)
    with open(str(tmp / "with_prj.prj"), "w") as f:
        f.write(QgsCoordinateReferenceSystem("EPSG:6372").toWkt())
    rs = plugin.open_path(src)
    assert rs.layer.crs().authid() == "EPSG:6372"


def test_probe_algorithm_transforms_points(env):
    plugin, iface, tmp = env
    import processing
    from qgis.core import QgsCoordinateTransform

    r2d = os.path.join(DATA, "demo_valley", "r2d_valley.slf")
    utm, geo = QgsCoordinateReferenceSystem("EPSG:32614"), QgsCoordinateReferenceSystem("EPSG:4326")
    tr = QgsCoordinateTransform(utm, geo, QgsProject.instance())
    p = tr.transform(QgsPointXY(512000 + 1700, 2155000 - 40))
    pts = QgsVectorLayer("Point?crs=EPSG:4326", "pts", "memory")
    f = QgsFeature()
    f.setGeometry(QgsGeometry.fromPointXY(p))
    pts.dataProvider().addFeatures([f])
    out = processing.run("telemactoolkit:probetimeseries", {"INPUT": r2d, "POINTS": pts, "VARIABLES": "WATER DEPTH",
                                                            "MESH_CRS": "EPSG:32614", "OUTPUT": "memory:"})["OUTPUT"]
    vals = [g["water_depth"] for g in out.getFeatures()]
    assert len(vals) == 31 and max(vals) > 3.5


def test_arrow_list_follows_the_map_vector(env):
    """The panel's arrow list: Off + every vector; choosing a vector on the map moves the arrows to it."""
    plugin, iface, tmp = env
    from telemac_toolkit.mesher.regular import regular_grid_mesh
    from telemac_toolkit.qgis_bridge import styles
    from telemac_toolkit.telemac_io import write_selafin

    x, y, ik = regular_grid_mesh(500000.0, 2000000.0, 10.0, 10.0, 10, 10, 0.0, "alternate")
    one = np.ones_like(x)
    p = str(tmp / "three_vectors.slf")
    write_selafin(p, x, y, ik, {"VELOCITY U": one, "VELOCITY V": 0 * one, "WATER DEPTH": 0.3 * one,
                                "FLOWRATE ALONG X": 0 * one, "FLOWRATE ALONG Y": 0.3 * one,
                                "WIND ALONG X": -3 * one, "WIND ALONG Y": -3 * one},
                  units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M", "FLOWRATE ALONG X": "M2/S",
                         "FLOWRATE ALONG Y": "M2/S", "WIND ALONG X": "M/S", "WIND ALONG Y": "M/S"})
    rs = plugin.open_path(p)
    d = plugin.dock
    assert d.rs is rs

    def show(v):  # QComboBox.findData does not match Python tuples
        d.map_var.setCurrentIndex(next(k for k in range(d.map_var.count()) if d.map_var.itemData(k) == v))

    labels = [d.arrows.itemText(k) for k in range(d.arrows.count())]
    assert labels == ["Off", "Velocity (m/s)", "Unit discharge (m2/s)", "Wind (m/s)"]
    assert d.arrows.currentText() == "Velocity (m/s)"
    assert "m/s" in d.ar_ref.text()
    show((3, 4))  # unit discharge on the map
    assert d.arrows.currentText() == "Unit discharge (m2/s)"
    assert styles.arrow_group(rs) == rs.group_of_pair[(3, 4)]
    assert "m2/s" in d.ar_ref.text()
    show(2)  # depth: no direction, arrows unchanged
    assert d.arrows.currentText() == "Unit discharge (m2/s)"
    d.arrows.setCurrentIndex(d.arrows.findText("Wind (m/s)"))
    assert styles.arrow_group(rs) == rs.group_of_pair[(5, 6)]
    d.arrows.setCurrentIndex(0)
    assert styles.arrow_group(rs) is None
    show((0, 1))  # arrows off stay off
    assert styles.arrow_group(rs) is None and d.arrows.currentIndex() == 0
