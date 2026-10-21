"""Render screenshots of the plugin on the demo valley flood (headless QGIS).

    QT_QPA_PLATFORM=offscreen python scripts/render_demo.py OUTDIR
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from qgis.core import Qgis, QgsApplication, QgsCoordinateReferenceSystem, QgsProject, QgsRectangle, QgsTemporalNavigationObject  # noqa: E402
from qgis.gui import QgsMapCanvas, QgsMessageBar  # noqa: E402
from qgis.PyQt.QtCore import QSize, Qt  # noqa: E402
from qgis.PyQt.QtWidgets import QMainWindow, QMenu  # noqa: E402


class FakeIface:
    def __init__(self):
        self.win = QMainWindow()
        self.win.resize(1500, 860)
        self.canvas = QgsMapCanvas(self.win)
        self.canvas.setCanvasColor(Qt.GlobalColor.white)
        self.nav = QgsTemporalNavigationObject(self.canvas)
        self.canvas.setTemporalController(self.nav)
        self.win.setCentralWidget(self.canvas)
        self.bar = QgsMessageBar(self.win)
        self.menu = QMenu()

    def mainWindow(self): return self.win
    def mapCanvas(self): return self.canvas
    def messageBar(self): return self.bar
    def addToolBar(self, name): return self.win.addToolBar(name)
    def addPluginToMenu(self, n, a): self.menu.addAction(a)
    def removePluginMenu(self, n, a): self.menu.removeAction(a)
    def addDockWidget(self, area, d): self.win.addDockWidget(area, d)
    def removeDockWidget(self, d): self.win.removeDockWidget(d)
    def registerCustomDropHandler(self, h): pass
    def unregisterCustomDropHandler(self, h): pass


def settle(app, canvas):
    for _ in range(5):
        app.processEvents()
        canvas.waitWhileRendering()
        app.processEvents()


def main(out):
    os.makedirs(out, exist_ok=True)
    app = QgsApplication([], True)
    app.initQgis()
    QgsProject.instance().setCrs(QgsCoordinateReferenceSystem("EPSG:32614"))
    iface = FakeIface()
    from telemac_toolkit import classFactory

    plugin = classFactory(iface)
    plugin.initGui()
    iface.win.show()
    canvas = iface.canvas
    canvas.setDestinationCrs(QgsProject.instance().crs())
    path = os.path.join(ROOT, "tests", "data", "demo_valley", "r2d_valley.slf")
    rs = plugin.open_path(path)
    dock = plugin.dock
    dock.setMinimumWidth(440)
    canvas.setLayers([rs.layer])
    e = rs.layer.extent()
    canvas.setExtent(QgsRectangle(e.xMinimum() - 20, e.yMinimum() - 60, e.xMaximum() + 20, e.yMaximum() + 60))
    # probes: channel upstream, floodplain, channel downstream
    ox, oy = 512000.0, 2155000.0
    for px, py in ((ox + 300, oy + 70 * 0.5), (ox + 1050, oy + 120), (ox + 1700, oy - 40)):
        dock.add_probe(px, py)
    for (px, py) in dock.probes:
        from qgis.core import QgsPointXY
        from qgis.gui import QgsVertexMarker
        from qgis.PyQt.QtGui import QColor
        from telemac_toolkit.gui.plot_widget import PALETTE

        m = QgsVertexMarker(canvas)
        m.setCenter(QgsPointXY(px, py))
        m.setIconType(QgsVertexMarker.IconType.ICON_CIRCLE)
        m.setIconSize(12)
        m.setPenWidth(3)
        m.setColor(QColor(PALETTE[len(plugin.probe_tool.markers) % len(PALETTE)]))
        plugin.probe_tool.markers.append(m)
    k = dock.ts_var.findData(rs.var_index("water_depth"))
    dock.ts_var.setCurrentIndex(k)
    # cross-valley profile at x = 1050 m
    dock.set_profile_line([(ox + 1050, oy - 290), (ox + 1050, oy + 290)])
    peak = int(rs.peak_frame(rs.var_index("water_depth")))
    tag = Qgis.version().split("-")[0]
    for fr, name in ((0, "start"), (peak, "peak")):
        dock.set_frame(fr)
        for tab in (1, 2):
            dock.tabs.setCurrentIndex(tab)
            settle(app, canvas)
            canvas.refresh()
            settle(app, canvas)
            iface.win.grab().save(os.path.join(out, "window_%s_tab%d_%s.png" % (name, tab, tag)))
    dock.tabs.setCurrentIndex(0)
    dock.update_cursor(ox + 1050, oy + 120)
    settle(app, canvas)
    iface.win.grab().save(os.path.join(out, "window_peak_cursor_%s.png" % tag))
    print("peak frame", peak, rs.time_label(peak))
    plugin.unload()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/shots")
