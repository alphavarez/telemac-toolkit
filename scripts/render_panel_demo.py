"""Screenshots of the Mesh workflow panel and of the flume demo (headless QGIS).

    QT_QPA_PLATFORM=offscreen python scripts/render_panel_demo.py VALLEY_CASEDIR FLUME_CASEDIR OUTDIR
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _p in (os.path.join(sys.prefix, "share", "qgis", "python", "plugins"),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.append(_p)

from qgis.core import (  # noqa: E402
    QgsApplication,
    QgsDateTimeRange,
    QgsMapRendererParallelJob,
    QgsMapSettings,
    QgsMeshLayer,
    QgsProject,
    QgsRasterLayer,
    QgsRectangle,
)
from qgis.gui import QgsMapCanvas  # noqa: E402
from qgis.PyQt.QtCore import QSize, Qt  # noqa: E402
from qgis.PyQt.QtGui import QColor  # noqa: E402
from qgis.PyQt.QtWidgets import QMainWindow  # noqa: E402


class Iface:
    def __init__(self):
        self.win = QMainWindow()
        self.win.resize(1500, 900)
        self.canvas = QgsMapCanvas(self.win)
        self.canvas.setCanvasColor(Qt.GlobalColor.white)
        self.win.setCentralWidget(self.canvas)

    def mapCanvas(self):
        return self.canvas

    def mainWindow(self):
        return self.win


def settle(app, canvas):
    for _ in range(6):
        app.processEvents()
        canvas.waitWhileRendering()
        app.processEvents()


def render(layers, extent, path, crs, size=(1500, 300), temporal=None):
    s = QgsMapSettings()
    s.setLayers(layers)
    s.setOutputSize(QSize(*size))
    s.setExtent(extent)
    s.setBackgroundColor(QColor("white"))
    s.setDestinationCrs(crs)
    s.setOutputDpi(110)
    if temporal is not None:
        s.setIsTemporal(True)
        s.setTemporalRange(temporal)
    job = QgsMapRendererParallelJob(s)
    job.start()
    job.waitForFinished()
    job.renderedImage().save(path)
    print("wrote", path)


def main(valley, flume, out):
    os.makedirs(out, exist_ok=True)
    app = QgsApplication([], True)
    app.initQgis()
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    from telemac_toolkit.gui.mesh_panel import MeshPanel
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.qgis_bridge.results import load_result
    from telemac_toolkit.qgis_bridge import styles

    # --- the panel on the valley project ------------------------------------------
    gpkg = os.path.join(valley, "valley_mesh.gpkg")
    lay = M.open_mesh_project(gpkg)
    M.add_project_to_map(lay, gpkg)
    crs = lay["domain"].crs()
    QgsProject.instance().setCrs(crs)
    dem = QgsRasterLayer(os.path.join(valley, "dem_valley.tif"), "DEM")
    QgsProject.instance().addMapLayer(dem)
    mesh = QgsMeshLayer(os.path.join(valley, "geo_valley_mesh.slf"), "geo_valley_mesh", "mdal")
    mesh.setCrs(crs)
    M.style_geometry_mesh(mesh)
    QgsProject.instance().addMapLayer(mesh)
    iface = Iface()
    panel = MeshPanel(iface, iface.win)
    iface.win.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, panel)
    panel.setMinimumWidth(430)
    panel.set_project(gpkg)
    for i in range(panel.rasters.count()):
        if panel.rasters.item(i).text() == "DEM":
            panel.rasters.item(i).setCheckState(Qt.CheckState.Checked)
    panel.size.setValue(30.0)
    panel.growth.setValue(1.15)
    panel.output.setFilePath(os.path.join(out, "panel_valley.slf"))
    panel.generate(wait=True)  # a real run: the log in the picture is the run's own
    canvas = iface.canvas
    canvas.setDestinationCrs(crs)
    canvas.setLayers([lay["boundaries"], lay["channels"], mesh, dem])
    e = mesh.extent()
    canvas.setExtent(QgsRectangle(e.xMinimum() + 100, e.yMinimum() + 60, e.xMaximum() - 900, e.yMaximum() - 60))
    iface.win.show()
    for page in (1, 3):
        panel.toolbox.setCurrentIndex(page)
        settle(app, canvas)
        canvas.refresh()
        settle(app, canvas)
        iface.win.grab().save(os.path.join(out, "panel_step%d.png" % (page + 1)))
        print("wrote panel_step%d.png" % (page + 1))

    # --- the flume, local coordinates -----------------------------------------------
    fmesh = QgsMeshLayer(os.path.join(flume, "geo_flume.slf"), "flume", "mdal")
    lcrs = M.local_crs()
    fmesh.setCrs(lcrs)
    M.style_geometry_mesh(fmesh)
    render([fmesh], QgsRectangle(-0.2, -0.15, 6.2, 1.15), os.path.join(out, "flume_mesh.png"), lcrs, size=(1500, 300))
    res = os.path.join(flume, "r2d_flume.slf")
    if os.path.exists(res):
        rs = load_result(res, crs=lcrs, add_to_project=True)
        k = rs.var_index("water_depth")
        styles.set_active_variable(rs, k)
        t = rs.reference_datetime().addSecs(int(rs.times[-1]))
        render([rs.layer], QgsRectangle(-0.3, -0.25, 20.3, 1.25), os.path.join(out, "flume_depth.png"), lcrs, size=(1500, 160),
               temporal=QgsDateTimeRange(t, t.addSecs(1)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "/tmp/shots")
