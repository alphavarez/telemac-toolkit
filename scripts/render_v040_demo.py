"""Pictures for 0.4.0 from the bend run (scripts/make_bend_demo.py + TELEMAC + scripts/check_bend_demo.py).

    QT_QPA_PLATFORM=offscreen python scripts/render_v040_demo.py BEND_CASEDIR OUTDIR

arrows_before.png / arrows_after.png  the same frame with the 0.3.0 arrows and the new defaults
results_panel_arrows.png              the Results panel with the arrow settings open
bend_vorticity.png                    vorticity from the Calculator
bend_mesh.png                         the structured bend from the channel template
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests", "qgis"))
for _p in (os.path.join(sys.prefix, "share", "qgis", "python", "plugins"),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.append(_p)

from qgis.core import (  # noqa: E402
    QgsApplication,
    QgsDateTimeRange,
    QgsMapRendererParallelJob,
    QgsMapSettings,
    QgsMeshLayer,
    QgsMeshRendererVectorArrowSettings,
    QgsProject,
    QgsRectangle,
)
from qgis.PyQt.QtCore import QSize, Qt  # noqa: E402
from qgis.PyQt.QtGui import QColor  # noqa: E402

EXTENT = QgsRectangle(-0.4, -0.9, 6.0, 5.0)


def render(layer, path, crs, t=None, size=(760, 700), extent=EXTENT):
    s = QgsMapSettings()
    s.setLayers([layer])
    s.setOutputSize(QSize(*size))
    s.setExtent(extent)
    s.setBackgroundColor(QColor("white"))
    s.setDestinationCrs(crs)
    s.setOutputDpi(96)
    if t is not None:
        s.setIsTemporal(True)
        s.setTemporalRange(QgsDateTimeRange(t, t.addSecs(1)))
    job = QgsMapRendererParallelJob(s)
    job.start()
    job.waitForFinished()
    job.renderedImage().save(path)
    print("wrote", path)


def main(case, out):
    os.makedirs(out, exist_ok=True)
    app = QgsApplication([], True)
    app.initQgis()
    from processing.core.Processing import Processing

    Processing.initialize()
    from test_plugin_qgis import FakeIface

    from telemac_toolkit import classFactory
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.qgis_bridge import styles

    QgsProject.instance().clear()
    iface = FakeIface()
    iface.win.resize(1500, 950)
    plugin = classFactory(iface)
    plugin.interactive = False
    plugin.initGui()
    rs = plugin.open_path(os.path.join(case, "r2d_bend.slf"))
    crs = rs.layer.crs()
    k = rs.var_index("scalar_velocity")
    pair = next(iter(rs.group_of_pair))
    styles.set_active_variable(rs, pair if k is None else k)
    t = rs.reference_datetime().addSecs(int(rs.times[-1]))
    g = rs.group_of_pair[pair]
    # 0.3.0 arrows: grid of 36 px with arrows up to 28.8 "mm" (the grid was in pixels, the length in mm)
    st = rs.layer.rendererSettings()
    vs = st.vectorSettings(g)
    vs.setOnUserDefinedGrid(True)
    vs.setUserGridCellWidth(36)
    vs.setUserGridCellHeight(36)
    a = vs.arrowSettings()
    a.setShaftLengthMethod(QgsMeshRendererVectorArrowSettings.ArrowScalingMethod.MinMax)
    a.setMinShaftLength(3)
    a.setMaxShaftLength(36 * 0.8)
    vs.setArrowsSettings(a)
    st.setVectorSettings(g, vs)
    rs.layer.setRendererSettings(st)
    render(rs.layer, os.path.join(out, "arrows_before.png"), crs, t)
    # 0.4.0 defaults
    styles.apply_arrow_style(rs, styles.ArrowStyle())
    render(rs.layer, os.path.join(out, "arrows_after.png"), crs, t)
    # panel, arrow settings open
    d = plugin.dock
    d.set_frame(rs.nframes - 1)
    want = pair if k is None else k
    for i in range(d.map_var.count()):
        if d.map_var.itemData(i) == want:
            d.map_var.setCurrentIndex(i)
    d.arrow_box.setCollapsed(False)
    d.ar_spacing.setValue(28)
    d._apply_arrows()
    canvas = iface.canvas
    canvas.setDestinationCrs(crs)
    canvas.setLayers([rs.layer])
    canvas.setExtent(EXTENT)
    iface.win.show()
    for _ in range(8):
        app.processEvents()
        canvas.waitWhileRendering()
    d.setMinimumWidth(460)
    for _ in range(4):
        app.processEvents()
        canvas.refresh()
        canvas.waitWhileRendering()
    iface.win.grab().save(os.path.join(out, "results_panel_arrows.png"))
    print("wrote results_panel_arrows.png")
    # vorticity from the calculator
    calc = os.path.join(case, "r2d_bend_calc.slf")
    if os.path.exists(calc):
        vr = QgsMeshLayer(calc, "calc", "mdal")
        vr.setCrs(crs)
        from telemac_toolkit.qgis_bridge.results import ResultSet

        crs_rs = ResultSet(calc, calc, vr)
        styles.rename_groups(crs_rs)
        gi = crs_rs.group_of_var[crs_rs.var_index("VORTICITY")]
        ss = vr.rendererSettings()
        sc = ss.scalarSettings(gi)
        from qgis.core import QgsColorRampShader, QgsGradientColorRamp, QgsGradientStop

        lim = 0.6
        ramp = QgsGradientColorRamp(QColor("#2166ac"), QColor("#b2182b"))
        ramp.setStops([QgsGradientStop(0.5, QColor("#f7f7f7"))])
        sc.setClassificationMinimumMaximum(-lim, lim)
        sh = QgsColorRampShader(-lim, lim)
        sh.setSourceColorRamp(ramp)
        sh.setColorRampItemList(styles._items_from_ramp(ramp, -lim, lim))
        sh.setClip(False)
        sc.setColorRampShader(sh)
        ss.setScalarSettings(gi, sc)
        ss.setActiveScalarDatasetGroup(gi)
        vr.setRendererSettings(ss)
        render(vr, os.path.join(out, "bend_vorticity.png"), crs)
    # the question asked when nothing tells what the coordinates are
    from telemac_toolkit.gui.crs_dialog import CrsDialog
    from telemac_toolkit.qgis_bridge.crs import local_crs
    from telemac_toolkit.qgis_bridge.results import file_extent

    dlg = CrsDialog(iface.win, os.path.join(case, "r2d_bend.slf"), file_extent(os.path.join(case, "r2d_bend.slf")), local_crs())
    dlg.show()
    for _ in range(5):
        app.processEvents()
    dlg.grab().save(os.path.join(out, "crs_dialog.png"))
    print("wrote crs_dialog.png")
    dlg.close()
    # the mesh of the template
    geo = QgsMeshLayer(os.path.join(case, "geo_bend.slf"), "geo", "mdal")
    geo.setCrs(crs)
    M.style_geometry_mesh(geo)
    render(geo, os.path.join(out, "bend_mesh.png"), crs, size=(760, 700),
           extent=QgsRectangle(2.2, -0.75, 5.8, 2.3))
    plugin.unload()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/tmp/shots")
