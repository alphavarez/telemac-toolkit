"""Render pictures of the mesher demo (run scripts/make_mesh_demo.py and TELEMAC first).

    QT_QPA_PLATFORM=offscreen python scripts/render_mesh_demo.py CASEDIR OUTDIR
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from qgis.core import (  # noqa: E402
    QgsApplication,
    QgsColorRampShader,
    QgsCoordinateReferenceSystem,
    QgsDateTimeRange,
    QgsMapRendererParallelJob,
    QgsMapSettings,
    QgsMeshLayer,
    QgsProject,
    QgsRasterLayer,
    QgsRasterShader,
    QgsRectangle,
    QgsSingleBandPseudoColorRenderer,
    QgsVectorLayer,
)
from qgis.PyQt.QtCore import QSize, Qt  # noqa: E402
from qgis.PyQt.QtGui import QColor  # noqa: E402

X0, Y0 = 512000.0, 2155000.0


def render(layers, extent, path, size=(1500, 560), crs=None, temporal=None):
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


def dem_style(lyr):
    from telemac_toolkit.qgis_bridge.styles import _terrain_ramp

    stats = lyr.dataProvider().bandStatistics(1)
    lo, hi = stats.minimumValue, stats.maximumValue
    ramp = _terrain_ramp()
    sh = QgsColorRampShader(lo, hi, ramp)
    sh.classifyColorRamp(9)
    rs = QgsRasterShader()
    rs.setRasterShaderFunction(sh)
    r = QgsSingleBandPseudoColorRenderer(lyr.dataProvider(), 1, rs)
    r.setClassificationMin(lo)
    r.setClassificationMax(hi)
    lyr.setRenderer(r)
    lyr.setOpacity(0.85)


def main(case, out):
    os.makedirs(out, exist_ok=True)
    app = QgsApplication([], True)
    app.initQgis()
    crs = QgsCoordinateReferenceSystem("EPSG:32614")
    QgsProject.instance().setCrs(crs)
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.qgis_bridge.results import load_result
    from telemac_toolkit.qgis_bridge import styles

    dem = QgsRasterLayer(os.path.join(case, "dem_valley.tif"), "DEM")
    dem_style(dem)
    gpkg = os.path.join(case, "valley_mesh.gpkg")
    vec = {}
    for name in ("domain", "channels", "boundaries"):
        v = QgsVectorLayer("%s|layername=%s" % (gpkg, name), name, "ogr")
        M.style_project_layer(v, name)
        vec[name] = v
    full = QgsRectangle(X0 - 60, Y0 - 330, X0 + 2060, Y0 + 330)
    render([vec["boundaries"], vec["channels"], vec["domain"], dem], full, os.path.join(out, "mesher_inputs.png"), crs=crs)

    mesh = QgsMeshLayer(os.path.join(case, "geo_valley_mesh.slf"), "mesh", "mdal")
    mesh.setCrs(crs)
    M.style_geometry_mesh(mesh)
    zoom = QgsRectangle(X0 + 20, Y0 - 150, X0 + 760, Y0 + 150)
    render([vec["channels"], mesh], zoom, os.path.join(out, "mesher_mesh.png"), size=(1500, 610), crs=crs)

    res_path = os.path.join(case, "r2d_valley_mesh.slf")
    if os.path.exists(res_path):
        rs = load_result(res_path, crs=crs, add_to_project=True)
        k = rs.var_index("water_depth")
        styles.set_active_variable(rs, k)
        peak = int(rs.peak_frame(k))
        t = rs.reference_datetime().addSecs(int(rs.times[peak]))
        settings = rs.layer.rendererSettings()
        ms = settings.nativeMeshSettings()
        ms.setEnabled(True)
        ms.setColor(QColor(60, 60, 70, 60))
        ms.setLineWidth(0.1)
        settings.setNativeMeshSettings(ms)
        rs.layer.setRendererSettings(settings)
        render([rs.layer, dem], full, os.path.join(out, "mesher_result_peak.png"), crs=crs,
               temporal=QgsDateTimeRange(t, t.addSecs(1)))
        print("peak", rs.time_label(peak))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/tmp/shots")
