"""Rebuild the demo valley case with the TELEMAC Toolkit mesher, from QGIS layers and a DEM.

Run with QGIS's Python (headless is fine):

    QT_QPA_PLATFORM=offscreen python scripts/make_mesh_demo.py OUTDIR

Writes a DEM GeoTIFF, a mesh project GeoPackage (domain, channel banks and
axis, boundary lines), then runs the Processing algorithm "Generate mesh"
(unstructured floodplain + structured channel, bottom from the DEM, .cli and .liq), and
finally a hotstart file and a TELEMAC-2D steering file so the case runs as is.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _p in (os.path.join(sys.prefix, "share", "qgis", "python", "plugins"),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.append(_p)

from osgeo import gdal, osr  # noqa: E402
from qgis.core import (  # noqa: E402
    QgsApplication,
    QgsCoordinateReferenceSystem,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsRasterLayer,
)

X0, Y0 = 512000.0, 2155000.0
L, W = 2000.0, 600.0
EPSG = 32614


def terrain(x, y):
    """Same valley as scripts/make_demo_case.py (relative coordinates)."""
    yc = 70.0 * np.sin(2 * np.pi * x / 900.0)
    d = np.abs(y - yc)
    floor = 10.0 - 0.002 * x + 0.006 * np.maximum(d - 25.0, 0.0) + 2.5e-5 * np.maximum(d - 25.0, 0.0) ** 2
    bank = np.clip((d - 18.0) / 8.0, 0.0, 1.0)
    return floor - 3.0 * (1.0 - bank)


def write_dem(path, res=2.0, margin=40.0):
    nx = int((L + 2 * margin) / res)
    ny = int((W + 2 * margin) / res)
    xs = -margin + (np.arange(nx) + 0.5) * res
    ys = W / 2 + margin - (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    Z = terrain(X, Y).astype(np.float32)
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, nx, ny, 1, gdal.GDT_Float32, options=["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform((X0 - margin, res, 0.0, Y0 + W / 2 + margin, 0.0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.WriteArray(Z)
    band.SetNoDataValue(-9999.0)
    ds = None
    return path


def add(layer, geom, **attrs):
    f = QgsFeature(layer.fields())
    f.setGeometry(geom)
    for k, v in attrs.items():
        f.setAttribute(k, v)
    layer.startEditing()
    layer.addFeature(f)
    if not layer.commitChanges():
        raise RuntimeError(layer.commitErrors())


def line(pts):
    return QgsGeometry.fromPolylineXY([QgsPointXY(X0 + a, Y0 + b) for a, b in pts])


def main(out, run_alg=True):
    os.makedirs(out, exist_ok=True)
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    import processing

    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile, read_cli, write_selafin
    from telemac_toolkit.telemac_io.cli import liquid_boundaries

    crs = QgsCoordinateReferenceSystem("EPSG:%d" % EPSG)
    dem = write_dem(os.path.join(out, "dem_valley.tif"))
    gpkg = os.path.join(out, "valley_mesh.gpkg")
    if os.path.exists(gpkg):
        os.remove(gpkg)
    lay = M.create_mesh_project(gpkg, crs, add_to_project=False)

    add(lay["domain"], QgsGeometry.fromPolygonXY([[QgsPointXY(X0 + a, Y0 + b) for a, b in
                                                     ((0, -W / 2), (L, -W / 2), (L, W / 2), (0, W / 2), (0, -W / 2))]]),
        name="valley")
    xs = np.linspace(0.0, L, 201)
    yc = 70.0 * np.sin(2 * np.pi * xs / 900.0)
    half = 24.0  # channel strip half width: bed (18 m) plus most of the banks
    add(lay["channels"], line(zip(xs, yc - half)), channel="river", role="left", nodes_across=9, spacing=8.0)
    add(lay["channels"], line(zip(xs, yc + half)), channel="river", role="right", nodes_across=9, spacing=8.0)
    add(lay["channels"], line(zip(xs, yc)), channel="river", role="axis", nodes_across=9, spacing=8.0)
    # finer floodplain next to the river, coarse far away (growth ratio does the rest)
    add(lay["boundaries"], line([(0, -280), (0, 280)]), name="upstream", type="Q", value=25.0)
    add(lay["boundaries"], line([(L, -280), (L, 280)]), name="downstream", type="H", value=4.4)

    geo = os.path.join(out, "geo_valley_mesh.slf")
    cli = os.path.join(out, "geo_valley_mesh.cli")  # written next to the geometry
    liq = os.path.join(out, "valley_mesh.liq")
    if run_alg:
        res = processing.run("telemactoolkit:generatemesh", {
            "DOMAIN": lay["domain"], "CHANNELS": lay["channels"], "BOUNDARIES": lay["boundaries"],
            "DEFAULT_SIZE": 30.0, "GROWTH": 1.15, "SMOOTHING": 8,
            "RASTERS": [QgsRasterLayer(dem, "dem")], "Z_FACTOR": 1.0,
            "OUTPUT": geo, "OUTPUT_LIQ": liq, "LOAD": False,
        })
        print(res)

    # hydrograph with TELEMAC's liquid boundary numbers, read back from the .cli
    from telemac_toolkit.mesher.boundary import attach_mesh

    with SelafinFile(geo) as g:
        table = attach_mesh(read_cli(cli), g.header.x, g.header.y, g.header.ikle)
    nums = {table.comments[b["start"]]: b["number"] for b in liquid_boundaries(table)}
    print("liquid boundaries", nums)
    hyd_t = [0, 1800, 3600, 7200, 10800, 14400, 18000]
    hyd_q = [25, 25, 60, 420, 160, 60, 30]
    cols = sorted([("Q(%d)" % nums["upstream"], "m3/s", "Q"), ("SL(%d)" % nums["downstream"], "m", "SL")],
                  key=lambda c: int(c[0].split("(")[1].rstrip(")")))
    with open(os.path.join(out, "hydrograph_mesh.liq"), "w", newline="\n") as f:
        f.write("# Flood hydrograph upstream, fixed stage downstream (mesh built in QGIS)\n")
        f.write("T " + " ".join(c[0] for c in cols) + "\n")
        f.write("s " + " ".join(c[1] for c in cols) + "\n")
        for tt, qq in zip(hyd_t, hyd_q):
            v = {"Q": "%.1f" % qq, "SL": "4.40"}
            f.write("%g " % tt + " ".join(v[c[2]] for c in cols) + "\n")

    # hotstart: base flow 0.8 m deep in the channel, still water downstream
    with SelafinFile(geo) as s:
        x, y, ik = s.x.copy(), s.y.copy(), s.ikle.copy()
        z = s.read("BOTTOM", 0)
    xr, yr = x - X0, y - Y0
    in_channel = np.abs(yr - 70.0 * np.sin(2 * np.pi * xr / 900.0)) < 22.0
    s0 = np.where(in_channel, np.maximum(4.4, z + 0.8), np.maximum(4.4, z))
    zeros = np.zeros_like(z)
    write_selafin(os.path.join(out, "init_valley_mesh.slf"), x, y, ik,
                  {"VELOCITY U": zeros, "VELOCITY V": zeros, "WATER DEPTH": s0 - z, "FREE SURFACE": s0, "BOTTOM": z},
                  units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M", "FREE SURFACE": "M", "BOTTOM": "M"},
                  title="Initial state - demo valley (QGIS mesh)")
    cas = """/ TELEMAC-2D demo case: mesh, bottom and boundary conditions built in QGIS (scripts/make_mesh_demo.py)
TITLE = 'Demo meandering valley flood - QGIS mesh'
GEOMETRY FILE = 'geo_valley_mesh.slf'
BOUNDARY CONDITIONS FILE = 'geo_valley_mesh.cli'
LIQUID BOUNDARIES FILE = 'hydrograph_mesh.liq'
RESULTS FILE = 'r2d_valley_mesh.slf'
VARIABLES FOR GRAPHIC PRINTOUTS = 'U,V,H,S,B,F,Q,M'
ORIGINAL DATE OF TIME = 2026;9;14
ORIGINAL HOUR OF TIME = 6;0;0
TIME STEP = 1.
DURATION = 18000.
GRAPHIC PRINTOUT PERIOD = 600
LISTING PRINTOUT PERIOD = 600
MASS-BALANCE = YES
PREVIOUS COMPUTATION FILE = 'init_valley_mesh.slf'
INITIAL TIME SET TO ZERO = YES
LAW OF BOTTOM FRICTION = 4
FRICTION COEFFICIENT = 0.035
TURBULENCE MODEL = 1
VELOCITY DIFFUSIVITY = 0.5
TIDAL FLATS = YES
OPTION FOR THE TREATMENT OF TIDAL FLATS = 1
TREATMENT OF NEGATIVE DEPTHS = 2
MASS-LUMPING ON H = 1.
CONTINUITY CORRECTION = YES
SUPG OPTION = 0;0
SCHEME FOR ADVECTION OF VELOCITIES = 1
FREE SURFACE GRADIENT COMPATIBILITY = 0.9
SOLVER = 1
SOLVER ACCURACY = 1.E-6
MAXIMUM NUMBER OF ITERATIONS FOR SOLVER = 200
INFORMATION ABOUT SOLVER = NO
"""
    with open(os.path.join(out, "t2d_valley_mesh.cas"), "w", newline="\n") as f:
        f.write(cas)
    print("nodes", len(x), "elements", len(ik))
    return out


if __name__ == "__main__":
    app = QgsApplication([], False)
    app.initQgis()
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "tests", "data", "demo_mesh"))
