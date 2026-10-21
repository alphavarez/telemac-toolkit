"""A laboratory flume built in local coordinates, from (0, 0), with the Processing tools only.

    QT_QPA_PLATFORM=offscreen python scripts/make_flume_demo.py OUTDIR

20 m x 1 m rectangular flume, bed slope 0.001, structured 5 cm grid, friction
written into the geometry file (BOTTOM FRICTION = 0.025, Manning) while the
steering file says 0.05: if TELEMAC reads the geometry value, the flow depth
in the middle of the flume is the normal depth for n = 0.025 (0.144 m), not
the one for n = 0.05 (0.218 m).
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _p in (os.path.join(sys.prefix, "share", "qgis", "python", "plugins"),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.append(_p)

from qgis.core import QgsApplication  # noqa: E402

L, W, S, Z0 = 20.0, 1.0, 0.001, 0.30
Q = 0.05  # m3/s
N = 0.025  # Manning, geometry file
h_n = (Q / W * N / math.sqrt(S)) ** 0.6  # normal depth, wide channel
Z_OUT = Z0 - S * L


def main(out):
    os.makedirs(out, exist_ok=True)
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    import processing

    gpkg = os.path.join(out, "flume.gpkg")
    if os.path.exists(gpkg):
        os.remove(gpkg)
    r = processing.run("telemactoolkit:channeltemplate", {
        "NEW_PROJECT": gpkg, "LENGTH": L, "WIDTH": W, "CELL": 0.05, "STRUCTURED": True, "BED_Z0": Z0, "SLOPE": S,
        "SECTION": 0, "Q": Q, "H": Z_OUT + h_n})
    print("bed formula:", r["BED_FORMULA"])
    geo = os.path.join(out, "geo_flume.slf")
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": gpkg, "FRICTION": N, "OUTPUT": geo,
                                                       "OUTPUT_LIQ": os.path.join(out, "flume.liq"), "LOAD": False})
    print(g)
    cas = """/ Laboratory flume in local coordinates, built with TELEMAC Toolkit (scripts/make_flume_demo.py)
TITLE = 'Flume 20 m x 1 m, slope 0.001'
GEOMETRY FILE = 'geo_flume.slf'
BOUNDARY CONDITIONS FILE = 'geo_flume.cli'
LIQUID BOUNDARIES FILE = 'flume.liq'
RESULTS FILE = 'r2d_flume.slf'
VARIABLES FOR GRAPHIC PRINTOUTS = 'U,V,H,S,B,W'
TIME STEP = 0.05
NUMBER OF TIME STEPS = 24000
GRAPHIC PRINTOUT PERIOD = 1200
LISTING PRINTOUT PERIOD = 1200
MASS-BALANCE = YES
INITIAL CONDITIONS = 'CONSTANT ELEVATION'
INITIAL ELEVATION = %.4f
LAW OF BOTTOM FRICTION = 4
FRICTION COEFFICIENT = 0.05
TURBULENCE MODEL = 1
VELOCITY DIFFUSIVITY = 1.E-3
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
""" % (Z_OUT + h_n)
    with open(os.path.join(out, "t2d_flume.cas"), "w", newline="\n") as f:
        f.write(cas)
    print("normal depth %.4f m, outlet level %.4f m" % (h_n, Z_OUT + h_n))
    return out


if __name__ == "__main__":
    app = QgsApplication([], False)
    app.initQgis()
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "tests", "data", "demo_flume"))
