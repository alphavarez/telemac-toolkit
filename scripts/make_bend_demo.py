"""A 180-degree laboratory bend built with the channel template, in local coordinates.

    QT_QPA_PLATFORM=offscreen python scripts/make_bend_demo.py OUTDIR

Rectangular flume 1 m wide: 3 m straight inlet, bend of radius 2 m on the axis
over 180 degrees, 3 m straight outlet; bed slope 0.001 along the axis (a
formula of the channel coordinate s), structured 5 cm cells normal to the
axis, Manning n = 0.025 written into the geometry file, Q = 0.05 m3/s.

Checks after the TELEMAC run (scripts/check_bend_demo.py): the discharge
through sections along the bend equals the inflow, and the water surface is
higher on the outer bank by about U^2 W / (g R) (Rozovskii's estimate).
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

W, R, ANGLE, L_IN, L_OUT = 1.0, 2.0, 180.0, 3.0, 3.0
S, Z0, Q, N, CELL = 0.001, 0.30, 0.05, 0.025, 0.05
LENGTH = L_IN + R * math.radians(ANGLE) + L_OUT
h_n = (Q / W * N / math.sqrt(S)) ** 0.6  # normal depth, wide channel
Z_OUT = Z0 - S * LENGTH


def main(out):
    os.makedirs(out, exist_ok=True)
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    import processing

    gpkg = os.path.join(out, "bend.gpkg")
    if os.path.exists(gpkg):
        os.remove(gpkg)
    r = processing.run("telemactoolkit:channeltemplate", {
        "NEW_PROJECT": gpkg, "PLANFORM": 1, "WIDTH": W, "INLET": L_IN, "OUTLET": L_OUT, "RADIUS": R, "BEND_ANGLE": ANGLE,
        "TURN": 0, "CELL": CELL, "STRUCTURED": True, "BED_Z0": Z0, "SLOPE": S, "SECTION": 0, "Q": Q, "H": Z_OUT + h_n})
    print("bed formula:", r["BED_FORMULA"])
    geo = os.path.join(out, "geo_bend.slf")
    processing.run("telemactoolkit:generatemesh", {"PROJECT": gpkg, "FRICTION": N, "OUTPUT": geo,
                                                   "OUTPUT_LIQ": os.path.join(out, "bend.liq"), "LOAD": False})
    cas = """/ 180-degree laboratory bend in local coordinates, built with TELEMAC Toolkit (scripts/make_bend_demo.py)
TITLE = 'Bend R = 2 m, 180 deg, W = 1 m'
GEOMETRY FILE = 'geo_bend.slf'
BOUNDARY CONDITIONS FILE = 'geo_bend.cli'
LIQUID BOUNDARIES FILE = 'bend.liq'
RESULTS FILE = 'r2d_bend.slf'
VARIABLES FOR GRAPHIC PRINTOUTS = 'U,V,H,S,B,W'
TIME STEP = 0.02
NUMBER OF TIME STEPS = 30000
GRAPHIC PRINTOUT PERIOD = 1500
LISTING PRINTOUT PERIOD = 1500
MASS-BALANCE = YES
INITIAL CONDITIONS = 'CONSTANT ELEVATION'
INITIAL ELEVATION = %.4f
LAW OF BOTTOM FRICTION = 4
FRICTION COEFFICIENT = 0.05
TURBULENCE MODEL = 1
VELOCITY DIFFUSIVITY = 1.E-4
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
    with open(os.path.join(out, "t2d_bend.cas"), "w", newline="\n") as f:
        f.write(cas)
    print("axis length %.4f m, normal depth %.4f m, outlet level %.4f m" % (LENGTH, h_n, Z_OUT + h_n))
    return out


if __name__ == "__main__":
    app = QgsApplication([], False)
    app.initQgis()
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "tests", "data", "demo_bend"))
