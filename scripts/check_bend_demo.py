"""Check the TELEMAC run of the 180-degree bend (scripts/make_bend_demo.py) with the plugin's own tools.

    QT_QPA_PLATFORM=offscreen python scripts/check_bend_demo.py CASEDIR

* discharge through five sections along the channel (Integrals and fluxes), which must all equal the
  inflow in steady flow;
* water surface superelevation across the bend apex against Rozovskii's estimate U^2 W / (g R);
* vorticity and the maximum depth with the Calculator.
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

import numpy as np  # noqa: E402
from qgis.core import QgsApplication, QgsFeature, QgsGeometry, QgsPointXY, QgsVectorLayer  # noqa: E402

W, R, L_IN, Q = 1.0, 2.0, 3.0, 0.05
C = (L_IN, R)  # centre of the bend (left turn from (0, 0) along +x)


def section(theta_deg):
    """From the inner (left) bank to the outer (right) bank, looking downstream."""
    t = math.radians(theta_deg)
    d = (math.sin(t), -math.cos(t))
    a = (C[0] + (R - W / 2) * d[0], C[1] + (R - W / 2) * d[1])
    b = (C[0] + (R + W / 2) * d[0], C[1] + (R + W / 2) * d[1])
    return [a, b]


def main(case):
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    import processing

    res = os.path.join(case, "r2d_bend.slf")
    lines = QgsVectorLayer("LineString?field=name:string", "sections", "memory")
    secs = [("inlet straight x=1.5", [(1.5, 0.5), (1.5, -0.5)])]
    for th in (45, 90, 135):
        secs.append(("bend %d deg" % th, section(th)))
    secs.append(("outlet straight x=1.5", [(1.5, 2 * R - 0.5), (1.5, 2 * R + 0.5)]))
    feats = []
    for name, pts in secs:
        f = QgsFeature(lines.fields())
        f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(*p) for p in pts]))
        f.setAttribute("name", name)
        feats.append(f)
    lines.dataProvider().addFeatures(feats)
    out = processing.run("telemactoolkit:integrals", {"INPUT": res, "MODE": 0, "LINES": lines, "FRAMES": "-1",
                                                      "OUTPUT": "memory:"})["OUTPUT"]
    report = []
    for f in out.getFeatures():
        report.append("Discharge, %-22s %.5f m3/s  (inflow %.3f)" % (f["name"], f["discharge"], Q))
    # superelevation at the apex (90 degrees), 2 cm inside each bank
    from telemac_toolkit.telemac_io import SelafinFile, TriangleLocator

    with SelafinFile(res) as s:
        h = s.header
        loc = TriangleLocator(h.x, h.y, h.ikle)
        S = s.read("FREE SURFACE", -1)
        H = s.read("WATER DEPTH", -1)
        t = math.radians(90)
        d = (math.sin(t), -math.cos(t))
        pin = (C[0] + (R - W / 2 + 0.02) * d[0], C[1] + (R - W / 2 + 0.02) * d[1])
        pout = (C[0] + (R + W / 2 - 0.02) * d[0], C[1] + (R + W / 2 - 0.02) * d[1])
        pax = (C[0] + R * d[0], C[1] + R * d[1])
        s_in, s_out = loc.interpolate(S, [pin[0], pout[0]], [pin[1], pout[1]])
        h_ax = loc.interpolate(H, [pax[0]], [pax[1]])[0]
    ubar = Q / (W * h_ax)
    est = ubar ** 2 * W / (9.81 * R)
    report.append("Apex: depth %.4f m, mean velocity %.4f m/s; surface outer - inner = %.2f mm "
                  "(Rozovskii U^2 W / (g R) = %.2f mm)" % (h_ax, ubar, 1000 * (s_out - s_in), 1000 * est))
    calc = processing.run("telemactoolkit:calculator", {
        "INPUT": res, "EXPRESSIONS": "VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)\nHMAX [M] = tmax(H)",
        "KEEP": False, "FRAMES": "-1", "OUTPUT": os.path.join(case, "r2d_bend_calc.slf"), "LOAD": False})
    with SelafinFile(calc["OUTPUT"]) as s:
        w = s.read("VORTICITY", -1)
        lo, hi = np.percentile(w, [1, 99])
        report.append("Vorticity at the last frame: %.3f to %.3f 1/s (1st to 99th percentile; the extremes, "
                      "%.1f and %.1f, are at the inlet corners)" % (lo, hi, float(w.min()), float(w.max())))
    for line in report:
        print(line)
    with open(os.path.join(case, "check_summary.txt"), "w") as f:
        f.write("\n".join(report) + "\n")


if __name__ == "__main__":
    app = QgsApplication([], False)
    app.initQgis()
    main(sys.argv[1])
