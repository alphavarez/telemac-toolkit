"""Build a small TELEMAC-2D flood case with telemac_io only.

A 2 km meandering river valley with floodplains; a flood hydrograph enters
upstream (prescribed Q) and the downstream stage is fixed (prescribed H).
Writes geometry (.slf), boundary conditions (.cli), liquid boundaries (.liq)
and a steering file (.cas). The geometry is placed at UTM-like coordinates to
exercise the float32 origin offset.

    python scripts/make_demo_case.py OUTDIR
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "telemac_toolkit"))

from telemac_io import BoundaryTable, SelafinFile, write_cli, write_selafin  # noqa: E402
from telemac_io.mesh import boundary_loops  # noqa: E402

X0, Y0 = 512000.0, 2155000.0  # UTM 14N-like origin (southern Mexico latitudes)
L, W = 2000.0, 600.0
DX = 12.5


def terrain(x, y):
    """Valley floor sloping 2 m/km, a meandering trapezoidal channel 3 m deep."""
    yc = 70.0 * np.sin(2 * np.pi * x / 900.0)
    d = np.abs(y - yc)
    floor = 10.0 - 0.002 * x + 0.006 * np.maximum(d - 25.0, 0.0) + 2.5e-5 * np.maximum(d - 25.0, 0.0) ** 2
    bank = np.clip((d - 18.0) / 8.0, 0.0, 1.0)  # 18 m flat bed, 8 m banks
    return floor - 3.0 * (1.0 - bank)


def build_mesh():
    nx = int(L / DX) + 1
    ny = int(W / DX) + 1
    xs = np.linspace(0, L, nx)
    ys = np.linspace(-W / 2, W / 2, ny)
    X, Y = np.meshgrid(xs, ys)  # (ny, nx)
    # gentle jitter on interior nodes so the mesh is not a perfect grid
    rng = np.random.default_rng(7)
    jit = rng.uniform(-0.18, 0.18, size=(2,) + X.shape) * DX
    jit[:, 0, :] = jit[:, -1, :] = 0
    jit[:, :, 0] = jit[:, :, -1] = 0
    X = X + jit[0]
    Y = Y + jit[1]
    x, y = X.ravel(), Y.ravel()
    tris = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            b, c, d = a + 1, a + nx + 1, a + nx
            if (i + j) % 2 == 0:
                tris += [[a, b, c], [a, c, d]]
            else:
                tris += [[a, b, d], [b, c, d]]
    return x, y, np.array(tris), nx, ny


def main(out):
    os.makedirs(out, exist_ok=True)
    x, y, ikle, nx, ny = build_mesh()
    z = terrain(x, y)
    geo = os.path.join(out, "geo_valley.slf")
    write_selafin(geo, x + X0, y + Y0, ikle, {"BOTTOM": z}, units={"BOTTOM": "M"},
                  title="Demo meandering valley - TELEMAC Toolkit")

    # initial state (hotstart): base flow 0.8 m deep in the channel, still water downstream
    yc = 70.0 * np.sin(2 * np.pi * x / 900.0)
    in_channel = np.abs(y - yc) < 22.0
    s0 = np.where(in_channel, np.maximum(4.4, z + 0.8), np.maximum(4.4, z))
    h0 = s0 - z
    zeros = np.zeros_like(z)
    write_selafin(os.path.join(out, "init_valley.slf"), x + X0, y + Y0, ikle,
                  {"VELOCITY U": zeros, "VELOCITY V": zeros, "WATER DEPTH": h0, "FREE SURFACE": s0, "BOTTOM": z},
                  units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M", "FREE SURFACE": "M", "BOTTOM": "M"},
                  title="Initial state - demo valley")

    # boundary conditions in boundary order (outer loop counter-clockwise)
    with SelafinFile(geo) as s:
        loop = boundary_loops(s.ikle, s.x, s.y)[0]
        bx = s.x[loop] - X0
    t = BoundaryTable.walls(loop)
    r = t.records
    eps = 1e-6
    upstream = np.abs(bx - 0.0) < eps
    downstream = np.abs(bx - L) < eps
    by = (np.asarray(SelafinFile(geo).y)[loop] - Y0)
    corners = np.abs(np.abs(by) - W / 2) < eps
    up = upstream & ~corners
    dn = downstream & ~corners
    r["lihbor"][up], r["liubor"][up], r["livbor"][up] = 4, 5, 5   # prescribed Q
    r["lihbor"][dn], r["liubor"][dn], r["livbor"][dn] = 5, 4, 4   # prescribed H
    for k in np.flatnonzero(up):
        t.comments[k] = "upstream"
    for k in np.flatnonzero(dn):
        t.comments[k] = "downstream"
    write_cli(os.path.join(out, "bc_valley.cli"), t)

    # liquid boundary numbering: TELEMAC numbers them in boundary order from K=1
    order = []
    for k in range(len(loop)):
        lab = "up" if up[k] else "dn" if dn[k] else None
        if lab and (not order or order[-1] != lab):
            order.append(lab)
    num = {lab: i + 1 for i, lab in enumerate(order)}
    hyd_t = np.array([0, 1800, 3600, 7200, 10800, 14400, 18000], dtype=float)
    hyd_q = np.array([25, 25, 60, 420, 160, 60, 30], dtype=float)
    stage = 4.4  # downstream water level (m)
    with open(os.path.join(out, "hydrograph.liq"), "w", newline="\n") as f:
        f.write("# Flood hydrograph upstream, fixed stage downstream\n")
        cols = sorted([("Q(%d)" % num["up"], "m3/s"), ("SL(%d)" % num["dn"], "m")], key=lambda c: int(c[0].split("(")[1].rstrip(")")))
        f.write("T " + " ".join(c[0] for c in cols) + "\n")
        f.write("s " + " ".join(c[1] for c in cols) + "\n")
        for tt, qq in zip(hyd_t, hyd_q):
            vals = {"Q": "%.1f" % qq, "SL": "%.2f" % stage}
            f.write("%g " % tt + " ".join(vals[c[0].split("(")[0]] for c in cols) + "\n")

    cas = """/ TELEMAC-2D demo case written by TELEMAC Toolkit (scripts/make_demo_case.py)
TITLE = 'Demo meandering valley flood'
GEOMETRY FILE = 'geo_valley.slf'
BOUNDARY CONDITIONS FILE = 'bc_valley.cli'
LIQUID BOUNDARIES FILE = 'hydrograph.liq'
RESULTS FILE = 'r2d_valley.slf'
VARIABLES FOR GRAPHIC PRINTOUTS = 'U,V,H,S,B,F,Q,M'
ORIGINAL DATE OF TIME = 2026;9;14
ORIGINAL HOUR OF TIME = 6;0;0
TIME STEP = 1.
DURATION = 18000.
GRAPHIC PRINTOUT PERIOD = 600
LISTING PRINTOUT PERIOD = 600
MASS-BALANCE = YES
PREVIOUS COMPUTATION FILE = 'init_valley.slf'
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
    with open(os.path.join(out, "t2d_valley.cas"), "w", newline="\n") as f:
        f.write(cas)
    print("nodes", x.size, "elements", ikle.shape[0], "liquid boundaries", num)
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "tests", "data", "demo_valley"))
