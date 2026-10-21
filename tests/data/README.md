# Test data

| File | Origin | Licence |
| --- | --- | --- |
| `r2d_tidal_flats.slf`, `r3d_tidal_flats.slf`, `r3d_bump.slf`, `r1d_tomsail.slf` | TELEMAC-MASCARET example outputs, as redistributed in the test suite of [xarray-selafin](https://github.com/seareport/xarray-selafin) | GPL-3.0 |
| `demo_valley/*` | Built with `scripts/make_demo_case.py` (this repository, `telemac_io` only) and run with openTELEMAC v9.1.1 from conda-forge | GPL-3.0 |

`demo_valley` is a 2 km meandering river valley with floodplains: a flood hydrograph
(25 → 420 → 30 m³/s) enters upstream and the downstream stage is fixed at 4.4 m.
TELEMAC-2D ran 5 h of simulated time in 81 s with a volume loss of 5 × 10⁻⁹ m³
(`run_summary.txt`). The geometry uses UTM-like coordinates (x ≈ 512 000, y ≈ 2 155 000)
written with an origin offset, so the file also tests float32 precision handling.
TELEMAC reported supercritical inflow at a few upstream boundary nodes during the
peak; that is a property of this synthetic case, not of the files.

Regenerate:

    python scripts/make_demo_case.py tests/data/demo_valley
    cd tests/data/demo_valley && telemac2d.py t2d_valley.cas

| `demo_mesh/*` | Steering file, `.cli`, `.liq` and run summary of the valley rebuilt in QGIS with `scripts/make_mesh_demo.py` (layers + DEM → Generate mesh) | GPL-3.0 |
| `demo_flume/*` | Steering file, `.cli`, `.liq` and run summary of a laboratory flume built in local coordinates with `scripts/make_flume_demo.py` | GPL-3.0 |

`demo_mesh` checks that a mesh, its bottom and its boundary conditions made only from QGIS
layers run in TELEMAC-2D. `demo_flume` checks local coordinates from (0, 0), the bed
formula and the friction written into the geometry file: TELEMAC reads `BOTTOM FRICTION`
from the geometry (n = 0.025) instead of the steering-file value (0.05), and the flow settles
at the normal depth for n = 0.025 (0.1448 m simulated, 0.1439 m by Manning's formula).
Both scripts regenerate the geometry; run TELEMAC in the output folder:

    QT_QPA_PLATFORM=offscreen <qgis-python> scripts/make_flume_demo.py OUTDIR
    cd OUTDIR && telemac2d.py t2d_flume.cas

| `demo_bend/*` | Steering file, `.cli`, `.liq`, run summary and checks of a 180° laboratory bend built with the channel template, `scripts/make_bend_demo.py` | GPL-3.0 |

`demo_bend` checks the curved channel template (structured cells normal to the axis, bed as a
formula of the distance along the axis), TELEMAC's numbering of liquid boundaries (the outlet
is boundary 1 here, because TELEMAC starts from the boundary point with the smallest x + y),
the discharge through sections (0.04999 to 0.05001 m³/s for 0.05 m³/s in) and the water
surface superelevation across the bend (6.3 mm, 6.0 mm by Rozovskii's U²W/(gR)):

    QT_QPA_PLATFORM=offscreen <qgis-python> scripts/make_bend_demo.py OUTDIR
    cd OUTDIR && telemac2d.py t2d_bend.cas
    QT_QPA_PLATFORM=offscreen <qgis-python> scripts/check_bend_demo.py OUTDIR
