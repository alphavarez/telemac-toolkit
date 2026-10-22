# 10. Scripting with Python

## Processing tools from the QGIS Python console

Every tool runs with `processing.run`; the keys are in the [reference](09-processing-reference.md).

```python
import processing

# a 180° bend flume in a new project, meshed for TELEMAC
processing.run("telemactoolkit:channeltemplate", {
    "NEW_PROJECT": "C:/work/bend.gpkg", "PLANFORM": 1, "WIDTH": 1.0, "INLET": 3, "OUTLET": 3,
    "RADIUS": 2.0, "BEND_ANGLE": 180, "TURN": 0, "CELL": 0.05, "STRUCTURED": True,
    "BED_Z0": 0.30, "SLOPE": 0.001, "Q": 0.05, "H": 0.43})
processing.run("telemactoolkit:generatemesh", {
    "PROJECT": "C:/work/bend.gpkg", "FRICTION": 0.025, "OUTPUT": "C:/work/geo_bend.slf",
    "OUTPUT_LIQ": "C:/work/bend.liq"})

# maxima of a result
processing.run("telemactoolkit:calculator", {
    "INPUT": "C:/work/r2d_bend.slf", "EXPRESSIONS": "HMAX [M] = tmax(H)\nVMAX [M/S] = tmax(mag(U, V))",
    "KEEP": False, "OUTPUT": "C:/work/r2d_bend_max.slf"})

# discharge through the lines of a layer, for every time step
res = processing.run("telemactoolkit:integrals", {
    "INPUT": "C:/work/r2d_bend.slf", "MODE": 0, "LINES": "C:/work/sections.gpkg", "OUTPUT": "memory:"})
for f in res["OUTPUT"].getFeatures():
    print(f["name"], f["time_s"], f["discharge"])
```

The full demo cases are scripts: `scripts/make_flume_demo.py`, `scripts/make_bend_demo.py` (+
`check_bend_demo.py`), `scripts/make_mesh_demo.py` (a valley from layers and a DEM).

## From the command line, without opening QGIS

`qgis_process` runs the same tools from a terminal or a batch file. It comes with QGIS (on Windows,
`qgis_process-qgis.bat` or `qgis_process-qgis-ltr.bat` in the `bin` folder of the QGIS installation) and
finds the plugin in your QGIS profile once it is enabled for it:

```bat
qgis_process plugins enable telemac_toolkit
qgis_process run telemactoolkit:calculator -- INPUT=C:/work/r2d_bend.slf "EXPRESSIONS=HMAX [M] = tmax(H)" OUTPUT=C:/work/r2d_bend_max.slf
qgis_process help telemactoolkit:generatemesh
```

With many parameters, pass them as JSON on standard input:

```bat
echo {"inputs": {"PROJECT": "C:/work/bend.gpkg", "FRICTION": 0.025, "OUTPUT": "C:/work/geo_bend.slf"}} | qgis_process run telemactoolkit:generatemesh -
```

This is the way to chain the plugin with a TELEMAC run in a script (build the mesh, run `telemac2d.py`,
compute maxima and discharges), for example with TELEMAC in a Docker container and the files in a shared
folder.

## The libraries without QGIS

`telemac_toolkit/telemac_io` (SELAFIN, `.cli`, calculator, integrals) and `telemac_toolkit/mesher`
(generator, planforms, writers) need only Python and NumPy (the mesher also SciPy, for the triangulation),
so they also run inside a TELEMAC container or on a cluster, next to the solver:

```python
import sys; sys.path.insert(0, "path/to/telemac_toolkit")
from telemac_io import SelafinFile, write_selafin
from telemac_io.calculator import run_calculator, MeshOperators
from telemac_io.integrals import Section

with SelafinFile("r2d.slf") as slf:
    print(slf)                                     # dimensions, variables, frames
    depth = slf.read("WATER DEPTH", frame=-1)      # last frame
    ts = slf.series("FREE SURFACE", nodes=[10])    # time series, read lazily
    h = slf.header
    ops = MeshOperators(h.x, h.y, h.ikle)
    w = ops.curl(slf.read("VELOCITY U", -1), slf.read("VELOCITY V", -1))
    sec = Section.build([(1.5, 0.5), (1.5, -0.5)], h.x, h.y, h.ikle)
    q = sec.discharge(depth, slf.read("VELOCITY U", -1), slf.read("VELOCITY V", -1))

run_calculator("r2d.slf", "r2d_calc.slf", "FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))")
write_selafin("geo.slf", x, y, ikle, {"BOTTOM": z}, units={"BOTTOM": "M"})
```

The reader takes single and double precision, big and little endian, 2D and 3D files, files with a date
record and incomplete files from running simulations.

[← Processing reference](09-processing-reference.md) · [Contents](README.md) · [Next: Troubleshooting →](11-troubleshooting.md)
