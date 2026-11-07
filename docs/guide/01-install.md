# 1. Install

## Requirements

- QGIS 3.44 LTR (Qt5) or QGIS 4.x (Qt6), on Windows, Linux or macOS.
- Nothing else. The plugin uses only Python and NumPy, which ship with QGIS. TELEMAC itself is not needed
  to build meshes or view results; you need it only to run simulations.

## From the QGIS plugin repository

*Plugins → Manage and Install Plugins → All*, search for **TELEMAC Toolkit**, *Install Plugin*. While the
plugin is marked experimental, tick *Show also experimental plugins* in the *Settings* tab of the plugin
manager first. QGIS tells you when a new version is available (*Upgradeable* tab).

## From a zip file

1. Download `telemac_toolkit-<version>.zip` (the installable zip, not the source zip).
2. *Plugins → Manage and Install Plugins → Install from ZIP*, choose the file, *Install Plugin*.
3. If an older version was installed, it is replaced. Restart QGIS so the panels and menus load cleanly.

## Optional: Gmsh

The Gmsh meshing methods need the **gmsh** Python package (GPL-2.0-or-later, like this plugin's licence) in
QGIS's own Python. Everything else works without it, and the panel greys out the Gmsh methods until it is
installed.

- **Windows**: open the *OSGeo4W Shell* that comes with QGIS and run `python -m pip install gmsh`.
- **macOS**: `/Applications/QGIS.app/Contents/MacOS/bin/python3 -m pip install gmsh` (the path of your QGIS
  application).
- **Linux**: `python3 -m pip install gmsh` with the Python QGIS uses; the wheel also needs the system library
  libGLU (`libglu1-mesa` on Debian and Ubuntu).

Restart QGIS afterwards. *Compare meshing methods* tells you whether Gmsh was found.

## Optional: MMG

*Optimise or adapt a mesh (MMG)* runs MMG's standalone program `mmg2d_O3` (LGPL); nothing is installed into
QGIS's Python. Get it from [mmgtools.org](https://www.mmgtools.org) (Windows, macOS and Linux builds) or from
conda-forge (`mmgsuite`), then give the path of `mmg2d_O3` the first time you run the tool; it is remembered.
Putting it on the PATH, or in an `MMG2D` environment variable, works too. MMG's Python package (mmgpy) is not
used: it needs a newer NumPy than some QGIS installations ship.

## What appears

- A **TELEMAC Toolkit** toolbar: open result, Results panel, live cursor, probe, profile, Mesh workflow
  panel, new mesh project.
- A **Plugins → TELEMAC Toolkit** menu with the same tools, the calculator, integrals and fluxes, SELAFIN
  information, 3D extraction and a **Mesh** submenu with every mesh tool.
- A **TELEMAC Toolkit** provider in the Processing toolbox with all the tools ([reference](09-processing-reference.md)).

## Updating and uninstalling

Updates come through the plugin manager (or by installing a newer zip). Mesh projects made with 0.2.0
and later keep working: older layer names (`outline`, `density`, `hard_lines`, `soft_lines`,
`hard_points`, `grids`) are read as their current equivalents. To uninstall: plugin manager → *Installed*
→ *Uninstall plugin*. Your files are not touched.

[← Contents](README.md) · [Next: First steps →](02-first-steps.md)
