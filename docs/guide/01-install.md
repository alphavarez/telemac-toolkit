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
