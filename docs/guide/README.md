# TELEMAC Toolkit for QGIS — User guide

Version 0.4.0 · QGIS 3.44 LTR and QGIS 4.x · Windows, Linux, macOS

TELEMAC Toolkit builds meshes and analyses results for TELEMAC inside QGIS, with nothing but QGIS and
Python. The mesher is made for TELEMAC first and writes the same mesh for other 2D programs (BASEMENT,
TUFLOW FV, Iber, Delft3D FM, MIKE 21 FM, ADCIRC/SCHISM, HEC-RAS).

![Results panel with velocity arrows on a 180° laboratory bend](../img/results_panel_arrows.png)

## Contents

1. [Install](01-install.md)
2. [First steps: a result in five minutes, a flume in ten](02-first-steps.md)
3. [Viewing results](03-results.md): opening files, coordinates, time, arrows, cursor, time series, profiles
4. [Calculator, integrals and fluxes](04-calculator-integrals.md)
5. [The mesh workflow](05-mesh-workflow.md): project, drawing, bottom, roughness, boundaries, meshing, quality
6. [Channel and flume templates](06-channel-templates.md): straight, bend, meanders, bed formulas in s and n
7. [Editing an existing mesh](07-edit-existing-mesh.md): bathymetry, dikes, dredging, friction
8. [Formats, coordinates and boundary conditions](08-formats-coordinates-boundaries.md)
9. [Processing tools reference](09-processing-reference.md) (every tool and parameter)
10. [Scripting with Python](10-python.md)
11. [Troubleshooting and FAQ](11-troubleshooting.md)

## What is where

| You want to | Use |
| --- | --- |
| Open and animate a TELEMAC result | Toolbar *Open TELEMAC result…* or drag the `.slf` onto QGIS; the **TELEMAC Results** panel |
| Read values, plot time series, draw profiles | Results panel tabs *Cursor*, *Time series*, *Profile* |
| Compute new variables (Froude, shear stress, vorticity, maxima) | Results panel → calculator button, or *Plugins → TELEMAC Toolkit → Calculator* |
| Discharge through a section, volumes | *Plugins → TELEMAC Toolkit → Integrals and fluxes* |
| Build a mesh | Toolbar *Mesh workflow panel* |
| Build a laboratory flume, a bend or meanders | Mesh workflow panel, step 1 → *Channel or flume template* |
| Change the bottom (or friction) of an existing mesh | *Plugins → TELEMAC Toolkit → Mesh → Edit node values of a mesh* |
| Convert a mesh to another program | *Mesh → Export or convert mesh* |
| Run anything in batch or in a model | Processing toolbox → **TELEMAC Toolkit** |

## Words used in this guide

- **Result**: a SELAFIN file (`.slf`, `.res`) written by TELEMAC, with one or more time steps (frames).
- **Geometry file**: the SELAFIN file TELEMAC reads as mesh, with the `BOTTOM` variable and possibly
  `BOTTOM FRICTION`.
- **Mesh project**: a GeoPackage made by the plugin with one layer per kind of object you draw (domain,
  refinement areas, channels, boundaries…) and the meshing settings.
- **Local coordinates**: a metric grid with (0, 0) at the origin and no place on Earth, for laboratory
  flumes and schematic channels.
- **Channel coordinates s and n**: distance along a channel axis line (s) and across it (n, positive on the
  left looking downstream).
