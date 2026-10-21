# TELEMAC Toolkit for QGIS

Build meshes and analyse results for TELEMAC inside QGIS, with nothing but QGIS and Python; the same
mesh can be written for BASEMENT, TUFLOW FV, Iber, Delft3D FM, MIKE 21 FM, ADCIRC/SCHISM and HEC-RAS.

- Results: SELAFIN results (2D and 3D) as styled mesh layers with the QGIS Temporal Controller, arrow
  controls, live cursor, probe time series, animated profiles, a calculator with gradients, vorticity and
  time statistics, discharges through sections and volumes.
- Mesh: a step-by-step workflow panel; unstructured triangles with refinement, breaklines, islands and
  fixed points, stitched to structured channels and grids; channel templates (straight flume, bend,
  meanders); bottom from rasters, points or formulas; roughness zones; boundary conditions numbered as
  TELEMAC numbers them; editing of existing meshes; eight output formats.

No dependency beyond NumPy, which ships with QGIS. QGIS 3.44 LTR and QGIS 4.x, Windows, Linux, macOS.

User guide: https://github.com/alphavarez/telemac-toolkit/blob/main/docs/guide/README.md
Source code and issues: https://github.com/alphavarez/telemac-toolkit

Licence: GPL-3.0-or-later (see LICENSE). Author: Kevin Douglas Alvarez Segales.
