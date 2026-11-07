# TELEMAC Toolkit for QGIS

Build meshes and analyse results for TELEMAC inside QGIS, with nothing but QGIS and Python; the same
mesh can be written for BASEMENT, TUFLOW FV, Iber, Delft3D FM, MIKE 21 FM, ADCIRC/SCHISM and HEC-RAS.

- Results: SELAFIN results (2D and 3D) as styled mesh layers with the QGIS Temporal Controller, arrow
  controls, live cursor, probe time series, animated profiles, a calculator with gradients, vorticity and
  time statistics, discharges through sections and volumes.
- Mesh: a step-by-step workflow panel; unstructured triangles with refinement, breaklines, islands and
  fixed points, stitched to structured channels and grids; a river reach helper (centreline, bank
  lines or banks found in a DEM; river networks with confluences; nodes across that follow the width); channel templates (straight flume, bend, meanders); a
  quality-guaranteed method (smallest angle up to 30°), optional Gmsh methods and a tool that compares meshing methods (and Blue Kenue meshes); terrain-driven sizes, breaklines from a DEM and MMG optimisation or adaptation; editing of generated meshes (QGIS Mesh Digitizing, re-meshing inside polygons) with their values and boundary conditions kept; bottom from rasters, points or formulas; roughness zones; boundary conditions numbered as
  TELEMAC numbers them; editing of existing meshes; eight output formats.

No dependency beyond NumPy, which ships with QGIS. QGIS 3.44 LTR and QGIS 4.x, Windows, Linux, macOS.

User guide: https://github.com/alphavarez/telemac-toolkit/blob/main/docs/guide/README.md
Guía de usuario en español: https://github.com/alphavarez/telemac-toolkit/blob/main/docs/guia/README.md
Source code and issues: https://github.com/alphavarez/telemac-toolkit

Licence: GPL-3.0-or-later (see LICENSE). Author: Kevin Douglas Alvarez Segales.
