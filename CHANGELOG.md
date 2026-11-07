# TELEMAC Toolkit — changelog

Newest first. Each version is tested in QGIS 3.44 LTR and QGIS 4; from 0.5.0 also with TELEMAC-2D v9.1 runs.
0.5.0 gathers the development builds 0.4.4 to 0.8.6, which were not published.

## 0.5.1 — 2026-11-06

- Security scan (Bandit) clean: errors the plugin works around are written to the TELEMAC Toolkit log instead of ignored; the MMG call runs only the mmg2d program you set (checked to be a program), with an argument list and no shell

## 0.5.0 — 2026-11-06

### Meshing

- Meshing methods: Balanced and Quality-guaranteed (smallest angle 20 to 30 degrees); with the optional gmsh package Gmsh Delaunay, Frontal-Delaunay, MeshAdapt and quadrilaterals; Compare meshing methods (also meshes made in Blue Kenue or elsewhere)
- Mesh type in step 4: Triangles, Rectangles where drawn (structured channels and grid zones) or Quadrilaterals everywhere (Gmsh)
- Sizes from the terrain (bed tolerance), element sizes as a raster, sizes suggested from the DEM's pixels; breaklines from a DEM (crests and toes)
- The mesh boundary follows the drawn outline within a set distance (Follow the outline within; 5% of the element size by default); digitizing zigzags no longer become nodes and slivers
- Node count estimated before meshing, with a node limit (2 million by default, adjustable); about 20% faster on large meshes
- Optimise or adapt a mesh with MMG (shapes, terrain sizes, or sizes from a TELEMAC result)
- Fixed: outlines digitized with many close vertices could hang Generate mesh; vertices of rings far from the origin could be lost; a domain edge in line with a channel end could loop forever; edge flips could worsen the worst angle or change structured parts

### Rivers and channels

- River reach helper: from a centreline (constant width or banks found in a DEM) or two bank lines; river networks with confluences and bifurcations; nodes across follow the width (at least 7); cross-sections turned only in tight bends
- Channel templates composed from segments: straight, bends, elbows, sine-generated and Kinoshita meanders, tapers and expansions, slopes and drops; cross-sections rectangular, trapezoidal, compound, V, parabolic, asymmetric or any formula
- Channels from a drawn centreline, or from TXT, CSV or DXF files (as a centreline with widths, or as outlines)
- Structured cells cut along one diagonal (the shorter one in channels), which stops TELEMAC's depths oscillating from node to node

### Drawing and editing

- Check the drawing measures the geometry against the element size (near-duplicate vertices, spikes, close corners, dense breaklines, narrow gaps, line ends just off the outline, invalid polygons), shows the problems on the map and fixes the safe ones; also as a Processing tool
- Hover tips on every project layer (size, area, perimeter, nodes, boundary codes), attribute forms with drop-down lists, snapping, tracing and Advanced Digitizing set when a project opens; object counts follow the drawing
- Edit a generated mesh in QGIS and rebuild it with every variable and its .cli; re-mesh inside polygons; boundary lines from a .cli; Blue Kenue lines (.i2s, .i3s) and meshes (.t3s) read
- Blue Kenue's boundary types and tracer codes (LITBOR) in the boundary forms; custom codes

### Run files and TELEMAC

- Run files (step 6): rating curve, liquid boundaries file (.liq) from CSV or Excel, initial conditions (also depth above the channel axes), steering file with a time step that divides the duration, NERD advection and velocity profiles on discharge boundaries; warnings for dry discharge boundaries
- SELAFIN coordinates kept exact in map coordinates (integer origin above 4 km, 0.01 mm); double precision (SERAFIND) as an option, declared in the steering file
- The whole chain tested in TELEMAC-2D v9.1 (meandering river, river network, Gmsh mesh, re-meshed and edited meshes, laboratory expansion)

### Fixes and other changes

- Generate mesh overwrites a mesh loaded in QGIS on Windows (before: OSError Errno 22); the local laboratory grid is a named QGIS CRS; a CRS without an EPSG code is matched to its EPSG equivalent

## 0.4.3 — 2026-10-22

- Arrows follow the vector shown on the map (velocity, unit discharge, wind); an Arrows list in the Results panel
- Unit discharge (FLOWRATE ALONG X/Y, DEBIT SUIVANT X/Y) and wind are recognised as vectors
- Calculator: NAME X / NAME Y outputs are written as NAME ALONG X/Y, so QGIS draws them with arrows

## 0.4.2 — 2026-10-21

- New plugin icon, legible in the plugin manager and on the plugin site
- Author name spelling
- Calculator from the command line (qgis_process): no error after the result is written

## 0.4.1 — 2026-10-21

- Code quality: bottom formulas evaluated by a whitelisted expression walker (no eval)
- Handled errors are logged to the TELEMAC Toolkit tab of the Log Messages panel instead of ignored
- Flake8 clean-up (unused imports, ambiguous names)

## 0.4.0 — 2026-10-21

- Arrow settings: grid spacing, proportional or fixed-scale length, hide below, colour by magnitude
- Results open in view: CRS from a .prj or the project, or asked (local or map CRS); zoom and CRS buttons
- Calculator: algebra, ddx, ddy, grad, div, curl, strain, Laplacian, tmax, tmin, tmean, tint, tpeak, ddt
- Integrals and fluxes: discharge through sections, line and area integrals, volumes
- Edit node values of a mesh: replace, add, raise, lower, inside polygons; .cli kept
- Channel templates: straight flume, circular bend, sine-generated meanders; bed formulas in s and n
- Liquid boundaries numbered exactly as TELEMAC numbers them
- Automatic default element size; legend without duplicates

## 0.3.0

- Mesh workflow panel; clearer project layers (0.2.0 projects still open)
- Local coordinates and straight flume template; files shifted to (0, 0) with optional rotation
- Roughness zones; output to SELAFIN, 2DM, Gmsh, GiD, UGRID, MIKE 21, ADCIRC and HEC-RAS layers
- Export or convert mesh; extract part of a mesh; find a node or element by number

## 0.2.0

- Mesher: unstructured triangles with refinement, constraint lines, fixed points and islands
- Structured channels and regular grids; interpolation from rasters and points
- Boundary conditions from lines; mesh quality

## 0.1.0

- SELAFIN read/write library; results as styled mesh layers, 2D and 3D
- Results panel with the Temporal Controller; live cursor, time series, animated profiles
- Processing tools: file information, 2D from 3D, subset, time series at points, mesh boundary
