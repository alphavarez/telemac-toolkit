# 9. Processing tools reference

Every tool of the plugin is in the Processing toolbox under **TELEMAC Toolkit**, so it can run in batch mode, in the Model Designer or from Python (`processing.run("telemactoolkit:<id>", {...})`). This page is generated from the tools themselves (`scripts/make_algorithm_reference.py`).

## Results (group SELAFIN)

### Calculator (formulas, gradients, vorticity, time maxima)

Id: `telemactoolkit:calculator`

New variables computed from the variables of a SELAFIN file, written to a new file (the original variables can be kept). One output per line: NAME [UNIT] = expression, in NumPy syntax. Each output can be used by its name in the following lines. Two outputs named 'NAME X' and 'NAME Y' (or U and V) are shown as a vector with arrows.

Variables: TELEMAC symbols (U, V, H water depth, S free surface, B bottom, F Froude, Q unit discharge, M velocity magnitude, ...), names with underscores (WATER_DEPTH) or canonical ids (water_depth); x, y, t (time in s), area (area around each node), g = 9.81.

Functions:

- ddx(f), ddy(f): spatial derivatives (P1 gradient, averaged to the nodes)
- grad(f): magnitude of the gradient, e.g. grad(S) = water surface slope
- div(u, v): divergence du/dx + dv/dy
- curl(u, v) or vorticity(u, v): vertical vorticity dv/dx - du/dy (1/s)
- strain(u, v): strain rate magnitude sqrt(2ux² + 2vy² + (uy+vx)²)
- okubo_weiss(u, v): sn² + ss² - w², negative inside eddies
- laplacian(f): Laplacian, div(grad f)
- smooth(f): one pass of nodal smoothing
- mag(u, v): vector magnitude
- bearing(u, v): direction towards which the vector points, degrees from north, clockwise
- tmax(f), tmin(f), tmean(f), tstd(f): maximum, minimum, mean, standard deviation over all frames
- tsum(f), tint(f): sum over frames; time integral (trapezoidal, in f x s)
- tpeak(f): time (s) of the maximum of f
- tfirst(f), tlast(f): f at the first or the last frame
- ddt(f): time derivative (central differences), e.g. ddt(S) = rate of rise
- wet(eps): 1 where the water depth is above eps, else 0
- where(c, a, b), clip, minimum, maximum, abs, sqrt, exp, log, sin, cos, arctan2, ...: NumPy functions

Spatial derivatives use the linear finite-element interpolation of TELEMAC: the gradient is constant in each triangle and averaged to the nodes (exact for linear fields; second order on regular meshes, first order on irregular ones). Near dry areas mask the result, e.g. where(H > 0.01, curl(U, V), 0).

Examples:

- FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))
- SHEAR [PA] = 1000 * g * 0.03**2 * mag(U, V)**2 / maximum(H, 0.01)**(1/3)
- HMAX [M] = tmax(H); TIME HMAX [S] = tpeak(H)
- SLOPE X = ddx(S); SLOPE Y = ddy(S)
- RISE [M/S] = ddt(S)
- DEPTH CHANGE [M] = tlast(H) - tfirst(H)

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | file |  |
| Outputs (NAME [UNIT] = expression, one per line) | `EXPRESSIONS` | text | (example text) |
| Keep the variables of the file | `KEEP` | yes/no | yes |
| Frames (empty = all; e.g. -1, 0:10, 0:100:5) *(optional)* | `FRAMES` | text |  |
| Depth below which wet() is 0 (m) *(advanced)* | `DRY` | number | 0.01 |
| CRS saved next to the output (.prj; empty = copy the input's) *(optional)* *(advanced)* | `MESH_CRS` | CRS |  |
| Output SELAFIN | `OUTPUT` | output file |  |
| Open the result when done | `LOAD` | yes/no | yes |

### Extract 2D from a 3D result

Id: `telemactoolkit:extract2d`

Write a 2D SELAFIN file from a TELEMAC-3D result: the depth average (ELEVATION Z becomes FREE SURFACE and the lowest plane becomes BOTTOM), the bottom plane, the surface plane or any plane by index (0 = bottom). QGIS can display the 2D file; it cannot read 3D SELAFIN files.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| 3D SELAFIN file | `INPUT` | file |  |
| Extract | `MODE` | choice: Depth average; Bottom plane; Surface plane; Plane by index | Depth average |
| Plane index (0 = bottom) | `PLANE` | number | 0 |
| 2D SELAFIN file | `OUTPUT` | output file |  |

### Extract frames and variables

Id: `telemactoolkit:subset`

Copy selected frames and variables to a new SELAFIN file. Frames: comma-separated indices (0 = first, -1 = last) or ranges such as 0:10 or 0:100:5. Empty = all. Variables: comma-separated names. Empty = all.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | file |  |
| Frames *(optional)* | `FRAMES` | text | -1 |
| Variables *(optional)* | `VARIABLES` | text |  |
| Output SELAFIN | `OUTPUT` | output file |  |

### Integrals and fluxes (discharge, volume, line and area integrals)

Id: `telemactoolkit:integrals`

Integrals of a SELAFIN result for every frame, as a table:

- Discharge through section lines: the integral of h u.n along each line, with depth and velocity interpolated as in the model (the line is cut at every mesh edge and the product integrated exactly). Positive from the right of the line to its left: draw a section from the left bank to the right bank (looking downstream) and the flow downstream is positive.
- Flux of an expression: the same with a weight, e.g. a tracer concentration for its mass flux.
- Line integral and mean of an expression along each line (e.g. mean velocity across a section).
- Area integral and mean over each polygon (elements whose centre is inside) or the whole mesh: H gives the water volume, 1 the area, where(H > 0.01, 1, 0) the wet area.

Expressions follow the Calculator syntax. Lines and polygons are transformed into the mesh CRS (a .prj next to the file or the CRS given here).

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN result | `INPUT` | file |  |
| Integral | `MODE` | choice: Discharge through section lines (m³/s); Flux of an expression through section lines (expression × unit discharge, e.g. a concentration); Line integral and mean of an expression along lines; Area integral and mean of an expression over polygons (or the whole mesh); H gives the volume | Discharge through section lines (m³/s) |
| Expression (flux weight, line or area integrand) *(optional)* | `EXPRESSION` | text | H |
| Section lines *(optional)* | `LINES` | vector layer |  |
| Polygons (empty = whole mesh) *(optional)* | `POLYGONS` | vector layer |  |
| Name field *(optional)* | `NAME_FIELD` | field | name |
| Frames (empty = all) *(optional)* | `FRAMES` | text |  |
| CRS of the mesh (empty = .prj next to the file, else the layers' CRS) *(optional)* | `MESH_CRS` | CRS |  |
| Integrals | `OUTPUT` | output layer |  |

### Mesh boundary as lines

Id: `telemactoolkit:meshboundary`

Extract the boundary of a SELAFIN mesh as closed lines: the outer boundary first, then each island, with node counts and boundary numbering ranges (IPOBO).

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | file |  |
| CRS of the mesh coordinates *(optional)* | `CRS` | CRS |  |
| Boundary | `OUTPUT` | output layer |  |

### SELAFIN file information

Id: `telemactoolkit:selafininfo`

Report the content of a TELEMAC SELAFIN file: dimensions, variables, frames, start date, precision and mesh quality. Warns when float32 coordinates lose precision.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | file |  |
| Report | `OUTPUT` | output file |  |

### Time series at points

Id: `telemactoolkit:probetimeseries`

Interpolate variables of a SELAFIN result at every point of a layer, for every frame. Output is a table with one row per point and frame (long format), ready for plotting or export to CSV. Points outside the mesh are skipped. Point coordinates must be in the same CRS as the mesh, or give the mesh CRS so the points are transformed (SELAFIN files carry no CRS).

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | file |  |
| Points | `POINTS` | vector layer |  |
| Variables (empty = all) *(optional)* | `VARIABLES` | text |  |
| CRS of the mesh coordinates (empty = same as the points) *(optional)* | `MESH_CRS` | CRS |  |
| Time series | `OUTPUT` | output layer |  |

## Mesh (group Mesh)

### Boundary conditions from lines (.cli)

Id: `telemactoolkit:boundaryconditions`

Write the TELEMAC boundary conditions file (.cli). Every boundary node is a wall (2 2 2) except the nodes near a boundary line, which take its type: wall, Q (4 5 5), H (5 4 4), QH (5 5 5), UV (4 6 6), UVH (5 6 6), free (4 4 4), waves (1 1 1) or a three-digit code. Nodes are listed in TELEMAC order (outer boundary counter-clockwise, then islands), so the liquid boundary numbers in the log are the ones the steering file and the .liq file use. Optionally writes a .liq template with one column per prescribed value.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh (SELAFIN) | `INPUT` | file |  |
| CRS of the mesh (empty = .prj next to the file, else project CRS) *(optional)* | `MESH_CRS` | CRS |  |
| Boundary lines | `LINES` | vector layer |  |
| Type field | `TYPE_FIELD` | field | type |
| Name field *(optional)* | `NAME_FIELD` | field | name |
| Value field (Q or H for the .liq template) *(optional)* | `VALUE_FIELD` | field | value |
| Snapping distance (0 = half the boundary edge length) | `TOLERANCE` | number | 0.0 |
| Boundary conditions file | `OUTPUT` | output file |  |
| Liquid boundaries file template *(optional)* | `OUTPUT_LIQ` | output file |  |

### Channel or flume template (straight, bend, meanders)

Id: `telemactoolkit:channeltemplate`

Draw a laboratory flume or a schematic channel in a mesh project, with a bed you define by numbers instead of a DEM. Uses: flume experiments and their numerical replicas, verification against analytical solutions (normal depth, backwater, bend flow), sensitivity tests of mesh size, friction or turbulence, and idealised reaches for teaching or sediment studies. Natural rivers are drawn on the map instead (domain, channels and a DEM).

Planform:

- Straight: a rectangle of the given length.
- Bend: straight inlet, circular bend (radius on the axis, angle such as 90 or 180 degrees, left or right), straight outlet: the classic curved-flume experiments.
- Meanders: a sine-generated curve (Langbein & Leopold), the shape of natural meanders: heading theta(s) = theta0 sin(2 pi s / wavelength), with the wavelength measured along the axis.

The channel starts at the origin ((0, 0) by default) heading along the direction angle. It adds the domain, an inlet line (Q) and an outlet line (H) and, with Structured, a structured strip whose cross-sections are normal to the axis, also in bends. The bed is saved as a formula in the project: of x and y for a straight channel, of the channel coordinates s (along the axis) and n (across it, left positive) otherwise, with a bed level at the inlet, a slope along the axis and a rectangular or trapezoidal section. Generate mesh applies it. Leave the project empty to create one in local coordinates.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh project (empty = create one) *(optional)* | `PROJECT` | file |  |
| New mesh project (when none is given) *(optional)* | `NEW_PROJECT` | output file |  |
| Planform | `PLANFORM` | choice: Straight; Bend (straight inlet, circular bend, straight outlet); Meanders (sine-generated) | Straight |
| Width (m) | `WIDTH` | number | 1.0 |
| Straight: length (m) | `LENGTH` | number | 20.0 |
| Bend and meanders: straight inlet length (m) | `INLET` | number | 5.0 |
| Bend and meanders: straight outlet length (m) | `OUTLET` | number | 5.0 |
| Bend: radius on the axis (m) | `RADIUS` | number | 2.0 |
| Bend: angle (degrees) | `BEND_ANGLE` | number | 180.0 |
| Bend: turning | `TURN` | choice: Left (counter-clockwise); Right (clockwise) | Left (counter-clockwise) |
| Meanders: wavelength along the axis (m) | `WAVELENGTH` | number | 10.0 |
| Meanders: largest deflection angle theta0 (degrees) | `DEFLECTION` | number | 60.0 |
| Meanders: number of wavelengths | `COUNT` | number | 2.0 |
| Origin (inlet, centre of the inlet for curved channels) *(optional)* | `ORIGIN` | point | 0,0 |
| Direction at the inlet (degrees, counter-clockwise from x) | `ANGLE` | number | 0.0 |
| Element size (m) | `CELL` | number | 0.05 |
| Structured mesh over the channel | `STRUCTURED` | yes/no | yes |
| Bed level at the inlet (m) | `BED_Z0` | number | 0.0 |
| Bed slope along the axis (m/m, positive downstream) | `SLOPE` | number | 0.0 |
| Cross-section | `SECTION` | choice: Rectangular (flat bed); Trapezoidal | Rectangular (flat bed) |
| Trapezoid: bottom width (m) | `BOTTOM_WIDTH` | number | 0.5 |
| Trapezoid: side slope (horizontal per 1 vertical) | `SIDE_SLOPE` | number | 1.0 |
| Trapezoid: bank height (m) | `BANK_HEIGHT` | number | 0.2 |
| Inflow (m³/s, optional) *(optional)* | `Q` | number |  |
| Outlet water level (m, optional) *(optional)* | `H` | number |  |

### Edit node values of a mesh (bottom, friction, …)

Id: `telemactoolkit:interpolatebottom`

Change a node variable of an existing SELAFIN geometry (BOTTOM, BOTTOM FRICTION or any other) without touching the mesh: same nodes, elements, boundary numbering, origin and date, so the boundary conditions file (.cli) still fits. Other variables are kept.

New values come, in priority order, from rasters (bilinear), survey points (linear on their triangulation), a formula, or a constant. The formula may use x and y, or the channel coordinates s and n along an axis line (e.g. carve a trapezoidal channel: 'z0 - 0.001*s + clip((abs(n) - 5)/2, 0, 3)').

Where: the whole mesh, or only the nodes inside the polygons, with an optional transition width along the polygon edge (the change fades from full inside to nothing at the edge).

How: replace, add (raise or lower by the value), raise to (keep the higher value: fills, dikes), lower to (keep the lower: dredging, cuts).

Typical edits: a new survey over a reach (rasters, polygon), a dredged channel (formula in s and n, lower to), a dike crest (constant, raise to, thin polygon), a friction zone (variable BOTTOM FRICTION, constant, replace).

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh (SELAFIN) | `INPUT` | file |  |
| CRS of the mesh (empty = .prj next to the file, else project CRS) *(optional)* | `MESH_CRS` | CRS |  |
| Variable | `VARIABLE` | text | BOTTOM |
| New values: rasters (priority order) *(optional)* | `RASTERS` | layers |  |
| New values: survey points *(optional)* | `POINTS` | vector layer |  |
| Survey value field (empty = point Z) *(optional)* | `POINT_FIELD` | field |  |
| New values: formula of x, y or s, n (lowest priority) *(optional)* | `FORMULA` | text |  |
| Channel axis line for s and n in the formula *(optional)* | `AXIS` | vector layer |  |
| New values: constant (when no other source) *(optional)* | `CONSTANT` | number |  |
| Operation | `OPERATION` | choice: Replace with the new values; Add the new values (raise or lower by); Raise to the new values (keep the higher: fills, dikes, embankments); Lower to the new values (keep the lower: dredging, cuts, channels) | Replace with the new values |
| Only inside these polygons (empty = whole mesh) *(optional)* | `AREA` | vector layer |  |
| Transition width inside the polygon edge (m, 0 = sharp) | `BLEND` | number | 0.0 |
| Unit of a new variable *(optional)* *(advanced)* | `UNIT` | text | M |
| Unit factor of rasters and points *(advanced)* | `Z_FACTOR` | number | 1.0 |
| Fill nodes no source covers from the nearest node *(advanced)* | `FILL` | yes/no | yes |
| Output geometry file | `OUTPUT` | output file |  |
| Open the mesh when done | `LOAD` | yes/no | yes |

### Export or convert mesh

Id: `telemactoolkit:exportmesh`

Write an existing mesh (SELAFIN, 2DM, UGRID or any mesh QGIS reads) for another program: TELEMAC, SMS 2DM (BASEMENT, TUFLOW FV), Gmsh, GiD (Iber and other GiD-based programs), UGRID (Delft3D FM), MIKE 21 FM or ADCIRC/SCHISM. Optional boundary lines and roughness zones add boundary strings, material ids and friction. The file can be written in local coordinates (a point moved to (0, 0), optional rotation). HEC-RAS builds its own cells, so it is exported from Generate mesh, not here.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh | `INPUT` | file |  |
| CRS of the mesh (empty = from the file or the project) *(optional)* | `MESH_CRS` | CRS |  |
| Boundary condition lines *(optional)* | `BOUNDARIES` | vector layer |  |
| Roughness zones *(optional)* | `ROUGHNESS_ZONES` | vector layer |  |
| Friction outside the roughness zones *(optional)* *(advanced)* | `FRICTION` | number |  |
| Output format | `FORMAT` | choice: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | SMS 2DM – BASEMENT, TUFLOW FV, SMS |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | yes/no | no |
| Coordinates in the file *(advanced)* | `FRAME` | choice: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(optional)* *(advanced)* | `FRAME_POINT` | point |  |
| Rotation of the local axes (degrees, counter-clockwise) *(advanced)* | `FRAME_ANGLE` | number | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(advanced)* | `STORE_OFFSET` | yes/no | yes |
| Mesh file | `OUTPUT` | output file |  |
| TELEMAC liquid boundaries template (.liq) *(optional)* | `OUTPUT_LIQ` | output file |  |
| Open the result when done | `LOAD` | yes/no | yes |

### Extract part of a mesh

Id: `telemactoolkit:cutmesh`

Keep the elements whose centre lies inside the polygons (or outside them) and write the smaller mesh, with its node values, in any output format. Useful to cut a study reach out of a regional mesh, or to remove an area before meshing it again.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh | `INPUT` | file |  |
| CRS of the mesh (empty = from the file or the project) *(optional)* | `MESH_CRS` | CRS |  |
| Polygons | `POLYGONS` | vector layer |  |
| Keep what is outside the polygons | `OUTSIDE` | yes/no | no |
| Output format | `FORMAT` | choice: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | yes/no | no |
| Coordinates in the file *(advanced)* | `FRAME` | choice: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(optional)* *(advanced)* | `FRAME_POINT` | point |  |
| Rotation of the local axes (degrees, counter-clockwise) *(advanced)* | `FRAME_ANGLE` | number | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(advanced)* | `STORE_OFFSET` | yes/no | yes |
| Mesh file | `OUTPUT` | output file |  |
| Open the result when done | `LOAD` | yes/no | yes |

### Generate mesh

Id: `telemactoolkit:generatemesh`

Build a mesh from a mesh project (or from individual layers, which take precedence) and write it for TELEMAC or another program.

Unstructured triangles fill the domain with an element size that grows away from refinement areas, points and lines at the growth rate; the mesh follows breaklines vertex by vertex and guide lines at the local size, keeps fixed points, and stitches in structured channel strips and regular grid zones node for node.

Bottom: rasters in priority order, then survey points, then a formula of x and y (for flumes); nodes nothing covers take the nearest value. Roughness zones give a material id per element and a friction value per node (BOTTOM FRICTION for TELEMAC). Boundary lines give the TELEMAC .cli (and a .liq template), 2DM nodestrings, Gmsh physical lines or ADCIRC/MIKE boundary lists.

Advanced: write the file in local coordinates (lower-left corner or a chosen point at (0, 0), optional rotation); TELEMAC files can keep the shift in their header.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh project *(optional)* | `PROJECT` | file |  |
| Domain (polygons; holes are islands) *(optional)* | `DOMAIN` | vector layer |  |
| Islands *(optional)* | `ISLANDS` | vector layer |  |
| Refinement areas (size) *(optional)* | `REFINEMENT_AREAS` | vector layer |  |
| Refinement points (size) *(optional)* | `REFINEMENT_POINTS` | vector layer |  |
| Constraint lines (mode: breakline or guide) *(optional)* | `CONSTRAINT_LINES` | vector layer |  |
| Fixed points *(optional)* | `FIXED_POINTS` | vector layer |  |
| Structured channels (role: left, right, axis) *(optional)* | `CHANNELS` | vector layer |  |
| Structured grid zones (dx, dy) *(optional)* | `GRID_ZONES` | vector layer |  |
| Roughness zones (material, friction) *(optional)* | `ROUGHNESS_ZONES` | vector layer |  |
| Boundary condition lines (type, name, value) *(optional)* | `BOUNDARIES` | vector layer |  |
| Default element size (0 = from the project) | `DEFAULT_SIZE` | number | 0.0 |
| Size growth rate | `GROWTH` | number | 1.2 |
| Smoothing iterations *(advanced)* | `SMOOTHING` | number | 8 |
| Bottom: rasters (priority order) *(optional)* | `RASTERS` | layers |  |
| Bottom: survey points *(optional)* | `POINTS` | vector layer |  |
| Survey value field (empty = point Z) *(optional)* | `POINT_FIELD` | field |  |
| Bottom: formula of x, y or channel coordinates s, n (empty = from the project) *(optional)* | `FORMULA` | text |  |
| Elevation unit factor *(advanced)* | `Z_FACTOR` | number | 1.0 |
| Friction outside the roughness zones *(optional)* *(advanced)* | `FRICTION` | number |  |
| Output format | `FORMAT` | choice: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | yes/no | no |
| Coordinates in the file *(advanced)* | `FRAME` | choice: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(optional)* *(advanced)* | `FRAME_POINT` | point |  |
| Rotation of the local axes (degrees, counter-clockwise) *(advanced)* | `FRAME_ANGLE` | number | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(advanced)* | `STORE_OFFSET` | yes/no | yes |
| Mesh file | `OUTPUT` | output file |  |
| TELEMAC liquid boundaries template (.liq) *(optional)* | `OUTPUT_LIQ` | output file |  |
| Open the result when done | `LOAD` | yes/no | yes |

### Mesh quality

Id: `telemactoolkit:meshquality`

Element quality of a mesh (SELAFIN, 2DM, UGRID, ...) as a polygon layer: area, smallest and largest angle, longest edge, and a flag for elements below the minimum angle, above the maximum angle or with all their nodes on the boundary (over-constrained elements, which block the flow in a corner). The log gives the overall statistics.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Mesh | `INPUT` | file |  |
| CRS of the mesh (empty = from the file or the project) *(optional)* | `MESH_CRS` | CRS |  |
| Flag angles below (degrees) | `MIN_ANGLE` | number | 20.0 |
| Flag angles above (degrees) | `MAX_ANGLE` | number | 130.0 |
| Only flagged elements | `ONLY_FLAGGED` | yes/no | no |
| Element quality | `OUTPUT` | output layer |  |

### New mesh project

Id: `telemactoolkit:createmeshproject`

Create a GeoPackage with one empty, styled layer per kind of object and add it to the map:

- domain – polygons to mesh; holes are islands; an optional size per polygon
- islands – polygons left out of the mesh
- refinement_areas, refinement_points – target element size
- constraint_lines – lines the mesh follows: mode 'breakline' keeps every vertex (dikes, roads, banks); mode 'guide' is resampled to the local size
- fixed_points – nodes that must exist (gauges, structures)
- channels – left bank, right bank and an optional axis (role field) of a channel, meshed as a structured strip with nodes_across and spacing
- grid_zones – rectangles meshed as regular grids (dx, dy)
- roughness_zones – material id and friction value
- boundaries – boundary lines with a type (wall, Q, H, QH, UV, UVH, free, waves or a code such as 455), a name and an optional value

Choose local coordinates for a laboratory flume or a schematic channel: the project then uses a metric grid with no geographic reference, so you can draw from (0, 0).

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| CRS (projected, in metres) | `CRS` | CRS | ProjectCrs |
| Local coordinates (laboratory or schematic channel, no map CRS) | `LOCAL` | yes/no | no |
| Default element size (m) *(optional)* | `DEFAULT_SIZE` | number | 0.0 |
| Mesh project | `OUTPUT` | output file |  |

### Regular grid mesh

Id: `telemactoolkit:regulargrid`

A structured mesh of nx by ny rectangular cells, starting at the origin (lower-left corner, (0, 0) by default) and rotated counter-clockwise by the angle. Cells are split into two triangles with alternating diagonals (no directional bias), or kept as quadrilaterals for formats that accept them. For grids inside a larger mesh, draw rectangles in grid_zones and use Generate mesh.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Origin (lower-left corner) | `ORIGIN` | point | 0,0 |
| Cell size along x | `DX` | number | 1.0 |
| Cell size along y | `DY` | number | 1.0 |
| Cells along x | `NX` | number | 100 |
| Cells along y | `NY` | number | 20 |
| Rotation (degrees, counter-clockwise) | `ANGLE` | number | 0.0 |
| Diagonals | `PATTERN` | choice: alternate; right; left | alternate |
| CRS | `CRS` | CRS | ProjectCrs |
| Bottom: formula of x and y (optional) *(optional)* | `FORMULA` | text |  |
| Output format | `FORMAT` | choice: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | yes/no | no |
| Coordinates in the file *(advanced)* | `FRAME` | choice: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(optional)* *(advanced)* | `FRAME_POINT` | point |  |
| Rotation of the local axes (degrees, counter-clockwise) *(advanced)* | `FRAME_ANGLE` | number | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(advanced)* | `STORE_OFFSET` | yes/no | yes |
| Mesh file | `OUTPUT` | output file |  |
| Open the mesh when done | `LOAD` | yes/no | yes |

### Structured channel mesh

Id: `telemactoolkit:channelmesh`

A structured mesh of a channel between a left and a right bank line, optionally split by an axis (thalweg) line. Lines are grouped by the 'channel' field and told apart by the 'role' field (left, right, axis). Nodes across and spacing along come from the nodes_across and spacing fields, or from the defaults below. Use it alone, or let Generate mesh stitch it into the floodplain mesh.

| Parameter | Key | Type | Default |
| --- | --- | --- | --- |
| Channel lines | `CHANNELS` | vector layer |  |
| Nodes across (default) | `CROSS_NODES` | number | 7 |
| Spacing along (default, m) | `INTERVAL` | number | 10.0 |
| Output format | `FORMAT` | choice: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | yes/no | no |
| Coordinates in the file *(advanced)* | `FRAME` | choice: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(optional)* *(advanced)* | `FRAME_POINT` | point |  |
| Rotation of the local axes (degrees, counter-clockwise) *(advanced)* | `FRAME_ANGLE` | number | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(advanced)* | `STORE_OFFSET` | yes/no | yes |
| Mesh file | `OUTPUT` | output file |  |
| Open the mesh when done | `LOAD` | yes/no | yes |

[← Formats, coordinates and boundary conditions](08-formats-coordinates-boundaries.md) · [Contents](README.md) · [Next: Scripting with Python →](10-python.md)
