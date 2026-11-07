# 9. Referencia de las herramientas de Procesos

Todas las herramientas del plugin están en la caja de herramientas de Procesos, bajo **TELEMAC Toolkit**, así que se pueden ejecutar por lotes, en el Diseñador de modelos o desde Python (`processing.run("telemactoolkit:<id>", {...})`). Esta página se genera a partir de las propias herramientas (`scripts/make_algorithm_reference.py`). **Los nombres de las herramientas, sus textos de ayuda y los nombres de los parámetros están en inglés, tal como los muestra QGIS**, para que los encuentres igual en pantalla; los títulos, las columnas y los tipos están en español. La explicación de cada función está en las páginas de esta guía.

## Resultados (grupo SELAFIN)

### Calculator (formulas, gradients, vorticity, time maxima)

Id: `telemactoolkit:calculator`

New variables computed from the variables of a SELAFIN file, written to a new file (the original variables can be kept). One output per line: NAME [UNIT] = expression, in NumPy syntax. Each output can be used by its name in the following lines. Two outputs named 'NAME X' and 'NAME Y' (or U and V) are shown as a vector with arrows (X and Y are written to the file as 'NAME ALONG X/Y', the names QGIS reads as a vector).

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

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | archivo |  |
| Outputs (NAME [UNIT] = expression, one per line) | `EXPRESSIONS` | texto | (texto de ejemplo) |
| Keep the variables of the file | `KEEP` | sí/no | sí |
| Frames (empty = all; e.g. -1, 0:10, 0:100:5) *(opcional)* | `FRAMES` | texto |  |
| Depth below which wet() is 0 (m) *(avanzado)* | `DRY` | número | 0.01 |
| CRS saved next to the output (.prj; empty = copy the input's) *(opcional)* *(avanzado)* | `MESH_CRS` | SRC |  |
| Output SELAFIN | `OUTPUT` | archivo de salida |  |
| Open the result when done | `LOAD` | sí/no | sí |

### Extract 2D from a 3D result

Id: `telemactoolkit:extract2d`

Write a 2D SELAFIN file from a TELEMAC-3D result: the depth average (ELEVATION Z becomes FREE SURFACE and the lowest plane becomes BOTTOM), the bottom plane, the surface plane or any plane by index (0 = bottom). QGIS can display the 2D file; it cannot read 3D SELAFIN files.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| 3D SELAFIN file | `INPUT` | archivo |  |
| Extract | `MODE` | opción: Depth average; Bottom plane; Surface plane; Plane by index | Depth average |
| Plane index (0 = bottom) | `PLANE` | número | 0 |
| 2D SELAFIN file | `OUTPUT` | archivo de salida |  |

### Extract frames and variables

Id: `telemactoolkit:subset`

Copy selected frames and variables to a new SELAFIN file. Frames: comma-separated indices (0 = first, -1 = last) or ranges such as 0:10 or 0:100:5. Empty = all. Variables: comma-separated names. Empty = all.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | archivo |  |
| Frames *(opcional)* | `FRAMES` | texto | -1 |
| Variables *(opcional)* | `VARIABLES` | texto |  |
| Output SELAFIN | `OUTPUT` | archivo de salida |  |

### Integrals and fluxes (discharge, volume, line and area integrals)

Id: `telemactoolkit:integrals`

Integrals of a SELAFIN result for every frame, as a table:

- Discharge through section lines: the integral of h u.n along each line, with depth and velocity interpolated as in the model (the line is cut at every mesh edge and the product integrated exactly). Positive from the right of the line to its left: draw a section from the left bank to the right bank (looking downstream) and the flow downstream is positive.
- Flux of an expression: the same with a weight, e.g. a tracer concentration for its mass flux.
- Line integral and mean of an expression along each line (e.g. mean velocity across a section).
- Area integral and mean over each polygon (elements whose centre is inside) or the whole mesh: H gives the water volume, 1 the area, where(H > 0.01, 1, 0) the wet area.

Expressions follow the Calculator syntax. Lines and polygons are transformed into the mesh CRS (a .prj next to the file or the CRS given here).

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN result | `INPUT` | archivo |  |
| Integral | `MODE` | opción: Discharge through section lines (m³/s); Flux of an expression through section lines (expression × unit discharge, e.g. a concentration); Line integral and mean of an expression along lines; Area integral and mean of an expression over polygons (or the whole mesh); H gives the volume | Discharge through section lines (m³/s) |
| Expression (flux weight, line or area integrand) *(opcional)* | `EXPRESSION` | texto | H |
| Section lines *(opcional)* | `LINES` | capa vectorial |  |
| Polygons (empty = whole mesh) *(opcional)* | `POLYGONS` | capa vectorial |  |
| Name field *(opcional)* | `NAME_FIELD` | campo | name |
| Frames (empty = all) *(opcional)* | `FRAMES` | texto |  |
| CRS of the mesh (empty = .prj next to the file, else the layers' CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Integrals | `OUTPUT` | capa de salida |  |

### Mesh boundary as lines

Id: `telemactoolkit:meshboundary`

Extract the boundary of a SELAFIN mesh as closed lines: the outer boundary first, then each island, with node counts and boundary numbering ranges (IPOBO).

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | archivo |  |
| CRS of the mesh coordinates *(opcional)* | `CRS` | SRC |  |
| Boundary | `OUTPUT` | capa de salida |  |

### SELAFIN file information

Id: `telemactoolkit:selafininfo`

Report the content of a TELEMAC SELAFIN file: dimensions, variables, frames, start date, precision and mesh quality. Warns when float32 coordinates lose precision.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | archivo |  |
| Report | `OUTPUT` | archivo de salida |  |

### Time series at points

Id: `telemactoolkit:probetimeseries`

Interpolate variables of a SELAFIN result at every point of a layer, for every frame. Output is a table with one row per point and frame (long format), ready for plotting or export to CSV. Points outside the mesh are skipped. Point coordinates must be in the same CRS as the mesh, or give the mesh CRS so the points are transformed (SELAFIN files carry no CRS).

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| SELAFIN file | `INPUT` | archivo |  |
| Points | `POINTS` | capa vectorial |  |
| Variables (empty = all) *(opcional)* | `VARIABLES` | texto |  |
| CRS of the mesh coordinates (empty = same as the points) *(opcional)* | `MESH_CRS` | SRC |  |
| Time series | `OUTPUT` | capa de salida |  |

## Malla (grupo Mesh)

### Boundary conditions from lines (.cli)

Id: `telemactoolkit:boundaryconditions`

Write the TELEMAC boundary conditions file (.cli). Every boundary node is a wall (2 2 2) except the nodes near a boundary line, which take its type: wall, Q (4 5 5), H (5 4 4), QH (5 5 5), UV (4 6 6), UVH (5 6 6), free (4 4 4), waves (1 1 1) or a three-digit code. Nodes are listed in TELEMAC order (outer boundary counter-clockwise, then islands), so the liquid boundary numbers in the log are the ones the steering file and the .liq file use. Optionally writes a .liq template with one column per prescribed value.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN) | `INPUT` | archivo |  |
| CRS of the mesh (empty = .prj next to the file, else project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Boundary lines | `LINES` | capa vectorial |  |
| Type field | `TYPE_FIELD` | campo | type |
| Name field *(opcional)* | `NAME_FIELD` | campo | name |
| Value field (Q or H for the .liq template) *(opcional)* | `VALUE_FIELD` | campo | value |
| Snapping distance (0 = half the boundary edge length) | `TOLERANCE` | número | 0.0 |
| Boundary conditions file | `OUTPUT` | archivo de salida |  |
| Liquid boundaries file template *(opcional)* | `OUTPUT_LIQ` | archivo de salida |  |

### Boundary lines from a .cli

Id: `telemactoolkit:clilines`

One line per liquid boundary of a TELEMAC mesh, with the name, type and value fields of a mesh project's boundary layer: the boundary conditions of an existing mesh as lines you can edit on the map (move their ends along the boundary, change the type, add a boundary), then write back with Boundary conditions from lines. Snap to the mesh boundary (Mesh boundary as lines) while editing. Values come from the first row of the liquid boundaries file when given.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| TELEMAC mesh (SELAFIN) | `MESH` | archivo |  |
| Boundary conditions (empty = .cli next to the mesh) *(opcional)* | `CLI` | archivo |  |
| Liquid boundaries file (for the values) *(opcional)* | `LIQ` | archivo |  |
| CRS of the mesh (empty = project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Boundary lines | `OUTPUT` | capa de salida |  |

### Breaklines from a DEM

Id: `telemactoolkit:dembreaklines`

Lines along the sharp bends of the terrain, for the mesh to follow: crests (levee and bank tops, road embankments: the ground bends down on both sides) and toes (bank feet, channel edges: it bends up). A line is drawn where an element of the given size would miss the bend by more than the tolerance (|curvature| > 8 tolerance / size²). The cells beyond that are thinned to one-cell-wide lines, traced, simplified to the cell size, and lines shorter than the minimum length dropped. With a mesh project the lines are clipped to its domain and can be added to its constraint lines.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| DEM | `DEM` | raster |  |
| Mesh project (extent, and where to add the lines) *(opcional)* | `PROJECT` | archivo |  |
| Extent (empty = the project's domain or the DEM) *(opcional)* | `EXTENT` | extent |  |
| Bed tolerance (m) | `TOLERANCE` | número | 0.05 |
| Element size the lines are for (m, 0 = the project's default size) | `SIZE` | número | 0.0 |
| Lines | `KINDS` | opción: Crests and toes; Crests (levee and bank tops, embankments); Toes (bank feet, channel edges) | Crests and toes |
| Cell size (m, 0 = the DEM's pixel) *(avanzado)* | `CELL` | número | 0.0 |
| Shortest line kept (m, 0 = three element sizes) *(avanzado)* | `MIN_LENGTH` | número | 0.0 |
| Add the lines to the project's constraint lines | `ADD` | sí/no | no |
| Added as | `MODE` | opción: Guide lines (resampled); Breaklines (every vertex kept) | Guide lines (resampled) |
| Breaklines | `OUTPUT` | capa de salida |  |

### Channel or flume template (straight, bend, meanders, composed, drawn, imported)

Id: `telemactoolkit:channeltemplate`

Draw a laboratory flume or a schematic channel in a mesh project, with a bed you define by numbers instead of a DEM. Uses: flume experiments and their numerical replicas, verification against analytical solutions (normal depth, backwater, bend flow), sensitivity tests of mesh size, friction or turbulence, and idealised reaches for teaching or sediment studies. Natural rivers are drawn on the map instead (domain, channels and a DEM).

Planform:

- Straight: a rectangle of the given length.
- Bend: straight inlet, circular bend (radius on the axis, angle such as 90 or 180 degrees, left or right), straight outlet: the classic curved-flume experiments.
- Meanders: a sine-generated curve (Langbein & Leopold), the shape of natural meanders: heading theta(s) = theta0 sin(2 pi s / wavelength), with the wavelength measured along the axis.

The channel starts at the origin ((0, 0) by default) heading along the direction angle. It adds the domain, an inlet line (Q) and an outlet line (H) and, with Structured, a structured strip whose cross-sections are normal to the axis, also in bends. The bed is saved as a formula in the project: of x and y for a straight channel, of the channel coordinates s (along the axis) and n (across it, left positive) otherwise, with a bed level at the inlet, a slope along the axis and a rectangular or trapezoidal section. Generate mesh applies it. Leave the project empty to create one in local coordinates.

- Composed from segments: any channel built from a text, one segment per line: width W, straight L, bend left|right R ANGLE, elbow left|right ANGLE (a sharp corner), meanders WAVELENGTH THETA0 COUNT, kinoshita WAVELENGTH THETA0 COUNT JS JF, taper WIDTH LENGTH (gradual width change), expand WIDTH (sudden expansion or contraction), slope S, drop D (a step in the bed), bed Z, and section KIND ... (the cross-section from there on: rect; trapezoid BOTTOM_WIDTH SIDE_SLOPE BANK_HEIGHT [OFFSET]; compound MAIN_WIDTH SIDE_SLOPE FLOODPLAIN_HEIGHT FLOODPLAIN_SLOPE [OFFSET]; v DEPTH [OFFSET]; parabolic DEPTH [OFFSET]; asym BOTTOM_WIDTH LEFT_SLOPE RIGHT_SLOPE BANK_HEIGHT [OFFSET]; custom FORMULA of n and b). Pieces between sudden width changes and sharp corners get their own structured strips; the corners and steps get triangles.
- Centreline from a line layer: a line you draw (in local coordinates for a flume), with a constant width or a width field; corners sharper than the corner angle become elbows.
- Imported file: TXT or CSV (x y [width] per line; a line without numbers starts a new polyline) or DXF, read as a centreline or as an outline (a basin or any shape: the largest closed line is the domain).

Cross-sections (n across the axis, left positive; b half the local width): rectangular; trapezoidal; compound (a main channel and floodplains at a height, optionally rising towards the walls); V-shaped and parabolic (depth at the axis); asymmetric trapezoid (different side slopes); or any formula of n and b. The offset moves the main channel to the left (+) or right (-) of the axis.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh project (empty = create one) *(opcional)* | `PROJECT` | archivo |  |
| New mesh project (when none is given) *(opcional)* | `NEW_PROJECT` | archivo de salida |  |
| Planform | `PLANFORM` | opción: Straight; Bend (straight inlet, circular bend, straight outlet); Meanders (sine-generated); Composed from segments (text below); Centreline from a line layer; Imported file (TXT, CSV or DXF) | Straight |
| Width (m) | `WIDTH` | número | 1.0 |
| Straight: length (m) | `LENGTH` | número | 20.0 |
| Bend and meanders: straight inlet length (m) | `INLET` | número | 5.0 |
| Bend and meanders: straight outlet length (m) | `OUTLET` | número | 5.0 |
| Bend: radius on the axis (m) | `RADIUS` | número | 2.0 |
| Bend: angle (degrees) | `BEND_ANGLE` | número | 180.0 |
| Bend: turning | `TURN` | opción: Left (counter-clockwise); Right (clockwise) | Left (counter-clockwise) |
| Meanders: wavelength along the axis (m) | `WAVELENGTH` | número | 10.0 |
| Meanders: largest deflection angle theta0 (degrees) | `DEFLECTION` | número | 60.0 |
| Meanders: number of wavelengths | `COUNT` | número | 2.0 |
| Origin (inlet, centre of the inlet for curved channels) *(opcional)* | `ORIGIN` | punto | 0,0 |
| Direction at the inlet (degrees, counter-clockwise from x) | `ANGLE` | número | 0.0 |
| Element size (m) | `CELL` | número | 0.05 |
| Structured mesh over the channel | `STRUCTURED` | sí/no | sí |
| Bed level at the inlet (m) | `BED_Z0` | número | 0.0 |
| Bed slope along the axis (m/m, positive downstream) | `SLOPE` | número | 0.0 |
| Cross-section | `SECTION` | opción: Rectangular (flat bed); Trapezoidal; Compound (main channel and floodplains); V-shaped; Parabolic; Asymmetric trapezoid; Custom formula of n and b | Rectangular (flat bed) |
| Trapezoid: bottom width (m) | `BOTTOM_WIDTH` | número | 0.5 |
| Trapezoid: side slope (horizontal per 1 vertical) | `SIDE_SLOPE` | número | 1.0 |
| Trapezoid: bank height (m) | `BANK_HEIGHT` | número | 0.2 |
| Composed: segments, one per line *(opcional)* | `SEGMENTS` | texto | (texto de ejemplo) |
| Centreline: line layer (local coordinates for a flume) *(opcional)* | `CENTRELINE` | capa vectorial |  |
| Centreline: width field (optional; one width per line) *(opcional)* | `WIDTH_FIELD` | campo |  |
| Imported file (TXT, CSV or DXF) *(opcional)* | `IMPORT_FILE` | archivo |  |
| Imported file: read as | `IMPORT_AS` | opción: Centreline (the longest line; a third column gives widths); Outline (the largest closed line is the domain, closed lines inside it islands, open lines breaklines) | Centreline (the longest line; a third column gives widths) |
| Drawn or imported: corners sharper than this are elbows (degrees) *(avanzado)* | `CORNER` | número | 30.0 |
| Compound: floodplain rise towards the walls (m/m) *(avanzado)* | `FLOODPLAIN_SLOPE` | número | 0.0 |
| Asymmetric: right side slope (horizontal per 1 vertical; the left one is the side slope) *(avanzado)* | `RIGHT_SLOPE` | número | 2.0 |
| Main channel offset from the axis (m, left positive) *(avanzado)* | `SECTION_OFFSET` | número | 0.0 |
| Custom section: formula of n and b (height above the axis bed) *(opcional)* *(avanzado)* | `SECTION_FORMULA` | texto | 0.1*(n/b)**2 |
| Inflow (m³/s, optional) *(opcional)* | `Q` | número |  |
| Outlet water level (m, optional) *(opcional)* | `H` | número |  |

### Check the drawing of a mesh project

Id: `telemactoolkit:checkdrawing`

Finds what QGIS does not flag but the mesher pays for, measured against the element size wanted at each place: near-duplicate vertices (closer than 1% of the size), spikes where a line turns back on itself, corners of an outline closer than a quarter of the size, breakline vertices closer than a tenth of it (every breakline vertex is a node), parts of the outline or islands within half an element of each other (slivers), line ends and grid-zone corners just off the outline, invalid polygons, islands crossing the outline and overlapping domain polygons. One point per problem. With Fix, the safe fixes are made and saved: near-duplicates and small spikes removed, invalid polygons repaired, line ends snapped onto the outline.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh project | `PROJECT` | archivo |  |
| Default element size (0 = from the project) | `DEFAULT_SIZE` | número | 0.0 |
| Make the safe fixes and save them | `FIX` | sí/no | no |
| Drawing problems | `OUTPUT` | capa de salida |  |

### Compare meshing methods

Id: `telemactoolkit:comparemeshing`

Meshes one mesh project with each chosen method and measures every result, and any mesh files you add (from Blue Kenue, Salomé or an earlier run), on equal terms. The table and the HTML report give, per mesh:

- Element shape: smallest angle, share of elements below 30° and 20°, largest angle, aspect ratio (1 = equilateral), edge lengths.
- Size changes: area ratio of neighbouring elements (1 = equal).
- Outline: how far the mesh boundary strays from the drawn domain.
- Terrain (with a DEM): how far TELEMAC's linear bed inside each element is from the DEM (RMS, 95th percentile, largest) and, with a water level, the error in the water volume under it.
- Time step at a Courant number of 1, for the water level (or depth) and velocity given: the largest stable step TELEMAC's explicit parts allow, set by the smallest elements in deep, fast water.
- Nodes, elements, density and meshing time.

The best value of each row is marked in the report. No single score is given: a finer mesh follows the DEM better and needs a smaller time step; the report shows that trade-off, the choice stays yours.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh project (empty = only the mesh files below) *(opcional)* | `PROJECT` | archivo |  |
| Methods *(opcional)* | `METHODS` | opción: Balanced (built-in: regular nodes, smoothing); Quality-guaranteed (built-in: refine until every angle reaches the target); Gmsh Delaunay (fastest on large meshes); Gmsh Frontal-Delaunay (highest element quality); Gmsh MeshAdapt (most robust); Gmsh quadrilaterals (for programs that read them; TELEMAC gets two triangles each) | [0, 1, 3] |
| Mesh files to compare as well *(opcional)* | `MESHES` | capas |  |
| DEM (bed error, volume, time step) *(opcional)* | `DEM` | raster |  |
| Water level for the volume and time step (m, optional) *(opcional)* | `LEVEL` | número |  |
| Or a depth for the time step (m) | `DEPTH` | número | 2.0 |
| Velocity for the time step (m/s) | `VELOCITY` | número | 1.0 |
| Default element size (0 = the project's) *(avanzado)* | `DEFAULT_SIZE` | número | 0.0 |
| Comparison table | `TABLE` | capa de salida |  |
| Report *(opcional)* | `REPORT` | archivo de salida |  |
| Folder for the generated meshes (SELAFIN) *(opcional)* | `MESH_FOLDER` | folderDestination |  |

### Edit node values of a mesh (bottom, friction, …)

Id: `telemactoolkit:interpolatebottom`

Change a node variable of an existing SELAFIN geometry (BOTTOM, BOTTOM FRICTION or any other) without touching the mesh: same nodes, elements, boundary numbering, origin and date, so the boundary conditions file (.cli) still fits. Other variables are kept.

New values come, in priority order, from rasters (bilinear), survey points (linear on their triangulation), a formula, or a constant. The formula may use x and y, or the channel coordinates s and n along an axis line (e.g. carve a trapezoidal channel: 'z0 - 0.001*s + clip((abs(n) - 5)/2, 0, 3)').

Where: the whole mesh, or only the nodes inside the polygons, with an optional transition width along the polygon edge (the change fades from full inside to nothing at the edge).

How: replace, add (raise or lower by the value), raise to (keep the higher value: fills, dikes), lower to (keep the lower: dredging, cuts).

Typical edits: a new survey over a reach (rasters, polygon), a dredged channel (formula in s and n, lower to), a dike crest (constant, raise to, thin polygon), a friction zone (variable BOTTOM FRICTION, constant, replace).

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN) | `INPUT` | archivo |  |
| CRS of the mesh (empty = .prj next to the file, else project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Variable | `VARIABLE` | texto | BOTTOM |
| New values: rasters (priority order) *(opcional)* | `RASTERS` | capas |  |
| New values: survey points *(opcional)* | `POINTS` | capa vectorial |  |
| Survey value field (empty = point Z) *(opcional)* | `POINT_FIELD` | campo |  |
| New values: formula of x, y or s, n (lowest priority) *(opcional)* | `FORMULA` | texto |  |
| Channel axis line for s and n in the formula *(opcional)* | `AXIS` | capa vectorial |  |
| New values: constant (when no other source) *(opcional)* | `CONSTANT` | número |  |
| Operation | `OPERATION` | opción: Replace with the new values; Add the new values (raise or lower by); Raise to the new values (keep the higher: fills, dikes, embankments); Lower to the new values (keep the lower: dredging, cuts, channels) | Replace with the new values |
| Only inside these polygons (empty = whole mesh) *(opcional)* | `AREA` | capa vectorial |  |
| Transition width inside the polygon edge (m, 0 = sharp) | `BLEND` | número | 0.0 |
| Unit of a new variable *(opcional)* *(avanzado)* | `UNIT` | texto | M |
| Unit factor of rasters and points *(avanzado)* | `Z_FACTOR` | número | 1.0 |
| Fill nodes no source covers from the nearest node *(avanzado)* | `FILL` | sí/no | sí |
| Output geometry file | `OUTPUT` | archivo de salida |  |
| Open the mesh when done | `LOAD` | sí/no | sí |

### Element sizes from terrain (raster)

Id: `telemactoolkit:terrainsizes`

The element sizes terrain refinement asks for, as a raster to look at before meshing: small where linear elements would miss the DEM by more than the tolerance (banks, levees, embankments), the largest size on flat ground, and never growing faster than the growth rate. TELEMAC interpolates the bed linearly in each element; over an element of size h that misses a curvature kappa by about h² kappa / 8, hence h = sqrt(8 tolerance / kappa). Generate mesh applies the same sizes with its terrain tolerance.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| DEM | `DEM` | raster |  |
| Mesh project (its domain sets the extent) *(opcional)* | `PROJECT` | archivo |  |
| Extent (empty = the project's domain or the DEM) *(opcional)* | `EXTENT` | extent |  |
| Bed tolerance (m) | `TOLERANCE` | número | 0.05 |
| Smallest size (m) | `MIN_SIZE` | número | 2.0 |
| Largest size (m) | `MAX_SIZE` | número | 30.0 |
| Growth rate | `GROWTH` | número | 1.2 |
| Cell size (m, 0 = half the smallest size) *(avanzado)* | `CELL` | número | 0.0 |
| Sizes | `OUTPUT` | rasterDestination |  |

### Export or convert mesh

Id: `telemactoolkit:exportmesh`

Write an existing mesh (SELAFIN, 2DM, UGRID or any mesh QGIS reads) for another program: TELEMAC, SMS 2DM (BASEMENT, TUFLOW FV), Gmsh, GiD (Iber and other GiD-based programs), UGRID (Delft3D FM), MIKE 21 FM or ADCIRC/SCHISM. Optional boundary lines and roughness zones add boundary strings, material ids and friction. The file can be written in local coordinates (a point moved to (0, 0), optional rotation). HEC-RAS builds its own cells, so it is exported from Generate mesh, not here.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh | `INPUT` | archivo |  |
| CRS of the mesh (empty = from the file or the project) *(opcional)* | `MESH_CRS` | SRC |  |
| Boundary condition lines *(opcional)* | `BOUNDARIES` | capa vectorial |  |
| Roughness zones *(opcional)* | `ROUGHNESS_ZONES` | capa vectorial |  |
| Friction outside the roughness zones *(opcional)* *(avanzado)* | `FRICTION` | número |  |
| Output format | `FORMAT` | opción: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | SMS 2DM – BASEMENT, TUFLOW FV, SMS |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | sí/no | no |
| Coordinates in the file *(avanzado)* | `FRAME` | opción: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(opcional)* *(avanzado)* | `FRAME_POINT` | punto |  |
| Rotation of the local axes (degrees, counter-clockwise) *(avanzado)* | `FRAME_ANGLE` | número | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(avanzado)* | `STORE_OFFSET` | sí/no | sí |
| TELEMAC: double precision (SERAFIND; the steering file then needs GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes) *(avanzado)* | `DOUBLE` | sí/no | no |
| Mesh file | `OUTPUT` | archivo de salida |  |
| TELEMAC liquid boundaries template (.liq) *(opcional)* | `OUTPUT_LIQ` | archivo de salida |  |
| Open the result when done | `LOAD` | sí/no | sí |

### Extract part of a mesh

Id: `telemactoolkit:cutmesh`

Keep the elements whose centre lies inside the polygons (or outside them) and write the smaller mesh, with its node values, in any output format. Useful to cut a study reach out of a regional mesh, or to remove an area before meshing it again.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh | `INPUT` | archivo |  |
| CRS of the mesh (empty = from the file or the project) *(opcional)* | `MESH_CRS` | SRC |  |
| Polygons | `POLYGONS` | capa vectorial |  |
| Keep what is outside the polygons | `OUTSIDE` | sí/no | no |
| Output format | `FORMAT` | opción: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | sí/no | no |
| Coordinates in the file *(avanzado)* | `FRAME` | opción: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(opcional)* *(avanzado)* | `FRAME_POINT` | punto |  |
| Rotation of the local axes (degrees, counter-clockwise) *(avanzado)* | `FRAME_ANGLE` | número | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(avanzado)* | `STORE_OFFSET` | sí/no | sí |
| TELEMAC: double precision (SERAFIND; the steering file then needs GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes) *(avanzado)* | `DOUBLE` | sí/no | no |
| Mesh file | `OUTPUT` | archivo de salida |  |
| Open the result when done | `LOAD` | sí/no | sí |

### Generate mesh

Id: `telemactoolkit:generatemesh`

Build a mesh from a mesh project (or from individual layers, which take precedence) and write it for TELEMAC or another program.

Unstructured triangles fill the domain with an element size that grows away from refinement areas, points and lines at the growth rate; the mesh follows breaklines vertex by vertex and guide lines at the local size, keeps fixed points, and stitches in structured channel strips and regular grid zones node for node.

Bottom: rasters in priority order, then survey points, then a formula of x and y (for flumes); nodes nothing covers take the nearest value. Roughness zones give a material id per element and a friction value per node (BOTTOM FRICTION for TELEMAC). Boundary lines give the TELEMAC .cli (and a .liq template), 2DM nodestrings, Gmsh physical lines or ADCIRC/MIKE boundary lists.

Advanced: write the file in local coordinates (lower-left corner or a chosen point at (0, 0), optional rotation); TELEMAC files can keep the shift in their header.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh project *(opcional)* | `PROJECT` | archivo |  |
| Domain (polygons; holes are islands) *(opcional)* | `DOMAIN` | capa vectorial |  |
| Islands *(opcional)* | `ISLANDS` | capa vectorial |  |
| Refinement areas (size) *(opcional)* | `REFINEMENT_AREAS` | capa vectorial |  |
| Refinement points (size) *(opcional)* | `REFINEMENT_POINTS` | capa vectorial |  |
| Constraint lines (mode: breakline or guide) *(opcional)* | `CONSTRAINT_LINES` | capa vectorial |  |
| Fixed points *(opcional)* | `FIXED_POINTS` | capa vectorial |  |
| Structured channels (role: left, right, axis) *(opcional)* | `CHANNELS` | capa vectorial |  |
| Structured grid zones (dx, dy) *(opcional)* | `GRID_ZONES` | capa vectorial |  |
| Roughness zones (material, friction) *(opcional)* | `ROUGHNESS_ZONES` | capa vectorial |  |
| Boundary condition lines (type, name, value) *(opcional)* | `BOUNDARIES` | capa vectorial |  |
| Default element size (0 = from the project) | `DEFAULT_SIZE` | número | 0.0 |
| Size growth rate | `GROWTH` | número | 1.2 |
| Smoothing iterations *(avanzado)* | `SMOOTHING` | número | 8 |
| Follow the outline and lines within (m, 0 = 5% of the element size) *(avanzado)* | `OUTLINE_TOLERANCE` | número | 0.0 |
| Method (Gmsh ones need the gmsh package) | `METHOD` | opción: Balanced (built-in: regular nodes, smoothing); Quality-guaranteed (built-in: refine until every angle reaches the target); Gmsh Delaunay (fastest on large meshes); Gmsh Frontal-Delaunay (highest element quality); Gmsh MeshAdapt (most robust); Gmsh quadrilaterals (for programs that read them; TELEMAC gets two triangles each) | Balanced (built-in: regular nodes, smoothing) |
| Quality-guaranteed: smallest angle (degrees, 20 to 30) | `MIN_ANGLE` | número | 30.0 |
| Terrain refinement: bed tolerance (m, 0 = off; uses the first bottom raster) | `TERRAIN_TOLERANCE` | número | 0.0 |
| Terrain refinement: smallest element size (m, 0 = default size / 15) *(avanzado)* | `TERRAIN_MIN_SIZE` | número | 0.0 |
| Node limit (0 = no limit) *(avanzado)* | `MAX_NODES` | número | 2000000 |
| Bottom: rasters (priority order) *(opcional)* | `RASTERS` | capas |  |
| Bottom: survey points *(opcional)* | `POINTS` | capa vectorial |  |
| Survey value field (empty = point Z) *(opcional)* | `POINT_FIELD` | campo |  |
| Bottom: formula of x, y or channel coordinates s, n (empty = from the project) *(opcional)* | `FORMULA` | texto |  |
| Elevation unit factor *(avanzado)* | `Z_FACTOR` | número | 1.0 |
| Friction outside the roughness zones *(opcional)* *(avanzado)* | `FRICTION` | número |  |
| Output format | `FORMAT` | opción: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | sí/no | no |
| Coordinates in the file *(avanzado)* | `FRAME` | opción: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(opcional)* *(avanzado)* | `FRAME_POINT` | punto |  |
| Rotation of the local axes (degrees, counter-clockwise) *(avanzado)* | `FRAME_ANGLE` | número | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(avanzado)* | `STORE_OFFSET` | sí/no | sí |
| TELEMAC: double precision (SERAFIND; the steering file then needs GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes) *(avanzado)* | `DOUBLE` | sí/no | no |
| Mesh file | `OUTPUT` | archivo de salida |  |
| TELEMAC liquid boundaries template (.liq) *(opcional)* | `OUTPUT_LIQ` | archivo de salida |  |
| Open the result when done | `LOAD` | sí/no | sí |

### Import Blue Kenue lines (.i2s, .i3s)

Id: `telemactoolkit:bluekenuelines`

Reads a Blue Kenue line set (.i2s, or .i3s with a z per vertex) into a line layer with each line's value, and can add the lines to a layer of a mesh project: constraint lines (breaklines or guide lines), the domain or islands (closed lines become polygons), or refinement areas with the line value as their size. With the same input lines a mesh from Blue Kenue and meshes from this plugin can be compared (Compare meshing methods, which also reads Blue Kenue .t3s meshes). Blue Kenue files carry no CRS: give it.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Blue Kenue line set | `INPUT` | archivo |  |
| CRS of the coordinates *(opcional)* | `CRS` | SRC |  |
| Mesh project to add the lines to (optional) *(opcional)* | `PROJECT` | archivo |  |
| Add them as | `TARGET` | opción: Nothing (only the output layer); Constraint lines, breakline; Constraint lines, guide line; Domain (closed lines); Islands (closed lines); Refinement areas (closed lines, the value as size) | Nothing (only the output layer) |
| Lines | `OUTPUT` | capa de salida |  |

### Mesh quality

Id: `telemactoolkit:meshquality`

Element quality of a mesh (SELAFIN, 2DM, UGRID, ...) as a polygon layer: area, smallest and largest angle, longest edge, and a flag for elements below the minimum angle, above the maximum angle or with all their nodes on the boundary (over-constrained elements, which block the flow in a corner). The log gives the overall statistics.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh | `INPUT` | archivo |  |
| CRS of the mesh (empty = from the file or the project) *(opcional)* | `MESH_CRS` | SRC |  |
| Flag angles below (degrees) | `MIN_ANGLE` | número | 20.0 |
| Flag angles above (degrees) | `MAX_ANGLE` | número | 130.0 |
| Only flagged elements | `ONLY_FLAGGED` | sí/no | no |
| Element quality | `OUTPUT` | capa de salida |  |

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

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| CRS (projected, in metres) | `CRS` | SRC | ProjectCrs |
| Local coordinates (laboratory or schematic channel, no map CRS) | `LOCAL` | sí/no | no |
| Default element size (m) *(opcional)* | `DEFAULT_SIZE` | número | 0.0 |
| Mesh project | `OUTPUT` | archivo de salida |  |

### Optimise or adapt a mesh (MMG)

Id: `telemactoolkit:mmgadapt`

Runs MMG's mmg2d on a TELEMAC mesh (MMG is not part of the plugin: install it from mmgtools.org or the conda-forge package mmgsuite, and give the path of mmg2d_O3 once; it is remembered).

- Optimise: better element shapes with the same sizes.
- Adapt to the terrain: sizes from the DEM, small where linear elements would miss it by more than the bed tolerance (see Element sizes from terrain).
- Adapt to a result: sizes from the curvature of a variable of a TELEMAC result on this mesh (free surface, velocity, depth), with a tolerance in its unit, for a second run with elements where the flow needs them (solution-adaptive refinement).

The boundary nodes do not move, so the boundary conditions carry over; elements inside the 'keep' polygons (structured channels, bridges) are left as they are, and sizes next to them are held to their edges. Every variable is carried over (unchanged nodes exactly, new ones interpolated); with a DEM the bottom is sampled from it. The .cli next to the mesh is rebuilt.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| TELEMAC mesh (SELAFIN, .cli next to it) | `INPUT` | archivo |  |
| Mode | `MODE` | opción: Optimise (better shapes, the same sizes); Adapt to the terrain (DEM and bed tolerance); Adapt to a TELEMAC result (a variable and its tolerance) | Optimise (better shapes, the same sizes) |
| DEM (terrain mode; also resamples the bottom) *(opcional)* | `DEM` | raster |  |
| Tolerance (m for the terrain; the variable's unit for a result) | `TOLERANCE` | número | 0.05 |
| TELEMAC result on this mesh (result mode) *(opcional)* | `RESULT` | archivo |  |
| Variable of the result *(opcional)* | `VARIABLE` | texto | FREE SURFACE |
| Frames of the result | `WHEN` | opción: Last frame; Largest value over all frames | Largest value over all frames |
| Smallest size (m) | `HMIN` | número | 1.0 |
| Largest size (m, 0 = the mesh's longest edge) | `HMAX` | número | 0.0 |
| Growth rate between neighbours | `HGRAD` | número | 1.2 |
| Areas to keep as they are (polygons) *(opcional)* | `KEEP` | capa vectorial |  |
| CRS of the mesh (empty = project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| mmg2d program (empty = remembered, MMG2D variable or PATH) *(opcional)* | `MMG` | archivo |  |
| New mesh | `OUTPUT` | archivo de salida |  |

### Re-mesh inside polygons

Id: `telemactoolkit:remeshregion`

Replace the elements of a TELEMAC mesh whose centre lies inside the polygons by a new mesh of the same area: to refine around a structure, fix a patch of bad elements, or coarsen an area. Where the area borders old elements its outline keeps its nodes, so the new elements join the old ones node for node; the size grows from the outline's edges to the size inside. Where it reaches the mesh boundary, nodes can be added there and take its boundary conditions.

Every variable of every frame is carried over: unchanged nodes keep their values, new nodes get values interpolated in the original elements. The .cli next to the mesh is rebuilt for the new numbering. Size: a number, a field of the polygons, or the median edge of the area's outline.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| TELEMAC mesh (SELAFIN, .cli next to it) | `INPUT` | archivo |  |
| Areas to re-mesh | `POLYGONS` | capa vectorial |  |
| Element size inside (m, 0 = the outline's edges) | `SIZE` | número | 0.0 |
| Or a size field of the polygons *(opcional)* | `SIZE_FIELD` | campo |  |
| Method | `METHOD` | opción: Balanced; Quality-guaranteed (smallest angle 30 degrees) | Balanced |
| CRS of the mesh (empty = project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Re-meshed mesh | `OUTPUT` | archivo de salida |  |

### Rebuild a TELEMAC mesh after editing it in QGIS

Id: `telemactoolkit:rebuildedited`

QGIS edits mesh layers with its Mesh Digitizing tools (add, remove and move vertices, flip edges, split faces, force the mesh along lines), but the saved file keeps only the geometry: the bottom, the friction and the other variables are lost, and the .cli no longer matches the nodes. This tool gives them back: nodes that did not move keep their values exactly, new or moved nodes get values interpolated in the original elements (nearest original node outside them), and the .cli is rebuilt in TELEMAC order. Boundary nodes that existed keep their conditions; new ones take the type of the original boundary edge they lie on when both of its ends share it, otherwise they are walls.

Edited mesh: the file saved by QGIS (SELAFIN or any mesh format QGIS reads). Original: the TELEMAC mesh before editing, with its .cli next to it.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Edited mesh (saved by QGIS) | `EDITED` | archivo |  |
| Original TELEMAC mesh (SELAFIN, .cli next to it) | `ORIGINAL` | archivo |  |
| Original boundary conditions (empty = .cli next to the original) *(opcional)* | `CLI` | archivo |  |
| Distance under which a node has not moved (m) *(avanzado)* | `TOLERANCE` | número | 0.0001 |
| Rebuilt mesh | `OUTPUT` | archivo de salida |  |

### Regular grid mesh

Id: `telemactoolkit:regulargrid`

A structured mesh of nx by ny rectangular cells, starting at the origin (lower-left corner, (0, 0) by default) and rotated counter-clockwise by the angle. Cells are split into two triangles with alternating diagonals (no directional bias), or kept as quadrilaterals for formats that accept them. For grids inside a larger mesh, draw rectangles in grid_zones and use Generate mesh.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Origin (lower-left corner) | `ORIGIN` | punto | 0,0 |
| Cell size along x | `DX` | número | 1.0 |
| Cell size along y | `DY` | número | 1.0 |
| Cells along x | `NX` | número | 100 |
| Cells along y | `NY` | número | 20 |
| Rotation (degrees, counter-clockwise) | `ANGLE` | número | 0.0 |
| Diagonals | `PATTERN` | opción: alternate; right; left | alternate |
| CRS | `CRS` | SRC | ProjectCrs |
| Bottom: formula of x and y (optional) *(opcional)* | `FORMULA` | texto |  |
| Output format | `FORMAT` | opción: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | sí/no | no |
| Coordinates in the file *(avanzado)* | `FRAME` | opción: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(opcional)* *(avanzado)* | `FRAME_POINT` | punto |  |
| Rotation of the local axes (degrees, counter-clockwise) *(avanzado)* | `FRAME_ANGLE` | número | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(avanzado)* | `STORE_OFFSET` | sí/no | sí |
| TELEMAC: double precision (SERAFIND; the steering file then needs GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes) *(avanzado)* | `DOUBLE` | sí/no | no |
| Mesh file | `OUTPUT` | archivo de salida |  |
| Open the mesh when done | `LOAD` | sí/no | sí |

### River reach or network from centrelines or bank lines

Id: `telemactoolkit:riverreach`

Draws a natural river in a mesh project, ready for Generate mesh: the left bank, right bank and axis of a structured channel strip with matched cross-sections normal to the flow, the domain (or the project's own), and inflow (Q) and outflow (H) lines where the river meets the domain outline.

Banks:

- Centreline and a constant width.
- Centreline and a DEM: at every cross-section the terrain across the river is read, the bed found near the centreline, and each bank put where the terrain first rises the bankfull depth above the bed (searched up to the largest width). The widths are smoothed along the river and the axis moved midway between the banks.
- Two bank lines (one layer with two lines, in any order and direction).

With a DEM the river is turned to run downstream, from the higher bed to the lower one; without one, draw the centreline downstream. Jagged lines are smoothed over the smoothing length (default: the width) and, in tight bends, the cross-sections are turned until they no longer cross; any that still cross are reported.

Boundaries: the inflow covers the channel and the outflow the whole downstream end of the domain by default (water on the floodplain can leave). The lines stop short of their ends, so the bank and corner nodes stay walls. Values are optional; set them later in the boundary layer or the .liq file.

River networks: give several centrelines, each tributary ending on the river it joins (or lines meeting end to end at a bifurcation). Every reach between junctions gets its own strip, stopped short of the junction; the channel area at the junction is meshed with triangles of the strips' size, with its banks as guide lines, so confluences keep their corners. Free upstream ends get inflow lines, downstream ends outflow lines (lines must run downstream, or be turned by the DEM). An optional width field gives each line its own constant width.

Nodes across follow the width (default): where the river widens or narrows the number of nodes across changes by one row of triangles at a time, so the spacing across stays near the spacing along. In tight bends only the cross-sections there are turned, until none cross and no inner-bank cell is shorter than the spacing over the largest ratio; the others stay normal to the flow.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh project (empty = create one) *(opcional)* | `PROJECT` | archivo |  |
| New mesh project (when none is given) *(opcional)* | `NEW_PROJECT` | archivo de salida |  |
| River given as | `SOURCE` | opción: Centreline; Two bank lines | Centreline |
| Centreline(s) (line layer; several lines make a network; selected features if any) *(opcional)* | `CENTRELINE` | capa vectorial |  |
| Width field of the centrelines (one constant width per line, optional) *(opcional)* | `WIDTH_FIELD` | campo |  |
| Bank lines (line layer with the two banks) *(opcional)* | `BANKS` | capa vectorial |  |
| Width (centreline) | `WIDTH_MODE` | opción: Constant width; From the DEM (bankfull depth) | Constant width |
| Constant width (m) | `WIDTH` | número | 20.0 |
| DEM (banks from the terrain, downstream direction) *(opcional)* | `DEM` | raster |  |
| Bankfull depth above the bed (m) | `BANKFULL_DEPTH` | número | 1.5 |
| Largest width to search in the DEM (m) | `MAX_WIDTH` | número | 100.0 |
| Smallest width (m) *(avanzado)* | `MIN_WIDTH` | número | 0.0 |
| Cross-section spacing along the river (m) | `SPACING` | número | 10.0 |
| Nodes across the channel (0 = width / spacing, odd) | `NODES_ACROSS` | número | 0 |
| Nodes across follow the width (when 0 above) | `VARY_ACROSS` | sí/no | sí |
| Bends: inner-bank cells at least the spacing divided by *(avanzado)* | `MAX_RATIO` | número | 4.0 |
| Smoothing length of the centreline (m, 0 = the width) *(avanzado)* | `SMOOTHING` | número | 0.0 |
| Domain | `DOMAIN_MODE` | opción: Create it: the reach plus a floodplain on each side; Keep the domain of the project | Create it: the reach plus a floodplain on each side |
| Floodplain width on each side (m) | `FLOODPLAIN` | número | 200.0 |
| Element size on the floodplain (m, 0 = default size) | `FLOODPLAIN_SIZE` | número | 0.0 |
| Add inflow and outflow lines | `BOUNDARIES` | sí/no | sí |
| Inflow (Q) covers | `INFLOW_SPAN` | opción: The channel only; The whole end of the domain | The channel only |
| Outflow (H) covers | `OUTFLOW_SPAN` | opción: The whole end of the domain; The channel only | The whole end of the domain |
| Inflow (m³/s, optional) *(opcional)* | `Q` | número |  |
| Outflow water level (m, optional) *(opcional)* | `H` | número |  |
| Channel name | `NAME` | texto | river |

### Structured channel mesh

Id: `telemactoolkit:channelmesh`

A structured mesh of a channel between a left and a right bank line, optionally split by an axis (thalweg) line. Lines are grouped by the 'channel' field and told apart by the 'role' field (left, right, axis). Nodes across and spacing along come from the nodes_across and spacing fields, or from the defaults below. Use it alone, or let Generate mesh stitch it into the floodplain mesh.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Channel lines | `CHANNELS` | capa vectorial |  |
| Nodes across (default) | `CROSS_NODES` | número | 7 |
| Spacing along (default, m) | `INTERVAL` | número | 10.0 |
| Output format | `FORMAT` | opción: TELEMAC – SELAFIN geometry + .cli + .liq; SMS 2DM – BASEMENT, TUFLOW FV, SMS; Gmsh MSH 2.2; GiD ASCII mesh – GiD-based programs such as Iber; UGRID NetCDF – Delft3D FM / D-HYDRO; MIKE 21 FM mesh; ADCIRC / SCHISM grid (fort.14, hgrid.gr3); HEC-RAS – perimeter, breaklines and refinement regions | TELEMAC – SELAFIN geometry + .cli + .liq |
| Keep structured cells as quadrilaterals (formats that accept them) | `QUADS` | sí/no | no |
| Coordinates in the file *(avanzado)* | `FRAME` | opción: Map coordinates as drawn; Shift the lower-left corner to (0, 0); Shift a chosen point to (0, 0) | Map coordinates as drawn |
| Point moved to (0, 0) *(opcional)* *(avanzado)* | `FRAME_POINT` | punto |  |
| Rotation of the local axes (degrees, counter-clockwise) *(avanzado)* | `FRAME_ANGLE` | número | 0.0 |
| TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place *(avanzado)* | `STORE_OFFSET` | sí/no | sí |
| TELEMAC: double precision (SERAFIND; the steering file then needs GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes) *(avanzado)* | `DOUBLE` | sí/no | no |
| Mesh file | `OUTPUT` | archivo de salida |  |
| Open the mesh when done | `LOAD` | sí/no | sí |

### Suggest element sizes from a DEM

Id: `telemactoolkit:suggestsizes`

Sizes for a mesh from what the DEM can show, the terrain, the rivers and a node budget:

- smallest size: two pixels (finer elements only interpolate the same pixels);
- terrain tolerance: about the vertical accuracy of a DEM with such pixels (1 m LiDAR 0.05 m, 5 m DEMs 0.2 m, 15 m DEMs 0.5 m);
- default size (floodplain): what the median terrain inside the domain asks for at that tolerance, at least four pixels, at most a twentieth of the domain, within half the node budget;
- rivers (the project's structured channels): a spacing across of a eighth of the width (7 nodes or more), 1.5 times that along;
- the number of nodes with the terrain refinement (the tolerance is relaxed if it is over the budget) and the time step at Courant number 1.

Without a project, or before the domain is drawn, the DEM's own extent is used. With 'Apply', the default size, terrain tolerance and smallest size go into the project's settings (step 1 of the panel).

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| DEM | `DEM` | raster |  |
| Mesh project (domain and channels; optional) *(opcional)* | `PROJECT` | archivo |  |
| Node budget | `NODE_BUDGET` | número | 300000 |
| Typical depth for the time step (m) *(avanzado)* | `DEPTH` | número | 2.0 |
| Typical velocity for the time step (m/s) *(avanzado)* | `VELOCITY` | número | 1.5 |
| Apply to the project's settings | `APPLY` | sí/no | no |

## Malla (grupo Mesh)

### Check boundary conditions

Id: `telemactoolkit:checkboundaries`

Lists the liquid boundaries of a mesh with the numbers TELEMAC gives them, and checks them against the liquid boundaries file, the stage-discharge curves and the steering file when given (or found through the steering file): a column for every prescribed value, levels that leave a boundary dry, curves on boundaries that are not prescribed elevations, value counts in the steering file, two liquid boundaries touching without a wall node between them, boundaries with fewer than 5 nodes. The output layer draws each liquid boundary labelled with its number.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN geometry with BOTTOM) | `MESH` | archivo |  |
| Boundary conditions file (empty = .cli next to the mesh) *(opcional)* | `CLI` | archivo |  |
| Steering file (optional; gives the other files) *(opcional)* | `STEERING` | archivo |  |
| Liquid boundaries file *(opcional)* | `LIQ` | archivo |  |
| Stage-discharge curves file *(opcional)* | `STAGE` | archivo |  |
| CRS of the mesh (empty = project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Liquid boundaries | `OUTPUT` | capa de salida |  |

### Initial conditions file

Id: `telemactoolkit:initialconditions`

Write a SELAFIN file with the initial state of a run (velocities zero, water depth, free surface, bottom), for PREVIOUS COMPUTATION FILE with INITIAL TIME SET TO ZERO = YES. The water level is constant, or varies linearly along a channel axis from the upstream to the downstream level (nodes beside the axis take the level of their position along it), or the depth above the bed is constant. Nodes where the level is below the bed are dry.

Depth above the bed of the channel axes: for rivers and river networks. Every node takes the level of the nearest point of any axis line: the bed there (the nearest mesh node, smoothed along the line) plus the depth. The axis layer may be the mesh project's channels layer: only its axis lines are used. The banks and the floodplain above that level stay dry, and every branch starts with water.

Start close to the expected flow: for a river, the downstream level from the rating curve at the first discharge and a level rising upstream with the slope. A state far from the first boundary values makes TELEMAC drain or fill the model at the start, which can dry a discharge boundary.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN geometry with BOTTOM) | `MESH` | archivo |  |
| Initial water | `MODE` | opción: Constant water level; Water level along a channel axis (upstream and downstream levels); Depth above the bed; Depth above the bed of the channel axes (rivers and networks) | Constant water level |
| Water level (m), or upstream level along the axis | `LEVEL` | número | 0.0 |
| Downstream level along the axis (m) | `LEVEL_DOWN` | número | 0.0 |
| Channel axis (line, drawn downstream) *(opcional)* | `AXIS` | capa vectorial |  |
| Depth above the bed (m) | `DEPTH` | número | 0.0 |
| CRS of the mesh (for the axis; empty = project CRS) *(opcional)* | `MESH_CRS` | SRC |  |
| Initial conditions file | `OUTPUT` | archivo de salida |  |

### Liquid boundaries file from a table (hydrographs)

Id: `telemactoolkit:liquidboundariesfile`

Write TELEMAC's LIQUID BOUNDARIES FILE (.liq) from a table: a CSV or Excel sheet loaded in QGIS, or any layer, with a time column and one column per prescribed value (for example HEC-HMS output). The mesh and its .cli tell which columns TELEMAC expects, in its numbering: Q(n) for a discharge boundary, SL(n) for a water level (not needed where a stage-discharge curve gives it), U(n) and V(n) for velocities.

Columns: list them as COLUMN=field, for example 'Q(1)=flow_m3s; SL(2)=level'; or leave it empty and choose the value fields in the order of the expected columns (the log lists them). Times: numbers in the chosen unit, or dates and times; they are written in seconds from the first row, and the first date is given for ORIGINAL DATE OF TIME. Rows are sorted by time; rows with an empty value are skipped.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN geometry with BOTTOM) | `MESH` | archivo |  |
| Boundary conditions file (empty = .cli next to the mesh) *(opcional)* | `CLI` | archivo |  |
| Table (CSV, Excel or any layer) | `TABLE` | capa vectorial |  |
| Time field | `TIME_FIELD` | campo |  |
| Time unit | `TIME_UNIT` | opción: seconds; minutes; hours; days; date and time | hours |
| Value fields, in the order of the expected columns *(opcional)* | `VALUE_FIELDS` | campo |  |
| Or the columns by name, e.g. Q(1)=flow; SL(2)=level *(opcional)* | `COLUMNS` | texto |  |
| Stage-discharge curves file (its boundaries need no SL column) *(opcional)* | `STAGE` | archivo |  |
| Liquid boundaries file | `OUTPUT` | archivo de salida |  |

### Stage-discharge curve for an outflow boundary

Id: `telemactoolkit:stagedischarge`

Normal-depth rating curve Z(Q) of a liquid boundary, computed on the boundary itself (its nodes and their bottom, as TELEMAC sees the section) with Manning's formula Q = sqrt(S) A R^(2/3) / n, the section divided where the roughness changes. Written in TELEMAC's STAGE-DISCHARGE CURVES FILE format; an existing file keeps its other curves. The boundary must be a prescribed elevation (5 4 4); the steering file then needs STAGE-DISCHARGE CURVES = 1 for it (the log gives the lines).

Slope: the energy slope, usually the bed or water-surface slope of the reach downstream. Roughness: one Manning's n, or the BOTTOM FRICTION variable of the mesh (Manning law). Above the lower end of the section, its ends act as vertical walls, like the closed boundaries of the mesh. A warning is given if the flow at the section is supercritical, where a rating curve does not control the level.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN geometry with BOTTOM) | `MESH` | archivo |  |
| Boundary conditions file (empty = .cli next to the mesh) *(opcional)* | `CLI` | archivo |  |
| Liquid boundary number (0 = the only prescribed elevation) | `BOUNDARY` | número | 0 |
| Energy slope (m/m) | `SLOPE` | número | 0.001 |
| Manning's n | `MANNING` | número | 0.03 |
| Use the BOTTOM FRICTION of the mesh (Manning) when present | `USE_MESH_FRICTION` | sí/no | sí |
| Highest level of the curve (m, empty = highest ground) *(opcional)* *(avanzado)* | `TOP` | número |  |
| Points of the curve *(avanzado)* | `POINTS` | número | 40 |
| Stage-discharge curves file | `OUTPUT` | archivo de salida |  |
| Curve table *(opcional)* | `TABLE` | capa de salida |  |

### TELEMAC-2D steering file

Id: `telemactoolkit:steeringfile`

Write a TELEMAC-2D steering file (.cas) for a mesh, its .cli and the other files: liquid boundaries, stage-discharge curves, initial conditions. The liquid boundaries come from the .cli, in TELEMAC's numbering, with PRESCRIBED FLOWRATES and ELEVATIONS from the first row of the liquid boundaries file, STAGE-DISCHARGE CURVES for the boundaries in the curves file, and VELOCITY PROFILES = 4 (velocity in sqrt(h)) on discharge boundaries, which keeps them stable where the depth varies across the section. The time step can be estimated for a Courant number of 1 from a typical depth and velocity.

The numerical settings are those of this plugin's demo cases that ran in TELEMAC-2D v9.1 (wetting and drying with positive depths). Review every value before running.

| Parámetro | Clave | Tipo | Valor por defecto |
| --- | --- | --- | --- |
| Mesh (SELAFIN geometry with BOTTOM) | `MESH` | archivo |  |
| Boundary conditions file (empty = .cli next to the mesh) *(opcional)* | `CLI` | archivo |  |
| Liquid boundaries file *(opcional)* | `LIQ` | archivo |  |
| Stage-discharge curves file *(opcional)* | `STAGE` | archivo |  |
| Initial conditions file (empty = start dry or at a level) *(opcional)* | `INITIAL` | archivo |  |
| Or a constant initial water level (m) *(opcional)* | `INITIAL_LEVEL` | número |  |
| Duration (s) | `DURATION` | número | 3600.0 |
| Time step (s, 0 = estimate for Courant 1) | `TIME_STEP` | número | 0.0 |
| Typical depth for the estimate (m) | `DEPTH_GUESS` | número | 2.0 |
| Typical velocity for the estimate (m/s) | `VELOCITY_GUESS` | número | 1.0 |
| Results every (s) | `OUTPUT_PERIOD` | número | 600.0 |
| Friction law | `FRICTION_LAW` | opción: Manning; Strickler; Chezy; Nikuradse | Manning |
| Friction coefficient (where the mesh has none) | `FRICTION` | número | 0.03 |
| Title | `TITLE` | texto | TELEMAC-2D model |
| Start date and time (optional) *(opcional)* | `START` | datetime |  |
| Steering file | `OUTPUT` | archivo de salida |  |

[← Formatos, coordenadas y condiciones de frontera](08-formatos-coordenadas-fronteras.md) · [Contenido](README.md) · [Siguiente: Scripts con Python →](10-python.md)
