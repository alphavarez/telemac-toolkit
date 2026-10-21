# 5. The mesh workflow

Toolbar → **Mesh workflow panel**. The panel has five steps. Each step runs a Processing tool, so the same
work can be done from the toolbox, a model or a script ([reference](09-processing-reference.md)).

![Mesh workflow panel, step 4](../img/mesh_panel_export.png)

## Step 1 · Project

- **New…** creates a *mesh project*: a GeoPackage with one empty, styled layer per kind of object, added to
  the map in a group, plus a settings table. Answer **Yes** for **local coordinates** (a laboratory flume
  or a schematic channel drawn from (0, 0)); **No** to choose a map CRS. Use a projected CRS in metres,
  never a geographic one: sizes would be in degrees.
- **Open…** opens an existing project (projects from 0.2.0 included). When a project's layers are already
  on the map, the panel finds them by itself.
- **Coordinates** shows the project CRS ("local (laboratory)" or an EPSG code).
- **Default size**: element size where no refinement applies, in metres. *auto* uses about 1/40 of the
  square root of the domain area (0.1 m for a 20 m × 1 m flume, 25 m for a 2 km² valley) and says so in
  the log.
- **Growth rate**: how fast the size grows away from refined areas (1.2 = 20 % per element). Lower is
  smoother and has more elements.
- **Smoothing**: iterations that move interior nodes to improve angles (8 is a good default).
- **Channel or flume template…** draws a flume, a bend or meanders with their bed formula — see
  [Channel and flume templates](06-channel-templates.md).

## Step 2 · Draw

The table lists the project layers and how many objects each has. Select a layer and **Draw in this
layer** (or use QGIS's own editing tools: snapping, tracing, copying from other layers all work).

| Layer | What it is for | Fields |
| --- | --- | --- |
| **Domain (holes are islands)** | Area to mesh. Holes in a polygon are left out. Several polygons are meshed together. | `name`, `size` (optional) |
| **Islands (not meshed)** | Areas left out: buildings, piers, dry islands | `name` |
| **Refinement areas** | Finer (or coarser) elements inside | `size` (m) |
| **Refinement points** | Element size at a point: structures, gauges | `size` |
| **Constraint lines** | Lines the mesh follows. `mode` = **breakline** (every vertex kept: dikes, road crests, banks, walls) or **guide** (resampled to the local size: thalwegs, contours). Empty mode = breakline. | `name`, `mode`, `size` (optional) |
| **Fixed points** | Nodes that must exist: gauges, outfalls, structure corners | `name`, `size` (optional) |
| **Structured channels** | A channel meshed as a structured strip: one **left** and one **right** bank line and an optional **axis**, told apart by `role`, grouped by `channel` | `channel`, `role`, `nodes_across`, `spacing` |
| **Structured grid zones** | Rectangles (any orientation) meshed as regular grids | `name`, `dx`, `dy` |
| **Roughness zones** | Material id per element and friction value per node; later zones win where they overlap | `name`, `material`, `friction` |
| **Boundary conditions** | Lines along the domain edge with a boundary type | `name`, `type`, `value` |

Tips:

- **Snap** lines to the domain (QGIS snapping toolbar) so they end exactly on its edge.
- **Structured channels**: banks run in the same direction (the plugin reverses one if needed);
  `nodes_across` with an axis must be odd; `spacing` is the distance between cross-sections. When the
  left bank, right bank and axis have the **same number of vertices and each vertex k forms a
  cross-section**, the strip uses those cross-sections as drawn (normal to the axis in bends); otherwise
  each bank is spaced by its own length.
- **Boundary types**: `wall`, `Q` (discharge, 4 5 5), `H` (water level, 5 4 4), `QH`, `UV`, `UVH`, `free`,
  `waves`, or a TELEMAC code such as `455`. The boundary line only has to run along the edge; every
  boundary node within half an edge length takes its type. Stop lines a little short of corners so the
  corner nodes stay walls.
- **Check the drawing** reports what the mesher will see — domains, islands, structured parts, lines,
  problems — before meshing.

## Step 3 · Bottom and roughness

- **Rasters**: tick the DEMs to use and drag to set their priority; the first raster that covers a node
  gives its value (bilinear). Rasters in another CRS are sampled after transforming the nodes.
- **Survey points**: a point layer and its value field (or the point Z), interpolated linearly on their
  triangulation.
- **Formula**: the bottom as a formula of `x` and `y` (`0.3 - 0.002*x`), or of the channel coordinates
  `s` and `n` along a channel axis ([templates](06-channel-templates.md#bed-formulas-in-channel-coordinates)).
  Used where no raster or point covers a node.
- Nodes nothing covers take the value of the nearest node that has one (and the log says how many).
- **Default friction**: friction value outside the roughness zones. With roughness zones, TELEMAC files
  get a `BOTTOM FRICTION` variable, which TELEMAC-2D reads from the geometry file (the listing says
  `FONSTR : FRICTION COEFFICIENTS READ IN THE GEOMETRY FILE`) and uses instead of the steering-file value;
  set `LAW OF BOTTOM FRICTION` to match your values (4 = Manning). Other formats get a material id per
  element.

## Step 4 · Mesh and export

- **Format**: TELEMAC (SELAFIN + `.cli`) or another program — see
  [Formats](08-formats-coordinates-boundaries.md#formats).
- **Keep structured cells as quadrilaterals**, for formats that accept them.
- **Coordinates**: map coordinates as drawn; lower-left corner at (0, 0); or a chosen point at (0, 0),
  with an optional rotation. For TELEMAC, *keep the shift in the file header* lets TELEMAC compute near
  (0, 0) while QGIS still shows the mesh in place.
- **Mesh file**, then **Generate mesh**. The log reports the mesh (nodes, elements, smallest angle),
  where the bottom came from, each liquid boundary with **the number TELEMAC will give it**, and the
  steering-file lines (`PRESCRIBED FLOWRATES`, `PRESCRIBED ELEVATIONS`). For TELEMAC the panel also
  writes a `.liq` template next to the mesh, with one column per prescribed value (`Q(1)`, `SL(2)`…) in
  TELEMAC's numbering: replace its values with your hydrographs.

How the mesh is built: the domain is filled with triangles whose size follows the refinement objects and
grows at the growth rate; breaklines keep their vertices, guide lines are resampled, crossings become
nodes; QGIS's GEOS library makes the Delaunay triangulation and missing constraint segments are split
until they are edges; oversized triangles get a node at their circumcentre; interior nodes are smoothed
without inverting any element; structured channels and grids are stitched in node for node; triangles
with three boundary nodes are fixed by an edge flip. A 30 000-node mesh takes about 8 seconds.

## Step 5 · Check

- **Mesh quality layer**: one polygon per element with its area, smallest and largest angle, longest edge,
  and a flag for angles below 20° or above 130° and for elements with all their nodes on the boundary
  (they block the flow in corners). The log gives the statistics.
- **Find** a node or an element by its number, as TELEMAC listings print them: the map zooms to it and the
  log gives its coordinates (and bottom for nodes).

[← Calculator](04-calculator-integrals.md) · [Contents](README.md) · [Next: Channel and flume templates →](06-channel-templates.md)
