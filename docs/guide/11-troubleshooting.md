# 11. Troubleshooting and FAQ

**Generate mesh stays at 6 % (or any low percentage) and never finishes.**
Before 0.5.0 this was almost always an element size far too small for the size of the domain: the
progress bar cannot move while the program places millions of nodes. The sizes (default size, refinement
areas, points and lines) are in map units, so 0.05 or 0.2 on a domain 100 km wide is impossible. From 0.5.0
the plugin estimates the number of nodes first and stops with a message giving the estimate, the domain
size and the smallest size asked for. Raise the sizes (a river reach of 5 km is usually fine with 10–25 m
in the channel and 50–100 m on the floodplain) or shrink the refined areas. If you really need more than
2 million nodes (a long reach of a large river at fine resolution), raise the **Node limit** in step 4. Check the unit too: coordinates in an unknown CRS are read as
metres.

**I opened a result and the map is blank.**
The result is outside the current view or in a CRS that does not match the project. Press the **zoom
button** of the Results panel (or right-click the layer → *Zoom to Layer(s)*). If it is still not visible,
check its coordinates with the **coordinates button**: a laboratory model must be in local coordinates, a
field model in its map CRS. Since 0.4.0 the plugin adjusts the map itself when it opens a result.

**The mesh appears in the wrong place, or the scale bar shows thousands of kilometres.**
The result was given the wrong CRS (for example metres read as degrees). Use the coordinates button and
tick *Save the choice next to the file*.

**The status bar says "Unknown CRS".**
QGIS's label for a CRS it has no name or code for. From 0.5.0 the plugin registers its local laboratory grid
as a QGIS user CRS the first time it uses it, so QGIS shows `USER:100000 - Local metric grid (laboratory)`
(the number may differ) and lists the layers with that code. A map CRS without an EPSG code is offered one
when the project opens (see step 1 in the [mesh workflow](05-mesh-workflow.md)). Either way the coordinates
are right.

**Generate mesh stops with "OSError: [Errno 22] Invalid argument" (or "cannot be written").**
On Windows a file that QGIS has loaded (a mesh layer of an earlier run) cannot be overwritten. From 0.5.0
the panel removes that layer before meshing and loads the new mesh in its place; the Processing tool does the
same, or writes the mesh beside it (`name_2.slf`) when it runs in the background, and says so. Before 0.5.0:
remove the old mesh layer from the Layers panel and run again.

**The mesh boundary does not match the domain outline exactly: a little inside here, a little outside there.**
The outline is resampled to the element size, and each boundary edge is a straight chord between two points
of the drawing; where the drawing bends gently between them (a vertex turning less than 30°), the edge cuts
across the bend: inside on a convex bend, outside on a concave one, mostly along the sides, away from the
corners. Before 0.5.0 the gap could reach about 10 % of the element size. From 0.5.0 such vertices stay
nodes where needed, so the gap stays below **Follow the outline within** in step 4 (*auto*: 5 % of the
element size, e.g. 1 m for 20 m elements); the log gives the largest one left. Set a distance in metres
(0.01 m keeps practically every vertex) for a tighter fit, or draw that stretch as a breakline, which
always keeps every vertex. 0.5.0 also fixes outlines digitized with many close vertices (a vertex every
metre or two, zigzagging): the mesher could hang at *Triangulating*.

**Re-mesh inside polygons says it "did not join the old elements cleanly".**
The new elements of an area must meet the old ones node for node; when they cannot (a polygon cutting very
thin elements, a size far from the edges along the area's outline), the tool says where and changes
nothing. Use a size closer to the edges of the outline (the log gives their median), or move the polygon
slightly.

**At a scale of 1:2 the mesh corners are a few millimetres off the drawn corners.**
SELAFIN stores coordinates in single precision (7 significant digits): at 125 000 m a coordinate rounds to
the nearest 7.8 mm, at 189 000 m to 15.6 mm, so before 0.5.0 a node could sit up to 8 mm from the drawing.
From 0.5.0 a mesh more than 4 km from the origin is written relative to an integer origin (IPARAM 3 and 4,
which TELEMAC and QGIS read), and nodes keep their position to 0.01 mm. Generate the mesh again.

**Arrows are huge, crowded or overlapping.**
*Arrow settings*: placement *Regular grid*, a spacing that suits the channel (15–25 px for narrow
channels), length *Proportional* with *Longest arrow* on auto. Raise *Hide below* to remove arrows in still
water.

**No arrows inside a narrow channel.**
The grid points fall outside the channel at this zoom: lower the spacing, zoom in, or use *Every mesh
node*.

**The arrows always show velocity, even when I choose another variable.**
Before 0.4.3 only velocity could be drawn as arrows, and *Map shows* did not move them. Since 0.4.3 the
*Arrows* list has every vector of the file (velocity, unit discharge, wind) and the arrows follow the vector
you choose in *Map shows*. A scalar (water depth, free surface) has no direction, so the arrows stay as they
were. A vector whose names QGIS does not read as a pair (for example `VENT X` / `VENT Y`) is shown as two
separate maps, without arrows.

**TELEMAC stops with "SL: MORE PRESCRIBED ELEVATIONS ARE REQUIRED" (or the same for flowrates).**
The liquid boundaries in the steering file or the `.liq` file are in the wrong order. Use the numbers and
the steering-file lines the plugin prints when it writes the `.cli`, and the `.liq` template it writes
(see [Liquid boundary numbers](08-formats-coordinates-boundaries.md#liquid-boundary-numbers)). Versions
before 0.4.0 started the numbering at the first boundary node instead of TELEMAC's start point, so some
meshes (a U-shaped bend, for example) got the order wrong: write the `.cli` and `.liq` again with 0.4.0.

**TELEMAC ignores my FRICTION COEFFICIENT.**
The geometry file has a `BOTTOM FRICTION` variable (roughness zones or a default friction), and TELEMAC-2D
reads it instead (the listing says `FONSTR : FRICTION COEFFICIENTS READ IN THE GEOMETRY FILE`). Remove the
zones and the default friction to use the steering-file value.

**"The formula uses s or n but no channel axis line is given."**
The project needs an axis line in its channels layer (role `axis`); with several axes, set `bottom_axis`
in the project settings to the channel name. In *Edit node values*, give the axis layer.

**The calculator says a name is unknown.**
Use the TELEMAC symbol (`H`, `U`…) or the variable name with `_` for spaces (`WATER_DEPTH`). The expression
box opened from the Results panel lists the names of the current file.

**The calculator says the output name is too long.**
SELAFIN variable names have at most 16 characters.

**Derivatives look noisy at the edge of the water.**
Velocity drops to zero at dry nodes. Mask: `where(H > 0.01, curl(U, V), 0)`.

**A 3D result opens as 2D.**
QGIS cannot read 3D SELAFIN files; the plugin shows a depth average or a plane written to
`<name>_2d_<choice>.slf`. Open the result again to choose another plane.

**"This is a sub-domain file of a parallel run".**
Merge the partitions with TELEMAC's `gretel` first.

**The mesher reports a constraint line that touches no boundary node / a boundary line that touches no
node.**
Snap the line to the domain edge (snapping toolbar) and check that it lies within half an element of the
boundary.

**Where do I report a problem?**
In the GitHub repository's *Issues*, with the QGIS version, the message from *View → Panels → Log Messages
→ TELEMAC Toolkit* and, if possible, a small file that shows it.

[← Python](10-python.md) · [Contents](README.md) · [Next: Boundary conditions and run files →](12-run-files.md)
