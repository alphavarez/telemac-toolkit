# 11. Troubleshooting and FAQ

**I opened a result and the map is blank.**
The result is outside the current view or in a CRS that does not match the project. Press the **zoom
button** of the Results panel (or right-click the layer → *Zoom to Layer(s)*). If it is still not visible,
check its coordinates with the **coordinates button**: a laboratory model must be in local coordinates, a
field model in its map CRS. Since 0.4.0 the plugin adjusts the map itself when it opens a result.

**The mesh appears in the wrong place, or the scale bar shows thousands of kilometres.**
The result was given the wrong CRS (for example metres read as degrees). Use the coordinates button and
tick *Save the choice next to the file*.

**The status bar says "Unknown CRS".**
That is QGIS's label for the local laboratory grid, as for any CRS without an EPSG code. It is not an
error.

**Arrows are huge, crowded or overlapping.**
*Arrow settings*: placement *Regular grid*, a spacing that suits the channel (15–25 px for narrow
channels), length *Proportional* with *Longest arrow* on auto. Raise *Hide below* to remove arrows in still
water.

**No arrows inside a narrow channel.**
The grid points fall outside the channel at this zoom: lower the spacing, zoom in, or use *Every mesh
node*.

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

[← Python](10-python.md) · [Contents](README.md)
