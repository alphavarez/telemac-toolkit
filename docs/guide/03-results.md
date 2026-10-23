# 3. Viewing results

## Opening a result

- *Open TELEMAC result…* (toolbar or menu), or drag `.slf`, `.res` or `.srf` files onto QGIS.
- **2D results** open as QGIS mesh layers with readable variable names ("Water depth (m)", "Velocity
  (m/s)") and a style per kind of variable: blues for depth, a terrain ramp for the bottom, Viridis for
  velocity, a critical-flow ramp for the Froude number.
- **3D results** (TELEMAC-3D): QGIS cannot read 3D SELAFIN files, so the plugin asks what to show — depth
  average, bottom plane, surface plane or any plane — writes a 2D file next to the result
  (`name_2d_average.slf`…) and opens that. The 2D file is reused while it is newer than the result.
- **French or English** variable names are both recognised.
- **Running simulations**: a result TELEMAC is still writing can be opened. The **refresh** button reads
  the new time steps and jumps to the last one.
- Results saved in a QGIS project are re-attached to the panel when the project is reopened.

The info line under the result list gives the number of nodes, elements and frames, the output interval,
the start date if the file has one, the coordinate system, and warnings (float32 precision of large
coordinates, an incomplete last frame, a sub-domain of a parallel run that must be merged first).

## Coordinates of a result

SELAFIN files store x and y but no coordinate reference system (CRS). The plugin decides, in this order:

1. **A `.prj` file next to the result** (`r2d.prj` or `r2d.slf.prj`, WKT or `EPSG:xxxx`).
2. **The project CRS**, when the file's coordinates lie inside that CRS's area of use. A 20 m flume with x
   from 0 to 20 does not fit UTM zone 14N; a mesh at x = 512 000, y = 2 155 000 does.
3. Otherwise it **asks**: local coordinates (laboratory flume or schematic channel) or a map CRS. The dialog
   shows the range of the coordinates and checks your choice ("✓ The coordinates fit this choice" or a
   warning). With *Save the choice next to the file* it writes the `.prj`, so the file opens directly next
   time.

![The coordinate question](../img/crs_dialog.png)

After opening, the map is **adjusted so the result is on screen**: the first layer of a project sets the
project CRS and the zoom; a result outside the current view is zoomed to; and a project in local
coordinates cannot show a result in a map CRS (or the reverse), so the project takes the CRS of the
result. To change the CRS later use the **coordinates button** of the Results panel; to re-centre the map
use the **zoom button** (or right-click the layer → *Zoom to Layer(s)*).

QGIS's status bar shows **"Unknown CRS"** for local coordinates. That is QGIS's label for any CRS without
an EPSG code; it is not an error. The Results and Mesh panels name it "local coordinates".

The origin offset some files store in their header (IPARAM 3 and 4) is applied, so coordinates are
absolute.

## Time

The Results panel, the QGIS Temporal Controller, the plots and the cursor share one clock.

- **First / previous / play-pause / next / last** buttons and the slider; **Speed** in frames per second;
  **Loop**.
- The time label shows the elapsed time ("8 h (28 800 s)"), the frame number and, when the file has a
  start date, the date and time.
- **Peak** (button at the end of the slider): jumps to the frame where the variable shown on the map is
  highest — at the first probe point if there is one, otherwise of its mean over the mesh.
- Click a time series plot to move the clock to that time.

## What the map shows

- **Map shows**: any variable of the file, plus vector magnitudes (velocity, unit discharge, wind) when the
  file has the components but not the magnitude.
- **Arrows**: the vector drawn as arrows (*Off*, *Velocity*, *Unit discharge*, *Wind* … — every vector of
  the file). Choosing a vector in *Map shows* moves the arrows to it, so the colours and the arrows show the
  same quantity. A scalar (water depth, free surface, bottom …) has no direction: it keeps the arrows that
  were shown, velocity at first, which is the usual view of depth with the flow direction on top.
  Vectors are the pairs `… U / … V`, `… ALONG X / … ALONG Y` and `… SUIVANT X / … SUIVANT Y`, the names QGIS
  reads as vectors; other `X / Y` pairs (`VENT X`, for example) are two separate maps.
- **Dry below**: water depth under this value is not drawn (transparent), on the map and in profiles.

Styles can be changed further in the layer properties (*Symbology*), as for any QGIS mesh layer.

### Arrow settings

Open **Arrow settings** in the Results panel. Changes apply at once, are remembered for the next results,
and are saved with the QGIS project.

| Setting | What it does |
| --- | --- |
| **Placement** | *Regular grid*: one arrow per cell of a grid on the screen, interpolated from the mesh, so fine parts of the mesh are not cluttered and the density stays the same at every zoom. *Every mesh node*: one arrow per node (dense meshes become crowded). |
| **Spacing** | Grid spacing in screen pixels. Smaller = more arrows. For narrow channels use 15–25 px so arrows fall inside the channel. |
| **Length** | *Proportional (fastest = longest)*: length proportional to the speed, the fastest value of the whole result gets the longest arrow. *Proportional, fixed scale*: a fixed number of mm per m/s, the same for every result — use it to compare scenarios. *Same length*: direction only (the colour of the map gives the speed). |
| **Longest arrow** | Length of the fastest arrow (or of every arrow with *Same length*). *auto* = 90 % of the spacing, so arrows never overlap. |
| **Scale** | mm per unit (m/s) for the fixed-scale mode. |
| **Hide below** | Arrows below this magnitude are not drawn (dry land, still water). |
| **Colour** | One colour, or *Coloured by magnitude* with the velocity ramp. |
| **Line width** | In mm. |

The grey line under the settings gives the **reference**: "0.491 m/s = 6.7 mm (fastest in the result:
0.491 m/s)". Put it in a map legend or a report caption.

![The same frame with 0.3.0 arrows (top) and 0.4.0 defaults (bottom)](../img/arrows_before.png)
![](../img/arrows_after.png)

## Cursor, time series and profiles

**Cursor** tab: turn on the live cursor and move the mouse over the mesh. Every variable at the current
time is interpolated inside the triangle under the mouse, plus the magnitude and the direction (bearing,
0° = north) of each vector.

**Time series** tab: turn on the probe and click points on the mesh (right-click clears them). Choose the
variable above the plot. A vertical line marks the current time; clicking the plot moves the clock.

**Profile** tab: turn on the profile tool, click the vertices of a line and right-click to finish. *Water
profile* draws the bed and the water surface (dry parts left blank); any other variable can be chosen.
The profile follows the clock, so playing the animation animates it.

**All plots**: right-click to copy the data, save it as CSV, or save the picture as PNG.

## Other result tools

In *Plugins → TELEMAC Toolkit* and the Processing toolbox: SELAFIN file information (dimensions, variables,
precision, mesh quality), Extract 2D from a 3D result, Extract frames and variables, Time series at points
(a table for every point of a layer), Mesh boundary as lines, and the [calculator and integrals](04-calculator-integrals.md).

[← First steps](02-first-steps.md) · [Contents](README.md) · [Next: Calculator →](04-calculator-integrals.md)
