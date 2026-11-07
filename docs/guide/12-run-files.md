# 12. Boundary conditions and run files

After *Generate mesh*, step **6 · Run files** of the Mesh workflow panel prepares the rest of a TELEMAC-2D
run. Each button opens a Processing tool with the generated mesh filled in. The tools work on any SELAFIN
mesh with its `.cli` and are also in the toolbox under *Boundary conditions and run*.

| Button | Writes |
| --- | --- |
| *Rating curve (outflow)…* | the stage-discharge curves file |
| *Hydrographs (.liq) from a table…* | the liquid boundaries file |
| *Initial state…* | an initial conditions file (SELAFIN) |
| *Steering file (.cas)…* | the TELEMAC-2D steering file |
| *Check boundary conditions…* | a layer of the liquid boundaries with their TELEMAC numbers, and a report |

The liquid boundaries are always numbered as TELEMAC numbers them, so the columns of the `.liq` file, the
curves and the lists of the steering file line up with what TELEMAC reads.

## Rating curve for the outflow

A prescribed water level at the downstream end is rarely known. A **stage-discharge curve** lets TELEMAC set
the level from the discharge that reaches the boundary. The tool computes the normal-depth curve on the
boundary itself, with its nodes and their bottom, exactly as TELEMAC sees the section:

- Manning's formula Q = √S · Σ Aᵢ Rᵢ^(2/3) / nᵢ, with the section divided where the roughness changes (the
  usual divided-channel method).
- *Energy slope*: usually the bed or water-surface slope of the reach downstream.
- *Manning's n*: one value, or the `BOTTOM FRICTION` variable of the mesh.
- Above the lower end of the section, its ends act as vertical walls, like the closed boundaries of the mesh.
- A warning is given if the flow at the section becomes supercritical, where a rating curve does not control
  the level.

The file follows TELEMAC's format, and an existing file keeps its other curves. The boundary must be a
prescribed elevation (`H`, 5 4 4). The log gives the steering-file lines: `STAGE-DISCHARGE CURVES = 1;0`
(one value per liquid boundary) and `STAGE-DISCHARGE CURVES FILE`.

## Hydrographs from a table

Load a CSV or Excel sheet in QGIS (for example HEC-HMS output), then choose:

- the **time field** and its unit: seconds, minutes, hours, days, or dates and times;
- the **value fields**, in the order of the columns TELEMAC expects (the log lists them: `Q(n)` for a discharge
  boundary, `SL(n)` for a level, `U(n)` and `V(n)` for velocities), or the columns by name, for example
  `Q(2)=flow_m3s; SL(1)=level`.

A boundary with a stage-discharge curve needs no `SL` column; give the curves file and it is left out. Times are
written in seconds from the first row. With dates, the first one is given for `ORIGINAL DATE OF TIME` and
`ORIGINAL HOUR OF TIME`.

## Initial state

TELEMAC starts more smoothly from a state close to the expected flow. The tool writes velocities of zero, the
water depth, the free surface and the bottom, from:

- a **constant water level**;
- a **level along a channel axis**, from an upstream to a downstream level (each node takes the level of its
  position along the axis), which suits a river;
- a constant **depth above the bed**;
- a **depth above the bed of the channel axes**, for rivers and river networks: every node takes the bed at
  the nearest point of any axis line plus the depth, so the banks and floodplain above that level stay dry and
  every branch starts with water. The axis layer can be the mesh project's channels layer (only its axis lines
  are read).

Nodes where the level is below the bed are dry. For a river, use the level of the rating curve at the first
discharge downstream, and that level plus slope × length upstream. A state far from the first boundary values
makes TELEMAC drain or fill the model at the start, which can leave a discharge boundary dry.

## Steering file

The steering file names the files above relative to itself, and contains:

- the liquid boundaries in TELEMAC's order, with `PRESCRIBED FLOWRATES` and `PRESCRIBED ELEVATIONS` from the first
  row of the `.liq` file;
- `STAGE-DISCHARGE CURVES` for the boundaries that have a curve;
- `VELOCITY PROFILES = 4` (velocity proportional to √h) on discharge boundaries;
- the time step: given, or estimated for a Courant number of 1 on the smallest element, from a typical depth and
  velocity;
- the friction law and coefficient; when the mesh has a `BOTTOM FRICTION` variable, TELEMAC reads it from there;
- numerical settings for wetting and drying with positive depths, those of this plugin's demo cases that ran in
  TELEMAC-2D v9.1, including `MAXIMUM NUMBER OF ITERATIONS FOR ADVECTION SCHEMES = 200` and, from 0.5.0, the
  NERD scheme for the velocities (`SCHEME FOR ADVECTION OF VELOCITIES = 14`).

Review every value before running; the comments at the top list the liquid boundaries.

### Tested in TELEMAC-2D

The whole chain was run in TELEMAC-2D v9.1.1 on a 2 km meandering river: river helper with the DEM, *Generate
mesh*, rating curve on the outflow, a 60 m³/s inflow from a CSV file, initial level along the axis, steering
file. After 3 hours the inflow was 59.9 m³/s, the outflow 60.0 m³/s, the volume error 10⁻¹⁵, and the level at
the outflow sat on the rating curve (4.212 m) within a millimetre. This run is part of the test suite.

Four findings from these runs are now built in:

- **Velocity profile on discharge boundaries.** With a constant velocity across the inflow (TELEMAC's default),
  the nodes next to the banks, where the water is shallower, emptied until the run stopped with
  `DEBIMP: PROBLEM ON BOUNDARY`. The √h profile keeps them wet.
- **Diagonals of structured cells.** Channel strips and grid zones used to cut their cells along alternating
  diagonals, so interior nodes had four and eight triangles in turn. TELEMAC's depths then oscillated from node
  to node. Cells are now cut along their shorter diagonal (one direction for grids), and every interior node has
  six triangles.
- **Advection of velocities (0.5.0).** With the method of characteristics, the same river on a finer mesh
  (time step 0.24 s instead of 0.55 s) stopped after 616 s: the inflow turned supercritical at a shallow node,
  dried out, and TELEMAC stopped with `DEBIMP: PROBLEM ON BOUNDARY`. The original mesh failed in the same way,
  so it was the numerics, not the mesh. Two other remedies were rejected: Thompson's method on the inflow ran
  but let 67 m³/s in instead of the 60 prescribed, and the virtual-depth velocity profile (5) failed in the
  same way. With the NERD scheme (14) both time steps ran for the 3 hours, with 60.000 m³/s in and out. It is
  about two to three times slower than characteristics; set `SCHEME FOR ADVECTION OF VELOCITIES = 1` back
  for speed once a model is known to be stable.
- **A time step that divides the duration (0.5.0).** TELEMAC runs round(duration / time step) steps, so a
  0.61 s step over 3 hours ended its last step at 10 800.05 s, past the end of a hydrograph that stopped at
  10 800 s, and the run stopped with `OUT OF RANGE OF THE FILE OF LIQUID BOUNDARIES` (TELEMAC allows 1 ms).
  The steering file now takes the largest step under the Courant limit that divides the duration exactly,
  with at most four decimals (0.61 s becomes 0.60 s, 18 000 steps). A time step you give is adjusted the same
  way, with a warning, and the tool warns when the liquid boundaries file ends before the duration.

## Check boundary conditions

The check lists every liquid boundary with its number, type, nodes, length and bed range, draws it on the map
labelled with its number, and warns about:

- a prescribed value without a column in the `.liq` file;
- a level that never rises above the bed of its boundary;
- columns for boundaries the mesh does not have;
- times that do not increase;
- a stage-discharge curve on a boundary that is not a prescribed elevation;
- lists of the steering file with the wrong number of values;
- two liquid boundaries of different types touching without a wall node between them;
- boundaries with fewer than five nodes.

Given the steering file, it finds the other files through it.

[← Troubleshooting](11-troubleshooting.md) · [Contents](README.md)
