# 4. Calculator, integrals and fluxes

## Calculator

Open it from the **calculator button** of the Results panel (the expression box then lists the variables
of the current result) or from *Plugins → TELEMAC Toolkit → Calculator*. It writes a **new SELAFIN file**
with the variables you define, plus the original ones if *Keep the variables of the file* is ticked, and
opens it in the Results panel, where colours, arrows, probes and profiles work on the new variables.

### Writing expressions

One output per line:

```
NAME [UNIT] = expression
```

- `NAME`: up to 16 characters (a SELAFIN limit), spaces allowed. `[UNIT]` is optional.
- Lines starting with `#` are comments. Several outputs can share a line separated by `;`.
- Each output can be used by name in the lines after it (spaces become `_`: `SLOPE X` → `SLOPE_X`).
- Two outputs named `… X` and `… Y` (or `… U` and `… V`) are shown as a **vector with arrows**.
- A line without `NAME =` is called `RESULT`.

### Variables

| Write | Means |
| --- | --- |
| `U`, `V`, `H`, `S`, `B`, `F`, `Q`, `M`, `W` … | TELEMAC symbols: velocity U and V, water depth, free surface, bottom, Froude, unit discharge, velocity magnitude, … (English or French files) |
| `WATER_DEPTH`, `VELOCITY_U`, `MY_TRACER` | Any variable by name, spaces as `_` |
| `water_depth`, `velocity_u` | Canonical names, the same in English and French files |
| `x`, `y` | Node coordinates |
| `t` | Time of the frame, in seconds |
| `area` | Area around each node (one third of the triangles touching it) |
| `g`, `pi`, `e` | 9.81, π, e |

### Operators and functions

**Arithmetic**: `+ - * / **` (power) `%` (remainder) and brackets. **Comparisons** `> < >= <= == !=`, `and`,
`or`, `where(condition, a, b)` (the IF of other calculators) and `a if condition else b`.

**Mathematics**: `abs sqrt exp log log10 sin cos tan arcsin arccos arctan arctan2 hypot min max minimum
maximum clip floor ceil sign degrees radians isnan`.

**Spatial calculus on the mesh**:

| Function | Result |
| --- | --- |
| `ddx(f)`, `ddy(f)` | ∂f/∂x, ∂f/∂y |
| `grad(f)` | magnitude of the gradient, e.g. `grad(S)` = water surface slope |
| `div(u, v)` | divergence ∂u/∂x + ∂v/∂y |
| `curl(u, v)` or `vorticity(u, v)` | vertical vorticity ∂v/∂x − ∂u/∂y (1/s) |
| `strain(u, v)` | strain rate √(2u<sub>x</sub>² + 2v<sub>y</sub>² + (u<sub>y</sub> + v<sub>x</sub>)²) |
| `okubo_weiss(u, v)` | s<sub>n</sub>² + s<sub>s</sub>² − ω²; negative inside eddies |
| `laplacian(f)` | ∇²f, as div(grad f) |
| `smooth(f)` | one pass of nodal smoothing |
| `mag(u, v)`, `bearing(u, v)` | magnitude; direction towards which the vector points, degrees from north, clockwise |

**Time**:

| Function | Result |
| --- | --- |
| `tmax(f)`, `tmin(f)`, `tmean(f)`, `tstd(f)` | maximum, minimum, mean, standard deviation over all frames |
| `tsum(f)`, `tint(f)` | sum over frames; time integral (trapezoidal, unit of f × s) |
| `tpeak(f)` | time (s) when f is highest |
| `tfirst(f)`, `tlast(f)` | f at the first or last frame |
| `ddt(f)` | rate of change in time (central differences) |
| `wet(eps)` | 1 where the water depth exceeds eps (default 0.01 m), else 0 |

An output that uses only time reductions (`HMAX = tmax(H)`) is the same in every frame; when the original
variables are not kept it is written once.

### Examples

```
# hydraulics
FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))
SHEAR [PA] = 1000 * g * 0.03**2 * (U**2 + V**2) / maximum(H, 0.01)**(1/3)
HV [M2/S] = H * mag(U, V)
SPECIFIC ENERGY [M] = H + (U**2 + V**2) / (2 * g)

# maps of extremes
HMAX [M] = tmax(H)
TIME HMAX [S] = tpeak(H)
VMAX [M/S] = tmax(mag(U, V))
WET TIME [S] = tint(wet(0.05))
DEPTH CHANGE [M] = tlast(H) - tfirst(H)

# flow structure
VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)
SLOPE X = ddx(S)
SLOPE Y = ddy(S)
RISE [M/S] = ddt(S)
```

`SHEAR` uses Manning's n = 0.03 (τ = ρ g n² |U|² / h<sup>1/3</sup>); change the 0.03 or use the file's
`BOTTOM FRICTION` variable (`N`) instead.

### Accuracy

Spatial derivatives use the linear interpolation on triangles that TELEMAC itself uses: the gradient is
constant in each triangle and averaged to the nodes, weighted by area. It is exact for linear fields,
second-order accurate on regular meshes and first-order on irregular ones. At wet/dry edges the velocity
jumps to zero, so mask derivatives with `where(H > 0.01, …, 0)`. On the boundary, derivatives are
one-sided. The operators are checked in the test suite against exact fields (linear fields, solid-body
rotation, pure strain, quadratic fields) and for convergence.

![Vorticity of the flow in a 180° bend](../img/bend_vorticity.png)

### Options

- **Frames**: empty = all; `-1` = last; `0:10`; `0:100:5` (every 5th).
- **Keep the variables of the file**: off gives a small file with only the new variables.
- *Advanced*: the depth used by `wet()`, and the CRS saved as `.prj` next to the output (by default the
  input's `.prj` is copied).

The calculator works on one file. Differences between two scenarios on the same mesh are planned.

## Integrals and fluxes

*Plugins → TELEMAC Toolkit → Integrals and fluxes*. The output is a table with one row per line or
polygon and per time step, which you can plot, join or export to CSV/Excel.

| Integral | What you give | Columns |
| --- | --- | --- |
| **Discharge through section lines** | a line layer | discharge (m³/s), length, wet length |
| **Flux of an expression** | lines and an expression, e.g. a tracer concentration `T1` | flux (expression × m³/s) |
| **Line integral and mean** | lines and an expression, e.g. `mag(U, V)` | integral, length, mean |
| **Area integral and mean** | polygons (or nothing = whole mesh) and an expression: `H` gives the volume, `1` the area, `wet(0.01)` the wet area | integral, area, mean, wet area |

**Sign of a discharge**: positive from the right of the line to its left, looking along the line as drawn.
Draw sections from the **left bank to the right bank, looking downstream** (as HEC-RAS cross-sections) and
downstream flow is positive.

**Accuracy**: the line is cut at every mesh edge it crosses; between two cuts depth and velocity are
linear and their product is integrated exactly. On a 180° bend run in TELEMAC-2D with 0.05 m³/s in, five
sections give 0.04999 to 0.05001 m³/s ([run summary](../../tests/data/demo_bend/run_summary.txt)). Area
integrals count the elements whose centre is inside the polygon, so a polygon edge cuts the mesh at
element resolution; the table gives the area actually counted.

Lines and polygons in another CRS are transformed to the mesh CRS (from the `.prj`, or the CRS you give).

[← Viewing results](03-results.md) · [Contents](README.md) · [Next: The mesh workflow →](05-mesh-workflow.md)
