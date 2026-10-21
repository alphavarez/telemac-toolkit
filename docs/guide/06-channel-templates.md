# 6. Channel and flume templates

Mesh workflow panel, step 1 → **Channel or flume template…** (or *Mesh → Channel or flume template* in the
menu and the toolbox). It draws a laboratory flume or a schematic channel in the current mesh project —
or in a new one in local coordinates — with its bed defined by numbers instead of a DEM.

## When to use it, and when not

Use it for:

- **flume experiments and their numerical replicas** (straight, curved or meandering flumes);
- **checks against analytical solutions**: normal depth, backwater curves, superelevation in bends;
- **sensitivity tests** of mesh size, friction law, turbulence model or time step on a clean geometry;
- **idealised reaches** for teaching, sediment transport or habitat studies.

A **natural river** with real topography is not drawn with the template: draw the domain, the channel
banks and axis and the boundaries on the map over a DEM ([mesh workflow](05-mesh-workflow.md)). The
template's ideas still help there: a curved structured channel follows the banks you draw, and a bed
formula in channel coordinates can carve an idealised channel into a real reach
([editing an existing mesh](07-edit-existing-mesh.md)).

## Planforms

| Planform | Parameters | Typical use |
| --- | --- | --- |
| **Straight** | length, width | flumes, normal depth, backwater |
| **Bend** | straight inlet length, bend radius (on the axis), angle (90°, 180°, …), turning left or right, straight outlet length, width | curved-flume experiments, secondary flow, superelevation |
| **Meanders** (sine-generated) | wavelength along the axis λ, largest deflection angle θ₀, number of wavelengths, straight inlet and outlet, width | meandering channels, bend sequences |

The meander shape is the sine-generated curve of Langbein and Leopold (1966): the heading of the axis is
θ(s) = θ₀ sin(2π s / λ). Its tightest bend has radius λ / (2π θ₀); the tool refuses shapes where that is
less than half the width (the inner bank would cross itself), and the same for a bend whose radius is less
than half the width.

Common parameters: **origin** (centre of the inlet for curved planforms; (0, 0) by default), **direction**
at the inlet (degrees, counter-clockwise from x), **element size**, **structured** mesh on or off, and
optionally the **inflow** and the **outlet water level** (written as the boundary values).

## What it adds to the project

- the **domain** polygon;
- an **inlet** line (type Q) and an **outlet** line (type H), stopping short of the corners so the corners
  stay walls;
- with *Structured*: for a straight flume a structured grid zone; for a bend or meanders the **left bank,
  right bank and axis** in the channels layer, with one vertex per cross-section, so the structured cells
  are normal to the axis even in the bends (narrower on the inner bank, wider on the outer one);
- the **bed formula** and the **default size** in the project settings.

![The 180° bend template: cross-sections normal to the axis](../img/bend_mesh.png)

## The bed

- **Bed level at the inlet** z₀ and **slope** S along the axis (positive downstream).
- **Cross-section**: rectangular (flat bed between vertical walls) or trapezoidal (bottom width, side slope
  as horizontal per 1 vertical, bank height; flat banks above).

### Bed formulas in channel coordinates

For a straight channel the formula is written in x and y. For bends and meanders it is written in the
**channel coordinates**:

- `s`: distance along the axis from its start (m);
- `n`: signed distance from the axis (m), positive on the left looking downstream.

A trapezoidal channel 0.6 m wide at the bottom, side slope 1:1, banks 0.2 m high, slope 0.001:

```
0.30 - 0.001*s + clip((abs(n) - 0.3)/1, 0, 0.2)
```

The formula is in the project settings (step 3 shows it) and can be edited: add a pool-riffle sequence
(`+ 0.02*sin(2*pi*s/5)`), a transverse bed slope in the bends, a scour hole, anything NumPy can write. Any
mesh project with an axis line in its channels layer can use `s` and `n`; when there are several axes the
setting `bottom_axis` names the channel to use.

## Example: the 180° bend

`scripts/make_bend_demo.py` builds a 1 m wide flume with a 3 m inlet, a 180° bend of radius 2 m and a 3 m
outlet, slope 0.001 along the axis, 5 cm cells, Manning n = 0.025 written into the geometry file and
0.05 m³/s. In TELEMAC-2D v9.1.1 the flow is steady after a few minutes; the discharge through five sections
is 0.04999 to 0.05001 m³/s, and the water surface is 6.3 mm higher on the outer bank at the apex,
against 6.0 mm by Rozovskii's estimate U²W/(gR) ([run summary](../../tests/data/demo_bend/run_summary.txt)).

Note that TELEMAC numbers this case's outlet as liquid boundary 1 (it starts from the boundary point with
the smallest x + y); the plugin numbers boundaries the same way, so its `.liq` file and steering-file lines
are right.

[← The mesh workflow](05-mesh-workflow.md) · [Contents](README.md) · [Next: Editing an existing mesh →](07-edit-existing-mesh.md)
