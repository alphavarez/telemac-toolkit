# 7. Editing an existing mesh

*Plugins → TELEMAC Toolkit → Mesh → Edit node values of a mesh* changes a node variable of an existing
SELAFIN geometry — the bottom, the friction or any other — without remeshing. The new file has the
**same nodes, elements, boundary numbering, origin and date**, so the boundary conditions file (`.cli`)
still fits; it is copied next to the new file. The other variables are kept.

## Where the new values come from

In priority order — the first source that covers a node gives its value:

1. **Rasters** (bilinear), in the order you list them;
2. **Survey points** (linear on their triangulation), with a value field or the point Z;
3. a **formula** of `x` and `y`, or of the channel coordinates `s` and `n` along an **axis line** you give;
4. a **constant**, when there is no other source (or, with *Fill nodes no source covers* turned off in the
   advanced options, for the nodes the other sources do not cover; with it on they take the nearest value).

*Unit factor* converts rasters and points in feet (0.3048) or centimetres (0.01).

## What happens to the old values

| Operation | New value | For |
| --- | --- | --- |
| **Replace** | the new value | a new survey, a corrected DEM, a friction map |
| **Add** | old + new value | raise or lower by a thickness (a negative constant lowers) |
| **Raise to** | the higher of old and new | fills, dikes, embankments, road crests |
| **Lower to** | the lower of old and new | dredging, cut channels, excavations |

## Where

- **The whole mesh**, or **only inside polygons**.
- **Transition width**: inside the polygon, the change fades from nothing at the edge to full at that
  distance from it, so a dredged area does not end in a vertical step.

## Recipes

**New survey over a reach.** Rasters = the new DEM; operation *Replace*; area = the polygon of the
survey; transition width 2–3 elements.

**Dredged channel with side slopes.** Draw the channel axis as a line. Formula (5 m bottom width, side
slopes 2:1, bottom 3 m below the formula's reference, slope 0.0005):
`20.0 - 0.0005*s - 3 + clip((abs(n) - 2.5)/2, 0, 3)`; operation *Lower to*; area = a polygon a little wider
than the cut.

**Dike or levee crest.** Constant = crest level; operation *Raise to*; area = a thin polygon along the
crest; transition width = the half-width of the dike base. For a sharp crest the mesh also needs a
breakline on it (mesh workflow, constraint lines).

**Friction zone.** Variable `BOTTOM FRICTION`; constant = the value; operation *Replace*; area = the
zone. If the file has no such variable it is added (0 elsewhere — give the whole mesh a value first, then
the zones).

## Other changes to an existing mesh

| Change | Tool |
| --- | --- |
| Boundary conditions (`.cli`, `.liq`) | *Mesh → Boundary conditions from lines* |
| Coordinates: shift to (0, 0), rotate | *Mesh → Export or convert mesh*, coordinates options ([details](08-formats-coordinates-boundaries.md#coordinates-in-the-file)) |
| Cut out a part | *Mesh → Extract part of a mesh* |
| Another format | *Mesh → Export or convert mesh* |
| New variables from formulas (e.g. friction from depth classes) | [Calculator](04-calculator-integrals.md) |
| Different mesh (new breakline, new refinement) | Mesh workflow: edit the project layers and generate again |

## For Blue Kenue users

In Blue Kenue the same tasks go through mapping a data set onto the mesh (*Map Object*), the calculator,
and building a new SELAFIN object from the mesh and its variables before saving it. Here the steps are
one tool that writes a ready geometry file, and every step can be repeated or batched from the Processing
toolbox or a model.

[← Channel templates](06-channel-templates.md) · [Contents](README.md) · [Next: Formats, coordinates and boundary conditions →](08-formats-coordinates-boundaries.md)
