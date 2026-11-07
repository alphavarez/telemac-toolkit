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

## Changing the mesh itself

The tools above change values at the nodes. These change the nodes and elements, and carry everything else
over: every variable of every frame (unchanged nodes keep their values exactly, new nodes get values
interpolated in the original elements, the nearest original node outside them), the `.cli` rebuilt for the
new numbering, and the `.prj`.

### Edit by hand in QGIS

1. Mesh workflow, step 5 → **Edit in QGIS (Mesh Digitizing)**. A copy of the mesh (`name_editing.slf`) opens
   in edit mode.
2. Use QGIS's **Mesh Digitizing** toolbar (*View → Toolbars*): add, remove or move vertices, flip edges, split
   faces, force the mesh along lines.
3. **Rebuild after editing** saves the edits and writes `name_edited.slf` and `name_edited.cli`, then adds the
   result to the map.

QGIS saves only the geometry of an edited mesh (the bottom and the other variables are dropped) and knows
nothing of the `.cli`; the rebuild puts both back. In the toolbox the same step is *Mesh → Rebuild a TELEMAC
mesh after editing it in QGIS*, which also accepts a mesh saved by QGIS in another format.

**Boundary conditions after editing**: boundary nodes that existed keep their conditions. A new boundary node
takes the type of the original boundary edge it lies on when both ends of that edge share it (a node added
in the middle of an inflow is an inflow node); otherwise it is a wall. The log reports how many nodes kept
their conditions and how many are new.

### Re-mesh inside polygons

*Mesh → Re-mesh inside polygons* (step 5 → **Re-mesh inside polygons…**) replaces the elements whose centre
lies inside the polygons by a new mesh of the same area: to refine around a bridge pier or a weir, repair a
patch of bad elements, or coarsen an area.

- **Size**: a number, a field of the polygons, or 0 for the median edge length of the area's outline. The
  size grows from the outline's edges to the size inside at the growth rate.
- **Method**: *Balanced* or *Quality-guaranteed* (smallest angle 30°), as in step 1.
- Where the area borders old elements, its outline keeps its nodes, so the new elements join the old ones
  node for node. Where it reaches the mesh boundary, new boundary nodes can appear; they take the boundary
  conditions as above.

### Boundary lines from a .cli

*Mesh → Boundary lines from a .cli* (step 5 → **Boundary lines from the .cli**) writes one line per liquid
boundary of an existing mesh, with the `name`, `type` and `value` fields of a mesh project's boundary layer
(the values from the first row of a `.liq` file when you give one). Edit the lines on the map — move their
ends along the boundary, change a type, add a boundary — and write them back with *Boundary conditions from
lines*. Snap to the mesh boundary while editing (*Mesh boundary as lines* gives it as a layer).

**Checked in TELEMAC-2D**: a river mesh re-meshed at its downstream end (new nodes on the outflow boundary)
and a river mesh edited in QGIS both run in TELEMAC-2D v9.1 with their rebuilt `.cli`, and the prescribed
discharge passes through (60 m³/s in and out).

## Other changes to an existing mesh

| Change | Tool |
| --- | --- |
| Boundary conditions (`.cli`, `.liq`) | *Mesh → Boundary conditions from lines*; the current ones as lines: *Boundary lines from a .cli* |
| Nodes and elements by hand | Step 5 → *Edit in QGIS*, then *Rebuild after editing* |
| Finer, coarser or better elements in an area | *Mesh → Re-mesh inside polygons* |
| Coordinates: shift to (0, 0), rotate | *Mesh → Export or convert mesh*, coordinates options ([details](08-formats-coordinates-boundaries.md#coordinates-in-the-file)) |
| Cut out a part | *Mesh → Extract part of a mesh* |
| Another format | *Mesh → Export or convert mesh* |
| New variables from formulas (e.g. friction from depth classes) | [Calculator](04-calculator-integrals.md) |
| Different mesh (new breakline, new refinement) | Mesh workflow: edit the project layers and generate again |

## For Blue Kenue users

In Blue Kenue the same tasks go through mapping a data set onto the mesh (*Map Object*), the calculator,
and building a new SELAFIN object from the mesh and its variables before saving it. Here the steps are
one tool that writes a ready geometry file, and every step can be repeated or batched from the Processing
toolbox or a model. Blue Kenue's node editing (moving, adding and deleting nodes, swapping edges) is QGIS's
Mesh Digitizing toolbar followed by *Rebuild after editing*, which also keeps the variables and the `.cli`.

**Bringing Blue Kenue work across**:

- *Mesh → Import Blue Kenue lines (.i2s, .i3s)* reads a line set into a line layer, with each line's value,
  and can add the lines straight to a mesh project: as breaklines or guide lines, as the domain or islands
  (closed lines become polygons), or as refinement areas with the line's value as their size. Blue Kenue files
  carry no CRS, so give it.
- Blue Kenue T3 meshes (`.t3s`) are read wherever the plugin reads a mesh: *Compare meshing methods*, *Export
  or convert mesh*, *Rebuild a TELEMAC mesh after editing it in QGIS*. With the same input lines, a Blue Kenue
  mesh and the plugin's methods can be compared measure for measure.

[← Channel templates](06-channel-templates.md) · [Contents](README.md) · [Next: Formats, coordinates and boundary conditions →](08-formats-coordinates-boundaries.md)
