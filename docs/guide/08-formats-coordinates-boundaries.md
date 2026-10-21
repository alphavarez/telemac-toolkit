# 8. Formats, coordinates and boundary conditions

## Formats

The mesher is built for TELEMAC; the same mesh can be written for other programs from step 4 of the mesh
workflow or with *Mesh → Export or convert mesh* (which also reads SELAFIN, 2DM, UGRID, MIKE and any mesh
QGIS opens).

| Format | For | Boundaries | Roughness | Quads |
| --- | --- | --- | --- | --- |
| SELAFIN `.slf` + `.cli` + `.liq` | TELEMAC-2D/3D and the TELEMAC suite | `.cli` in TELEMAC order, `.liq` template | `BOTTOM FRICTION` node variable | no |
| SMS 2DM `.2dm` | BASEMENT, TUFLOW FV, SMS, QGIS | nodestrings (+ `_nodestrings.csv` with their names) | `MAT` ids | yes |
| Gmsh `.msh` (2.2 ASCII) | Gmsh and solvers that read it | physical lines per boundary | physical surfaces per material | yes |
| GiD mesh `_gid.msh` | GiD and GiD-based programs such as Iber | — | material per element | yes |
| UGRID `.nc` | Delft3D FM / D-HYDRO | — | — | yes |
| MIKE 21 FM `.mesh` | MIKE 21 FM | boundary code per node | — | yes |
| ADCIRC `fort.14` (`.14`) | ADCIRC, SCHISM (`hgrid.gr3`) | open and land boundary node lists | — | no |
| HEC-RAS `.gpkg` | HEC-RAS (RAS Mapper) | boundary condition lines | — | — |

- **HEC-RAS** builds its own cells from a perimeter, breaklines and cell spacings, so it gets the drawing
  (perimeter, breaklines, refinement regions, boundary lines in a GeoPackage for RAS Mapper), not the
  triangles. Use *Generate mesh* for it, not *Export*.
- **Quadrilaterals**: structured channels and grids can stay quadrilaterals in formats that accept them
  (*Keep structured cells as quadrilaterals*).
- **Elevation**: written as node elevation (ADCIRC: depth = −elevation, its convention).
- **Checked**: SELAFIN, 2DM, UGRID and MIKE files are read back by QGIS in the tests, Gmsh files by
  meshio, and TELEMAC files run in TELEMAC-2D. GiD and ADCIRC files follow their published layouts but have
  not been loaded in Iber or ADCIRC yet: check the first one you make.

## Coordinates in the file

| Option | The file holds |
| --- | --- |
| **Map coordinates as drawn** | the coordinates of the project CRS |
| **Lower-left corner at (0, 0)** | coordinates shifted so the lower-left corner of the mesh is the origin |
| **A chosen point at (0, 0)** | coordinates shifted so that point is the origin |
| **Rotation** | axes rotated counter-clockwise by the angle (e.g. x along a flume drawn at an angle) |

For TELEMAC, a shift without rotation can be **kept in the file header** (IPARAM 3 and 4): TELEMAC computes
near (0, 0), with full single-precision accuracy, and QGIS still draws the mesh in place. Otherwise a
`<name>_frame.txt` note next to the file gives the shift and rotation to go back to map coordinates, and
the mesh opens in local coordinates.

Large coordinates (UTM, national grids) in single-precision SELAFIN files are stored with an integer origin
offset automatically, so float32 keeps sub-centimetre resolution.

### Local coordinates

The plugin's local grid ("Local metric grid (laboratory)") is a plane with (0, 0) at the origin and
lengths in metres. It behaves as a normal projected CRS in QGIS — measuring, snapping, labels, every
format — but it is not a place on Earth: it cannot be shown together with map layers in a real CRS. QGIS's
status bar calls it "Unknown CRS". A `.prj` file with its definition is written next to each mesh.

## Boundary conditions

*Generate mesh* (from the boundaries layer) and *Boundary conditions from lines* (for an existing mesh)
write the TELEMAC `.cli`:

- Every boundary node is a **wall** (2 2 2) except the nodes near a boundary line, which take its type:
  `Q` (4 5 5), `H` (5 4 4), `QH` (5 5 5), `UV` (4 6 6), `UVH` (5 6 6), `free` (4 4 4), `waves` (1 1 1),
  or a three-digit code.
- "Near" = within half the boundary edge length (or the snapping distance you give).
- Nodes are listed in TELEMAC order: outer boundary counter-clockwise, then islands.

### Liquid boundary numbers

TELEMAC numbers liquid boundaries in its own way: it starts from the **boundary point with the smallest
x + y** (the lowest-left corner) and walks counter-clockwise, numbering liquid boundaries as it meets
them; two adjacent liquid stretches with different types are two boundaries. The numbers matter: the
steering file (`PRESCRIBED FLOWRATES = 0.05;0.`) and the `.liq` columns (`Q(1)`, `SL(2)`) refer to them.

The plugin numbers boundaries exactly as TELEMAC does and reports them in the log ("Liquid boundary 1:
outlet, codes 5 4 4…"), with the steering-file lines to paste. The `.liq` template it writes has the
right columns. In a flume running along +x from (0, 0), the **outlet** is boundary 1, because the walk
starts at (0, 0) and runs along the right bank to the outlet first.

[← Editing an existing mesh](07-edit-existing-mesh.md) · [Contents](README.md) · [Next: Processing tools reference →](09-processing-reference.md)
