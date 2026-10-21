"""Mesh writers: TELEMAC first, other programs from the same mesh.

Pure Python + NumPy. Formats that need QGIS (UGRID through MDAL, HEC-RAS
geometry layers) live in ``qgis_bridge.meshing``.

| id      | program(s)                                  | file                    |
| ------- | ------------------------------------------- | ----------------------- |
| telemac | TELEMAC-2D/3D, TOMAWAC, SISYPHE/GAIA        | SELAFIN + .cli + .liq   |
| 2dm     | SMS, BASEMENT, TUFLOW FV, and QGIS           | .2dm (nodestrings, MAT) |
| gmsh    | Gmsh and the solvers that read it            | .msh (format 2.2)       |
| gid     | GiD and GiD-based programs (Iber, ...)       | .msh (GiD ASCII mesh)   |
| mike21  | MIKE 21 FM                                   | .mesh                   |
| adcirc  | ADCIRC, SCHISM                               | fort.14 / hgrid.gr3     |
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional

import numpy as np

from ..telemac_io import write_cli, write_selafin
from ..telemac_io.cli import BoundaryTable, liquid_boundaries
from .boundary import liq_template
from .meshdata import BoundaryString, Frame, Mesh

__all__ = ["FORMATS", "Format", "boundary_strings", "write_mesh", "write_telemac", "write_2dm", "write_gmsh",
           "write_gid", "write_adcirc", "write_mike21", "write_frame_note"]


@dataclass
class Format:
    key: str
    label: str
    ext: str
    quads: bool  # accepts quadrilaterals
    needs_qgis: bool = False
    note: str = ""
    suffix: str = ""  # added to the file name when two formats share an extension


FORMATS: Dict[str, Format] = {f.key: f for f in [
    Format("telemac", "TELEMAC – SELAFIN geometry + .cli + .liq", ".slf", False,
           note="Triangles; bottom and friction as node variables; boundary conditions file in TELEMAC order."),
    Format("2dm", "SMS 2DM – BASEMENT, TUFLOW FV, SMS", ".2dm", True,
           note="Elevation at nodes, material ID per element, one nodestring per open boundary."),
    Format("gmsh", "Gmsh MSH 2.2", ".msh", True,
           note="Physical groups: one per boundary string and one per material."),
    Format("gid", "GiD ASCII mesh – GiD-based programs such as Iber", ".msh", True,
           note="Elevation as z, material number per element. Import in GiD with Files > Import > GiD mesh.",
           suffix="_gid"),
    Format("ugrid", "UGRID NetCDF – Delft3D FM / D-HYDRO", ".nc", True, needs_qgis=True,
           note="Written by QGIS's MDAL library; elevation at nodes."),
    Format("mike21", "MIKE 21 FM mesh", ".mesh", True,
           note="Elevation at nodes, boundary code per node: 1 land, 2 and up open boundaries."),
    Format("adcirc", "ADCIRC / SCHISM grid (fort.14, hgrid.gr3)", ".14", False,
           note="Depth positive downwards; open boundaries and land boundaries listed."),
    Format("hecras", "HEC-RAS – perimeter, breaklines and refinement regions", ".gpkg", True, needs_qgis=True,
           note="HEC-RAS builds its own cells: the drawing, not the mesh, is exported, ready to import in RAS Mapper."),
]}


# ---------------------------------------------------------------------------
# boundary strings
# ---------------------------------------------------------------------------
def boundary_strings(mesh: Mesh, table: BoundaryTable, values: Optional[Dict[str, float]] = None) -> List[BoundaryString]:
    """Runs of boundary nodes with the same codes and name, loop by loop, in boundary order."""
    values = values or {}
    loops = mesh.boundary_loops()
    nodes = table.nodes
    codes = table.codes
    pos = {int(n): i for i, n in enumerate(nodes)}
    numbers = {}
    for b in liquid_boundaries(table):
        numbers[b["start"]] = b["number"]
    out: List[BoundaryString] = []
    for loop in loops:
        idx = [pos[int(n)] for n in loop if int(n) in pos]
        if not idx:
            continue

        def key(i):
            return tuple(int(c) for c in codes[i]), table.comments[i] if i < len(table.comments) else ""

        # start the walk at a change of segment so no run is cut in two at the loop start
        start = 0
        for j in range(len(idx)):
            if key(idx[j]) != key(idx[j - 1]):
                start = j
                break
        order = idx[start:] + idx[:start]
        run = [order[0]]
        for i in order[1:] + [None]:
            if i is not None and key(i) == key(run[-1]):
                run.append(i)
                continue
            c, name = key(run[0])
            num = 0
            for k in run:
                if k in numbers:
                    num = numbers[k]
            seq = list(run)
            out.append(BoundaryString(name or ("wall" if c == (2, 2, 2) else "boundary"), c,
                                      np.array([nodes[k] for k in seq], dtype=np.int64), values.get(name), num))
            if i is not None:
                run = [i]
    return out


def _liquid(strings: List[BoundaryString]) -> List[BoundaryString]:
    return [s for s in strings if not s.is_wall]


# ---------------------------------------------------------------------------
# TELEMAC
# ---------------------------------------------------------------------------
def write_telemac(path: str, mesh: Mesh, table: Optional[BoundaryTable] = None, cli_path: Optional[str] = None,
                  liq_path: Optional[str] = None, liq_values: Optional[Dict[str, float]] = None,
                  store_offset: bool = True, double: bool = False) -> List[str]:
    """SELAFIN geometry (BOTTOM + node variables), plus .cli and .liq when a boundary table is given."""
    files = []
    tri = mesh.all_triangles()
    variables = {"BOTTOM": mesh.bottom()}
    units = {"BOTTOM": "M"}
    for name, (unit, vals) in mesh.node_vars.items():
        variables[name[:16]] = vals
        units[name[:16]] = unit
    x, y, origin = mesh.x, mesh.y, None
    if not mesh.frame.is_identity and store_offset and mesh.frame.angle_deg == 0.0:
        # local coordinates in the file, the offset in IPARAM(3, 4): QGIS shows the mesh in place
        ox, oy = mesh.frame.origin
        x, y, origin = mesh.x + ox, mesh.y + oy, (ox, oy)
    elif mesh.frame.is_identity:
        origin = "auto"
    write_selafin(path, x, y, tri, variables, units=units, title=mesh.title, origin=origin, double=double)
    files.append(path)
    if table is not None and cli_path:
        write_cli(cli_path, table)
        files.append(cli_path)
        if liq_path:
            with open(liq_path, "w", encoding="utf-8", newline="\n") as f:
                f.write(liq_template(table, values=liq_values or {}))
            files.append(liq_path)
    return files


# ---------------------------------------------------------------------------
# SMS 2DM
# ---------------------------------------------------------------------------
def write_2dm(path: str, mesh: Mesh, strings: Optional[List[BoundaryString]] = None) -> List[str]:
    """2DM with E3T/E4Q (counter-clockwise), ND x y z, NS nodestrings and MAT names."""
    z = mesh.bottom()
    mat = mesh.material if mesh.material is not None else np.ones(mesh.nelem, dtype=np.int64)
    nt = len(mesh.triangles)
    lines = ["MESH2D", 'MESHNAME "%s"' % mesh.title.replace('"', "'"), "NUM_MATERIALS_PER_ELEM 1"]
    for k, t in enumerate(mesh.triangles):
        lines.append("E3T %d %d %d %d %d" % (k + 1, t[0] + 1, t[1] + 1, t[2] + 1, mat[k]))
    for k, q in enumerate(mesh.quads):
        lines.append("E4Q %d %d %d %d %d %d" % (nt + k + 1, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1, mat[nt + k]))
    for i in range(mesh.npoin):
        lines.append("ND %d %.10g %.10g %.10g" % (i + 1, mesh.x[i], mesh.y[i], z[i]))
    names_file = []
    for k, s in enumerate(_liquid(strings or [])):
        ids = [int(n) + 1 for n in s.nodes]
        ids[-1] = -ids[-1]
        for j in range(0, len(ids), 10):
            lines.append("NS " + " ".join(str(v) for v in ids[j:j + 10]))
        names_file.append("%d,%s,%d%d%d,%s" % (k + 1, s.name, *s.codes, "" if s.value is None else "%g" % s.value))
    for mid, name in sorted((mesh.material_names or {}).items()):
        if mid > 0:
            lines.append('MAT %d "%s"' % (mid, name.replace('"', "'")))
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    files = [path]
    if names_file:
        side = os.path.splitext(path)[0] + "_nodestrings.csv"
        with open(side, "w", encoding="utf-8", newline="\n") as f:
            f.write("nodestring,name,telemac_codes,value\n" + "\n".join(names_file) + "\n")
        files.append(side)
    return files


# ---------------------------------------------------------------------------
# Gmsh
# ---------------------------------------------------------------------------
def write_gmsh(path: str, mesh: Mesh, strings: Optional[List[BoundaryString]] = None) -> List[str]:
    """Gmsh MSH 2.2 ASCII: nodes with z, triangles/quads tagged by material, boundary lines by string."""
    z = mesh.bottom()
    mat = mesh.material if mesh.material is not None else np.ones(mesh.nelem, dtype=np.int64)
    strings = strings or []
    phys = []  # (dim, tag, name)
    line_tag = {}
    for k, s in enumerate(strings):
        key = s.name or "boundary"
        if key not in line_tag:
            line_tag[key] = 1000 + len(line_tag) + 1
            phys.append((1, line_tag[key], key))
    for mid in sorted(set(int(m) for m in mat)):
        phys.append((2, mid, (mesh.material_names or {}).get(mid, "material %d" % mid)))
    out = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$PhysicalNames", str(len(phys))]
    out += ['%d %d "%s"' % (d, t, n.replace('"', "'")) for d, t, n in phys]
    out += ["$EndPhysicalNames", "$Nodes", str(mesh.npoin)]
    out += ["%d %.10g %.10g %.10g" % (i + 1, mesh.x[i], mesh.y[i], z[i]) for i in range(mesh.npoin)]
    out += ["$EndNodes", "$Elements"]
    elems = []
    eid = 0
    for k, s in enumerate(strings):
        tag = line_tag[s.name or "boundary"]
        n = s.nodes
        for a, b in zip(n[:-1], n[1:]):
            eid += 1
            elems.append("%d 1 2 %d %d %d %d" % (eid, tag, k + 1, a + 1, b + 1))
    nt = len(mesh.triangles)
    for k, t in enumerate(mesh.triangles):
        eid += 1
        elems.append("%d 2 2 %d %d %d %d %d" % (eid, mat[k], mat[k], t[0] + 1, t[1] + 1, t[2] + 1))
    for k, q in enumerate(mesh.quads):
        eid += 1
        m = mat[nt + k]
        elems.append("%d 3 2 %d %d %d %d %d %d" % (eid, m, m, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1))
    out += [str(len(elems))] + elems + ["$EndElements"]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    return [path]


# ---------------------------------------------------------------------------
# GiD
# ---------------------------------------------------------------------------
def write_gid(path: str, mesh: Mesh) -> List[str]:
    """GiD ASCII mesh: one block per element type, coordinates in the first block."""
    z = mesh.bottom()
    mat = mesh.material if mesh.material is not None else np.ones(mesh.nelem, dtype=np.int64)
    out = ["mesh dimension 3 elemtype triangle nnode 3", "coordinates"]
    out += ["%d %.10g %.10g %.10g" % (i + 1, mesh.x[i], mesh.y[i], z[i]) for i in range(mesh.npoin)]
    out += ["end coordinates", "elements"]
    nt = len(mesh.triangles)
    out += ["%d %d %d %d %d" % (k + 1, t[0] + 1, t[1] + 1, t[2] + 1, mat[k]) for k, t in enumerate(mesh.triangles)]
    out += ["end elements"]
    if mesh.has_quads:
        out += ["mesh dimension 3 elemtype quadrilateral nnode 4", "coordinates", "end coordinates", "elements"]
        out += ["%d %d %d %d %d %d" % (nt + k + 1, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1, mat[nt + k])
                for k, q in enumerate(mesh.quads)]
        out += ["end elements"]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    return [path]


# ---------------------------------------------------------------------------
# ADCIRC / SCHISM
# ---------------------------------------------------------------------------
def write_adcirc(path: str, mesh: Mesh, strings: Optional[List[BoundaryString]] = None) -> List[str]:
    """fort.14 / hgrid.gr3: triangles, depth = -bottom, open and land boundary node lists."""
    tri = mesh.all_triangles()
    depth = -mesh.bottom()
    out = [mesh.title[:80], "%d %d" % (len(tri), mesh.npoin)]
    out += ["%d %.10f %.10f %.6f" % (i + 1, mesh.x[i], mesh.y[i], depth[i]) for i in range(mesh.npoin)]
    out += ["%d 3 %d %d %d" % (k + 1, t[0] + 1, t[1] + 1, t[2] + 1) for k, t in enumerate(tri)]
    strings = strings or []
    open_b = [s for s in strings if not s.is_wall]
    land = [s for s in strings if s.is_wall]
    loops = mesh.boundary_loops()
    outer = set(int(n) for n in loops[0]) if loops else set()
    out.append("%d = Number of open boundaries" % len(open_b))
    out.append("%d = Total number of open boundary nodes" % sum(len(s.nodes) for s in open_b))
    for k, s in enumerate(open_b):
        out.append("%d = Number of nodes for open boundary %d (%s)" % (len(s.nodes), k + 1, s.name))
        out += [str(int(n) + 1) for n in s.nodes]
    out.append("%d = Number of land boundaries" % len(land))
    out.append("%d = Total number of land boundary nodes" % sum(len(s.nodes) for s in land))
    for k, s in enumerate(land):
        flag = 0 if int(s.nodes[0]) in outer else 1  # 0 mainland, 1 island
        out.append("%d %d = Number of nodes for land boundary %d" % (len(s.nodes), flag, k + 1))
        out += [str(int(n) + 1) for n in s.nodes]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    return [path]


# ---------------------------------------------------------------------------
# MIKE 21 FM
# ---------------------------------------------------------------------------
def write_mike21(path: str, mesh: Mesh, strings: Optional[List[BoundaryString]] = None, projection: str = "") -> List[str]:
    """MIKE .mesh: header (item code, nodes, projection), nodes with boundary codes, elements."""
    z = mesh.bottom()
    code = np.zeros(mesh.npoin, dtype=np.int64)
    loops = mesh.boundary_loops()
    for lp in loops:
        code[lp] = 1  # land by default
    k = 1
    for s in strings or []:
        if not s.is_wall:
            k += 1
            code[s.nodes] = k
    proj = projection.strip() or "NON-UTM"
    out = ["100079  1000  %d  %s" % (mesh.npoin, proj)]  # item type (bathymetry), unit (metre), nodes, projection
    out += ["%d %.10g %.10g %.10g %d" % (i + 1, mesh.x[i], mesh.y[i], z[i], code[i]) for i in range(mesh.npoin)]
    if mesh.has_quads:
        out.append("%d 4 25" % mesh.nelem)
        out += ["%d %d %d %d 0" % (e + 1, t[0] + 1, t[1] + 1, t[2] + 1) for e, t in enumerate(mesh.triangles)]
        nt = len(mesh.triangles)
        out += ["%d %d %d %d %d" % (nt + e + 1, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1) for e, q in enumerate(mesh.quads)]
    else:
        out.append("%d 3 21" % len(mesh.triangles))
        out += ["%d %d %d %d" % (e + 1, t[0] + 1, t[1] + 1, t[2] + 1) for e, t in enumerate(mesh.triangles)]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    return [path]


def write_frame_note(path: str, frame: Frame, crs_name: str = "") -> Optional[str]:
    """A text note next to a mesh written in a local frame, to go back to map coordinates."""
    if frame.is_identity:
        return None
    note = os.path.splitext(path)[0] + "_frame.txt"
    a = frame.angle_deg
    with open(note, "w", encoding="utf-8", newline="\n") as f:
        f.write("Local coordinates written by TELEMAC Toolkit\n")
        f.write("Map CRS of the drawing: %s\n" % (crs_name or "unknown"))
        f.write("Origin X0 = %.6f\nOrigin Y0 = %.6f\nRotation = %.6f degrees (counter-clockwise)\n" % (frame.origin[0], frame.origin[1], a))
        f.write("x_local = (X - X0) cos(r) + (Y - Y0) sin(r)\ny_local = -(X - X0) sin(r) + (Y - Y0) cos(r)\n")
        f.write("X = X0 + x_local cos(r) - y_local sin(r)\nY = Y0 + x_local sin(r) + y_local cos(r)\n")
    return note


WRITERS: Dict[str, Callable] = {"2dm": write_2dm, "gmsh": write_gmsh, "adcirc": write_adcirc}


def write_mesh(fmt: str, path: str, mesh: Mesh, table: Optional[BoundaryTable] = None,
               values: Optional[Dict[str, float]] = None, **kw) -> List[str]:
    """Write ``mesh`` in one of the pure-Python formats; returns the files written."""
    if fmt not in FORMATS or FORMATS[fmt].needs_qgis:
        raise ValueError("format %r needs QGIS or is unknown" % fmt)
    if not FORMATS[fmt].quads and mesh.has_quads:
        mesh = Mesh(mesh.x, mesh.y, mesh.all_triangles(), z=mesh.z, node_vars=mesh.node_vars,
                    material=None if mesh.material is None else np.concatenate(
                        [mesh.material[:len(mesh.triangles)], np.tile(mesh.material[len(mesh.triangles):], 2)]),
                    material_names=mesh.material_names, crs_wkt=mesh.crs_wkt, frame=mesh.frame, title=mesh.title)
    strings = boundary_strings(mesh, table, values) if table is not None else []
    if fmt == "telemac":
        files = write_telemac(path, mesh, table, kw.get("cli_path"), kw.get("liq_path"), values,
                              kw.get("store_offset", True), kw.get("double", False))
    elif fmt == "gid":
        files = write_gid(path, mesh)
    elif fmt == "mike21":
        files = write_mike21(path, mesh, strings, kw.get("projection", ""))
    else:
        files = WRITERS[fmt](path, mesh, strings)
    return files
