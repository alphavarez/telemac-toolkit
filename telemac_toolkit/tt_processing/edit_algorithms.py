"""Processing tools for editing existing TELEMAC meshes: rebuild after editing in QGIS, re-mesh inside
polygons, boundary lines from a .cli."""

from __future__ import annotations

import os
import shutil

import numpy as np
from qgis.core import (
    Qgis,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingException,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsWkbTypes,
)

from ..mesher.backends import default_triangulator
from ..mesher.boundary import attach_mesh
from ..mesher.edit import match_nodes, rebuild_boundary, region_elements, remesh_region, transfer
from ..telemac_io import SelafinFile, read_cli, write_cli
from ..telemac_io.cli import liquid_boundaries
from ..telemac_io.selafin import SelafinWriter
from .algorithms import SLF_FILTER, _Base, _field
from .mesh_algorithms import CLI_FILTER, LIQ_FILTER, _advanced

_DOUBLE = Qgis.ProcessingNumberParameterType.Double
MESH_FILTER = "Meshes (*.slf *.2dm *.nc *.ply *.mesh *.mesh);;All files (*.*)"
KIND_NAMES = {(4, 5, 5): "Q", (5, 4, 4): "H", (5, 5, 5): "QH", (4, 6, 6): "UV", (5, 6, 6): "UVH", (4, 4, 4): "free", (1, 1, 1): "waves"}


def _read_slf(path):
    if not path or not os.path.exists(path):
        raise QgsProcessingException("File not found: %s" % path)
    with SelafinFile(path) as s:
        h = s.header
        if h.is_3d or h.ndp != 3:
            raise QgsProcessingException("A 2D triangular SELAFIN file is required.")
        frames = np.array([s.read_frame(k) for k in range(s.nframes)]) if s.nframes else np.zeros((0, h.nvar, h.npoin))
        return {"x": h.x.copy(), "y": h.y.copy(), "ikle": h.ikle.copy(), "names": [n.strip() for n in h.var_names],
                "units": [u.strip() for u in h.var_units], "times": np.asarray(s.times, dtype=np.float64), "frames": frames,
                "title": h.title.strip() if isinstance(h.title, str) else "", "date": h.date}


def _write_slf(path, src, x, y, ikle, frames):
    with SelafinWriter(path, src["names"], src["units"], x, y, ikle, title=src["title"] or "TELEMAC Toolkit",
                       date=src["date"]) as w:
        for k in range(frames.shape[0]):
            w.write_frame(float(src["times"][k]), [frames[k, i] for i in range(frames.shape[1])])


def _copy_prj(src, dst):
    for cand in (os.path.splitext(src)[0] + ".prj", src + ".prj"):
        if os.path.exists(cand):
            shutil.copyfile(cand, os.path.splitext(dst)[0] + ".prj")
            return


class _EditBase(_Base):
    def group(self):
        return self.tr("Mesh")

    def groupId(self):
        return "mesh"

    def _carry(self, feedback, src, nx, ny, nik, cli_in, out, match):
        """Values of every variable and frame, and the .cli, onto the new mesh; writes both."""
        nf, nv = src["frames"].shape[0], src["frames"].shape[1]
        flat = src["frames"].reshape(nf * nv, -1) if nf * nv else np.zeros((0, len(src["x"])))
        vals, copied, interp, near = transfer(src["x"], src["y"], src["ikle"], flat, nx, ny, match)
        frames = vals.reshape(nf, nv, -1) if nf * nv else np.zeros((0, 0, len(nx)))
        feedback.pushInfo(self.tr("Node values: %d copied from unchanged nodes, %d interpolated in the original elements, "
                                  "%d from the nearest original node") % (copied, interp, near))
        _write_slf(out, src, nx, ny, nik, frames)
        res = {"OUTPUT": out}
        if cli_in and os.path.exists(cli_in):
            table = read_cli(cli_in)
            attach_mesh(table, src["x"], src["y"], src["ikle"])
            t2, rep = rebuild_boundary(table, src["x"], src["y"], nx, ny, nik, match)
            cli_out = os.path.splitext(out)[0] + ".cli"
            write_cli(cli_out, t2)
            for line in rep:
                feedback.pushInfo(line)
            for b in liquid_boundaries(t2):
                feedback.pushInfo(self.tr("Liquid boundary %d: %s, %d nodes") % (b["number"], t2.comments[b["start"]] or "", len(b["nodes"])))
            res["OUTPUT_CLI"] = cli_out
        else:
            feedback.pushInfo(self.tr("No .cli next to the original mesh: only the mesh was written."))
        return res


# --------------------------------------------------------------------------
class RebuildEditedAlgorithm(_EditBase):
    def name(self):
        return "rebuildedited"

    def displayName(self):
        return self.tr("Rebuild a TELEMAC mesh after editing it in QGIS")

    def shortHelpString(self):
        return self.tr(
            "QGIS edits mesh layers with its Mesh Digitizing tools (add, remove and move vertices, flip edges, split "
            "faces, force the mesh along lines), but the saved file keeps only the geometry: the bottom, the friction "
            "and the other variables are lost, and the .cli no longer matches the nodes. This tool gives them back: "
            "nodes that did not move keep their values exactly, new or moved nodes get values interpolated in the "
            "original elements (nearest original node outside them), and the .cli is rebuilt in TELEMAC order. "
            "Boundary nodes that existed keep their conditions; new ones take the type of the original boundary edge "
            "they lie on when both of its ends share it, otherwise they are walls.\n\n"
            "Edited mesh: the file saved by QGIS (SELAFIN or any mesh format QGIS reads). Original: the TELEMAC mesh "
            "before editing, with its .cli next to it.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("EDITED", self.tr("Edited mesh (saved by QGIS)"), fileFilter=MESH_FILTER))
        self.addParameter(QgsProcessingParameterFile("ORIGINAL", self.tr("Original TELEMAC mesh (SELAFIN, .cli next to it)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFile("CLI", self.tr("Original boundary conditions (empty = .cli next to the original)"),
                                                     fileFilter=CLI_FILTER, optional=True))
        self.addParameter(_advanced(QgsProcessingParameterNumber("TOLERANCE", self.tr("Distance under which a node has not moved (m)"),
                                                                 _DOUBLE, 1e-4, minValue=0.0)))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Rebuilt mesh"), SLF_FILTER))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge import meshing as M

        orig = self.parameterAsFile(parameters, "ORIGINAL", context)
        src = _read_slf(orig)
        edited = self.parameterAsFile(parameters, "EDITED", context)
        try:
            mesh, _crs = M.read_mesh_file(edited)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        ik = np.asarray(mesh.all_triangles(), dtype=np.int64)
        nx, ny = np.asarray(mesh.x, float), np.asarray(mesh.y, float)
        used = np.zeros(len(nx), dtype=bool)
        used[ik.ravel()] = True
        if not used.all():
            remap = -np.ones(len(nx), dtype=np.int64)
            remap[used] = np.arange(int(used.sum()))
            nx, ny, ik = nx[used], ny[used], remap[ik]
            feedback.pushInfo(self.tr("%d nodes not used by any element dropped") % int((~used).sum()))
        tol = self.parameterAsDouble(parameters, "TOLERANCE", context) or 1e-6
        match = match_nodes(src["x"], src["y"], nx, ny, tol)
        feedback.pushInfo(self.tr("Edited mesh: %d nodes, %d elements (original: %d nodes, %d elements)")
                          % (len(nx), len(ik), len(src["x"]), len(src["ikle"])))
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        cli = self.parameterAsFile(parameters, "CLI", context) or os.path.splitext(orig)[0] + ".cli"
        res = self._carry(feedback, src, nx, ny, ik, cli, out, match)
        _copy_prj(orig, out)
        return res


# --------------------------------------------------------------------------
class RemeshRegionAlgorithm(_EditBase):
    METHODS = ["Balanced", "Quality-guaranteed (smallest angle 30 degrees)"]

    def name(self):
        return "remeshregion"

    def displayName(self):
        return self.tr("Re-mesh inside polygons")

    def shortHelpString(self):
        return self.tr(
            "Replace the elements of a TELEMAC mesh whose centre lies inside the polygons by a new mesh of the same "
            "area: to refine around a structure, fix a patch of bad elements, or coarsen an area. Where the area "
            "borders old elements its outline keeps its nodes, so the new elements join the old ones node for node; "
            "the size grows from the outline's edges to the size inside. Where it reaches the mesh boundary, nodes can "
            "be added there and take its boundary conditions.\n\n"
            "Every variable of every frame is carried over: unchanged nodes keep their values, new nodes get values "
            "interpolated in the original elements. The .cli next to the mesh is rebuilt for the new numbering. "
            "Size: a number, a field of the polygons, or the median edge of the area's outline.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("TELEMAC mesh (SELAFIN, .cli next to it)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFeatureSource("POLYGONS", self.tr("Areas to re-mesh"), [Qgis.ProcessingSourceType.VectorPolygon]))
        self.addParameter(QgsProcessingParameterNumber("SIZE", self.tr("Element size inside (m, 0 = the outline's edges)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterField("SIZE_FIELD", self.tr("Or a size field of the polygons"), parentLayerParameterName="POLYGONS",
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterEnum("METHOD", self.tr("Method"), [self.tr(m) for m in self.METHODS], defaultValue=0))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Re-meshed mesh"), SLF_FILTER))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge import meshing as M

        path = self.parameterAsFile(parameters, "INPUT", context)
        src = _read_slf(path)
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid():
            from ..qgis_bridge.crs import sidecar_crs

            crs = sidecar_crs(path) or (context.project().crs() if context.project() else crs)
        layer = self.parameterAsSource(parameters, "POLYGONS", context)
        sfield = self.parameterAsString(parameters, "SIZE_FIELD", context)
        size0 = self.parameterAsDouble(parameters, "SIZE", context)
        min_angle = 30.0 if self.parameterAsEnum(parameters, "METHOD", context) == 1 else None
        x, y, ik = src["x"], src["y"], src["ikle"]
        tri = default_triangulator()
        n_areas = 0
        for f, g in M._features(layer, crs, context.project()):
            for part in M.polygon_parts(g):
                inside = region_elements(x, y, ik, part)
                if not inside.any():
                    feedback.pushWarning(self.tr("Polygon %d holds no element centre; skipped") % f.id())
                    continue
                s = size0
                if sfield:
                    try:
                        s = float(f[sfield]) or size0
                    except (TypeError, ValueError):
                        pass
                try:
                    nx, ny, nik, _old, rep = remesh_region(x, y, ik, inside, s or None, tri, min_angle=min_angle)
                except ValueError as exc:
                    raise QgsProcessingException(str(exc))
                for line in rep:
                    feedback.pushInfo(line)
                x, y, ik = nx, ny, nik
                n_areas += 1
        if not n_areas:
            raise QgsProcessingException(self.tr("No polygon covers an element of the mesh (check the CRS)."))
        match = match_nodes(src["x"], src["y"], x, y, max(1e-9 * max(np.ptp(x), np.ptp(y)), 1e-9))
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        res = self._carry(feedback, src, x, y, ik, os.path.splitext(path)[0] + ".cli", out, match)
        feedback.pushInfo(self.tr("Mesh: %d nodes, %d elements (before: %d, %d)") % (len(x), len(ik), len(src["x"]), len(src["ikle"])))
        _copy_prj(path, out)
        return res


# --------------------------------------------------------------------------
class CliLinesAlgorithm(_EditBase):
    def name(self):
        return "clilines"

    def displayName(self):
        return self.tr("Boundary lines from a .cli")

    def shortHelpString(self):
        return self.tr(
            "One line per liquid boundary of a TELEMAC mesh, with the name, type and value fields of a mesh project's "
            "boundary layer: the boundary conditions of an existing mesh as lines you can edit on the map (move their "
            "ends along the boundary, change the type, add a boundary), then write back with Boundary conditions from "
            "lines. Snap to the mesh boundary (Mesh boundary as lines) while editing. Values come from the first row "
            "of the liquid boundaries file when given.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("MESH", self.tr("TELEMAC mesh (SELAFIN)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFile("CLI", self.tr("Boundary conditions (empty = .cli next to the mesh)"), fileFilter=CLI_FILTER,
                                                     optional=True))
        self.addParameter(QgsProcessingParameterFile("LIQ", self.tr("Liquid boundaries file (for the values)"), fileFilter=LIQ_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Boundary lines"), Qgis.ProcessingSourceType.VectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        from ..telemac_io.liq import read_liq
        from .bc_algorithms import load_mesh

        x, y, _ik, _z, _fr, table, _cli = load_mesh(self.parameterAsFile(parameters, "MESH", context),
                                                    self.parameterAsFile(parameters, "CLI", context) or None)
        liq = {}
        lp = self.parameterAsFile(parameters, "LIQ", context)
        if lp:
            liq = read_liq(lp)
        fields = QgsFields()
        for n, k in (("name", "string"), ("type", "string"), ("value", "double"), ("number", "int"), ("tracer", "string")):
            fields.append(_field(n, k))
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid() and context.project() is not None:
            crs = context.project().crs()
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, QgsWkbTypes.Type.LineString, crs)
        count = 0
        for b in liquid_boundaries(table):
            codes = tuple(int(c) for c in b["codes"])
            kind = KIND_NAMES.get(codes, "%d%d%d" % codes)
            nodes = np.asarray(b["nodes"])
            if nodes.size < 2:
                continue
            col = ("Q(%d)" if codes[1] == 5 else "SL(%d)") % b["number"]
            val = float(liq[col][0]) if col in liq else None
            ft = QgsFeature(fields)
            ft.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(float(x[n]), float(y[n])) for n in nodes]))
            lit = int(table.records["litbor"][b["start"]])
            ft.setAttributes([table.comments[b["start"]] or "boundary %d" % b["number"], kind, val, b["number"],
                              {2: "wall", 5: "prescribed", 4: "free"}.get(lit)])
            sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
            count += 1
        feedback.pushInfo(self.tr("%d liquid boundaries as lines; every other boundary node is a wall") % count)
        return {"OUTPUT": dest}


EDIT_ALGORITHMS = [RebuildEditedAlgorithm, RemeshRegionAlgorithm, CliLinesAlgorithm]
