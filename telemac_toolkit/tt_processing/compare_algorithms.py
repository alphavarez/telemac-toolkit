"""Compare meshing methods on one mesh project, and meshes made elsewhere (Blue Kenue, Salomé, ...);
import Blue Kenue line sets."""

from __future__ import annotations

import html
import os
import time

import numpy as np
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsFeature,
    QgsGeometry,
    QgsPoint,
    QgsPointXY,
    QgsFeatureSink,
    QgsFields,
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterNumber,
    QgsProcessingParameterCrs,
    QgsProcessingParameterRasterLayer,
    QgsWkbTypes,
)

from ..mesher import gmsh_engine
from ..mesher.metrics import METRIC_INFO, best_is_low, mesh_metrics
from .algorithms import _Base, _field
from .mesh_algorithms import _DOUBLE, METHOD_KEYS, METHOD_LABELS, _advanced

REPORT_FILTER = "HTML report (*.html)"


def method_choices():
    """(key, label) of every method, built-in first; Gmsh ones only when Gmsh can be imported."""
    out = []
    for k, label in zip(METHOD_KEYS, METHOD_LABELS):
        if k.startswith("gmsh") and not gmsh_engine.available():
            continue
        out.append((k, label))
    return out


class CompareMeshingAlgorithm(_Base):
    def group(self):
        return self.tr("Mesh")

    def groupId(self):
        return "mesh"

    def name(self):
        return "comparemeshing"

    def displayName(self):
        return self.tr("Compare meshing methods")

    def shortHelpString(self):
        return self.tr(
            "Meshes one mesh project with each chosen method and measures every result, and any mesh files you add "
            "(from Blue Kenue, Salomé or an earlier run), on equal terms. The table and the HTML report give, per mesh:\n\n"
            "• Element shape: smallest angle, share of elements below 30° and 20°, largest angle, aspect ratio "
            "(1 = equilateral), edge lengths.\n"
            "• Size changes: area ratio of neighbouring elements (1 = equal).\n"
            "• Outline: how far the mesh boundary strays from the drawn domain.\n"
            "• Terrain (with a DEM): how far TELEMAC's linear bed inside each element is from the DEM (RMS, 95th "
            "percentile, largest) and, with a water level, the error in the water volume under it.\n"
            "• Time step at a Courant number of 1, for the water level (or depth) and velocity given: the largest "
            "stable step TELEMAC's explicit parts allow, set by the smallest elements in deep, fast water.\n"
            "• Nodes, elements, density and meshing time.\n\n"
            "The best value of each row is marked in the report. No single score is given: a finer mesh follows the "
            "DEM better and needs a smaller time step; the report shows that trade-off, the choice stays yours.")

    def initAlgorithm(self, config=None):
        choices = method_choices()
        self._choices = choices
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (empty = only the mesh files below)"),
                                                     extension="gpkg", optional=True))
        default = [i for i, (k, _l) in enumerate(choices) if k in ("balanced", "quality", "gmsh_frontal")]
        self.addParameter(QgsProcessingParameterEnum("METHODS", self.tr("Methods"), [lab for _k, lab in choices],
                                                     allowMultiple=True, defaultValue=default, optional=True))
        self.addParameter(QgsProcessingParameterMultipleLayers("MESHES", self.tr("Mesh files to compare as well"),
                                                               QgsProcessing.SourceType.TypeFile, optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM (bed error, volume, time step)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("LEVEL", self.tr("Water level for the volume and time step (m, optional)"),
                                                       _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterNumber("DEPTH", self.tr("Or a depth for the time step (m)"), _DOUBLE, 2.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("VELOCITY", self.tr("Velocity for the time step (m/s)"), _DOUBLE, 1.0, minValue=0.0))
        self.addParameter(_advanced(QgsProcessingParameterNumber("DEFAULT_SIZE", self.tr("Default element size (0 = the project's)"),
                                                                 _DOUBLE, 0.0, minValue=0.0)))
        self.addParameter(QgsProcessingParameterFeatureSink("TABLE", self.tr("Comparison table"), QgsProcessing.SourceType.TypeVector))
        self.addParameter(QgsProcessingParameterFileDestination("REPORT", self.tr("Report"), REPORT_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFolderDestination("MESH_FOLDER", self.tr("Folder for the generated meshes (SELAFIN)"),
                                                                  optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge import meshing as M

        choices = method_choices()
        idx = self.parameterAsEnums(parameters, "METHODS", context)
        methods = [choices[i] for i in idx if i < len(choices)]
        path = self.parameterAsFile(parameters, "PROJECT", context)
        files = self.parameterAsFileList(parameters, "MESHES", context) if hasattr(self, "parameterAsFileList") else []
        if not files:
            files = [str(v) for v in (parameters.get("MESHES") or [])]
        if not (path and methods) and not files:
            raise QgsProcessingException(self.tr("Give a mesh project and methods, or mesh files."))
        crs = None
        dem_layer = self.parameterAsRasterLayer(parameters, "DEM", context)
        level = parameters.get("LEVEL")
        level = None if level in (None, "") else float(level)
        depth = self.parameterAsDouble(parameters, "DEPTH", context)
        vel = self.parameterAsDouble(parameters, "VELOCITY", context)
        rows = []
        outline = None
        folder = self.parameterAsString(parameters, "MESH_FOLDER", context)
        if path and methods:
            layers = M.open_mesh_project(path)
            st = M.read_settings(path)
            crs = layers["domain"].crs()
            size = self.parameterAsDouble(parameters, "DEFAULT_SIZE", context)
            if not size:
                try:
                    size = float(st.get("default_size") or 0)
                except ValueError:
                    size = 0.0
            size = size or M.auto_default_size(layers, context.project())
            if not size:
                raise QgsProcessingException(self.tr("The project has no default element size."))
            growth = float(st.get("growth") or 1.2)
            smoothing = int(float(st.get("smoothing") or 8))
            try:
                outline_tol = float(st.get("outline_tolerance") or 0) or None
            except ValueError:
                outline_tol = None
            outline = M.domain_rings(layers, crs, context.project())
            dem = M.RasterSampler(dem_layer, 1, crs, project=context.project()).sample if dem_layer is not None else None
            for k, (key, label) in enumerate(methods):
                if feedback.isCanceled():
                    break
                feedback.setProgressText(self.tr("Meshing with %s") % label)
                t0 = time.perf_counter()
                try:
                    build = M.build_mesh_input(layers, crs, size, growth, smoothing, context.project(), max_nodes=None,
                                               min_angle=30.0 if key == "quality" else None, outline_tolerance=outline_tol)
                    res = M.run_mesher(build, None, method=key)
                except (ValueError, RuntimeError) as exc:
                    feedback.reportError(self.tr("%s: %s") % (label, exc), fatalError=False)
                    continue
                secs = time.perf_counter() - t0
                m = mesh_metrics(res.x, res.y, res.ikle, outline, dem, None, level, None if level is not None else depth, vel)
                m["seconds"] = secs
                rows.append((label, m))
                feedback.pushInfo(self.tr("%s: %d nodes, %d elements in %.1f s, smallest angle %.1f°")
                                  % (label, len(res.x), len(res.ikle), secs, m["angle_min"]))
                if folder:
                    os.makedirs(folder, exist_ok=True)
                    from ..telemac_io.selafin import SelafinWriter

                    out = os.path.join(folder, "%s.slf" % key)
                    names, units, vals = [], [], np.zeros((0, len(res.x)))
                    if dem is not None:  # the bed from the DEM, as Generate mesh would give it
                        names, units, vals = ["BOTTOM"], ["M"], np.asarray(dem(res.x, res.y), dtype=np.float64)[None, :]
                    with SelafinWriter(out, names, units, res.x, res.y, res.ikle, title="TELEMAC Toolkit - %s" % label) as w:
                        w.write_frame(0.0, vals)
                feedback.setProgress(100.0 * (k + 1) / (len(methods) + len(files)))
        for j, f in enumerate(files):
            try:
                mesh, fcrs = M.read_mesh_file(f)
            except ValueError as exc:
                feedback.reportError("%s: %s" % (os.path.basename(f), exc), fatalError=False)
                continue
            ik = np.asarray(mesh.all_triangles(), dtype=np.int64)
            dem = None
            if dem_layer is not None:
                dem = M.RasterSampler(dem_layer, 1, crs or fcrs or dem_layer.crs(), project=context.project()).sample
            z = np.asarray(mesh.z, dtype=np.float64) if getattr(mesh, "z", None) is not None and len(mesh.z) == len(mesh.x) else None
            m = mesh_metrics(mesh.x, mesh.y, ik, outline, dem, z, level, None if level is not None else depth, vel)
            rows.append((os.path.basename(f), m))
            feedback.pushInfo(self.tr("%s: %d nodes, %d elements, smallest angle %.1f°") % (os.path.basename(f), len(mesh.x), len(ik), m["angle_min"]))
            feedback.setProgress(100.0 * (len(methods) + j + 1) / (len(methods) + len(files)))
        if not rows:
            raise QgsProcessingException(self.tr("No mesh could be measured."))
        keys = [k for k in METRIC_INFO if any(k in m for _n, m in rows)]
        fields = QgsFields()
        fields.append(_field("mesh", "string"))
        for k in keys:
            fields.append(_field(k, "double"))
        sink, dest = self.parameterAsSink(parameters, "TABLE", context, fields, QgsWkbTypes.Type.NoGeometry)
        for name, m in rows:
            ft = QgsFeature(fields)
            ft.setAttributes([name] + [None if k not in m else float(m[k]) * (100.0 if METRIC_INFO[k][1] == "%" else 1.0) for k in keys])
            sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
        result = {"TABLE": dest}
        report = self.parameterAsFileOutput(parameters, "REPORT", context)
        if report:
            with open(report, "w", encoding="utf-8") as fh:
                fh.write(html_report(rows, keys, os.path.basename(path) if path else "", level, depth, vel))
            result["REPORT"] = report
            feedback.pushInfo(self.tr("Report: %s") % report)
        return result


def _fmt(v, unit):
    if v is None:
        return "–"
    if unit == "%":
        return "%.2f" % (100 * v)
    a = abs(v)
    if a >= 1000:
        return "{:,.0f}".format(v)
    if a >= 10:
        return "%.1f" % v
    if a >= 0.01:
        return "%.3f" % v
    return "%.2g" % v


def html_report(rows, keys, project="", level=None, depth=None, velocity=None) -> str:
    """A self-contained HTML table: one column per mesh, one row per measure, the best value marked."""
    names = [n for n, _m in rows]
    out = ["<!doctype html><html><head><meta charset='utf-8'><title>Meshing comparison</title><style>",
           "body{font-family:system-ui,sans-serif;margin:24px;color:#1b1b1b;background:#fff}"
           "table{border-collapse:collapse;font-size:14px}th,td{border:1px solid #ccc;padding:6px 10px;text-align:right}"
           "th:first-child,td:first-child{text-align:left}thead th{background:#f2f2f2}"
           "td.best{background:#dff3e4;font-weight:600}td.unit{color:#666;text-align:left}"
           "p.note{color:#444;max-width:60em}</style></head><body>",
           "<h1>Meshing comparison</h1>"]
    sub = []
    if project:
        sub.append("project %s" % html.escape(project))
    if level is not None:
        sub.append("water level %g m" % level)
    elif depth is not None:
        sub.append("depth %g m" % depth)
    if velocity is not None:
        sub.append("velocity %g m/s" % velocity)
    if sub:
        out.append("<p>%s</p>" % ", ".join(sub))
    out.append("<table><thead><tr><th>Measure</th><th>Unit</th>%s</tr></thead><tbody>"
               % "".join("<th>%s</th>" % html.escape(n) for n in names))
    for k in keys:
        label, unit, _low = METRIC_INFO[k]
        vals = [m.get(k) for _n, m in rows]
        good = best_is_low(k)
        best = None
        # errors that can be negative (the volume) are judged by their size
        score = [None if v is None else (abs(v) if k == "volume_error" else v) for v in vals]
        avail = [v for v in score if v is not None]
        if good is not None and len(avail) > 1 and not np.allclose(avail, avail[0], rtol=1e-3, atol=1e-12):
            best = min(avail) if good else max(avail)  # rows where every mesh is the same are not marked
        cells = []
        for v, sc in zip(vals, score):
            cls = " class='best'" if best is not None and sc is not None and np.isclose(sc, best, rtol=1e-3) else ""
            cells.append("<td%s>%s</td>" % (cls, _fmt(v, unit)))
        out.append("<tr><td>%s</td><td class='unit'>%s</td>%s</tr>" % (html.escape(label), html.escape(unit), "".join(cells)))
    out.append("</tbody></table>")
    out.append("<p class='note'>Marked: the best value of the row (rows where the meshes do not differ are not marked; "
               "the volume error by its size). Bed error: TELEMAC's linear bed inside each element "
               "against the DEM, at seven points per element. Time step: the largest stable step at a Courant number of 1. "
               "A finer mesh follows the terrain better and needs a smaller time step; the choice stays with the "
               "modeller.</p></body></html>")
    return "\n".join(out)


class BlueKenueLinesAlgorithm(_Base):
    TARGETS = ["Nothing (only the output layer)", "Constraint lines, breakline", "Constraint lines, guide line",
               "Domain (closed lines)", "Islands (closed lines)", "Refinement areas (closed lines, the value as size)"]
    TARGET_KEYS = [None, ("constraint_lines", "breakline"), ("constraint_lines", "guide"), ("domain", None),
                   ("islands", None), ("refinement_areas", None)]

    def group(self):
        return self.tr("Mesh")

    def groupId(self):
        return "mesh"

    def name(self):
        return "bluekenuelines"

    def displayName(self):
        return self.tr("Import Blue Kenue lines (.i2s, .i3s)")

    def shortHelpString(self):
        return self.tr(
            "Reads a Blue Kenue line set (.i2s, or .i3s with a z per vertex) into a line layer with each line's value, "
            "and can add the lines to a layer of a mesh project: constraint lines (breaklines or guide lines), the "
            "domain or islands (closed lines become polygons), or refinement areas with the line value as their size. "
            "With the same input lines a mesh from Blue Kenue and meshes from this plugin can be compared (Compare "
            "meshing methods, which also reads Blue Kenue .t3s meshes). Blue Kenue files carry no CRS: give it.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Blue Kenue line set"),
                                                     fileFilter="Blue Kenue lines (*.i2s *.i3s);;All files (*.*)"))
        self.addParameter(QgsProcessingParameterCrs("CRS", self.tr("CRS of the coordinates"), optional=True))
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project to add the lines to (optional)"),
                                                     extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterEnum("TARGET", self.tr("Add them as"), [self.tr(t) for t in self.TARGETS], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Lines"), QgsProcessing.SourceType.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        from ..telemac_io.bluekenue import read_lines

        path = self.parameterAsFile(parameters, "INPUT", context)
        try:
            lines = read_lines(path)
        except (OSError, ValueError) as exc:
            raise QgsProcessingException(str(exc))
        crs = self.parameterAsCrs(parameters, "CRS", context)
        if not crs.isValid():
            crs = context.project().crs() if context.project() is not None else QgsCoordinateReferenceSystem()
        has_z = any(c.shape[1] > 2 for c, _v in lines)
        fields = QgsFields()
        fields.append(_field("value", "double"))
        fields.append(_field("closed", "int"))
        wkb = QgsWkbTypes.Type.LineStringZ if has_z else QgsWkbTypes.Type.LineString
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, wkb, crs)
        n_closed = 0
        for c, v in lines:
            closed = len(c) > 3 and np.hypot(*(c[0, :2] - c[-1, :2])) <= 1e-9 * max(1.0, float(np.ptp(c[:, :2])))
            n_closed += int(closed)
            ft = QgsFeature(fields)
            if has_z:
                pts = [QgsPoint(float(p[0]), float(p[1]), float(p[2]) if len(p) > 2 else 0.0) for p in c]
                ft.setGeometry(QgsGeometry.fromPolyline(pts))
            else:
                ft.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(float(p[0]), float(p[1])) for p in c]))
            ft.setAttributes([float(v), int(closed)])
            sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
        feedback.pushInfo(self.tr("%d lines (%d closed) from %s") % (len(lines), n_closed, os.path.basename(path)))
        project = self.parameterAsFile(parameters, "PROJECT", context)
        target = self.TARGET_KEYS[self.parameterAsEnum(parameters, "TARGET", context)]
        if project and target:
            from ..qgis_bridge import meshing as M

            layers = M.open_mesh_project(project)
            key, mode = target
            lyr = layers[key]
            if lyr.crs().isValid() and crs.isValid() and lyr.crs() != crs:
                raise QgsProcessingException(self.tr("The project is in %s, the lines in %s: give the lines' CRS or reproject.")
                                             % (lyr.crs().authid(), crs.authid()))
            polygon = key in ("domain", "islands", "refinement_areas")
            lyr.startEditing()
            added = 0
            for c, v in lines:
                f = QgsFeature(lyr.fields())
                xy = [QgsPointXY(float(p[0]), float(p[1])) for p in c]
                if polygon:
                    if len(xy) < 3:
                        continue
                    if xy[0] != xy[-1]:
                        xy.append(xy[0])
                    f.setGeometry(QgsGeometry.fromPolygonXY([xy]))
                    if key == "refinement_areas" and v > 0 and lyr.fields().indexOf("size") >= 0:
                        f.setAttribute("size", float(v))
                else:
                    f.setGeometry(QgsGeometry.fromPolylineXY(xy))
                    if lyr.fields().indexOf("mode") >= 0:
                        f.setAttribute("mode", mode)
                if lyr.fields().indexOf("name") >= 0:
                    f.setAttribute("name", os.path.splitext(os.path.basename(path))[0])
                lyr.addFeature(f)
                added += 1
            if not lyr.commitChanges():
                raise QgsProcessingException("; ".join(lyr.commitErrors()))
            feedback.pushInfo(self.tr("%d added to the project layer %s") % (added, key))
            M.refresh_project_layers(project, context.project())
        return {"OUTPUT": dest}


COMPARE_ALGORITHMS = [CompareMeshingAlgorithm, BlueKenueLinesAlgorithm]
