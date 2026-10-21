"""Processing algorithms (usable in batch mode, the Model Designer and scripts)."""

from __future__ import annotations

import html
import os

import numpy as np
from qgis.core import (
    Qgis,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
)
from qgis.PyQt.QtCore import QCoreApplication, QMetaType

from ..telemac_io import SelafinFile, TriangleLocator, boundary_loops, extract_2d, mesh_statistics, subset_frames
from ..telemac_io.selafin import precision_warning
from ..telemac_io.variables import display_name, vector_pairs

SLF_FILTER = "SELAFIN (*.slf *.res *.srf *.geo *.SLF);;All files (*.*)"


def _field(name, kind):
    """QgsField that works on QGIS 3.44 and 4.x."""
    t = {"int": QMetaType.Type.Int, "double": QMetaType.Type.Double, "string": QMetaType.Type.QString}[kind]
    return QgsField(name, t)


class _Base(QgsProcessingAlgorithm):
    def tr(self, s):
        return QCoreApplication.translate("TelemacToolkit", s)

    def createInstance(self):
        return type(self)()

    def group(self):
        return self.tr("SELAFIN")

    def groupId(self):
        return "selafin"

    def _open(self, parameters, context, key="INPUT"):
        path = self.parameterAsFile(parameters, key, context)
        if not path or not os.path.exists(path):
            raise QgsProcessingException(self.tr("File not found: %s") % path)
        try:
            return SelafinFile(path)
        except Exception as exc:
            raise QgsProcessingException(str(exc))

    @staticmethod
    def _parse_list(text):
        return [t.strip() for t in (text or "").replace(";", ",").split(",") if t.strip()]


# --------------------------------------------------------------------------
class SelafinInfoAlgorithm(_Base):
    def name(self):
        return "selafininfo"

    def displayName(self):
        return self.tr("SELAFIN file information")

    def shortHelpString(self):
        return self.tr("Report the content of a TELEMAC SELAFIN file: dimensions, variables, frames, start date, "
                       "precision and mesh quality. Warns when float32 coordinates lose precision.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Report"), "HTML (*.html)"))

    def processAlgorithm(self, parameters, context, feedback):
        slf = self._open(parameters, context)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        with slf:
            d = slf.info()
            h = slf.header
            st = mesh_statistics(h.x2, h.y2, h.ikle2) if h.ndp in (3, 6) and h.nelem else {}
            warn = precision_warning(h)
            rows = [
                ("Title", d["title"]), ("Dimension", d["dimension"] + (" (%d planes)" % d["planes"] if h.is_3d else "")),
                ("Precision", "%s, %s-endian" % (d["precision"], d["endian"])), ("Nodes", d["nodes"]),
                ("Elements", d["elements"]), ("Nodes per element", d["nodes_per_element"]), ("Frames", d["frames"]),
                ("First / last time (s)", "%s / %s" % (d["first_time_s"], d["last_time_s"])),
                ("Output step (s)", d["time_step_s"]), ("Start date", d["start_date"] or "—"),
                ("Origin (IPARAM 3, 4)", "%g, %g" % d["origin"]), ("Extent", "%.3f, %.3f, %.3f, %.3f" % d["extent"] if "extent" in d else "—"),
                ("Partitioned sub-domain", "yes" if d["partitioned"] else "no"),
            ]
            for k, v in st.items():
                rows.append(("Mesh: " + k.replace("_", " "), "%.4g" % v if isinstance(v, float) else v))
            body = "".join("<tr><th align='left'>%s</th><td>%s</td></tr>" % (html.escape(str(a)), html.escape(str(b))) for a, b in rows)
            vars_ = "".join("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
                i + 1, html.escape(v["name"]), html.escape(v["unit"]), html.escape(display_name(v["name"], v["unit"])))
                for i, v in enumerate(d["variables"]))
            pairs = ", ".join("%s + %s" % (h.var_names[i], h.var_names[j]) for i, j, _ in vector_pairs(h.var_names)) or "—"
            doc = ("<html><head><meta charset='utf-8'><title>SELAFIN information</title></head><body style='font-family:sans-serif'>"
                   "<h2>%s</h2>%s<table>%s</table><h3>Variables</h3><table><tr><th>#</th><th>Name</th><th>Unit</th><th>Label</th></tr>%s</table>"
                   "<p>Vector pairs: %s</p></body></html>") % (
                html.escape(os.path.basename(slf.path)),
                ("<p style='color:#b35900'>⚠ %s</p>" % html.escape(warn)) if warn else "", body, vars_, html.escape(pairs))
            for a, b in rows:
                feedback.pushInfo("%s: %s" % (a, b))
            if warn:
                feedback.reportError(warn, fatalError=False)
        with open(out, "w", encoding="utf-8") as f:
            f.write(doc)
        return {"OUTPUT": out}


# --------------------------------------------------------------------------
class Extract2DAlgorithm(_Base):
    MODES = ["average", "bottom", "surface", "index"]

    def name(self):
        return "extract2d"

    def displayName(self):
        return self.tr("Extract 2D from a 3D result")

    def shortHelpString(self):
        return self.tr("Write a 2D SELAFIN file from a TELEMAC-3D result: the depth average (ELEVATION Z becomes FREE "
                       "SURFACE and the lowest plane becomes BOTTOM), the bottom plane, the surface plane or any plane "
                       "by index (0 = bottom). QGIS can display the 2D file; it cannot read 3D SELAFIN files.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("3D SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterEnum("MODE", self.tr("Extract"), [
            self.tr("Depth average"), self.tr("Bottom plane"), self.tr("Surface plane"), self.tr("Plane by index")], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("PLANE", self.tr("Plane index (0 = bottom)"),
                                                       Qgis.ProcessingNumberParameterType.Integer,
                                                       defaultValue=0))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("2D SELAFIN file"), "SELAFIN (*.slf)"))

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsFile(parameters, "INPUT", context)
        mode = self.MODES[self.parameterAsEnum(parameters, "MODE", context)]
        plane = self.parameterAsInt(parameters, "PLANE", context) if mode == "index" else mode
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)

        def prog(frac):
            feedback.setProgress(100 * frac)
            return feedback.isCanceled()

        try:
            extract_2d(src, out, plane, progress=prog)
        except Exception as exc:
            raise QgsProcessingException(str(exc))
        return {"OUTPUT": out}


# --------------------------------------------------------------------------
class SubsetAlgorithm(_Base):
    def name(self):
        return "subset"

    def displayName(self):
        return self.tr("Extract frames and variables")

    def shortHelpString(self):
        return self.tr("Copy selected frames and variables to a new SELAFIN file. Frames: comma-separated indices "
                       "(0 = first, -1 = last) or ranges such as 0:10 or 0:100:5. Empty = all. Variables: comma-separated "
                       "names. Empty = all.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterString("FRAMES", self.tr("Frames"), defaultValue="-1", optional=True))
        self.addParameter(QgsProcessingParameterString("VARIABLES", self.tr("Variables"), defaultValue="", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Output SELAFIN"), "SELAFIN (*.slf)"))

    @staticmethod
    def parse_frames(text, n):
        text = (text or "").strip()
        if not text:
            return list(range(n))
        out = []
        for tok in text.replace(";", ",").split(","):
            tok = tok.strip()
            if not tok:
                continue
            if ":" in tok:
                parts = [int(p) if p.strip() else None for p in tok.split(":")]
                out.extend(range(n)[slice(*parts)])
            else:
                i = int(tok)
                out.append(i + n if i < 0 else i)
        return out

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsFile(parameters, "INPUT", context)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        with SelafinFile(src) as s:
            n = s.nframes
            try:
                frames = self.parse_frames(self.parameterAsString(parameters, "FRAMES", context), n)
            except ValueError as exc:
                raise QgsProcessingException(self.tr("Could not read the frame list: %s") % exc)
        vars_ = self._parse_list(self.parameterAsString(parameters, "VARIABLES", context)) or None
        try:
            subset_frames(src, out, frames, vars_)
        except (KeyError, IndexError) as exc:
            raise QgsProcessingException(str(exc))
        feedback.pushInfo(self.tr("%d frames written") % len(frames))
        return {"OUTPUT": out}


# --------------------------------------------------------------------------
class ProbeAlgorithm(_Base):
    def name(self):
        return "probetimeseries"

    def displayName(self):
        return self.tr("Time series at points")

    def shortHelpString(self):
        return self.tr("Interpolate variables of a SELAFIN result at every point of a layer, for every frame. "
                       "Output is a table with one row per point and frame (long format), ready for plotting or "
                       "export to CSV. Points outside the mesh are skipped. Point coordinates must be in the same "
                       "CRS as the mesh, or give the mesh CRS so the points are transformed (SELAFIN files carry no CRS).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFeatureSource("POINTS", self.tr("Points"), [Qgis.ProcessingSourceType.VectorPoint]))
        self.addParameter(QgsProcessingParameterString("VARIABLES", self.tr("Variables (empty = all)"), defaultValue="", optional=True))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh coordinates (empty = same as the points)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Time series"), Qgis.ProcessingSourceType.Vector))

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsSource(parameters, "POINTS", context)
        with self._open(parameters, context) as slf:
            h = slf.header
            if h.is_3d:
                raise QgsProcessingException(self.tr("3D file: extract a 2D file first (Extract 2D from a 3D result)."))
            wanted = self._parse_list(self.parameterAsString(parameters, "VARIABLES", context))
            try:
                vidx = [h.var_index(v) for v in wanted] if wanted else list(range(h.nvar))
            except KeyError as exc:
                raise QgsProcessingException(str(exc))
            fields = QgsFields()
            fields.append(_field("point_id", "int"))
            fields.append(_field("frame", "int"))
            fields.append(_field("time_s", "double"))
            if h.has_calendar_date:
                fields.append(_field("datetime", "string"))
            names = []
            for i in vidx:
                nm = h.var_names[i].strip().lower().replace(" ", "_").replace("'", "")[:30]
                names.append(nm)
                fields.append(_field(nm, "double"))
            sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, Qgis.WkbType.NoGeometry)
            loc = TriangleLocator(h.x, h.y, h.ikle)
            mesh_crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
            xform = None
            if mesh_crs.isValid() and src.sourceCrs().isValid() and mesh_crs != src.sourceCrs():
                xform = QgsCoordinateTransform(src.sourceCrs(), mesh_crs, context.transformContext())
            times = slf.times
            dts = slf.datetimes if h.has_calendar_date else None
            total = max(src.featureCount(), 1)
            for k, feat in enumerate(src.getFeatures()):
                if feedback.isCanceled():
                    break
                geom = feat.geometry()
                if geom.isEmpty():
                    continue
                pt = geom.asPoint() if not geom.isMultipart() else geom.asMultiPoint()[0]
                if xform is not None:
                    pt = xform.transform(pt)
                r = loc.locate(pt.x(), pt.y())
                if r is None:
                    feedback.pushWarning(self.tr("Point %s is outside the mesh") % feat.id())
                    continue
                series = [slf.series_weighted(i, r[1], r[2]) for i in vidx]
                for fr in range(slf.nframes):
                    f = QgsFeature(fields)
                    attrs = [int(feat.id()), fr, float(times[fr])]
                    if dts is not None:
                        attrs.append(dts[fr].isoformat())
                    attrs += [float(s[fr]) for s in series]
                    f.setAttributes(attrs)
                    sink.addFeature(f, QgsFeatureSink.Flag.FastInsert)
                feedback.setProgress(100.0 * (k + 1) / total)
        return {"OUTPUT": dest}


# --------------------------------------------------------------------------
class BoundaryAlgorithm(_Base):
    def name(self):
        return "meshboundary"

    def displayName(self):
        return self.tr("Mesh boundary as lines")

    def shortHelpString(self):
        return self.tr("Extract the boundary of a SELAFIN mesh as closed lines: the outer boundary first, then each "
                       "island, with node counts and boundary numbering ranges (IPOBO).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterCrs("CRS", self.tr("CRS of the mesh coordinates"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Boundary"), Qgis.ProcessingSourceType.VectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        crs = self.parameterAsCrs(parameters, "CRS", context)
        with self._open(parameters, context) as slf:
            h = slf.header
            if h.ndp not in (3, 6):
                raise QgsProcessingException(self.tr("The file has no triangular mesh."))
            x, y = h.x2, h.y2
            loops = boundary_loops(h.ikle2, x, y)
            ipobo = h.ipobo[: h.npoin2]
        fields = QgsFields()
        fields.append(_field("loop", "int"))
        fields.append(_field("kind", "string"))
        fields.append(_field("nodes", "int"))
        fields.append(_field("length_m", "double"))
        fields.append(_field("first_node", "int"))
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields,
                                          Qgis.WkbType.LineString, crs)
        for k, lp in enumerate(loops):
            pts = [QgsPointXY(float(x[n]), float(y[n])) for n in lp] + [QgsPointXY(float(x[lp[0]]), float(y[lp[0]]))]
            g = QgsGeometry.fromPolylineXY(pts)
            f = QgsFeature(fields)
            f.setGeometry(g)
            f.setAttributes([k + 1, "outer" if k == 0 else "island", int(len(lp)), float(g.length()), int(lp[0]) + 1])
            sink.addFeature(f, QgsFeatureSink.Flag.FastInsert)
        feedback.pushInfo(self.tr("%d boundary loops, %d boundary nodes") % (len(loops), int((ipobo > 0).sum())))
        return {"OUTPUT": dest}


ALGORITHMS = [SelafinInfoAlgorithm, Extract2DAlgorithm, SubsetAlgorithm, ProbeAlgorithm, BoundaryAlgorithm]
