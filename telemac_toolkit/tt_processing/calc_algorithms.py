"""Processing algorithms of the calculator: new variables from formulas, integrals and fluxes."""

from __future__ import annotations

import os
import shutil

import numpy as np
from qgis.core import (
    Qgis,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
)

from ..telemac_io import SelafinFile
from ..telemac_io.calculator import (
    FUNCTIONS_HELP, Calculator, CalculatorError, Output, describe_variables, run_calculator)
from ..telemac_io.integrals import Section, element_mask_in_rings
from ..telemac_io.variables import canonical
from .algorithms import SLF_FILTER, SubsetAlgorithm, _Base, _field

_DOUBLE = Qgis.ProcessingNumberParameterType.Double

EXAMPLES = """# One output per line:  NAME [UNIT] = expression     (lines starting with # are ignored)
# Variables: TELEMAC symbols (U, V, H, S, B, F, Q, M ...), x, y, t (s), area, g
VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)
HMAX [M] = tmax(H)
"""


def template_for(path: str) -> str:
    """The expression box pre-filled for a file: its variables as comments, then two examples."""
    try:
        with SelafinFile(path) as s:
            pairs = describe_variables(s.var_names)
            names = [n for _sym, n in pairs]
    except Exception:
        return EXAMPLES
    lines = ["# One output per line:  NAME [UNIT] = expression     (lines starting with # are ignored)",
             "# Variables of this file: " + ", ".join("%s = %s" % (sym, n) for sym, n in pairs),
             "# Also: x, y, t (time, s), area (around each node), g = 9.81"]
    cids = {canonical(n) for n in names}
    if {"velocity_u", "velocity_v"} <= cids:
        lines.append("VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)" if "water_depth" in cids else
                     "VORTICITY [1/S] = curl(U, V)")
    if "water_depth" in cids:
        lines.append("HMAX [M] = tmax(H)")
    if len(lines) == 3:
        sym = pairs[0][0] if pairs else "B"
        lines.append("RESULT = %s" % sym)
    return "\n".join(lines) + "\n"


def _plugin():
    try:
        from qgis.utils import plugins

        return plugins.get(__package__.split(".")[0])
    except Exception:
        return None


class CalculatorAlgorithm(_Base):
    def name(self):
        return "calculator"

    def displayName(self):
        return self.tr("Calculator (formulas, gradients, vorticity, time maxima)")

    def shortHelpString(self):
        funcs = "\n".join("• %s: %s" % (a, b) for a, b in FUNCTIONS_HELP)
        return self.tr(
            "New variables computed from the variables of a SELAFIN file, written to a new file (the original "
            "variables can be kept). One output per line: NAME [UNIT] = expression, in NumPy syntax. Each output "
            "can be used by its name in the following lines. Two outputs named 'NAME X' and 'NAME Y' (or U and V) "
            "are shown as a vector with arrows (X and Y are written to the file as 'NAME ALONG X/Y', the names "
            "QGIS reads as a vector).\n\n"
            "Variables: TELEMAC symbols (U, V, H water depth, S free surface, B bottom, F Froude, Q unit "
            "discharge, M velocity magnitude, ...), names with underscores (WATER_DEPTH) or canonical ids "
            "(water_depth); x, y, t (time in s), area (area around each node), g = 9.81.\n\n"
            "Functions:\n%s\n\n"
            "Spatial derivatives use the linear finite-element interpolation of TELEMAC: the gradient is "
            "constant in each triangle and averaged to the nodes (exact for linear fields; second order on "
            "regular meshes, first order on irregular ones). Near dry areas mask the result, e.g. "
            "where(H > 0.01, curl(U, V), 0).\n\n"
            "Examples:\n"
            "• FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))\n"
            "• SHEAR [PA] = 1000 * g * 0.03**2 * mag(U, V)**2 / maximum(H, 0.01)**(1/3)\n"
            "• HMAX [M] = tmax(H); TIME HMAX [S] = tpeak(H)\n"
            "• SLOPE X = ddx(S); SLOPE Y = ddy(S)\n"
            "• RISE [M/S] = ddt(S)\n"
            "• DEPTH CHANGE [M] = tlast(H) - tfirst(H)") % funcs

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN file"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterString("EXPRESSIONS", self.tr("Outputs (NAME [UNIT] = expression, one per line)"),
                                                       EXAMPLES, multiLine=True))
        self.addParameter(QgsProcessingParameterBoolean("KEEP", self.tr("Keep the variables of the file"), True))
        self.addParameter(QgsProcessingParameterString("FRAMES", self.tr("Frames (empty = all; e.g. -1, 0:10, 0:100:5)"), "", optional=True))
        p = QgsProcessingParameterNumber("DRY", self.tr("Depth below which wet() is 0 (m)"), _DOUBLE, 0.01, minValue=0.0)
        p.setFlags(p.flags() | Qgis.ProcessingParameterFlag.Advanced)
        self.addParameter(p)
        p = QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS saved next to the output (.prj; empty = copy the input's)"), optional=True)
        p.setFlags(p.flags() | Qgis.ProcessingParameterFlag.Advanced)
        self.addParameter(p)
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Output SELAFIN"), "SELAFIN (*.slf)"))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the result when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsFile(parameters, "INPUT", context)
        if not src or not os.path.exists(src):
            raise QgsProcessingException(self.tr("File not found: %s") % src)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        if os.path.abspath(out) == os.path.abspath(src):
            raise QgsProcessingException(self.tr("Write to a new file, not over the input."))
        text = self.parameterAsString(parameters, "EXPRESSIONS", context)
        frames = None
        ftxt = self.parameterAsString(parameters, "FRAMES", context)
        if ftxt and ftxt.strip():
            with SelafinFile(src) as s:
                try:
                    frames = SubsetAlgorithm.parse_frames(ftxt, s.nframes)
                except ValueError as exc:
                    raise QgsProcessingException(self.tr("Could not read the frame list: %s") % exc)

        def prog(frac):
            feedback.setProgress(100.0 * frac)
            return feedback.isCanceled()

        try:
            rep = run_calculator(src, out, text, self.parameterAsBoolean(parameters, "KEEP", context), frames, prog,
                                 self.parameterAsDouble(parameters, "DRY", context))
        except CalculatorError as exc:
            raise QgsProcessingException(str(exc))
        for o in rep.outputs:
            lo, hi = rep.ranges[o.name]
            feedback.pushInfo("%s%s = %s   [%s%s]" % (o.name, " (%s)" % o.unit if o.unit else "", o.expr,
                                                       ("%.4g to %.4g" % (lo, hi)) if np.isfinite(lo) else "no finite value",
                                                       "" if o.time_dependent else ", same in every frame"))
        for a, b in rep.renamed.items():
            feedback.pushInfo(self.tr("%s written as %s, so that QGIS draws the X/Y pair as a vector with arrows")
                              % (a, b))
        feedback.pushInfo(self.tr("%d frame(s) written to %s") % (len(rep.frames), out))
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        prj_in = os.path.splitext(src)[0] + ".prj"
        if crs.isValid():
            from ..qgis_bridge.crs import write_prj

            write_prj(out, crs)
        elif os.path.exists(prj_in):
            shutil.copyfile(prj_in, os.path.splitext(out)[0] + ".prj")
        self._load = out if self.parameterAsBoolean(parameters, "LOAD", context) else None
        return {"OUTPUT": out}

    def postProcessAlgorithm(self, context, feedback):
        out = getattr(self, "_load", None)
        self._load = None
        if out:
            pl = _plugin()
            # qgis_process and standalone scripts load the plugin without a QGIS window (iface is None)
            if pl is not None and hasattr(pl, "open_path") and getattr(pl, "iface", None) is not None:
                pl.open_path(out)
            else:
                from qgis.core import QgsProcessingContext, QgsProcessingUtils

                name = os.path.splitext(os.path.basename(out))[0]
                context.addLayerToLoadOnCompletion(out, QgsProcessingContext.LayerDetails(
                    name, context.project(), name, QgsProcessingUtils.LayerHint.Mesh))
        return {}


class IntegralsAlgorithm(_Base):
    MODES = ["Discharge through section lines (m³/s)",
             "Flux of an expression through section lines (expression × unit discharge, e.g. a concentration)",
             "Line integral and mean of an expression along lines",
             "Area integral and mean of an expression over polygons (or the whole mesh); H gives the volume"]

    def name(self):
        return "integrals"

    def displayName(self):
        return self.tr("Integrals and fluxes (discharge, volume, line and area integrals)")

    def shortHelpString(self):
        return self.tr(
            "Integrals of a SELAFIN result for every frame, as a table:\n"
            "• Discharge through section lines: the integral of h u.n along each line, with depth and velocity "
            "interpolated as in the model (the line is cut at every mesh edge and the product integrated exactly). "
            "Positive from the right of the line to its left: draw a section from the left bank to the right bank "
            "(looking downstream) and the flow downstream is positive.\n"
            "• Flux of an expression: the same with a weight, e.g. a tracer concentration for its mass flux.\n"
            "• Line integral and mean of an expression along each line (e.g. mean velocity across a section).\n"
            "• Area integral and mean over each polygon (elements whose centre is inside) or the whole mesh: "
            "H gives the water volume, 1 the area, where(H > 0.01, 1, 0) the wet area.\n"
            "Expressions follow the Calculator syntax. Lines and polygons are transformed into the mesh CRS "
            "(a .prj next to the file or the CRS given here).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("SELAFIN result"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterEnum("MODE", self.tr("Integral"), [self.tr(m) for m in self.MODES], defaultValue=0))
        self.addParameter(QgsProcessingParameterString("EXPRESSION", self.tr("Expression (flux weight, line or area integrand)"), "H", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("LINES", self.tr("Section lines"), [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("POLYGONS", self.tr("Polygons (empty = whole mesh)"),
                                                              [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterField("NAME_FIELD", self.tr("Name field"), "name", "LINES", optional=True))
        self.addParameter(QgsProcessingParameterString("FRAMES", self.tr("Frames (empty = all)"), "", optional=True))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = .prj next to the file, else the layers' CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Integrals"), Qgis.ProcessingSourceType.Vector))

    def _shapes(self, parameters, context, key, mesh_crs, polygon):
        src = self.parameterAsSource(parameters, key, context)
        if src is None:
            return []
        xf = None
        if mesh_crs is not None and mesh_crs.isValid() and src.sourceCrs().isValid() and mesh_crs != src.sourceCrs():
            xf = QgsCoordinateTransform(src.sourceCrs(), mesh_crs, context.transformContext())
        nf = self.parameterAsString(parameters, "NAME_FIELD", context) or "name"
        out = []
        for f in src.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            if xf is not None:
                g.transform(xf)
            idx = f.fields().indexOf(nf)
            name = str(f.attribute(idx)) if idx >= 0 and f.attribute(idx) not in (None, "") else "%s %d" % (
                "polygon" if polygon else "line", f.id())
            if polygon:
                rings = []
                for poly in (g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]):
                    rings += [np.array([(p.x(), p.y()) for p in r]) for r in poly]
                out.append((f.id(), name, rings))
            else:
                parts = g.asMultiPolyline() if g.isMultipart() else [g.asPolyline()]
                for k, part in enumerate(parts):
                    out.append((f.id(), name if len(parts) == 1 else "%s/%d" % (name, k + 1),
                                np.array([(p.x(), p.y()) for p in part])))
        return out

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge.crs import sidecar_crs

        path = self.parameterAsFile(parameters, "INPUT", context)
        mode = self.parameterAsEnum(parameters, "MODE", context)
        expr = (self.parameterAsString(parameters, "EXPRESSION", context) or "").strip()
        mesh_crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not mesh_crs.isValid():
            mesh_crs = sidecar_crs(path)
        with self._open(parameters, context) as slf:
            h = slf.header
            if h.is_3d or h.ndp != 3:
                raise QgsProcessingException(self.tr("A 2D triangular result is required (extract 2D from 3D first)."))
            cids = [canonical(n) for n in h.var_names]
            calc = None
            if mode in (1, 2, 3):
                if not expr:
                    raise QgsProcessingException(self.tr("Give an expression."))
                try:
                    calc = Calculator(slf, [Output("E", "", expr)])
                except CalculatorError as exc:
                    raise QgsProcessingException(str(exc))
            ftxt = self.parameterAsString(parameters, "FRAMES", context)
            try:
                frames = SubsetAlgorithm.parse_frames(ftxt, slf.nframes)
            except ValueError as exc:
                raise QgsProcessingException(str(exc))
            fields = QgsFields()
            for n, k in (("feature", "int"), ("name", "string"), ("frame", "int"), ("time_s", "double")):
                fields.append(_field(n, k))
            if h.has_calendar_date:
                fields.append(_field("datetime", "string"))
            for n in (("discharge", "flux", "integral", "integral")[mode], "length" if mode < 3 else "area", "mean",
                      "wet_length" if mode < 3 else "wet_area"):
                fields.append(_field(n, "double"))
            sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, Qgis.WkbType.NoGeometry)
            from ..telemac_io import TriangleLocator

            hidx = cids.index("water_depth") if "water_depth" in cids else None
            if mode < 3:
                shapes = self._shapes(parameters, context, "LINES", mesh_crs, False)
                if not shapes:
                    raise QgsProcessingException(self.tr("Give section lines."))
                loc = TriangleLocator(h.x, h.y, h.ikle)
                secs = [(fid, name, Section.build(line, h.x, h.y, h.ikle, loc)) for fid, name, line in shapes]
                for fid, name, sec in secs:
                    if not sec.inside.any():
                        feedback.pushWarning(self.tr("%s does not cross the mesh (check the CRS)") % name)
                if mode in (0, 1):
                    if {"water_depth", "velocity_u", "velocity_v"} <= set(cids):
                        iq = (cids.index("water_depth"), cids.index("velocity_u"), cids.index("velocity_v"))
                        how = "depth x velocity"
                    elif {"flowrate_x", "flowrate_y"} <= set(cids):
                        iq = (None, cids.index("flowrate_x"), cids.index("flowrate_y"))
                        how = "FLOWRATE ALONG X/Y"
                    else:
                        raise QgsProcessingException(self.tr("The file needs WATER DEPTH and VELOCITY U/V, or FLOWRATE "
                                                             "ALONG X/Y, for a discharge."))
                    feedback.pushInfo(self.tr("Discharge from %s") % how)
            else:
                shapes = self._shapes(parameters, context, "POLYGONS", mesh_crs, True) or [(0, "whole mesh", None)]
                masks = [(fid, name, None if rings is None else element_mask_in_rings(h.x, h.y, h.ikle, rings))
                         for fid, name, rings in shapes]
                ops = calc.ops
                for fid, name, m in masks:
                    if m is not None and not m.any():
                        feedback.pushWarning(self.tr("%s contains no element centre (check the CRS)") % name)
            dts = slf.datetimes if h.has_calendar_date else None
            last = {}
            for k, fi in enumerate(frames):
                if feedback.isCanceled():
                    break
                data = slf.read_frame(fi)
                wet = (np.asarray(data[hidx], dtype=np.float64) > 0.01).astype(np.float64) if hidx is not None else None
                val = calc.evaluate(calc.outputs[0], fi) if calc is not None else None
                rows = []
                if mode < 3:
                    for fid, name, sec in secs:
                        if mode == 2:
                            v = sec.integral(val)
                        elif iq[0] is not None:
                            v = sec.discharge(data[iq[0]], data[iq[1]], data[iq[2]], val if mode == 1 else None)
                        else:
                            v = sec.flux(data[iq[1]], data[iq[2]], val if mode == 1 else None)
                        meas = sec.covered_length()
                        mean = v / meas if (mode == 2 and meas > 0) else None
                        rows.append((fid, name, v, meas, mean, sec.covered_length(wet) if wet is not None else None))
                else:
                    for fid, name, m in masks:
                        v = ops.integral(val, m)
                        area = float(ops.area[m].sum()) if m is not None else float(ops.area.sum())
                        rows.append((fid, name, v, area, v / area if area > 0 else None,
                                     ops.integral(wet, m) if wet is not None else None))
                for fid, name, v, meas, mean, wm in rows:
                    f = QgsFeature(fields)
                    attrs = [int(fid), name, int(fi), float(slf.times[fi])]
                    if dts is not None:
                        attrs.append(dts[fi].isoformat())
                    attrs += [float(v), float(meas), None if mean is None else float(mean), None if wm is None else float(wm)]
                    f.setAttributes(attrs)
                    sink.addFeature(f, QgsFeatureSink.Flag.FastInsert)
                    last[name] = (fi, v, mean)
                feedback.setProgress(100.0 * (k + 1) / max(len(frames), 1))
            label = ("discharge", "flux", "integral", "integral")[mode]
            for name, (fi, v, mean) in last.items():
                feedback.pushInfo(self.tr("%s, frame %d (t = %g s): %s = %.6g%s") % (
                    name, fi, float(slf.times[fi]), label, v, "" if mean is None else ", mean %.6g" % mean))
        return {"OUTPUT": dest}


CALC_ALGORITHMS = [CalculatorAlgorithm, IntegralsAlgorithm]
