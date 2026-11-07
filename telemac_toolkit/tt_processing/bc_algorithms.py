"""Processing tools for boundary conditions and run files: rating curves, hydrographs, checks,
initial conditions and the steering file."""

from __future__ import annotations

import datetime as _dt
import os
from typing import Dict, Optional

import numpy as np
from qgis.core import (
    Qgis,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterDateTime,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterString,
    QgsWkbTypes,
)

from ..mesher.boundary import attach_mesh
from ..mesher.hydraulics import boundary_report, critical_depth_check, initial_state, rating_curve, times_to_seconds
from ..telemac_io import SelafinFile, read_cli, write_selafin
from ..telemac_io.cli import liquid_boundaries
from ..telemac_io.liq import read_liq, read_stage_discharge, read_steering, write_liq, write_stage_discharge
from ..telemac_io.steering import FRICTION_LAWS, courant_time_step, fitting_time_step, steering_template
from ..telemac_io.variables import canonical
from .algorithms import SLF_FILTER, _Base, _field
from .mesh_algorithms import CLI_FILTER, LIQ_FILTER, _advanced

_DOUBLE = Qgis.ProcessingNumberParameterType.Double
_INT = Qgis.ProcessingNumberParameterType.Integer
TXT_FILTER = "Text files (*.txt *.qz);;All files (*.*)"
CAS_FILTER = "Steering files (*.cas);;All files (*.*)"


def _cli_next_to(path: str) -> str:
    return os.path.splitext(path)[0] + ".cli"


def _var(slf: SelafinFile, cid: str) -> Optional[str]:
    for n in slf.header.var_names:
        if canonical(n) == cid:
            return n
    return None


def load_mesh(mesh_path: str, cli_path: Optional[str] = None):
    """(x, y, ikle, bottom, friction or None, BoundaryTable numbered as TELEMAC numbers it, cli path)."""
    if not mesh_path or not os.path.exists(mesh_path):
        raise QgsProcessingException("File not found: %s" % mesh_path)
    with SelafinFile(mesh_path) as s:
        h = s.header
        if h.is_3d or h.ndp != 3:
            raise QgsProcessingException("A 2D triangular mesh is required.")
        x, y, ik = h.x.copy(), h.y.copy(), h.ikle.copy()
        zb = _var(s, "bottom")
        if zb is None:
            raise QgsProcessingException("%s has no BOTTOM (FOND) variable." % os.path.basename(mesh_path))
        z = s.read(zb, 0)
        fn = _var(s, "bottom_friction")
        fr = s.read(fn, 0) if fn else None
    cli = cli_path or _cli_next_to(mesh_path)
    if not os.path.exists(cli):
        raise QgsProcessingException("Boundary conditions file not found: %s" % cli)
    table = read_cli(cli)
    if len(table) and int(table.records["n"].max()) > len(x):
        raise QgsProcessingException("%s refers to node %d but the mesh has %d nodes: they do not belong together."
                                     % (os.path.basename(cli), int(table.records["n"].max()), len(x)))
    attach_mesh(table, x, y, ik)
    return x, y, ik, z, fr, table, cli


class _BCBase(_Base):
    def group(self):
        return self.tr("Boundary conditions and run")

    def groupId(self):
        return "boundaries"

    def _mesh_params(self):
        self.addParameter(QgsProcessingParameterFile("MESH", self.tr("Mesh (SELAFIN geometry with BOTTOM)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterFile("CLI", self.tr("Boundary conditions file (empty = .cli next to the mesh)"),
                                                     fileFilter=CLI_FILTER, optional=True))

    def _load(self, parameters, context):
        return load_mesh(self.parameterAsFile(parameters, "MESH", context), self.parameterAsFile(parameters, "CLI", context) or None)


# --------------------------------------------------------------------------
class StageDischargeAlgorithm(_BCBase):
    def name(self):
        return "stagedischarge"

    def displayName(self):
        return self.tr("Stage-discharge curve for an outflow boundary")

    def shortHelpString(self):
        return self.tr(
            "Normal-depth rating curve Z(Q) of a liquid boundary, computed on the boundary itself (its nodes and "
            "their bottom, as TELEMAC sees the section) with Manning's formula Q = sqrt(S) A R^(2/3) / n, the section "
            "divided where the roughness changes. Written in TELEMAC's STAGE-DISCHARGE CURVES FILE format; an existing "
            "file keeps its other curves. The boundary must be a prescribed elevation (5 4 4); the steering file then "
            "needs STAGE-DISCHARGE CURVES = 1 for it (the log gives the lines).\n\n"
            "Slope: the energy slope, usually the bed or water-surface slope of the reach downstream. Roughness: one "
            "Manning's n, or the BOTTOM FRICTION variable of the mesh (Manning law). Above the lower end of the "
            "section, its ends act as vertical walls, like the closed boundaries of the mesh. A warning is given if "
            "the flow at the section is supercritical, where a rating curve does not control the level.")

    def initAlgorithm(self, config=None):
        self._mesh_params()
        self.addParameter(QgsProcessingParameterNumber("BOUNDARY", self.tr("Liquid boundary number (0 = the only prescribed elevation)"),
                                                       _INT, 0, minValue=0))
        self.addParameter(QgsProcessingParameterNumber("SLOPE", self.tr("Energy slope (m/m)"), _DOUBLE, 0.001, minValue=1e-7))
        self.addParameter(QgsProcessingParameterNumber("MANNING", self.tr("Manning's n"), _DOUBLE, 0.03, minValue=1e-4))
        self.addParameter(QgsProcessingParameterBoolean(
            "USE_MESH_FRICTION", self.tr("Use the BOTTOM FRICTION of the mesh (Manning) when present"), True))
        self.addParameter(_advanced(QgsProcessingParameterNumber("TOP", self.tr("Highest level of the curve (m, empty = highest ground)"),
                                                                 _DOUBLE, optional=True)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("POINTS", self.tr("Points of the curve"), _INT, 40, minValue=5, maxValue=500)))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Stage-discharge curves file"), TXT_FILTER))
        self.addParameter(QgsProcessingParameterFeatureSink("TABLE", self.tr("Curve table"), Qgis.ProcessingSourceType.Vector,
                                                            optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        x, y, ik, z, fr, table, _cli = self._load(parameters, context)
        lb = liquid_boundaries(table)
        num = self.parameterAsInt(parameters, "BOUNDARY", context)
        if num == 0:
            hs = [b for b in lb if tuple(b["codes"]) == (5, 4, 4)]
            if len(hs) != 1:
                raise QgsProcessingException(self.tr("Give the boundary number: the mesh has %d prescribed elevation boundaries.") % len(hs))
            b = hs[0]
        else:
            b = next((bb for bb in lb if bb["number"] == num), None)
            if b is None:
                raise QgsProcessingException(self.tr("The mesh has no liquid boundary %d (it has %d).") % (num, len(lb)))
        num = b["number"]
        if tuple(b["codes"]) != (5, 4, 4):
            feedback.pushWarning(self.tr("Boundary %d is %d %d %d; a stage-discharge curve needs a prescribed elevation (5 4 4).")
                                 % ((num,) + tuple(b["codes"])))
        nodes = np.asarray(b["nodes"], dtype=np.int64)
        st = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x[nodes]), np.diff(y[nodes])))])
        n = self.parameterAsDouble(parameters, "MANNING", context)
        nn = n
        if fr is not None and self.parameterAsBoolean(parameters, "USE_MESH_FRICTION", context):
            nn = fr[nodes]
            feedback.pushInfo(self.tr("Manning's n from the mesh: %.4g to %.4g") % (float(nn.min()), float(nn.max())))
        top = parameters.get("TOP")
        top = None if top in (None, "") else float(top)
        slope = self.parameterAsDouble(parameters, "SLOPE", context)
        try:
            Z, Q, notes = rating_curve(st, z[nodes], nn, slope, count=self.parameterAsInt(parameters, "POINTS", context), top=top)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        for m in notes:
            feedback.pushInfo(m)
        fru = critical_depth_check(Z, Q, st, z[nodes])
        if np.any(fru[1:] > 1.0):
            feedback.pushWarning(self.tr("The flow is supercritical at the section for levels above %.3f m (Froude up to %.2f): "
                                         "a rating curve does not control the level there.")
                                 % (float(Z[1:][fru[1:] > 1][0]), float(fru.max())))
        path = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        curves = {}
        if os.path.exists(path):
            try:
                curves = read_stage_discharge(path)
                if curves:
                    feedback.pushInfo(self.tr("Kept the other curves of %s: %s") % (os.path.basename(path),
                                                                                    ", ".join(str(k) for k in curves if k != num)))
            except ValueError:
                curves = {}
        curves[num] = (Q, Z)
        name = table.comments[b["start"]] if b["start"] < len(table.comments) else ""
        write_stage_discharge(path, curves, {num: "normal depth, S = %g, %s" % (slope, ("n from the mesh" if np.ndim(nn) else "n = %g" % n))
                                             + (" (%s)" % name if name else "")})
        nl = len(lb)
        feedback.pushInfo(self.tr("Steering file: STAGE-DISCHARGE CURVES = %s") % ";".join(
            "1" if bb["number"] in curves else "0" for bb in lb))
        feedback.pushInfo(self.tr("Steering file: STAGE-DISCHARGE CURVES FILE = '%s'") % os.path.basename(path))
        feedback.pushInfo(self.tr("Boundary %d: %d nodes, bed %.3f to %.3f m; %d liquid boundaries in all")
                          % (num, len(nodes), float(z[nodes].min()), float(z[nodes].max()), nl))
        res = {"OUTPUT": path}
        if parameters.get("TABLE"):
            fields = QgsFields()
            for f in ("level", "discharge", "froude"):
                fields.append(_field(f, "double"))
            sink, dest = self.parameterAsSink(parameters, "TABLE", context, fields, QgsWkbTypes.Type.NoGeometry)
            for a, q, f in zip(Z, Q, fru):
                ft = QgsFeature(fields)
                ft.setAttributes([float(a), float(q), float(f)])
                sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
            res["TABLE"] = dest
        return res


# --------------------------------------------------------------------------
class LiquidBoundariesFileAlgorithm(_BCBase):
    UNITS = ["seconds", "minutes", "hours", "days", "date and time"]
    _UNIT_KEYS = ["s", "min", "h", "d", "datetime"]

    def name(self):
        return "liquidboundariesfile"

    def displayName(self):
        return self.tr("Liquid boundaries file from a table (hydrographs)")

    def shortHelpString(self):
        return self.tr(
            "Write TELEMAC's LIQUID BOUNDARIES FILE (.liq) from a table: a CSV or Excel sheet loaded in QGIS, or any "
            "layer, with a time column and one column per prescribed value (for example HEC-HMS output). The mesh and "
            "its .cli tell which columns TELEMAC expects, in its numbering: Q(n) for a discharge boundary, SL(n) for "
            "a water level (not needed where a stage-discharge curve gives it), U(n) and V(n) for velocities.\n\n"
            "Columns: list them as COLUMN=field, for example 'Q(1)=flow_m3s; SL(2)=level'; or leave it empty and "
            "choose the value fields in the order of the expected columns (the log lists them). Times: numbers in "
            "the chosen unit, or dates and times; they are written in seconds from the first row, and the first date "
            "is given for ORIGINAL DATE OF TIME. Rows are sorted by time; rows with an empty value are skipped.")

    def initAlgorithm(self, config=None):
        self._mesh_params()
        self.addParameter(QgsProcessingParameterFeatureSource("TABLE", self.tr("Table (CSV, Excel or any layer)"),
                                                              [Qgis.ProcessingSourceType.Vector]))
        self.addParameter(QgsProcessingParameterField("TIME_FIELD", self.tr("Time field"), parentLayerParameterName="TABLE"))
        self.addParameter(QgsProcessingParameterEnum("TIME_UNIT", self.tr("Time unit"), [self.tr(u) for u in self.UNITS], defaultValue=2))
        self.addParameter(QgsProcessingParameterField("VALUE_FIELDS", self.tr("Value fields, in the order of the expected columns"),
                                                      parentLayerParameterName="TABLE", allowMultiple=True, optional=True,
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric))
        self.addParameter(QgsProcessingParameterString("COLUMNS", self.tr("Or the columns by name, e.g. Q(1)=flow; SL(2)=level"),
                                                       optional=True))
        self.addParameter(QgsProcessingParameterFile("STAGE", self.tr("Stage-discharge curves file (its boundaries need no SL column)"),
                                                     fileFilter=TXT_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Liquid boundaries file"), LIQ_FILTER))

    @staticmethod
    def expected_columns(table, stage_numbers=()):
        out = []
        for b in liquid_boundaries(table):
            h, u, _v = (int(c) for c in b["codes"])
            n = b["number"]
            if u == 5:
                out.append(("Q(%d)" % n, "m3/s"))
            if h == 5 and n not in stage_numbers:
                out.append(("SL(%d)" % n, "m"))
            if u == 6:
                out += [("U(%d)" % n, "m/s"), ("V(%d)" % n, "m/s")]
        return out

    def processAlgorithm(self, parameters, context, feedback):
        _x, _y, _ik, _z, _fr, table, _cli = self._load(parameters, context)
        stage_path = self.parameterAsFile(parameters, "STAGE", context)
        stage_nums = set(read_stage_discharge(stage_path)) if stage_path else set()
        expected = self.expected_columns(table, stage_nums)
        if not expected:
            raise QgsProcessingException(self.tr("The mesh has no liquid boundary with a prescribed value."))
        feedback.pushInfo(self.tr("Columns TELEMAC expects: %s") % ", ".join(c for c, _u in expected))
        src = self.parameterAsSource(parameters, "TABLE", context)
        tfield = self.parameterAsString(parameters, "TIME_FIELD", context)
        mapping: Dict[str, str] = {}
        text = self.parameterAsString(parameters, "COLUMNS", context).strip()
        if text:
            for part in text.replace(",", ";").split(";"):
                if not part.strip():
                    continue
                if "=" not in part:
                    raise QgsProcessingException(self.tr("Write each column as NAME=field, e.g. Q(1)=flow: %r") % part)
                k, v = (t.strip() for t in part.split("=", 1))
                mapping[k.upper()] = v
        else:
            fields = self.parameterAsStrings(parameters, "VALUE_FIELDS", context)
            if len(fields) != len(expected):
                raise QgsProcessingException(self.tr("Choose %d value fields, one per expected column (%s), or write the columns by name.")
                                             % (len(expected), ", ".join(c for c, _u in expected)))
            mapping = {c: f for (c, _u), f in zip(expected, fields)}
        names = src.fields().names()
        for c, f in mapping.items():
            if f not in names:
                raise QgsProcessingException(self.tr("Field %r not found in the table") % f)
        missing = [c for c, _u in expected if c not in mapping]
        if missing:
            feedback.pushWarning(self.tr("No column for %s: TELEMAC will take the steering file values there.") % ", ".join(missing))
        unit_names = dict(expected)
        unit = self._UNIT_KEYS[self.parameterAsEnum(parameters, "TIME_UNIT", context)]
        rows = []
        for ft in src.getFeatures():
            tv = ft[tfield]
            if tv is None or (hasattr(tv, "isNull") and tv.isNull()):
                continue
            vals = []
            ok = True
            for c in mapping:
                v = ft[mapping[c]]
                try:
                    vals.append(float(v))
                except (TypeError, ValueError):
                    ok = False
                    break
            if not ok:
                continue
            if unit == "datetime":
                if hasattr(tv, "toPyDateTime"):
                    tv = tv.toPyDateTime()
                elif hasattr(tv, "toPyDate"):
                    tv = _dt.datetime.combine(tv.toPyDate(), _dt.time())
                elif isinstance(tv, str):
                    try:
                        tv = _dt.datetime.fromisoformat(tv.strip().replace("/", "-"))
                    except ValueError:
                        raise QgsProcessingException(self.tr("Cannot read %r as a date and time (use YYYY-MM-DD HH:MM)") % tv)
            else:
                try:
                    tv = float(tv)
                except (TypeError, ValueError):
                    raise QgsProcessingException(self.tr("Cannot read %r as a number of %s") % (tv, self.UNITS[self._UNIT_KEYS.index(unit)]))
            rows.append((tv, vals))
        if len(rows) < 2:
            raise QgsProcessingException(self.tr("The table needs at least two rows with a time and every value."))
        rows.sort(key=lambda r: r[0])
        secs, first = times_to_seconds([r[0] for r in rows], unit)
        keep = np.concatenate([[True], np.diff(secs) > 0])
        if not keep.all():
            feedback.pushWarning(self.tr("%d rows with a repeated time skipped") % int((~keep).sum()))
        secs = secs[keep]
        data = np.array([r[1] for r in rows], dtype=np.float64)[keep]
        cols = [(c, unit_names.get(c, ""), data[:, i]) for i, c in enumerate(mapping)]
        path = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        try:
            write_liq(path, secs, cols, comment="Liquid boundaries file written by TELEMAC Toolkit from %s" % src.sourceName())
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        feedback.pushInfo(self.tr("%d times from 0 to %g s; columns %s") % (len(secs), secs[-1], ", ".join(mapping)))
        if first is not None:
            feedback.pushInfo(self.tr("Steering file: ORIGINAL DATE OF TIME = %d;%d;%d") % (first.year, first.month, first.day))
            feedback.pushInfo(self.tr("Steering file: ORIGINAL HOUR OF TIME = %d;%d;%d") % (first.hour, first.minute, first.second))
        feedback.pushInfo(self.tr("Steering file: LIQUID BOUNDARIES FILE = '%s'") % os.path.basename(path))
        return {"OUTPUT": path}


# --------------------------------------------------------------------------
class CheckBoundariesAlgorithm(_BCBase):
    def name(self):
        return "checkboundaries"

    def displayName(self):
        return self.tr("Check boundary conditions")

    def shortHelpString(self):
        return self.tr(
            "Lists the liquid boundaries of a mesh with the numbers TELEMAC gives them, and checks them against the "
            "liquid boundaries file, the stage-discharge curves and the steering file when given (or found through "
            "the steering file): a column for every prescribed value, levels that leave a boundary dry, curves on "
            "boundaries that are not prescribed elevations, value counts in the steering file, two liquid boundaries "
            "touching without a wall node between them, boundaries with fewer than 5 nodes. The output layer draws "
            "each liquid boundary labelled with its number.")

    def initAlgorithm(self, config=None):
        self._mesh_params()
        self.addParameter(QgsProcessingParameterFile("STEERING", self.tr("Steering file (optional; gives the other files)"),
                                                     fileFilter=CAS_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFile("LIQ", self.tr("Liquid boundaries file"), fileFilter=LIQ_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFile("STAGE", self.tr("Stage-discharge curves file"), fileFilter=TXT_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Liquid boundaries"), Qgis.ProcessingSourceType.VectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        x, y, _ik, z, _fr, table, _cli = self._load(parameters, context)
        cas = None
        liq_path = self.parameterAsFile(parameters, "LIQ", context)
        stage_path = self.parameterAsFile(parameters, "STAGE", context)
        cas_path = self.parameterAsFile(parameters, "STEERING", context)
        if cas_path:
            cas = read_steering(cas_path)
            base = os.path.dirname(cas_path)
            if not liq_path and cas.get("LIQUID BOUNDARIES FILE"):
                liq_path = os.path.join(base, cas["LIQUID BOUNDARIES FILE"])
            if not stage_path and cas.get("STAGE-DISCHARGE CURVES FILE"):
                stage_path = os.path.join(base, cas["STAGE-DISCHARGE CURVES FILE"])
        liq = stage = None
        if liq_path:
            if not os.path.exists(liq_path):
                raise QgsProcessingException(self.tr("Liquid boundaries file not found: %s") % liq_path)
            liq = read_liq(liq_path)
        if stage_path:
            if not os.path.exists(stage_path):
                raise QgsProcessingException(self.tr("Stage-discharge curves file not found: %s") % stage_path)
            stage = read_stage_discharge(stage_path)
        bnds, msgs = boundary_report(x, y, z, table, liq, stage, cas)
        n_warn = 0
        for m in msgs:
            if "Warning" in m:
                n_warn += 1
                feedback.pushWarning(m.strip())
            else:
                feedback.pushInfo(m)
        feedback.pushInfo(self.tr("%d liquid boundaries, %d warnings") % (len(bnds), n_warn))
        fields = QgsFields()
        for name, kind in (("number", "int"), ("name", "string"), ("kind", "string"), ("codes", "string"), ("nodes", "int"),
                           ("length", "double"), ("bed_min", "double"), ("bed_max", "double"), ("label", "string")):
            fields.append(_field(name, kind))
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid() and context.project() is not None:
            crs = context.project().crs()
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, QgsWkbTypes.Type.LineString, crs)
        for b in bnds:
            ft = QgsFeature(fields)
            pts = [QgsPointXY(float(a), float(c)) for a, c in zip(b["x"], b["y"])]
            ft.setGeometry(QgsGeometry.fromPolylineXY(pts) if len(pts) > 1 else QgsGeometry.fromPointXY(pts[0]).buffer(0.01, 4))
            ft.setAttributes([b["number"], b["name"], b["kind"], "%d %d %d" % b["codes"], int(len(b["nodes"])), b["length"],
                              b["zmin"], b["zmax"], "%d %s (%s)" % (b["number"], b["name"], b["kind"])])
            sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
        self._dest = dest
        return {"OUTPUT": dest, "WARNINGS": n_warn}

    def postProcessAlgorithm(self, context, feedback):
        dest = getattr(self, "_dest", None)
        lyr = context.getMapLayer(dest) if dest else None
        if lyr is not None:
            from qgis.core import QgsPalLayerSettings, QgsVectorLayerSimpleLabeling

            pal = QgsPalLayerSettings()
            pal.fieldName = "label"
            pal.enabled = True
            lyr.setLabeling(QgsVectorLayerSimpleLabeling(pal))
            lyr.setLabelsEnabled(True)
        return {}


# --------------------------------------------------------------------------
def _axis_lines(src, crs, project=None):
    """Lines of an axis layer; from a mesh project's channels layer only the axis lines (role axis)."""
    from ..qgis_bridge import meshing as M

    names = src.fields().names()
    out = []
    for f, g in M._features(src, crs, project):
        if "role" in names and str(f["role"]).strip().lower() not in ("axis", "a", "thalweg", "centreline", "centerline"):
            continue
        out += [c for c in M.line_parts(g) if len(c) >= 2]
    return out


def axis_bed_level(x, y, z, lines, depth: float, chunk: int = 2048) -> np.ndarray:
    """Water level at every node: the bed at the nearest point of the axis lines (nearest mesh node, running median
    over 5 points along each line) plus ``depth``."""
    from ..mesher.river import resample

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    pts, beds = [], []
    xy = np.stack([x, y], axis=1)
    for c in lines:
        L = float(np.hypot(*np.diff(c, axis=0).T).sum())
        if L <= 0:
            continue
        a = resample(c, max(L / 2000.0, 1e-6))
        near = np.empty(len(a), dtype=np.int64)
        for i in range(0, len(a), chunk):
            dd = (a[i:i + chunk, None, 0] - x[None, :]) ** 2 + (a[i:i + chunk, None, 1] - y[None, :]) ** 2
            near[i:i + chunk] = np.argmin(dd, axis=1)
        b = z[near]
        pad = np.concatenate([np.repeat(b[:1], 2), b, np.repeat(b[-1:], 2)])
        b = np.median(np.lib.stride_tricks.sliding_window_view(pad, 5), axis=1)
        pts.append(a)
        beds.append(b)
    P = np.vstack(pts)
    B = np.concatenate(beds)
    level = np.empty(len(x))
    for i in range(0, len(x), chunk):
        dd = (xy[i:i + chunk, None, 0] - P[None, :, 0]) ** 2 + (xy[i:i + chunk, None, 1] - P[None, :, 1]) ** 2
        level[i:i + chunk] = B[np.argmin(dd, axis=1)] + float(depth)
    return level


class InitialConditionsAlgorithm(_BCBase):
    MODES = ["Constant water level", "Water level along a channel axis (upstream and downstream levels)",
             "Depth above the bed", "Depth above the bed of the channel axes (rivers and networks)"]

    def name(self):
        return "initialconditions"

    def displayName(self):
        return self.tr("Initial conditions file")

    def shortHelpString(self):
        return self.tr(
            "Write a SELAFIN file with the initial state of a run (velocities zero, water depth, free surface, "
            "bottom), for PREVIOUS COMPUTATION FILE with INITIAL TIME SET TO ZERO = YES. The water level is constant, "
            "or varies linearly along a channel axis from the upstream to the downstream level (nodes beside the axis "
            "take the level of their position along it), or the depth above the bed is constant. Nodes where the level "
            "is below the bed are dry.\n\n"
            "Depth above the bed of the channel axes: for rivers and river networks. Every node takes the level of the "
            "nearest point of any axis line: the bed there (the nearest mesh node, smoothed along the line) plus the "
            "depth. The axis layer may be the mesh project's channels layer: only its axis lines are used. The banks "
            "and the floodplain above that level stay dry, and every branch starts with water.\n\n"
            "Start close to the expected flow: for a river, the downstream level from the rating curve at the first "
            "discharge and a level rising upstream with the slope. A state far from the first boundary values makes "
            "TELEMAC drain or fill the model at the start, which can dry a discharge boundary.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("MESH", self.tr("Mesh (SELAFIN geometry with BOTTOM)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterEnum("MODE", self.tr("Initial water"), [self.tr(m) for m in self.MODES], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("LEVEL", self.tr("Water level (m), or upstream level along the axis"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterNumber("LEVEL_DOWN", self.tr("Downstream level along the axis (m)"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterFeatureSource("AXIS", self.tr("Channel axis (line, drawn downstream)"),
                                                              [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber("DEPTH", self.tr("Depth above the bed (m)"), _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (for the axis; empty = project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Initial conditions file"), SLF_FILTER))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.planform import axis_coordinates
        from ..qgis_bridge import meshing as M

        path = self.parameterAsFile(parameters, "MESH", context)
        if not path or not os.path.exists(path):
            raise QgsProcessingException(self.tr("File not found: %s") % path)
        with SelafinFile(path) as s:
            h = s.header
            x, y, ik = h.x.copy(), h.y.copy(), h.ikle.copy()
            geo_double = bool(h.is_double)
            zb = _var(s, "bottom")
            if zb is None:
                raise QgsProcessingException(self.tr("The mesh has no BOTTOM variable."))
            z = s.read(zb, 0)
        mode = self.parameterAsEnum(parameters, "MODE", context)
        lv = self.parameterAsDouble(parameters, "LEVEL", context)
        if mode == 0:
            hh, ss = initial_state(z, level=lv)
            desc = "constant level %g m" % lv
        elif mode == 1:
            src = self.parameterAsSource(parameters, "AXIS", context)
            if src is None:
                raise QgsProcessingException(self.tr("Give the channel axis."))
            crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
            if not crs.isValid() and context.project() is not None:
                crs = context.project().crs()
            lines = M.line_arrays(src, crs, context.project())
            if len(lines) != 1:
                raise QgsProcessingException(self.tr("The axis layer must hold one line (it has %d).") % len(lines))
            ax = lines[0]
            s_ax, _n = axis_coordinates(ax, x, y)
            L = float(np.hypot(*np.diff(ax, axis=0).T).sum())
            ld = self.parameterAsDouble(parameters, "LEVEL_DOWN", context)
            level = lv + (ld - lv) * np.clip(s_ax / L, 0.0, 1.0)
            hh, ss = initial_state(z, level=level)
            desc = "level from %g m upstream to %g m downstream along %.4g m of axis" % (lv, ld, L)
        elif mode == 2:
            d = self.parameterAsDouble(parameters, "DEPTH", context)
            hh, ss = initial_state(z, depth=d)
            desc = "depth %g m above the bed" % d
        else:
            src = self.parameterAsSource(parameters, "AXIS", context)
            if src is None:
                raise QgsProcessingException(self.tr("Give the channel axes (or the channels layer of the mesh project)."))
            crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
            if not crs.isValid() and context.project() is not None:
                crs = context.project().crs()
            d = self.parameterAsDouble(parameters, "DEPTH", context)
            lines = _axis_lines(src, crs, context.project())
            if not lines:
                raise QgsProcessingException(self.tr("The axis layer has no line."))
            level = axis_bed_level(x, y, z, lines, d)
            hh, ss = initial_state(z, level=level)
            desc = "depth %g m above the bed of %d channel axis line(s)" % (d, len(lines))
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        zero = np.zeros_like(z)
        write_selafin(out, x, y, ik, {"VELOCITY U": zero, "VELOCITY V": zero, "WATER DEPTH": hh, "FREE SURFACE": ss, "BOTTOM": z},
                      units={"VELOCITY U": "M/S", "VELOCITY V": "M/S", "WATER DEPTH": "M", "FREE SURFACE": "M", "BOTTOM": "M"},
                      title="Initial conditions - TELEMAC Toolkit", double=geo_double)
        wet = hh > 0
        feedback.pushInfo(self.tr("Initial state: %s; %d of %d nodes wet, depth up to %.3f m")
                          % (desc, int(wet.sum()), len(z), float(hh.max()) if hh.size else 0.0))
        if not wet.any():
            feedback.pushWarning(self.tr("Every node is dry: the level is below the whole bed."))
        feedback.pushInfo(self.tr("Steering file: PREVIOUS COMPUTATION FILE = '%s'") % os.path.basename(out))
        feedback.pushInfo(self.tr("Steering file: INITIAL TIME SET TO ZERO = YES"))
        prj = os.path.splitext(path)[0] + ".prj"
        if os.path.exists(prj):
            import shutil

            shutil.copyfile(prj, os.path.splitext(out)[0] + ".prj")
        return {"OUTPUT": out}


# --------------------------------------------------------------------------
class SteeringFileAlgorithm(_BCBase):
    LAWS = ["Manning", "Strickler", "Chezy", "Nikuradse"]

    def name(self):
        return "steeringfile"

    def displayName(self):
        return self.tr("TELEMAC-2D steering file")

    def shortHelpString(self):
        return self.tr(
            "Write a TELEMAC-2D steering file (.cas) for a mesh, its .cli and the other files: liquid boundaries, "
            "stage-discharge curves, initial conditions. The liquid boundaries come from the .cli, in TELEMAC's "
            "numbering, with PRESCRIBED FLOWRATES and ELEVATIONS from the first row of the liquid boundaries file, "
            "STAGE-DISCHARGE CURVES for the boundaries in the curves file, and VELOCITY PROFILES = 4 (velocity in "
            "sqrt(h)) on discharge boundaries, which keeps them stable where the depth varies across the section. "
            "The time step can be estimated for a Courant number of 1 from a typical depth and velocity.\n\n"
            "The numerical settings are those of this plugin's demo cases that ran in TELEMAC-2D v9.1 (wetting and "
            "drying with positive depths). Review every value before running.")

    def initAlgorithm(self, config=None):
        self._mesh_params()
        self.addParameter(QgsProcessingParameterFile("LIQ", self.tr("Liquid boundaries file"), fileFilter=LIQ_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFile("STAGE", self.tr("Stage-discharge curves file"), fileFilter=TXT_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterFile("INITIAL", self.tr("Initial conditions file (empty = start dry or at a level)"),
                                                     fileFilter=SLF_FILTER, optional=True))
        self.addParameter(QgsProcessingParameterNumber("INITIAL_LEVEL", self.tr("Or a constant initial water level (m)"), _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterNumber("DURATION", self.tr("Duration (s)"), _DOUBLE, 3600.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("TIME_STEP", self.tr("Time step (s, 0 = estimate for Courant 1)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("DEPTH_GUESS", self.tr("Typical depth for the estimate (m)"), _DOUBLE, 2.0, minValue=0.01))
        self.addParameter(QgsProcessingParameterNumber("VELOCITY_GUESS", self.tr("Typical velocity for the estimate (m/s)"),
                                                       _DOUBLE, 1.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("OUTPUT_PERIOD", self.tr("Results every (s)"), _DOUBLE, 600.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterEnum("FRICTION_LAW", self.tr("Friction law"), self.LAWS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("FRICTION", self.tr("Friction coefficient (where the mesh has none)"),
                                                       _DOUBLE, 0.03, minValue=0.0))
        self.addParameter(QgsProcessingParameterString("TITLE", self.tr("Title"), "TELEMAC-2D model"))
        self.addParameter(QgsProcessingParameterDateTime("START", self.tr("Start date and time (optional)"), optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Steering file"), CAS_FILTER))

    def processAlgorithm(self, parameters, context, feedback):
        x, y, ik, _z, fr, table, cli = self._load(parameters, context)
        mesh = self.parameterAsFile(parameters, "MESH", context)
        lb = liquid_boundaries(table)
        kinds = [{(4, 5, 5): "Q", (5, 4, 4): "H", (5, 5, 5): "QH"}.get(tuple(int(c) for c in b["codes"]), "other") for b in lb]
        liq_path = self.parameterAsFile(parameters, "LIQ", context) or None
        stage_path = self.parameterAsFile(parameters, "STAGE", context) or None
        liq = read_liq(liq_path) if liq_path else {}
        stage = read_stage_discharge(stage_path) if stage_path else {}
        q0 = [float(liq["Q(%d)" % b["number"]][0]) if "Q(%d)" % b["number"] in liq else 0.0 for b in lb]
        h0 = [float(liq["SL(%d)" % b["number"]][0]) if "SL(%d)" % b["number"] in liq else 0.0 for b in lb]
        dt = self.parameterAsDouble(parameters, "TIME_STEP", context)
        duration = self.parameterAsDouble(parameters, "DURATION", context)
        if dt <= 0:
            cfl = courant_time_step(x, y, ik, self.parameterAsDouble(parameters, "DEPTH_GUESS", context),
                                    self.parameterAsDouble(parameters, "VELOCITY_GUESS", context))
            dt = fitting_time_step(cfl, duration)
            feedback.pushInfo(self.tr("Time step %g s: Courant number 1 on the smallest element for the typical depth and "
                                      "velocity is %.4g s; %g s divides the duration exactly") % (dt, cfl, dt))
        else:
            fit = fitting_time_step(dt, duration)
            if abs(fit - dt) > 1e-12:
                feedback.pushWarning(self.tr("Time step %g s does not divide the duration %g s: %g s is used, so the last step "
                                             "ends on the duration") % (dt, duration, fit))
                dt = fit
        if liq and "T" in liq and len(liq["T"]) and float(np.max(liq["T"])) < duration - 1e-3:
            feedback.pushWarning(self.tr("The liquid boundaries file ends at %g s, before the duration (%g s): TELEMAC stops "
                                         "with OUT OF RANGE when it gets there. Extend the hydrograph.") % (float(np.max(liq["T"])), duration))
        level = parameters.get("INITIAL_LEVEL")
        level = None if level in (None, "") else float(level)
        init = self.parameterAsFile(parameters, "INITIAL", context) or None
        depth0 = None
        if init and os.path.exists(init):
            try:
                with SelafinFile(init) as si:
                    v = _var(si, "water_depth")
                    if v is not None and si.header.npoin == len(x):
                        depth0 = si.read(v, 0)
            except Exception:  # noqa: BLE001 - only a check
                depth0 = None
        elif level is not None and _z is not None:
            depth0 = np.maximum(level - np.asarray(_z, dtype=np.float64), 0.0)
        if depth0 is not None:
            for b, k in zip(lb, kinds):
                nodes = np.asarray(b["nodes"], dtype=np.int64)
                if k == "Q" and len(nodes) and float(np.max(depth0[nodes])) <= 1e-6:
                    feedback.pushWarning(self.tr("Discharge boundary %d (%s) is dry in the initial state: TELEMAC stops at the "
                                                 "first step (DEBIMP: problem on boundary number %d). Start with water there: "
                                                 "a higher initial level, or the initial conditions tool's 'depth above the "
                                                 "bed of the channel axes'.")
                                         % (b["number"], table.comments[b["start"]] or "", b["number"]))
        start = self.parameterAsDateTime(parameters, "START", context)
        st = None
        if start is not None and start.isValid():
            p = start.toPyDateTime()
            st = (p.year, p.month, p.day, p.hour, p.minute, p.second)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        law = self.LAWS[self.parameterAsEnum(parameters, "FRICTION_LAW", context)].lower()
        text = steering_template(
            out, mesh, cli, len(lb), liquid_file=liq_path, stage_file=stage_path, stage_boundaries=sorted(stage),
            prescribed_flowrates=q0, prescribed_elevations=h0, initial_file=self.parameterAsFile(parameters, "INITIAL", context) or None,
            initial_level=level, duration=self.parameterAsDouble(parameters, "DURATION", context), time_step=dt,
            output_period=self.parameterAsDouble(parameters, "OUTPUT_PERIOD", context), friction_law=law,
            friction=self.parameterAsDouble(parameters, "FRICTION", context), friction_in_geometry=fr is not None,
            title=self.parameterAsString(parameters, "TITLE", context), start=st, kinds=kinds,
            notes=["Liquid boundaries in TELEMAC order: " + ", ".join(
                "%d %s (%s)" % (b["number"], table.comments[b["start"]] or "", k) for b, k in zip(lb, kinds))])
        with open(out, "w", encoding="ascii", errors="replace", newline="\n") as f:
            f.write(text)
        if fr is not None:
            feedback.pushInfo(self.tr("The mesh has a BOTTOM FRICTION variable: TELEMAC reads the coefficients from it"))
        feedback.pushInfo(self.tr("Run it with: telemac2d.py %s") % os.path.basename(out))
        return {"OUTPUT": out}


BC_ALGORITHMS = [StageDischargeAlgorithm, LiquidBoundariesFileAlgorithm, CheckBoundariesAlgorithm, InitialConditionsAlgorithm,
                 SteeringFileAlgorithm]
__all__ = ["BC_ALGORITHMS", "load_mesh", "FRICTION_LAWS"]
