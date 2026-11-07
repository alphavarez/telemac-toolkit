"""Processing algorithms of the mesher (group "Mesh")."""

from __future__ import annotations

import math
import os

import numpy as np
from qgis.core import (
    QgsProcessingOutputNumber,
    Qgis,
    QgsCoordinateReferenceSystem,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsGeometry,
    QgsMeshLayer,
    QgsPointXY,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterNumber,
    QgsProcessingParameterPoint,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsProject,
    QgsVectorLayer,
)

from ..mesher.meshdata import Mesh
from ..mesher.regular import regular_grid_mesh
from ..mesher.writers import FORMATS
from ..mesher.interp import channel_formula, formula_names
from ..qgis_bridge import meshing as M
from ..log import log_exception
from ..qgis_bridge.crs import show_layer, sidecar_crs
from ..telemac_io import SelafinError, mesh_statistics
from ..telemac_io.variables import canonical
from .algorithms import SLF_FILTER, _Base, _field

GEO_FILTER = "SELAFIN geometry (*.slf);;All files (*.*)"
CLI_FILTER = "Boundary conditions (*.cli);;All files (*.*)"
LIQ_FILTER = "Liquid boundaries (*.liq);;All files (*.*)"
GPKG_FILTER = "GeoPackage (*.gpkg)"
MESH_FILTER = "Meshes (*.slf *.2dm *.nc *.mesh *.msh *.geo *.res);;All files (*.*)"
FORMAT_KEYS = list(FORMATS)
OUT_FILTER = ";;".join("%s (*%s)" % (f.label, f.ext) for f in FORMATS.values()) + ";;All files (*.*)"
MDAL_READABLE = (".slf", ".2dm", ".nc", ".mesh")
FRAME_LABELS = ["Map coordinates as drawn", "Shift the lower-left corner to (0, 0)", "Shift a chosen point to (0, 0)"]

_DOUBLE = Qgis.ProcessingNumberParameterType.Double
#: meshing methods, in the order of the Method list (M.METHOD_KEYS)
METHOD_KEYS = list(M.METHOD_KEYS)
METHOD_LABELS = ["Balanced (built-in: regular nodes, smoothing)",
                 "Quality-guaranteed (built-in: refine until every angle reaches the target)",
                 "Gmsh Delaunay (fastest on large meshes)",
                 "Gmsh Frontal-Delaunay (highest element quality)",
                 "Gmsh MeshAdapt (most robust)",
                 "Gmsh quadrilaterals (for programs that read them; TELEMAC gets two triangles each)"]
_INT = Qgis.ProcessingNumberParameterType.Integer


def _advanced(p):
    p.setFlags(p.flags() | Qgis.ProcessingParameterFlag.Advanced)
    return p


class _MeshBase(_Base):
    """Common helpers: load what was written into the project when the algorithm ends."""

    def group(self):
        return self.tr("Mesh")

    def groupId(self):
        return "mesh"

    def __init__(self):
        super().__init__()
        self._to_load = []

    def _queue_load(self, path, crs, name=None):
        self._to_load.append((path, crs, name or os.path.splitext(os.path.basename(path))[0]))

    def postProcessAlgorithm(self, context, feedback):
        project = context.project() or QgsProject.instance()
        for path, crs, name in self._to_load:
            ext = os.path.splitext(path)[1].lower()
            if ext == ".gpkg":
                group = project.layerTreeRoot().insertGroup(0, name)
                for lname in ("perimeter", "breaklines", "refinement_regions", "bc_lines"):
                    v = QgsVectorLayer("%s|layername=%s" % (path, lname), lname, "ogr")
                    if v.isValid():
                        project.addMapLayer(v, False)
                        group.addLayer(v)
                continue
            if ext not in MDAL_READABLE:
                continue
            lyr = QgsMeshLayer(path, name, "mdal")
            if not lyr.isValid():
                feedback.reportError(self.tr("Could not load %s as a mesh layer") % path, fatalError=False)
                continue
            if crs is not None and crs.isValid():
                lyr.setCrs(crs)
            M.style_geometry_mesh(lyr)
            first = len(project.mapLayers()) == 0
            project.addMapLayer(lyr)
            try:
                from qgis.utils import iface

                if iface is not None and project is QgsProject.instance():
                    show_layer(lyr, iface.mapCanvas(), project, first=first)
            except (ImportError, RuntimeError) as exc:  # no QGIS window (scripts, tests)
                log_exception("Zoom to %s" % name, exc)
        self._to_load = []
        return {}

    def _mesh_crs(self, parameters, context, path, key="MESH_CRS"):
        crs = self.parameterAsCrs(parameters, key, context)
        if crs is not None and crs.isValid():
            return crs
        crs = sidecar_crs(path)
        if crs is not None and crs.isValid():
            return crs
        try:
            from ..qgis_bridge.crs import suggest_crs
            from ..qgis_bridge.results import file_extent

            choice = suggest_crs(path, file_extent(path), context.project() or QgsProject.instance())
            if choice.crs is not None and choice.crs.isValid():
                return choice.crs
        except (OSError, ValueError, SelafinError) as exc:  # not a SELAFIN file we can read: fall back to the project CRS
            log_exception("CRS of %s" % os.path.basename(path), exc)
        return context.project().crs() if context.project() else QgsCoordinateReferenceSystem()

    def _progress(self, feedback, lo=0.0, hi=100.0):
        def cb(frac, msg=""):
            feedback.setProgress(lo + (hi - lo) * float(frac))
            if msg:
                feedback.setProgressText(msg)
            return feedback.isCanceled()

        return cb

    # export parameters shared by every algorithm that writes a mesh ------------
    def _add_export_params(self, default_format=0, quads=True):
        self.addParameter(QgsProcessingParameterEnum("FORMAT", self.tr("Output format"),
                                                     [FORMATS[k].label for k in FORMAT_KEYS], defaultValue=default_format))
        if quads:
            self.addParameter(QgsProcessingParameterBoolean(
                "QUADS", self.tr("Keep structured cells as quadrilaterals (formats that accept them)"), False))
        self.addParameter(_advanced(QgsProcessingParameterEnum("FRAME", self.tr("Coordinates in the file"),
                                                               [self.tr(t) for t in FRAME_LABELS], defaultValue=0)))
        self.addParameter(_advanced(QgsProcessingParameterPoint("FRAME_POINT", self.tr("Point moved to (0, 0)"), optional=True)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("FRAME_ANGLE", self.tr("Rotation of the local axes (degrees, counter-clockwise)"),
                                                                 _DOUBLE, 0.0)))
        self.addParameter(_advanced(QgsProcessingParameterBoolean(
            "STORE_OFFSET", self.tr("TELEMAC: store the shift in the file header, so QGIS still shows the mesh in place"), True)))
        self.addParameter(_advanced(QgsProcessingParameterBoolean(
            "DOUBLE", self.tr("TELEMAC: double precision (SERAFIND; the steering file then needs "
                              "GEOMETRY FILE FORMAT = 'SERAFIND', which the steering file tool writes)"), False)))

    def _export_options(self, parameters, context, crs):
        fmt = FORMAT_KEYS[self.parameterAsEnum(parameters, "FORMAT", context)]
        mode = M.FRAME_MODES[self.parameterAsEnum(parameters, "FRAME", context)]
        pt = (0.0, 0.0)
        if mode == "point" and parameters.get("FRAME_POINT") not in (None, ""):
            p = self.parameterAsPoint(parameters, "FRAME_POINT", context, crs)
            pt = (p.x(), p.y())
        quads = self.parameterAsBoolean(parameters, "QUADS", context) if self.parameterDefinition("QUADS") else False
        return M.ExportOptions(fmt=fmt, quads=quads, frame=mode, frame_point=pt,
                               angle=self.parameterAsDouble(parameters, "FRAME_ANGLE", context),
                               store_offset=self.parameterAsBoolean(parameters, "STORE_OFFSET", context),
                               double=bool(self.parameterDefinition("DOUBLE") and self.parameterAsBoolean(parameters, "DOUBLE", context)))


# --------------------------------------------------------------------------
class CreateMeshProjectAlgorithm(_MeshBase):
    def name(self):
        return "createmeshproject"

    def displayName(self):
        return self.tr("New mesh project")

    def shortHelpString(self):
        return self.tr(
            "Create a GeoPackage with one empty, styled layer per kind of object and add it to the map:\n"
            "• domain – polygons to mesh; holes are islands; an optional size per polygon\n"
            "• islands – polygons left out of the mesh\n"
            "• refinement_areas, refinement_points – target element size\n"
            "• constraint_lines – lines the mesh follows: mode 'breakline' keeps every vertex (dikes, roads, "
            "banks); mode 'guide' is resampled to the local size\n"
            "• fixed_points – nodes that must exist (gauges, structures)\n"
            "• channels – left bank, right bank and an optional axis (role field) of a channel, meshed as a "
            "structured strip with nodes_across and spacing\n"
            "• grid_zones – rectangles meshed as regular grids (dx, dy)\n"
            "• roughness_zones – material id and friction value\n"
            "• boundaries – boundary lines with a type (wall, Q, H, QH, UV, UVH, free, waves or a code such as "
            "455), a name and an optional value\n"
            "Choose local coordinates for a laboratory flume or a schematic channel: the project then uses a "
            "metric grid with no geographic reference, so you can draw from (0, 0).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterCrs("CRS", self.tr("CRS (projected, in metres)"), defaultValue="ProjectCrs"))
        self.addParameter(QgsProcessingParameterBoolean("LOCAL", self.tr("Local coordinates (laboratory or schematic channel, no map CRS)"), False))
        self.addParameter(QgsProcessingParameterNumber("DEFAULT_SIZE", self.tr("Default element size (m)"), _DOUBLE, 0.0, optional=True, minValue=0.0))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh project"), GPKG_FILTER))

    def flags(self):
        return super().flags() | Qgis.ProcessingAlgorithmFlag.NoThreading

    def processAlgorithm(self, parameters, context, feedback):
        local = self.parameterAsBoolean(parameters, "LOCAL", context)
        crs = M.local_crs() if local else self.parameterAsCrs(parameters, "CRS", context)
        path = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        if crs.isGeographic():
            feedback.pushWarning(self.tr("The CRS is geographic: sizes would be in degrees. Use a projected CRS or local coordinates."))
        if os.path.exists(path):
            raise QgsProcessingException(self.tr("%s already exists; choose a new file name.") % path)
        size = self.parameterAsDouble(parameters, "DEFAULT_SIZE", context)
        try:
            M.create_mesh_project(path, crs, project=context.project(), settings={"default_size": "%g" % size if size else ""})
        except Exception as exc:
            raise QgsProcessingException(str(exc))
        if local and context.project() is not None and not context.project().crs().isValid():
            context.project().setCrs(crs)
        return {"OUTPUT": path}


# --------------------------------------------------------------------------
_EXAMPLE_SEGMENTS = """# an S-shaped flume widening downstream (one segment per line)
straight 3
bend left 2 90
straight 1
bend right 2 90
taper 1.5 2
straight 4"""


class ChannelTemplateAlgorithm(_MeshBase):
    SECTIONS = ["Rectangular (flat bed)", "Trapezoidal", "Compound (main channel and floodplains)", "V-shaped", "Parabolic",
                "Asymmetric trapezoid", "Custom formula of n and b"]
    PLANFORMS = ["Straight", "Bend (straight inlet, circular bend, straight outlet)", "Meanders (sine-generated)",
                 "Composed from segments (text below)", "Centreline from a line layer", "Imported file (TXT, CSV or DXF)"]
    TURNS = ["Left (counter-clockwise)", "Right (clockwise)"]
    IMPORT_AS = ["Centreline (the longest line; a third column gives widths)",
                 "Outline (the largest closed line is the domain, closed lines inside it islands, open lines breaklines)"]

    def name(self):
        return "channeltemplate"

    def displayName(self):
        return self.tr("Channel or flume template (straight, bend, meanders, composed, drawn, imported)")

    def shortHelpString(self):
        return self.tr(
            "Draw a laboratory flume or a schematic channel in a mesh project, with a bed you define by numbers "
            "instead of a DEM. Uses: flume experiments and their numerical replicas, verification against analytical "
            "solutions (normal depth, backwater, bend flow), sensitivity tests of mesh size, friction or turbulence, "
            "and idealised reaches for teaching or sediment studies. Natural rivers are drawn on the map instead "
            "(domain, channels and a DEM).\n\n"
            "Planform:\n"
            "• Straight: a rectangle of the given length.\n"
            "• Bend: straight inlet, circular bend (radius on the axis, angle such as 90 or 180 degrees, left or "
            "right), straight outlet: the classic curved-flume experiments.\n"
            "• Meanders: a sine-generated curve (Langbein & Leopold), the shape of natural meanders: heading "
            "theta(s) = theta0 sin(2 pi s / wavelength), with the wavelength measured along the axis.\n\n"
            "The channel starts at the origin ((0, 0) by default) heading along the direction angle. It adds the "
            "domain, an inlet line (Q) and an outlet line (H) and, with Structured, a structured strip whose "
            "cross-sections are normal to the axis, also in bends. The bed is saved as a formula in the project: "
            "of x and y for a straight channel, of the channel coordinates s (along the axis) and n (across it, "
            "left positive) otherwise, with a bed level at the inlet, a slope along the axis and a rectangular or "
            "trapezoidal section. Generate mesh applies it. Leave the project empty to create one in local "
            "coordinates.\n\n"
            "• Composed from segments: any channel built from a text, one segment per line: width W, straight L, "
            "bend left|right R ANGLE, elbow left|right ANGLE (a sharp corner), meanders WAVELENGTH THETA0 COUNT, "
            "kinoshita WAVELENGTH THETA0 COUNT JS JF, taper WIDTH LENGTH (gradual width change), expand WIDTH (sudden "
            "expansion or contraction), slope S, drop D (a step in the bed), bed Z, and section KIND ... (the "
            "cross-section from there on: rect; trapezoid BOTTOM_WIDTH SIDE_SLOPE BANK_HEIGHT [OFFSET]; compound "
            "MAIN_WIDTH SIDE_SLOPE FLOODPLAIN_HEIGHT FLOODPLAIN_SLOPE [OFFSET]; v DEPTH [OFFSET]; parabolic DEPTH "
            "[OFFSET]; asym BOTTOM_WIDTH LEFT_SLOPE RIGHT_SLOPE BANK_HEIGHT [OFFSET]; custom FORMULA of n and b). "
            "Pieces between sudden width changes and sharp corners get their own structured strips; the corners and "
            "steps get triangles.\n"
            "• Centreline from a line layer: a line you draw (in local coordinates for a flume), with a constant "
            "width or a width field; corners sharper than the corner angle become elbows.\n"
            "• Imported file: TXT or CSV (x y [width] per line; a line without numbers starts a new polyline) or "
            "DXF, read as a centreline or as an outline (a basin or any shape: the largest closed line is the domain).\n\n"
            "Cross-sections (n across the axis, left positive; b half the local width): rectangular; trapezoidal; "
            "compound (a main channel and floodplains at a height, optionally rising towards the walls); V-shaped and "
            "parabolic (depth at the axis); asymmetric trapezoid (different side slopes); or any formula of n and b. "
            "The offset moves the main channel to the left (+) or right (-) of the axis.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (empty = create one)"), extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("NEW_PROJECT", self.tr("New mesh project (when none is given)"),
                                                                GPKG_FILTER, optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterEnum("PLANFORM", self.tr("Planform"), [self.tr(t) for t in self.PLANFORMS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("WIDTH", self.tr("Width (m)"), _DOUBLE, 1.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("LENGTH", self.tr("Straight: length (m)"), _DOUBLE, 20.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("INLET", self.tr("Bend and meanders: straight inlet length (m)"), _DOUBLE, 5.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("OUTLET", self.tr("Bend and meanders: straight outlet length (m)"), _DOUBLE, 5.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("RADIUS", self.tr("Bend: radius on the axis (m)"), _DOUBLE, 2.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("BEND_ANGLE", self.tr("Bend: angle (degrees)"), _DOUBLE, 180.0, minValue=1.0, maxValue=360.0))
        self.addParameter(QgsProcessingParameterEnum("TURN", self.tr("Bend: turning"), [self.tr(t) for t in self.TURNS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("WAVELENGTH", self.tr("Meanders: wavelength along the axis (m)"), _DOUBLE, 10.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("DEFLECTION", self.tr("Meanders: largest deflection angle theta0 (degrees)"), _DOUBLE, 60.0,
                                                       minValue=1.0, maxValue=125.0))
        self.addParameter(QgsProcessingParameterNumber("COUNT", self.tr("Meanders: number of wavelengths"), _DOUBLE, 2.0, minValue=0.25))
        self.addParameter(QgsProcessingParameterPoint("ORIGIN", self.tr("Origin (inlet, centre of the inlet for curved channels)"),
                                                      defaultValue="0,0", optional=True))
        self.addParameter(QgsProcessingParameterNumber("ANGLE", self.tr("Direction at the inlet (degrees, counter-clockwise from x)"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterNumber("CELL", self.tr("Element size (m)"), _DOUBLE, 0.05, minValue=1e-6))
        self.addParameter(QgsProcessingParameterBoolean("STRUCTURED", self.tr("Structured mesh over the channel"), True))
        self.addParameter(QgsProcessingParameterNumber("BED_Z0", self.tr("Bed level at the inlet (m)"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterNumber("SLOPE", self.tr("Bed slope along the axis (m/m, positive downstream)"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterEnum("SECTION", self.tr("Cross-section"), self.SECTIONS, defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("BOTTOM_WIDTH", self.tr("Trapezoid: bottom width (m)"), _DOUBLE, 0.5, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("SIDE_SLOPE", self.tr("Trapezoid: side slope (horizontal per 1 vertical)"), _DOUBLE, 1.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("BANK_HEIGHT", self.tr("Trapezoid: bank height (m)"), _DOUBLE, 0.2, minValue=0.0))
        self.addParameter(QgsProcessingParameterString("SEGMENTS", self.tr("Composed: segments, one per line"), _EXAMPLE_SEGMENTS,
                                                       multiLine=True, optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("CENTRELINE", self.tr("Centreline: line layer (local coordinates "
                                                                                    "for a flume)"),
                                                              [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterField("WIDTH_FIELD", self.tr("Centreline: width field (optional; one width per "
                                                                             "line)"), parentLayerParameterName="CENTRELINE",
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterFile("IMPORT_FILE", self.tr("Imported file (TXT, CSV or DXF)"), optional=True,
                                                     fileFilter="Shapes (*.txt *.csv *.xyz *.dxf);;All files (*)"))
        self.addParameter(QgsProcessingParameterEnum("IMPORT_AS", self.tr("Imported file: read as"), [self.tr(t) for t in self.IMPORT_AS],
                                                     defaultValue=0))
        self.addParameter(_advanced(QgsProcessingParameterNumber("CORNER", self.tr("Drawn or imported: corners sharper than this are "
                                                                                    "elbows (degrees)"), _DOUBLE, 30.0, minValue=1.0,
                                                                 maxValue=179.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("FLOODPLAIN_SLOPE", self.tr("Compound: floodplain rise towards the "
                                                                                              "walls (m/m)"), _DOUBLE, 0.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("RIGHT_SLOPE", self.tr("Asymmetric: right side slope (horizontal "
                                                                                         "per 1 vertical; the left one is the "
                                                                                         "side slope)"), _DOUBLE, 2.0, minValue=0.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("SECTION_OFFSET", self.tr("Main channel offset from the axis (m, "
                                                                                            "left positive)"), _DOUBLE, 0.0)))
        self.addParameter(_advanced(QgsProcessingParameterString("SECTION_FORMULA", self.tr("Custom section: formula of n and b "
                                                                                             "(height above the axis bed)"),
                                                                 "0.1*(n/b)**2", optional=True)))
        self.addParameter(QgsProcessingParameterNumber("Q", self.tr("Inflow (m³/s, optional)"), _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterNumber("H", self.tr("Outlet water level (m, optional)"), _DOUBLE, optional=True))

    def flags(self):
        return super().flags() | Qgis.ProcessingAlgorithmFlag.NoThreading

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher import planform as P
        from ..mesher.interp import channel_formula_sn

        path = self.parameterAsFile(parameters, "PROJECT", context)
        if not path:
            path = self.parameterAsFileOutput(parameters, "NEW_PROJECT", context)
            if not path:
                raise QgsProcessingException(self.tr("Give a mesh project or a file name for a new one."))
            if os.path.exists(path):
                raise QgsProcessingException(self.tr("%s already exists.") % path)
            new_crs = M.local_crs()
            src = self.parameterAsSource(parameters, "CENTRELINE", context) \
                if self.parameterAsEnum(parameters, "PLANFORM", context) == 4 else None
            if src is not None and src.sourceCrs().isValid():
                new_crs = src.sourceCrs()  # a drawn centreline stays where it is drawn
            layers = M.create_mesh_project(path, new_crs, project=context.project())
            if context.project() is not None and not context.project().crs().isValid():
                context.project().setCrs(new_crs)
        else:
            layers = M.open_mesh_project(path)
        crs = layers["domain"].crs()
        kind = self.parameterAsEnum(parameters, "PLANFORM", context)
        W = self.parameterAsDouble(parameters, "WIDTH", context)
        ang = self.parameterAsDouble(parameters, "ANGLE", context)
        cell = self.parameterAsDouble(parameters, "CELL", context)
        o = (0.0, 0.0)
        if parameters.get("ORIGIN") not in (None, ""):
            p = self.parameterAsPoint(parameters, "ORIGIN", context, crs)
            o = (p.x(), p.y())
        q = parameters.get("Q")
        h = parameters.get("H")
        q = None if q in (None, "") else float(q)
        h = None if h in (None, "") else float(h)
        structured = self.parameterAsBoolean(parameters, "STRUCTURED", context)
        sec_kind = self.parameterAsEnum(parameters, "SECTION", context)
        trap = sec_kind == 1
        z0 = self.parameterAsDouble(parameters, "BED_Z0", context)
        S = self.parameterAsDouble(parameters, "SLOPE", context)
        bw = self.parameterAsDouble(parameters, "BOTTOM_WIDTH", context) if trap else 0.0
        ss = self.parameterAsDouble(parameters, "SIDE_SLOPE", context) if trap else 0.0
        hb = self.parameterAsDouble(parameters, "BANK_HEIGHT", context) if trap else 0.0
        settings = {"default_size": "%g" % cell}
        if kind >= 3 or sec_kind >= 2:
            return self._composed(parameters, context, feedback, path, layers, crs, kind, W, o, ang, cell, structured, q, h, z0, S,
                                  settings)
        if kind == 0:
            L = self.parameterAsDouble(parameters, "LENGTH", context)
            rep = M.add_flume(layers, L, W, o, ang, cell, structured, q, h)
            a = math.radians(ang)
            ax0, ay0 = o[0] - 0.5 * W * math.sin(a), o[1] + 0.5 * W * math.cos(a)  # axis starts mid-inlet
            formula = channel_formula(z0, S, bw, ss, hb, ax0, ay0, ang)
            settings.update(bottom_formula=formula, bottom_axis="")
        else:
            inlet = self.parameterAsDouble(parameters, "INLET", context)
            outlet = self.parameterAsDouble(parameters, "OUTLET", context)
            try:
                if kind == 1:
                    pf = P.bend(inlet, self.parameterAsDouble(parameters, "RADIUS", context),
                                self.parameterAsDouble(parameters, "BEND_ANGLE", context), outlet, W, cell,
                                "left" if self.parameterAsEnum(parameters, "TURN", context) == 0 else "right")
                else:
                    pf = P.meanders(self.parameterAsDouble(parameters, "WAVELENGTH", context),
                                    self.parameterAsDouble(parameters, "DEFLECTION", context),
                                    self.parameterAsDouble(parameters, "COUNT", context), W, cell, inlet, outlet)
            except ValueError as exc:
                raise QgsProcessingException(str(exc))
            pf = pf.transformed(o, ang)
            rep = M.add_planform(layers, pf, cell, structured, q, h, name="flume")
            if not structured:
                # the bed formula needs the axis: keep it in the channels layer as a lone axis line
                lyr = layers["channels"]
                f = QgsFeature(lyr.fields())
                f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(float(x), float(y)) for x, y in pf.axis]))
                f.setAttribute("channel", "flume_axis")
                f.setAttribute("role", "axis")
                lyr.startEditing()
                lyr.addFeature(f)
                lyr.commitChanges()
                settings["bottom_axis"] = "flume_axis"
            else:
                settings["bottom_axis"] = "flume"
            formula = channel_formula_sn(z0, S, bw, ss, hb)
            settings["bottom_formula"] = formula
            if kind == 2:
                feedback.pushInfo(self.tr("Tightest bend radius %.3g m (wavelength / (2 pi theta0))") % P.min_radius(
                    self.parameterAsDouble(parameters, "WAVELENGTH", context), self.parameterAsDouble(parameters, "DEFLECTION", context)))
            z_out = z0 - S * pf.length
            feedback.pushInfo(self.tr("Axis length %.4g m; bed at the outlet %.4g m") % (pf.length, z_out))
        M.write_settings(path, settings)
        for r in rep:
            feedback.pushInfo(r)
        feedback.pushInfo(self.tr("Bed formula saved in the project: z = %s") % formula)
        M.refresh_project_layers(path, context.project())
        return {"OUTPUT": path, "BED_FORMULA": formula}

    def _section(self, parameters, context):
        from ..mesher.composer import Section

        k = self.parameterAsEnum(parameters, "SECTION", context)
        num = lambda key: self.parameterAsDouble(parameters, key, context)  # noqa: E731
        off = num("SECTION_OFFSET")
        if k == 1:
            return Section("trapezoid", (num("BOTTOM_WIDTH"), num("SIDE_SLOPE"), num("BANK_HEIGHT"), off))
        if k == 2:
            return Section("compound", (num("BOTTOM_WIDTH"), num("SIDE_SLOPE"), num("BANK_HEIGHT"), num("FLOODPLAIN_SLOPE"), off))
        if k == 3:
            return Section("v", (num("BANK_HEIGHT"), off))
        if k == 4:
            return Section("parabolic", (num("BANK_HEIGHT"), off))
        if k == 5:
            return Section("asym", (num("BOTTOM_WIDTH"), num("SIDE_SLOPE"), num("RIGHT_SLOPE"), num("BANK_HEIGHT"), off))
        if k == 6:
            from ..mesher.interp import _check_formula

            f = self.parameterAsString(parameters, "SECTION_FORMULA", context).strip() or "0"
            try:
                _check_formula(f)
            except (ValueError, SyntaxError) as exc:
                raise QgsProcessingException(self.tr("Custom section: %s") % exc)
            return Section("custom", (), f)
        return Section()

    def _composed(self, parameters, context, feedback, path, layers, crs, kind, W, o, ang, cell, structured, q, h, z0, S, settings):
        """Composed, drawn and imported channels, and the planforms of the first three kinds with the newer sections."""
        from ..mesher import composer as C

        sec = self._section(parameters, context)
        corner = self.parameterAsDouble(parameters, "CORNER", context)
        try:
            if kind == 0:
                comp = C.compose("straight %r" % self.parameterAsDouble(parameters, "LENGTH", context), W, cell, z0, S, sec)
            elif kind == 1:
                turn = "left" if self.parameterAsEnum(parameters, "TURN", context) == 0 else "right"
                comp = C.compose("straight %r; bend %s %r %r; straight %r" % (
                    self.parameterAsDouble(parameters, "INLET", context), turn, self.parameterAsDouble(parameters, "RADIUS", context),
                    self.parameterAsDouble(parameters, "BEND_ANGLE", context), self.parameterAsDouble(parameters, "OUTLET", context)),
                    W, cell, z0, S, sec)
            elif kind == 2:
                comp = C.compose("straight %r; meanders %r %r %r; straight %r" % (
                    self.parameterAsDouble(parameters, "INLET", context), self.parameterAsDouble(parameters, "WAVELENGTH", context),
                    self.parameterAsDouble(parameters, "DEFLECTION", context), self.parameterAsDouble(parameters, "COUNT", context),
                    self.parameterAsDouble(parameters, "OUTLET", context)), W, cell, z0, S, sec)
            elif kind == 3:
                comp = C.compose(self.parameterAsString(parameters, "SEGMENTS", context), W, cell, z0, S, sec)
            elif kind == 4:
                src = self.parameterAsSource(parameters, "CENTRELINE", context)
                if src is None:
                    raise QgsProcessingException(self.tr("Give the centreline layer."))
                items = M.line_arrays_with(src, crs, context.project(), self.parameterAsString(parameters, "WIDTH_FIELD", context) or None)
                if not items:
                    raise QgsProcessingException(self.tr("The centreline layer has no line."))
                line, wv = max(items, key=lambda it: len(it[0]))
                comp = C.from_centreline(line, cell, wv if wv and wv > 0 else W, None, corner, z0, S, sec)
                o, ang = (0.0, 0.0), 0.0  # drawn where it is
            else:
                fpath = self.parameterAsFile(parameters, "IMPORT_FILE", context)
                if not fpath or not os.path.exists(fpath):
                    raise QgsProcessingException(self.tr("Give the file to import."))
                closed, open_ = M.read_shape_file(fpath)
                if self.parameterAsEnum(parameters, "IMPORT_AS", context) == 1:
                    rep = M.add_outline(layers, [c[:, :2] for c in closed], [c[:, :2] for c in open_], cell)
                    xmin = min(float(c[:, 0].min()) for c in closed)
                    formula = "%g - %g*(x - %g)" % (z0, S, xmin) if S else "%g" % z0
                    settings.update(bottom_formula=formula, bottom_axis="", bottom_halfwidth="")
                    return self._finish(path, settings, rep, formula, feedback, context)
                lines = open_ + closed
                line = max(lines, key=len)
                widths = line[:, 2] if line.shape[1] > 2 and np.all(np.isfinite(line[:, 2])) else None
                comp = C.from_centreline(line[:, :2], cell, W, widths, corner, z0, S, sec)
                o, ang = (0.0, 0.0), 0.0
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        if kind == 0:  # a straight flume's origin is its inlet corner; the axis starts mid-inlet
            a = math.radians(ang)
            o = (o[0] - 0.5 * W * math.sin(a), o[1] + 0.5 * W * math.cos(a))
        if kind <= 3:
            comp = comp.transformed(o, ang)
        rep = M.add_composition(layers, comp, cell, structured, q, h, name="flume")
        settings.update(bottom_formula=comp.formula, bottom_axis="flume_axis",
                        bottom_halfwidth=M.profile_text(*comp.half_width))
        feedback.pushInfo(self.tr("Axis length %.4g m") % comp.length)
        return self._finish(path, settings, rep, comp.formula, feedback, context)

    def _finish(self, path, settings, rep, formula, feedback, context):
        M.write_settings(path, settings)
        for r in rep:
            if r.startswith("Warning: "):
                feedback.reportError(r[9:], fatalError=False)
            else:
                feedback.pushInfo(r)
        feedback.pushInfo(self.tr("Bed formula saved in the project: z = %s") % formula)
        M.refresh_project_layers(path, context.project())
        return {"OUTPUT": path, "BED_FORMULA": formula}


class RiverReachAlgorithm(_MeshBase):
    SOURCES = ["Centreline", "Two bank lines"]
    WIDTHS = ["Constant width", "From the DEM (bankfull depth)"]
    DOMAINS = ["Create it: the reach plus a floodplain on each side", "Keep the domain of the project"]
    SPANS_IN = ["The channel only", "The whole end of the domain"]
    SPANS_OUT = ["The whole end of the domain", "The channel only"]

    def name(self):
        return "riverreach"

    def displayName(self):
        return self.tr("River reach or network from centrelines or bank lines")

    def shortHelpString(self):
        return self.tr(
            "Draws a natural river in a mesh project, ready for Generate mesh: the left bank, right bank and axis "
            "of a structured channel strip with matched cross-sections normal to the flow, the domain (or the "
            "project's own), and inflow (Q) and outflow (H) lines where the river meets the domain outline.\n\n"
            "Banks:\n"
            "• Centreline and a constant width.\n"
            "• Centreline and a DEM: at every cross-section the terrain across the river is read, the bed found "
            "near the centreline, and each bank put where the terrain first rises the bankfull depth above the bed "
            "(searched up to the largest width). The widths are smoothed along the river and the axis moved midway "
            "between the banks.\n"
            "• Two bank lines (one layer with two lines, in any order and direction).\n\n"
            "With a DEM the river is turned to run downstream, from the higher bed to the lower one; without one, "
            "draw the centreline downstream. Jagged lines are smoothed over the smoothing length (default: the "
            "width) and, in tight bends, the cross-sections are turned until they no longer cross; any that still "
            "cross are reported.\n\n"
            "Boundaries: the inflow covers the channel and the outflow the whole downstream end of the domain by "
            "default (water on the floodplain can leave). The lines stop short of their ends, so the bank and "
            "corner nodes stay walls. Values are optional; set them later in the boundary layer or the .liq file.\n\n"
            "River networks: give several centrelines, each tributary ending on the river it joins (or lines meeting "
            "end to end at a bifurcation). Every reach between junctions gets its own strip, stopped short of the "
            "junction; the channel area at the junction is meshed with triangles of the strips' size, with its banks "
            "as guide lines, so confluences keep their corners. Free upstream ends get inflow lines, downstream ends "
            "outflow lines (lines must run downstream, or be turned by the DEM). An optional width field gives each "
            "line its own constant width.\n\n"
            "Nodes across follow the width (default): where the river widens or narrows the number of nodes across "
            "changes by one row of triangles at a time, so the spacing across stays near the spacing along. In tight "
            "bends only the cross-sections there are turned, until none cross and no inner-bank cell is shorter than "
            "the spacing over the largest ratio; the others stay normal to the flow.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (empty = create one)"), extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("NEW_PROJECT", self.tr("New mesh project (when none is given)"),
                                                                GPKG_FILTER, optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterEnum("SOURCE", self.tr("River given as"), [self.tr(t) for t in self.SOURCES], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSource("CENTRELINE", self.tr("Centreline(s) (line layer; several lines make a "
                                                                                    "network; selected features if any)"),
                                                              [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterField("WIDTH_FIELD", self.tr("Width field of the centrelines (one constant width "
                                                                             "per line, optional)"),
                                                      parentLayerParameterName="CENTRELINE", type=Qgis.ProcessingFieldParameterDataType.Numeric,
                                                      optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("BANKS", self.tr("Bank lines (line layer with the two banks)"),
                                                              [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterEnum("WIDTH_MODE", self.tr("Width (centreline)"), [self.tr(t) for t in self.WIDTHS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("WIDTH", self.tr("Constant width (m)"), _DOUBLE, 20.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM (banks from the terrain, downstream direction)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("BANKFULL_DEPTH", self.tr("Bankfull depth above the bed (m)"), _DOUBLE, 1.5, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("MAX_WIDTH", self.tr("Largest width to search in the DEM (m)"), _DOUBLE, 100.0, minValue=1e-6))
        self.addParameter(_advanced(QgsProcessingParameterNumber("MIN_WIDTH", self.tr("Smallest width (m)"), _DOUBLE, 0.0, minValue=0.0)))
        self.addParameter(QgsProcessingParameterNumber("SPACING", self.tr("Cross-section spacing along the river (m)"), _DOUBLE, 10.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("NODES_ACROSS", self.tr("Nodes across the channel (0 = width / spacing, odd)"),
                                                       _INT, 0, minValue=0))
        self.addParameter(QgsProcessingParameterBoolean("VARY_ACROSS", self.tr("Nodes across follow the width (when 0 above)"), True))
        self.addParameter(_advanced(QgsProcessingParameterNumber(
            "MAX_RATIO", self.tr("Bends: inner-bank cells at least the spacing divided by"), _DOUBLE, 4.0, minValue=1.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("SMOOTHING", self.tr("Smoothing length of the centreline (m, 0 = the width)"),
                                                                 _DOUBLE, 0.0, minValue=0.0)))
        self.addParameter(QgsProcessingParameterEnum("DOMAIN_MODE", self.tr("Domain"), [self.tr(t) for t in self.DOMAINS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("FLOODPLAIN", self.tr("Floodplain width on each side (m)"), _DOUBLE, 200.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("FLOODPLAIN_SIZE", self.tr("Element size on the floodplain (m, 0 = default size)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean("BOUNDARIES", self.tr("Add inflow and outflow lines"), True))
        self.addParameter(QgsProcessingParameterEnum("INFLOW_SPAN", self.tr("Inflow (Q) covers"), [self.tr(t) for t in self.SPANS_IN], defaultValue=0))
        self.addParameter(QgsProcessingParameterEnum("OUTFLOW_SPAN", self.tr("Outflow (H) covers"), [self.tr(t) for t in self.SPANS_OUT], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("Q", self.tr("Inflow (m³/s, optional)"), _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterNumber("H", self.tr("Outflow water level (m, optional)"), _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterString("NAME", self.tr("Channel name"), "river"))

    def flags(self):
        return super().flags() | Qgis.ProcessingAlgorithmFlag.NoThreading

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher import river as RV

        from_banks = self.parameterAsEnum(parameters, "SOURCE", context) == 1
        src = self.parameterAsSource(parameters, "BANKS" if from_banks else "CENTRELINE", context)
        if src is None:
            raise QgsProcessingException(self.tr("Give the bank lines.") if from_banks else self.tr("Give the centreline."))
        path = self.parameterAsFile(parameters, "PROJECT", context)
        if not path:
            path = self.parameterAsFileOutput(parameters, "NEW_PROJECT", context)
            if not path:
                raise QgsProcessingException(self.tr("Give a mesh project or a file name for a new one."))
            if os.path.exists(path):
                raise QgsProcessingException(self.tr("%s already exists.") % path)
            layers = M.create_mesh_project(path, src.sourceCrs(), project=context.project())
        else:
            layers = M.open_mesh_project(path)
        crs = layers["domain"].crs()
        wfield = "" if from_banks else self.parameterAsString(parameters, "WIDTH_FIELD", context)
        items = M.line_arrays_with(src, crs, context.project(), wfield or None)
        lines = [c for c, _v in items]
        dem_layer = self.parameterAsRasterLayer(parameters, "DEM", context)
        dem = None
        if dem_layer is not None:
            sampler = M.RasterSampler(dem_layer, 1, crs, project=context.project())
            dem = sampler.sample
        step = self.parameterAsDouble(parameters, "SPACING", context)
        sm = self.parameterAsDouble(parameters, "SMOOTHING", context) or None
        ratio = self.parameterAsDouble(parameters, "MAX_RATIO", context) or 4.0
        vary = self.parameterAsBoolean(parameters, "VARY_ACROSS", context)
        q = parameters.get("Q")
        h = parameters.get("H")
        common = dict(
            boundaries=self.parameterAsBoolean(parameters, "BOUNDARIES", context),
            inflow_span=("channel", "valley")[self.parameterAsEnum(parameters, "INFLOW_SPAN", context)],
            outflow_span=("valley", "channel")[self.parameterAsEnum(parameters, "OUTFLOW_SPAN", context)],
            inflow=None if q in (None, "") else float(q), outflow_level=None if h in (None, "") else float(h),
            name=self.parameterAsString(parameters, "NAME", context).strip() or "river", vary_across=vary)
        nodes = self.parameterAsInt(parameters, "NODES_ACROSS", context)
        fp = self.parameterAsDouble(parameters, "FLOODPLAIN", context)
        fp_size = self.parameterAsDouble(parameters, "FLOODPLAIN_SIZE", context) or None
        try:
            if from_banks:
                if len(lines) != 2:
                    raise QgsProcessingException(self.tr("The bank layer must hold exactly two lines (it has %d); select "
                                                         "the two banks.") % len(lines))
                reach = RV.reach_from_banks(lines[0], lines[1], step, dem=dem, smoothing=sm, max_ratio=ratio)
            else:
                if not lines:
                    raise QgsProcessingException(self.tr("The centreline layer has no line."))
                by_dem = self.parameterAsEnum(parameters, "WIDTH_MODE", context) == 1
                if by_dem and dem is None:
                    raise QgsProcessingException(self.tr("Banks from the DEM need a DEM."))
                width = self.parameterAsDouble(parameters, "WIDTH", context)
                widths = [v if v and v > 0 else (None if by_dem else width) for _c, v in items]
                bankfull = self.parameterAsDouble(parameters, "BANKFULL_DEPTH", context)
                max_w = self.parameterAsDouble(parameters, "MAX_WIDTH", context)
                min_w = self.parameterAsDouble(parameters, "MIN_WIDTH", context)
                if len(lines) > 1:
                    if self.parameterAsEnum(parameters, "DOMAIN_MODE", context) == 1:
                        raise QgsProcessingException(self.tr("A river network makes its own domain: choose 'Create it' for "
                                                             "the domain, or draw the reaches one at a time."))
                    from ..mesher import network as NW
                    from ..mesher.river import _downstream

                    if dem is not None:  # each line runs downstream
                        for i, c in enumerate(lines):
                            rev, z0, z1 = _downstream(RV.resample(c, step), dem)
                            if rev:
                                lines[i] = c[::-1].copy()
                                feedback.pushInfo(self.tr("Line %d reversed to run downstream (bed %.3f m to %.3f m)") % (i + 1, z0, z1))
                    net = NW.river_network(lines, step, widths, dem=dem, bankfull_depth=bankfull if by_dem else None,
                                           max_width=max_w, min_width=min_w, smoothing=sm, max_ratio=ratio)
                    rep = M.add_river_network(layers, net, step, nodes, fp, fp_size, **common)
                    reach = None
                else:
                    reach = RV.reach_from_centreline(
                        lines[0], step, width=widths[0], dem=dem,
                        bankfull_depth=bankfull if widths[0] is None else None, max_width=max_w, min_width=min_w,
                        smoothing=sm, max_ratio=ratio)
            if reach is not None or from_banks:
                rep = M.add_river(layers, reach, step, nodes, fp, fp_size,
                                  use_domain=self.parameterAsEnum(parameters, "DOMAIN_MODE", context) == 1,
                                  project=context.project(), **common)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        for r in rep:
            if r.startswith("Warning: "):
                feedback.reportError(r[9:], fatalError=False)
            else:
                feedback.pushInfo(r)
        st = M.read_settings(path)
        if not st.get("default_size"):
            fp = self.parameterAsDouble(parameters, "FLOODPLAIN_SIZE", context)
            M.write_settings(path, {"default_size": "%g" % (fp or 4 * step)})
            feedback.pushInfo(self.tr("Default element size set to %g m (change it in step 1 of the panel)") % (fp or 4 * step))
        if dem_layer is not None:
            feedback.pushInfo(self.tr("Tick the DEM %s in step 3 of the panel so Generate mesh takes the bottom from it.") % dem_layer.name())
        M.refresh_project_layers(path, context.project())
        return {"OUTPUT": path}


# --------------------------------------------------------------------------
_INPUTS = [  # parameter, project layer key, label, geometry
    ("DOMAIN", "domain", "Domain (polygons; holes are islands)", "polygon"),
    ("ISLANDS", "islands", "Islands", "polygon"),
    ("REFINEMENT_AREAS", "refinement_areas", "Refinement areas (size)", "polygon"),
    ("REFINEMENT_POINTS", "refinement_points", "Refinement points (size)", "point"),
    ("CONSTRAINT_LINES", "constraint_lines", "Constraint lines (mode: breakline or guide)", "line"),
    ("FIXED_POINTS", "fixed_points", "Fixed points", "point"),
    ("CHANNELS", "channels", "Structured channels (role: left, right, axis)", "line"),
    ("GRID_ZONES", "grid_zones", "Structured grid zones (dx, dy)", "polygon"),
    ("ROUGHNESS_ZONES", "roughness_zones", "Roughness zones (material, friction)", "polygon"),
    ("BOUNDARIES", "boundaries", "Boundary condition lines (type, name, value)", "line"),
]
_GEOM = {"polygon": Qgis.ProcessingSourceType.VectorPolygon, "line": Qgis.ProcessingSourceType.VectorLine,
         "point": Qgis.ProcessingSourceType.VectorPoint}


class GenerateMeshAlgorithm(_MeshBase):
    def name(self):
        return "generatemesh"

    def _free_output(self, out, context, feedback):
        """A mesh file loaded in QGIS cannot be overwritten on Windows: remove its layer (on the main thread; the new
        mesh is loaded instead), or write beside it."""
        from qgis.PyQt.QtCore import QCoreApplication, QThread

        if not M.layers_using([out], context.project()):
            return out
        app = QCoreApplication.instance()
        if app is not None and QThread.currentThread() == app.thread():
            names = M.free_output([out], context.project())
            feedback.pushInfo(self.tr("Removed the layer(s) %s, which read %s, so it can be overwritten")
                              % (", ".join(names), os.path.basename(out)))
            return out
        alt = M.free_name(out)
        feedback.reportError(self.tr("%s is loaded in QGIS and cannot be overwritten: the mesh is written to %s. Remove "
                                     "that layer to overwrite it next time.") % (os.path.basename(out), os.path.basename(alt)),
                             fatalError=False)
        return alt

    def displayName(self):
        return self.tr("Generate mesh")

    def shortHelpString(self):
        return self.tr(
            "Build a mesh from a mesh project (or from individual layers, which take precedence) and write it for "
            "TELEMAC or another program.\n\n"
            "Unstructured triangles fill the domain with an element size that grows away from refinement areas, "
            "points and lines at the growth rate; the mesh follows breaklines vertex by vertex and guide lines at "
            "the local size, keeps fixed points, and stitches in structured channel strips and regular grid zones "
            "node for node.\n\n"
            "Bottom: rasters in priority order, then survey points, then a formula of x and y (for flumes); nodes "
            "nothing covers take the nearest value. Roughness zones give a material id per element and a friction "
            "value per node (BOTTOM FRICTION for TELEMAC). Boundary lines give the TELEMAC .cli (and a .liq "
            "template), 2DM nodestrings, Gmsh physical lines or ADCIRC/MIKE boundary lists.\n\n"
            "Advanced: write the file in local coordinates (lower-left corner or a chosen point at (0, 0), "
            "optional rotation); TELEMAC files can keep the shift in their header.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project"), extension="gpkg", optional=True))
        for key, _k, label, geom in _INPUTS:
            self.addParameter(QgsProcessingParameterFeatureSource(key, self.tr(label), [_GEOM[geom]], optional=True))
        self.addParameter(QgsProcessingParameterNumber("DEFAULT_SIZE", self.tr("Default element size (0 = from the project)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("GROWTH", self.tr("Size growth rate"), _DOUBLE, 1.2, minValue=1.0, maxValue=3.0))
        self.addParameter(_advanced(QgsProcessingParameterNumber("SMOOTHING", self.tr("Smoothing iterations"), _INT, 8, minValue=0, maxValue=60)))
        self.addParameter(_advanced(QgsProcessingParameterNumber(
            "OUTLINE_TOLERANCE", self.tr("Follow the outline and lines within (m, 0 = 5% of the element size)"), _DOUBLE, 0.0,
            minValue=0.0)))
        self.addParameter(QgsProcessingParameterEnum("METHOD", self.tr("Method (Gmsh ones need the gmsh package)"),
                                                     [self.tr(t) for t in METHOD_LABELS], defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber("MIN_ANGLE", self.tr("Quality-guaranteed: smallest angle (degrees, 20 to 30)"),
                                                       _DOUBLE, 30.0, minValue=20.0, maxValue=30.0))
        self.addParameter(QgsProcessingParameterNumber(
            "TERRAIN_TOLERANCE", self.tr("Terrain refinement: bed tolerance (m, 0 = off; uses the first bottom raster)"),
            _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(_advanced(QgsProcessingParameterNumber(
            "TERRAIN_MIN_SIZE", self.tr("Terrain refinement: smallest element size (m, 0 = default size / 15)"), _DOUBLE, 0.0,
            minValue=0.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber(
            "MAX_NODES", self.tr("Node limit (0 = no limit)"), _INT, 2000000, minValue=0, maxValue=2000000000)))
        self.addParameter(QgsProcessingParameterMultipleLayers("RASTERS", self.tr("Bottom: rasters (priority order)"),
                                                               Qgis.ProcessingSourceType.Raster, optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("POINTS", self.tr("Bottom: survey points"),
                                                              [Qgis.ProcessingSourceType.VectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterField("POINT_FIELD", self.tr("Survey value field (empty = point Z)"), parentLayerParameterName="POINTS",
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterString("FORMULA", self.tr("Bottom: formula of x, y or channel coordinates s, n (empty = from the project)"), "", optional=True))
        self.addParameter(_advanced(QgsProcessingParameterNumber("Z_FACTOR", self.tr("Elevation unit factor"), _DOUBLE, 1.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("FRICTION", self.tr("Friction outside the roughness zones"), _DOUBLE, optional=True)))
        self._add_export_params()
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh file"), OUT_FILTER))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT_LIQ", self.tr("TELEMAC liquid boundaries template (.liq)"), LIQ_FILTER,
                                                                optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the result when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        project_path = self.parameterAsFile(parameters, "PROJECT", context)
        layers, settings = {}, dict(M.SETTINGS_DEFAULTS)
        if project_path:
            layers = M.open_mesh_project(project_path)
            settings = M.read_settings(project_path)
        for key, lk, _label, _g in _INPUTS:
            src = self.parameterAsSource(parameters, key, context)
            if src is not None:
                layers[lk] = src
        if M.normalize_layers(layers).get("domain") is None:
            raise QgsProcessingException(self.tr("Give a mesh project or a domain layer."))
        dom = M.normalize_layers(layers)["domain"]
        crs = dom.sourceCrs() if hasattr(dom, "sourceCrs") else dom.crs()
        if crs.isValid() and crs.isGeographic():
            feedback.pushWarning(self.tr("The domain CRS is geographic: sizes are in degrees. Reproject to a projected CRS."))
        size = self.parameterAsDouble(parameters, "DEFAULT_SIZE", context)
        if not size:
            try:
                size = float(settings.get("default_size") or 0)
            except ValueError:
                size = 0.0
        if not size:
            size = M.auto_default_size(layers, context.project())
            if not size:
                raise QgsProcessingException(self.tr("Give a default element size."))
            feedback.pushInfo(self.tr("No default element size set: using %g (about 1/40 of the square root of the "
                                      "domain area). Set it in the project to change it.") % size)
        formula = self.parameterAsString(parameters, "FORMULA", context) or settings.get("bottom_formula", "")
        zf = self.parameterAsDouble(parameters, "Z_FACTOR", context)
        try:
            axis = half = None
            if formula_names(formula) & {"s", "n", "b"}:
                axis, half = M.formula_axis(layers, settings, crs, context.project())
            sources = M.bottom_sources_from(self.parameterAsLayerList(parameters, "RASTERS", context) or [],
                                            self.parameterAsSource(parameters, "POINTS", context),
                                            self.parameterAsString(parameters, "POINT_FIELD", context), formula, zf, crs,
                                            context.project(), axis, half)
        except (ValueError, SyntaxError) as exc:
            raise QgsProcessingException(self.tr("Bottom: %s") % exc)
        fr = parameters.get("FRICTION")
        method = METHOD_KEYS[self.parameterAsEnum(parameters, "METHOD", context)]
        terrain = None
        tol = self.parameterAsDouble(parameters, "TERRAIN_TOLERANCE", context)
        if tol > 0:
            rasters = self.parameterAsLayerList(parameters, "RASTERS", context) or []
            if not rasters:
                raise QgsProcessingException(self.tr("Terrain refinement needs a bottom raster (the first one of the list is used)."))
            terrain = (M.RasterSampler(rasters[0], 1, crs, project=context.project()).sample, tol,
                       self.parameterAsDouble(parameters, "TERRAIN_MIN_SIZE", context))
        if method.startswith("gmsh_"):
            from ..mesher import gmsh_engine

            if not gmsh_engine.available():
                raise QgsProcessingException(M.gmsh_missing_message())
        opts = self._export_options(parameters, context, crs)
        opts.liq_path = self.parameterAsFileOutput(parameters, "OUTPUT_LIQ", context) or None
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        out = self._free_output(out, context, feedback)
        try:
            res = M.generate_and_export(layers, out, opts, crs, size, self.parameterAsDouble(parameters, "GROWTH", context),
                                        self.parameterAsInt(parameters, "SMOOTHING", context), sources, True,
                                        None if fr in (None, "") else float(fr), None, self._progress(feedback), context.project(),
                                        self.parameterAsInt(parameters, "MAX_NODES", context),
                                        self.parameterAsDouble(parameters, "MIN_ANGLE", context) if method == "quality" else None,
                                        method, terrain, self.parameterAsDouble(parameters, "OUTLINE_TOLERANCE", context) or None)
        except RuntimeError as exc:
            if str(exc) == "cancelled":
                return {}
            raise QgsProcessingException(str(exc))
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        for w in res.warnings:
            feedback.pushWarning(w)
        for line in res.report:
            (feedback.pushWarning if "Warning" in line else feedback.pushInfo)(line.strip())
        for line in res.cas:
            feedback.pushInfo(self.tr("Steering file: ") + line)
        if not sources:
            feedback.pushWarning(self.tr("No bottom source: the bottom is 0 everywhere."))
        st = res.stats
        feedback.pushInfo(self.tr("Mesh: %s; edges %.3g to %.3g; angles %.1f to %.1f deg") % (
            res.mesh.summary(), st["edge_min"], st["edge_max"], st["angle_min"], st["angle_max"]))
        for f in res.files:
            feedback.pushInfo(self.tr("Written: %s") % f)
        if self.parameterAsBoolean(parameters, "LOAD", context) and res.files:
            local = opts.frame != "keep" and not (opts.fmt == "telemac" and opts.store_offset and not opts.angle)
            self._queue_load(res.path, M.local_crs() if local else res.crs)
        out = {"OUTPUT": res.path}
        cli = [f for f in res.files if f.endswith(".cli")]
        if cli:
            out["OUTPUT_CLI"] = cli[0]
        if opts.liq_path:
            out["OUTPUT_LIQ"] = opts.liq_path
        return out


# --------------------------------------------------------------------------
class ExportMeshAlgorithm(_MeshBase):
    def name(self):
        return "exportmesh"

    def displayName(self):
        return self.tr("Export or convert mesh")

    def shortHelpString(self):
        return self.tr(
            "Write an existing mesh (SELAFIN, 2DM, UGRID or any mesh QGIS reads) for another program: TELEMAC, "
            "SMS 2DM (BASEMENT, TUFLOW FV), Gmsh, GiD (Iber and other GiD-based programs), UGRID (Delft3D FM), "
            "MIKE 21 FM or ADCIRC/SCHISM. Optional boundary lines and roughness zones add boundary strings, "
            "material ids and friction. The file can be written in local coordinates (a point moved to (0, 0), "
            "optional rotation). HEC-RAS builds its own cells, so it is exported from Generate mesh, not here.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Mesh"), fileFilter=MESH_FILTER))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = from the file or the project)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("BOUNDARIES", self.tr("Boundary condition lines"), [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("ROUGHNESS_ZONES", self.tr("Roughness zones"), [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(_advanced(QgsProcessingParameterNumber("FRICTION", self.tr("Friction outside the roughness zones"), _DOUBLE, optional=True)))
        self._add_export_params(default_format=1)
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh file"), OUT_FILTER))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT_LIQ", self.tr("TELEMAC liquid boundaries template (.liq)"), LIQ_FILTER,
                                                                optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the result when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        path = self.parameterAsFile(parameters, "INPUT", context)
        try:
            mesh, fcrs = M.read_mesh_file(path)
        except Exception as exc:
            raise QgsProcessingException(str(exc))
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid():
            crs = fcrs if fcrs is not None else (context.project().crs() if context.project() else QgsCoordinateReferenceSystem())
        opts = self._export_options(parameters, context, crs)
        if opts.fmt == "hecras":
            raise QgsProcessingException(self.tr("HEC-RAS builds its own cells: export the drawing with Generate mesh."))
        opts.quads = True if self.parameterAsBoolean(parameters, "QUADS", context) else False
        opts.liq_path = self.parameterAsFileOutput(parameters, "OUTPUT_LIQ", context) or None
        zsrc = self.parameterAsSource(parameters, "ROUGHNESS_ZONES", context)
        fr = parameters.get("FRICTION")
        if zsrc is not None or fr not in (None, ""):
            zones = M.zones_from_layer(zsrc, crs, context.project()) if zsrc is not None else []
            for line in mesh.assign_zones(zones, default_friction=None if fr in (None, "") else float(fr)):
                feedback.pushInfo(line)
        lines, values = None, {}
        bsrc = self.parameterAsSource(parameters, "BOUNDARIES", context)
        if bsrc is not None:
            lines, values, problems = M.boundary_lines_from_layer(bsrc, crs, project=context.project())
            for p in problems:
                feedback.pushWarning(p)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        files, report, cas = M.export_mesh(out, mesh, crs, opts, lines, values, project=context.project())
        for line in report:
            feedback.pushInfo(line.strip())
        for line in cas:
            feedback.pushInfo(self.tr("Steering file: ") + line)
        feedback.pushInfo(self.tr("Mesh: %s") % mesh.summary())
        if self.parameterAsBoolean(parameters, "LOAD", context) and files:
            local = opts.frame != "keep" and not (opts.fmt == "telemac" and opts.store_offset and not opts.angle)
            self._queue_load(files[0], M.local_crs() if local else crs)
        return {"OUTPUT": files[0] if files else out}


# --------------------------------------------------------------------------
class CutMeshAlgorithm(_MeshBase):
    def name(self):
        return "cutmesh"

    def displayName(self):
        return self.tr("Extract part of a mesh")

    def shortHelpString(self):
        return self.tr("Keep the elements whose centre lies inside the polygons (or outside them) and write the "
                       "smaller mesh, with its node values, in any output format. Useful to cut a study reach out of "
                       "a regional mesh, or to remove an area before meshing it again.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Mesh"), fileFilter=MESH_FILTER))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = from the file or the project)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("POLYGONS", self.tr("Polygons"), [Qgis.ProcessingSourceType.VectorPolygon]))
        self.addParameter(QgsProcessingParameterBoolean("OUTSIDE", self.tr("Keep what is outside the polygons"), False))
        self._add_export_params(default_format=0)
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh file"), OUT_FILTER))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the result when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.geometry import points_in_rings

        path = self.parameterAsFile(parameters, "INPUT", context)
        mesh, fcrs = M.read_mesh_file(path)
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid():
            crs = fcrs if fcrs is not None else context.project().crs()
        rings = [r for _f, g in M._features(self.parameterAsSource(parameters, "POLYGONS", context), crs, context.project())
                 for part in M.polygon_parts(g) for r in part]
        cx, cy = mesh.element_centroids()
        inside = points_in_rings(cx, cy, rings)
        keep = ~inside if self.parameterAsBoolean(parameters, "OUTSIDE", context) else inside
        nt = len(mesh.triangles)
        tri, quads = mesh.triangles[keep[:nt]], mesh.quads[keep[nt:]]
        used = np.zeros(mesh.npoin, dtype=bool)
        used[tri.ravel()] = True
        used[quads.ravel()] = True
        remap = -np.ones(mesh.npoin, dtype=np.int64)
        remap[used] = np.arange(int(used.sum()))
        sub = Mesh(mesh.x[used], mesh.y[used], remap[tri], remap[quads], None if mesh.z is None else mesh.z[used],
                   {k: (u, v[used]) for k, (u, v) in mesh.node_vars.items()}, crs_wkt=mesh.crs_wkt, title=mesh.title)
        if sub.nelem == 0:
            raise QgsProcessingException(self.tr("No element selected."))
        opts = self._export_options(parameters, context, crs)
        opts.quads = True
        files, report, _cas = M.export_mesh(self.parameterAsFileOutput(parameters, "OUTPUT", context), sub, crs, opts,
                                            project=context.project())
        for line in report:
            feedback.pushInfo(line)
        feedback.pushInfo(self.tr("Kept %s of %d elements") % (sub.summary(), mesh.nelem))
        if self.parameterAsBoolean(parameters, "LOAD", context) and files:
            self._queue_load(files[0], crs if opts.frame == "keep" else M.local_crs())
        return {"OUTPUT": files[0]}


# --------------------------------------------------------------------------
class ChannelMeshAlgorithm(_MeshBase):
    def name(self):
        return "channelmesh"

    def displayName(self):
        return self.tr("Structured channel mesh")

    def shortHelpString(self):
        return self.tr(
            "A structured mesh of a channel between a left and a right bank line, optionally split by an axis "
            "(thalweg) line. Lines are grouped by the 'channel' field and told apart by the 'role' field (left, "
            "right, axis). Nodes across and spacing along come from the nodes_across and spacing fields, or from "
            "the defaults below. Use it alone, or let Generate mesh stitch it into the floodplain mesh.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource("CHANNELS", self.tr("Channel lines"), [Qgis.ProcessingSourceType.VectorLine]))
        self.addParameter(QgsProcessingParameterNumber("CROSS_NODES", self.tr("Nodes across (default)"), _INT, 7, minValue=2))
        self.addParameter(QgsProcessingParameterNumber("INTERVAL", self.tr("Spacing along (default, m)"), _DOUBLE, 10.0, minValue=1e-6))
        self._add_export_params(default_format=0)
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh file"), OUT_FILTER))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the mesh when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsSource(parameters, "CHANNELS", context)
        crs = src.sourceCrs()
        n = self.parameterAsInt(parameters, "CROSS_NODES", context)
        iv = self.parameterAsDouble(parameters, "INTERVAL", context)
        warnings, info = [], []
        subs = M.channel_submeshes(src, crs, iv * 2.0, warnings, info, cross_default=n, interval_default=iv)
        for w in warnings:
            feedback.pushWarning(w)
        if not subs:
            raise QgsProcessingException(self.tr("No channel with a left and a right bank line."))
        for line in info:
            feedback.pushInfo(line)
        xs, ys, iks, qs, off = [], [], [], [], 0
        for sm in subs:
            xs.append(sm.x)
            ys.append(sm.y)
            iks.append(sm.ikle + off)
            qs.append(sm.quads + off)
            off += len(sm.x)
        x, y = np.concatenate(xs), np.concatenate(ys)
        mesh = Mesh(x, y, np.vstack(iks), crs_wkt=crs.toWkt())
        opts = self._export_options(parameters, context, crs)
        if opts.quads:
            mesh = Mesh(x, y, np.zeros((0, 3), dtype=np.int64), np.vstack(qs), crs_wkt=crs.toWkt())
        files, report, _cas = M.export_mesh(self.parameterAsFileOutput(parameters, "OUTPUT", context), mesh, crs, opts,
                                            project=context.project())
        for line in report:
            feedback.pushInfo(line)
        if self.parameterAsBoolean(parameters, "LOAD", context) and files:
            self._queue_load(files[0], crs if opts.frame == "keep" else M.local_crs())
        return {"OUTPUT": files[0]}


# --------------------------------------------------------------------------
class RegularGridAlgorithm(_MeshBase):
    PATTERNS = ["alternate", "right", "left"]

    def name(self):
        return "regulargrid"

    def displayName(self):
        return self.tr("Regular grid mesh")

    def shortHelpString(self):
        return self.tr(
            "A structured mesh of nx by ny rectangular cells, starting at the origin (lower-left corner, (0, 0) by "
            "default) and rotated counter-clockwise by the angle. Cells are split into two triangles with "
            "alternating diagonals (no directional bias), or kept as quadrilaterals for formats that accept them. "
            "For grids inside a larger mesh, draw rectangles in grid_zones and use Generate mesh.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterPoint("ORIGIN", self.tr("Origin (lower-left corner)"), defaultValue="0,0"))
        self.addParameter(QgsProcessingParameterNumber("DX", self.tr("Cell size along x"), _DOUBLE, 1.0, minValue=1e-9))
        self.addParameter(QgsProcessingParameterNumber("DY", self.tr("Cell size along y"), _DOUBLE, 1.0, minValue=1e-9))
        self.addParameter(QgsProcessingParameterNumber("NX", self.tr("Cells along x"), _INT, 100, minValue=1))
        self.addParameter(QgsProcessingParameterNumber("NY", self.tr("Cells along y"), _INT, 20, minValue=1))
        self.addParameter(QgsProcessingParameterNumber("ANGLE", self.tr("Rotation (degrees, counter-clockwise)"), _DOUBLE, 0.0))
        self.addParameter(QgsProcessingParameterEnum("PATTERN", self.tr("Diagonals"), self.PATTERNS, defaultValue=0))
        self.addParameter(QgsProcessingParameterCrs("CRS", self.tr("CRS"), defaultValue="ProjectCrs"))
        self.addParameter(QgsProcessingParameterString("FORMULA", self.tr("Bottom: formula of x and y (optional)"), "", optional=True))
        self._add_export_params(default_format=0)
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Mesh file"), OUT_FILTER))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the mesh when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.interp import FormulaSource
        from ..mesher.meshdata import pairs_to_quads

        crs = self.parameterAsCrs(parameters, "CRS", context)
        p = self.parameterAsPoint(parameters, "ORIGIN", context, crs)
        x, y, ik = regular_grid_mesh(p.x(), p.y(), self.parameterAsDouble(parameters, "DX", context),
                                     self.parameterAsDouble(parameters, "DY", context),
                                     self.parameterAsInt(parameters, "NX", context), self.parameterAsInt(parameters, "NY", context),
                                     self.parameterAsDouble(parameters, "ANGLE", context),
                                     self.PATTERNS[self.parameterAsEnum(parameters, "PATTERN", context)])
        opts = self._export_options(parameters, context, crs)
        mesh = Mesh(x, y, ik, crs_wkt=crs.toWkt())
        if opts.quads:
            mesh = Mesh(x, y, np.zeros((0, 3), dtype=np.int64), pairs_to_quads(x, y, ik), crs_wkt=crs.toWkt())
        formula = self.parameterAsString(parameters, "FORMULA", context)
        if formula:
            try:
                mesh.z = FormulaSource(formula).sample(x, y)
            except (ValueError, SyntaxError) as exc:
                raise QgsProcessingException(str(exc))
        files, report, _cas = M.export_mesh(self.parameterAsFileOutput(parameters, "OUTPUT", context), mesh, crs, opts,
                                            project=context.project())
        for line in report:
            feedback.pushInfo(line)
        feedback.pushInfo(mesh.summary())
        if self.parameterAsBoolean(parameters, "LOAD", context) and files:
            self._queue_load(files[0], crs if opts.frame == "keep" else M.local_crs())
        return {"OUTPUT": files[0]}


# --------------------------------------------------------------------------
class InterpolateAlgorithm(_MeshBase):
    OPERATIONS = ["Replace with the new values", "Add the new values (raise or lower by)",
                  "Raise to the new values (keep the higher: fills, dikes, embankments)",
                  "Lower to the new values (keep the lower: dredging, cuts, channels)"]

    def name(self):
        return "interpolatebottom"

    def displayName(self):
        return self.tr("Edit node values of a mesh (bottom, friction, …)")

    def shortHelpString(self):
        return self.tr(
            "Change a node variable of an existing SELAFIN geometry (BOTTOM, BOTTOM FRICTION or any other) "
            "without touching the mesh: same nodes, elements, boundary numbering, origin and date, so the "
            "boundary conditions file (.cli) still fits. Other variables are kept.\n\n"
            "New values come, in priority order, from rasters (bilinear), survey points (linear on their "
            "triangulation), a formula, or a constant. The formula may use x and y, or the channel "
            "coordinates s and n along an axis line (e.g. carve a trapezoidal channel: 'z0 - 0.001*s + "
            "clip((abs(n) - 5)/2, 0, 3)').\n\n"
            "Where: the whole mesh, or only the nodes inside the polygons, with an optional transition width "
            "along the polygon edge (the change fades from full inside to nothing at the edge).\n"
            "How: replace, add (raise or lower by the value), raise to (keep the higher value: fills, dikes), "
            "lower to (keep the lower: dredging, cuts).\n\n"
            "Typical edits: a new survey over a reach (rasters, polygon), a dredged channel (formula in s and n, "
            "lower to), a dike crest (constant, raise to, thin polygon), a friction zone (variable BOTTOM "
            "FRICTION, constant, replace).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Mesh (SELAFIN)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = .prj next to the file, else project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterString("VARIABLE", self.tr("Variable"), "BOTTOM"))
        self.addParameter(QgsProcessingParameterMultipleLayers("RASTERS", self.tr("New values: rasters (priority order)"),
                                                               Qgis.ProcessingSourceType.Raster, optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("POINTS", self.tr("New values: survey points"), [Qgis.ProcessingSourceType.VectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterField("POINT_FIELD", self.tr("Survey value field (empty = point Z)"), parentLayerParameterName="POINTS",
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterString("FORMULA", self.tr("New values: formula of x, y or s, n (lowest priority)"), "", optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("AXIS", self.tr("Channel axis line for s and n in the formula"),
                                                              [Qgis.ProcessingSourceType.VectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber("CONSTANT", self.tr("New values: constant (when no other source)"), _DOUBLE, optional=True))
        self.addParameter(QgsProcessingParameterEnum("OPERATION", self.tr("Operation"), [self.tr(o) for o in self.OPERATIONS], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSource("AREA", self.tr("Only inside these polygons (empty = whole mesh)"),
                                                              [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterNumber("BLEND", self.tr("Transition width inside the polygon edge (m, 0 = sharp)"), _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(_advanced(QgsProcessingParameterString("UNIT", self.tr("Unit of a new variable"), "M", optional=True)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("Z_FACTOR", self.tr("Unit factor of rasters and points"), _DOUBLE, 1.0)))
        self.addParameter(_advanced(QgsProcessingParameterBoolean("FILL", self.tr("Fill nodes no source covers from the nearest node"), True)))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Output geometry file"), GEO_FILTER))
        self.addParameter(QgsProcessingParameterBoolean("LOAD", self.tr("Open the mesh when done"), True))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.interp import formula_names
        from ..telemac_io.convert import writer_like
        from ..telemac_io.integrals import element_mask_in_rings  # noqa: F401  (same even-odd rule)

        path = self.parameterAsFile(parameters, "INPUT", context)
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        if os.path.abspath(out) == os.path.abspath(path):
            raise QgsProcessingException(self.tr("Write to a new file, not over the input."))
        crs = self._mesh_crs(parameters, context, path)
        with self._open(parameters, context) as slf:
            h = slf.header
            if h.is_3d or h.ndp != 3:
                raise QgsProcessingException(self.tr("A 2D triangular mesh is required."))
            if slf.nframes > 1:
                feedback.pushWarning(self.tr("The file has %d frames: the first one is edited and written (a geometry "
                                             "file has one).") % slf.nframes)
            x, y = h.x.copy(), h.y.copy()
            names = [n.strip() for n in h.var_names]
            units = [u.strip() for u in h.var_units]
            data = [np.asarray(r, dtype=np.float64) for r in slf.read_frame(0)] if slf.nframes else [np.zeros(h.npoin) for _ in names]
            t0 = float(slf.times[0]) if slf.nframes else 0.0
            formula = self.parameterAsString(parameters, "FORMULA", context)
            axis = None
            if formula and formula_names(formula) & {"s", "n"}:
                asrc = self.parameterAsSource(parameters, "AXIS", context)
                if asrc is None:
                    raise QgsProcessingException(self.tr("The formula uses s or n: give the channel axis line."))
                lines = [c for _f, g in M._features(asrc, crs, context.project()) for c in M.line_parts(g)]
                if not lines:
                    raise QgsProcessingException(self.tr("The axis layer has no line."))
                axis = np.asarray(max(lines, key=len), dtype=np.float64)[:, :2]
            try:
                sources = M.bottom_sources_from(self.parameterAsLayerList(parameters, "RASTERS", context) or [],
                                                self.parameterAsSource(parameters, "POINTS", context),
                                                self.parameterAsString(parameters, "POINT_FIELD", context), formula,
                                                self.parameterAsDouble(parameters, "Z_FACTOR", context), crs, context.project(), axis)
            except (ValueError, SyntaxError) as exc:
                raise QgsProcessingException(str(exc))
            const = parameters.get("CONSTANT")
            const = None if const in (None, "") else self.parameterAsDouble(parameters, "CONSTANT", context)
            # where
            weight = np.ones(h.npoin)
            area = self.parameterAsSource(parameters, "AREA", context)
            if area is not None:
                rings = [np.asarray(r, dtype=np.float64)[:, :2] for _f, g in M._features(area, crs, context.project())
                         for part in M.polygon_parts(g) for r in part]
                from ..mesher.geometry import points_in_rings

                inside = points_in_rings(x, y, rings)
                blend = self.parameterAsDouble(parameters, "BLEND", context)
                weight = inside.astype(np.float64)
                if blend > 0 and inside.any():
                    d = _distance_to_rings(x[inside], y[inside], rings)
                    weight[inside] = np.clip(d / blend, 0.0, 1.0)
                feedback.pushInfo(self.tr("%d nodes inside the polygons") % int(inside.sum()))
            todo = weight > 0
            if not todo.any():
                raise QgsProcessingException(self.tr("No node to change (check the polygons and the CRS)."))
            # new values, only where needed
            if sources:
                vals = np.full(h.npoin, np.nan)
                sub, _src, rep = M.interpolate_values(x[todo], y[todo], sources, self.parameterAsBoolean(parameters, "FILL", context))
                vals[todo] = sub
                for line in rep:
                    feedback.pushInfo(line)
                if const is not None:
                    vals[todo & ~np.isfinite(vals)] = const
            elif const is not None:
                vals = np.full(h.npoin, float(const))
            else:
                raise QgsProcessingException(self.tr("Give a raster, survey points, a formula or a constant."))
            var = self.parameterAsString(parameters, "VARIABLE", context).strip().upper() or "BOTTOM"
            key = next((i for i, n in enumerate(names) if n.upper() == var or canonical(n) == canonical(var)), None)
            if key is None:
                names.append(var)
                units.append((self.parameterAsString(parameters, "UNIT", context) or "").strip().upper())
                data.append(np.zeros(h.npoin))
                key = len(names) - 1
                feedback.pushInfo(self.tr("New variable %s (0 where it is not set)") % var)
            old = data[key]
            op = self.parameterAsEnum(parameters, "OPERATION", context)
            new = {0: vals, 1: old + vals, 2: np.fmax(old, vals), 3: np.fmin(old, vals)}[op]
            new = np.where(np.isfinite(new), new, old)
            res = old + weight * (new - old)
            changed = np.abs(res - old) > 1e-12
            data[key] = res
            if changed.any():
                dz = (res - old)[changed]
                feedback.pushInfo(self.tr("%s changed at %d nodes: by %.4g to %.4g (mean %.4g)") % (
                    names[key], int(changed.sum()), float(dz.min()), float(dz.max()), float(dz.mean())))
            else:
                feedback.pushWarning(self.tr("No node changed."))
            with writer_like(h, out, names, units) as w:
                w.write_frame(t0, data)
        M.write_prj(out, crs)
        if os.path.exists(os.path.splitext(path)[0] + ".cli") and not os.path.exists(os.path.splitext(out)[0] + ".cli"):
            import shutil

            shutil.copyfile(os.path.splitext(path)[0] + ".cli", os.path.splitext(out)[0] + ".cli")
            feedback.pushInfo(self.tr("Boundary conditions copied unchanged next to the new file (same nodes)."))
        if self.parameterAsBoolean(parameters, "LOAD", context):
            self._queue_load(out, crs)
        return {"OUTPUT": out}


def _distance_to_rings(px, py, rings, chunk=4000):
    """Distance from points to the nearest ring edge."""
    segs = []
    for r in rings:
        r = np.asarray(r, dtype=np.float64)[:, :2]
        if len(r) > 1 and not np.all(np.abs(r[0] - r[-1]) <= 1e-12 * max(1.0, float(np.abs(r).max()))):
            r = np.vstack([r, r[:1]])
        segs.append(np.stack([r[:-1], r[1:]], axis=1))
    S = np.concatenate(segs) if segs else np.zeros((0, 2, 2))
    a, d = S[:, 0], S[:, 1] - S[:, 0]
    L2 = np.maximum(np.einsum("ij,ij->i", d, d), 1e-300)
    out = np.empty(px.size)
    for k in range(0, px.size, chunk):
        X = px[k:k + chunk, None] - a[None, :, 0]
        Y = py[k:k + chunk, None] - a[None, :, 1]
        t = np.clip((X * d[None, :, 0] + Y * d[None, :, 1]) / L2[None, :], 0.0, 1.0)
        out[k:k + chunk] = np.sqrt(np.min((X - t * d[None, :, 0]) ** 2 + (Y - t * d[None, :, 1]) ** 2, axis=1))
    return out


# --------------------------------------------------------------------------
class BoundaryConditionsAlgorithm(_MeshBase):
    def name(self):
        return "boundaryconditions"

    def displayName(self):
        return self.tr("Boundary conditions from lines (.cli)")

    def shortHelpString(self):
        return self.tr(
            "Write the TELEMAC boundary conditions file (.cli). Every boundary node is a wall (2 2 2) except the "
            "nodes near a boundary line, which take its type: wall, Q (4 5 5), H (5 4 4), QH (5 5 5), UV (4 6 6), "
            "UVH (5 6 6), free (4 4 4), waves (1 1 1) or a three-digit code. Nodes are listed in TELEMAC order "
            "(outer boundary counter-clockwise, then islands), so the liquid boundary numbers in the log are the "
            "ones the steering file and the .liq file use. Optionally writes a .liq template with one column per "
            "prescribed value.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Mesh (SELAFIN)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = .prj next to the file, else project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFeatureSource("LINES", self.tr("Boundary lines"), [Qgis.ProcessingSourceType.VectorLine]))
        self.addParameter(QgsProcessingParameterField("TYPE_FIELD", self.tr("Type field"), "type", "LINES"))
        self.addParameter(QgsProcessingParameterField("NAME_FIELD", self.tr("Name field"), "name", "LINES", optional=True))
        self.addParameter(QgsProcessingParameterField("VALUE_FIELD", self.tr("Value field (Q or H for the .liq template)"), "value", "LINES",
                                                      type=Qgis.ProcessingFieldParameterDataType.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterNumber("TOLERANCE", self.tr("Snapping distance (0 = half the boundary edge length)"), _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("Boundary conditions file"), CLI_FILTER))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT_LIQ", self.tr("Liquid boundaries file template"), LIQ_FILTER,
                                                                optional=True, createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        path = self.parameterAsFile(parameters, "INPUT", context)
        crs = self._mesh_crs(parameters, context, path)
        with self._open(parameters, context) as slf:
            h = slf.header
            if h.is_3d or h.ndp != 3:
                raise QgsProcessingException(self.tr("A 2D triangular mesh is required."))
            x, y, ik = h.x.copy(), h.y.copy(), h.ikle.copy()
        src = self.parameterAsSource(parameters, "LINES", context)
        tf = self.parameterAsString(parameters, "TYPE_FIELD", context) or "type"
        nf = self.parameterAsString(parameters, "NAME_FIELD", context) or "name"
        vf = self.parameterAsString(parameters, "VALUE_FIELD", context) or None
        lines, values, problems = M.boundary_lines_from_layer(src, crs, tf, nf, vf, project=context.project())
        for p in problems:
            feedback.pushWarning(p)
        tol = self.parameterAsDouble(parameters, "TOLERANCE", context) or None
        cli = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        liq = self.parameterAsFileOutput(parameters, "OUTPUT_LIQ", context)
        table, report, cas = M.make_boundary_files(x, y, ik, lines, values, cli, liq or None, tol)
        for line in report:
            (feedback.pushWarning if "Warning" in line or "touches no" in line else feedback.pushInfo)(line.strip())
        for line in cas:
            feedback.pushInfo(self.tr("Steering file: ") + line)
        res = {"OUTPUT": cli}
        if liq:
            res["OUTPUT_LIQ"] = liq
        return res


# --------------------------------------------------------------------------
class MeshQualityAlgorithm(_MeshBase):
    def name(self):
        return "meshquality"

    def displayName(self):
        return self.tr("Mesh quality")

    def shortHelpString(self):
        return self.tr(
            "Element quality of a mesh (SELAFIN, 2DM, UGRID, ...) as a polygon layer: area, smallest and largest "
            "angle, longest edge, and a flag for elements below the minimum angle, above the maximum angle or with "
            "all their nodes on the boundary (over-constrained elements, which block the flow in a corner). The log "
            "gives the overall statistics.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("Mesh"), fileFilter=MESH_FILTER))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = from the file or the project)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("MIN_ANGLE", self.tr("Flag angles below (degrees)"), _DOUBLE, 20.0, minValue=0.0, maxValue=60.0))
        self.addParameter(QgsProcessingParameterNumber("MAX_ANGLE", self.tr("Flag angles above (degrees)"), _DOUBLE, 130.0, minValue=60.0, maxValue=180.0))
        self.addParameter(QgsProcessingParameterBoolean("ONLY_FLAGGED", self.tr("Only flagged elements"), False))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Element quality"), Qgis.ProcessingSourceType.VectorPolygon))

    def processAlgorithm(self, parameters, context, feedback):
        path = self.parameterAsFile(parameters, "INPUT", context)
        try:
            mesh, fcrs = M.read_mesh_file(path)
        except Exception as exc:
            raise QgsProcessingException(str(exc))
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid():
            crs = fcrs if fcrs is not None else (context.project().crs() if context.project() else QgsCoordinateReferenceSystem())
        x, y = mesh.x, mesh.y
        st = mesh_statistics(x, y, mesh.all_triangles())
        for k, v in st.items():
            feedback.pushInfo("%s: %s" % (k.replace("_", " "), ("%.4g" % v) if isinstance(v, float) else v))
        elems = [list(t) for t in mesh.triangles.tolist()] + [list(q) for q in mesh.quads.tolist()]
        lo = self.parameterAsDouble(parameters, "MIN_ANGLE", context)
        hi = self.parameterAsDouble(parameters, "MAX_ANGLE", context)
        allt = mesh.all_triangles()
        edges = np.sort(np.concatenate([allt[:, [0, 1]], allt[:, [1, 2]], allt[:, [2, 0]]]), axis=1)
        ue, cnt = np.unique(edges, axis=0, return_counts=True)
        onb = np.zeros(len(x), dtype=bool)
        onb[ue[cnt == 1].ravel()] = True
        fields = QgsFields()
        for n, k in (("element", "int"), ("nodes", "int"), ("area", "double"), ("min_angle", "double"), ("max_angle", "double"),
                     ("max_edge", "double"), ("over_constrained", "int"), ("flagged", "int")):
            fields.append(_field(n, k))
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, Qgis.WkbType.Polygon, crs)
        only = self.parameterAsBoolean(parameters, "ONLY_FLAGGED", context)
        nflag = nmin = nmax = nover = 0
        total = max(len(elems), 1)
        for t, e in enumerate(elems):
            if feedback.isCanceled():
                break
            px, py = x[e], y[e]
            k = len(e)
            ang = []
            for i in range(k):
                ax, ay = px[i - 1] - px[i], py[i - 1] - py[i]
                bx, by = px[(i + 1) % k] - px[i], py[(i + 1) % k] - py[i]
                c = (ax * bx + ay * by) / (math.hypot(ax, ay) * math.hypot(bx, by) + 1e-300)
                ang.append(math.degrees(math.acos(max(-1.0, min(1.0, c)))))
            area = 0.5 * abs(sum(px[i] * py[(i + 1) % k] - px[(i + 1) % k] * py[i] for i in range(k)))
            emax = max(math.hypot(px[(i + 1) % k] - px[i], py[(i + 1) % k] - py[i]) for i in range(k))
            over = bool(onb[e].all()) and k == 3
            fl = min(ang) < lo or max(ang) > hi or over
            nmin += min(ang) < lo
            nmax += max(ang) > hi
            nover += over
            nflag += fl
            if only and not fl:
                continue
            pts = [QgsPointXY(float(a), float(b)) for a, b in zip(px, py)]
            f = QgsFeature(fields)
            f.setGeometry(QgsGeometry.fromPolygonXY([pts + pts[:1]]))
            f.setAttributes([t + 1, k, float(area), float(min(ang)), float(max(ang)), float(emax), int(over), int(fl)])
            sink.addFeature(f, QgsFeatureSink.Flag.FastInsert)
            if t % 2000 == 0:
                feedback.setProgress(100.0 * t / total)
        feedback.pushInfo(self.tr("Flagged elements: %d of %d (angle < %g: %d, angle > %g: %d, over-constrained: %d)") % (
            nflag, len(elems), lo, nmin, hi, nmax, nover))
        return {"OUTPUT": dest}


# --------------------------------------------------------------------------
class CheckDrawingAlgorithm(_MeshBase):
    def name(self):
        return "checkdrawing"

    def displayName(self):
        return self.tr("Check the drawing of a mesh project")

    def shortHelpString(self):
        return self.tr(
            "Finds what QGIS does not flag but the mesher pays for, measured against the element size wanted at each "
            "place: near-duplicate vertices (closer than 1% of the size), spikes where a line turns back on itself, "
            "corners of an outline closer than a quarter of the size, breakline vertices closer than a tenth of it "
            "(every breakline vertex is a node), parts of the outline or islands within half an element of each other "
            "(slivers), line ends and grid-zone corners just off the outline, invalid polygons, islands crossing the "
            "outline and overlapping domain polygons. One point per problem. With Fix, the safe fixes are made and "
            "saved: near-duplicates and small spikes removed, invalid polygons repaired, line ends snapped onto the "
            "outline.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project"), extension="gpkg"))
        self.addParameter(QgsProcessingParameterNumber("DEFAULT_SIZE", self.tr("Default element size (0 = from the project)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterBoolean("FIX", self.tr("Make the safe fixes and save them"), False))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Drawing problems"), Qgis.ProcessingSourceType.VectorPoint))
        self.addOutput(QgsProcessingOutputNumber("PROBLEMS", self.tr("Problems found")))
        self.addOutput(QgsProcessingOutputNumber("FIXABLE", self.tr("Problems with a safe fix")))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge import drawing as D

        path = self.parameterAsFile(parameters, "PROJECT", context)
        layers = M.open_mesh_project(path)
        size = self.parameterAsDouble(parameters, "DEFAULT_SIZE", context)
        if not size:
            try:
                size = float(M.read_settings(path).get("default_size") or 0)
            except ValueError:
                size = 0.0
        size = size or M.auto_default_size(layers, context.project())
        if not size:
            raise QgsProcessingException(self.tr("The project has no domain to check."))
        try:
            build = M.build_mesh_input(layers, None, size, project=context.project(), max_nodes=None)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        chk = D.check_drawing(layers, build.inp.size, build.crs, context.project())
        if self.parameterAsBoolean(parameters, "FIX", context) and chk.fixable:
            for line in D.fix_drawing(layers, chk, build.inp.size, context.project(), commit=True):
                feedback.pushInfo(line)
            chk = D.check_drawing(layers, build.inp.size, build.crs, context.project())
        for line in chk.summary():
            feedback.pushInfo(line)
        for i in chk.issues:
            (feedback.pushWarning if not i.problem.fixable else feedback.pushInfo)(i.text)
        mem = D.issues_layer(chk)
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, mem.fields(), Qgis.WkbType.Point, build.crs)
        for f in mem.getFeatures():
            sink.addFeature(f)
        return {"OUTPUT": dest, "PROBLEMS": len(chk.issues), "FIXABLE": len(chk.fixable)}


MESH_ALGORITHMS = [CreateMeshProjectAlgorithm, ChannelTemplateAlgorithm, RiverReachAlgorithm, GenerateMeshAlgorithm, ExportMeshAlgorithm,
                   CutMeshAlgorithm, ChannelMeshAlgorithm, RegularGridAlgorithm, InterpolateAlgorithm,
                   BoundaryConditionsAlgorithm, MeshQualityAlgorithm, CheckDrawingAlgorithm]
