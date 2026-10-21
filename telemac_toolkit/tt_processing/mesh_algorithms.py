"""Processing algorithms of the mesher (group "Mesh")."""

from __future__ import annotations

import math
import os

import numpy as np
from qgis.core import (
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
                               store_offset=self.parameterAsBoolean(parameters, "STORE_OFFSET", context))


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
class ChannelTemplateAlgorithm(_MeshBase):
    SECTIONS = ["Rectangular (flat bed)", "Trapezoidal"]
    PLANFORMS = ["Straight", "Bend (straight inlet, circular bend, straight outlet)", "Meanders (sine-generated)"]
    TURNS = ["Left (counter-clockwise)", "Right (clockwise)"]

    def name(self):
        return "channeltemplate"

    def displayName(self):
        return self.tr("Channel or flume template (straight, bend, meanders)")

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
            "coordinates.")

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
            layers = M.create_mesh_project(path, M.local_crs(), project=context.project())
            if context.project() is not None and not context.project().crs().isValid():
                context.project().setCrs(M.local_crs())
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
        trap = self.parameterAsEnum(parameters, "SECTION", context) == 1
        z0 = self.parameterAsDouble(parameters, "BED_Z0", context)
        S = self.parameterAsDouble(parameters, "SLOPE", context)
        bw = self.parameterAsDouble(parameters, "BOTTOM_WIDTH", context) if trap else 0.0
        ss = self.parameterAsDouble(parameters, "SIDE_SLOPE", context) if trap else 0.0
        hb = self.parameterAsDouble(parameters, "BANK_HEIGHT", context) if trap else 0.0
        settings = {"default_size": "%g" % cell}
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
            axis = None
            if formula_names(formula) & {"s", "n"}:
                axis = M.channel_axis(layers, settings.get("bottom_axis", ""), crs, context.project())
            sources = M.bottom_sources_from(self.parameterAsLayerList(parameters, "RASTERS", context) or [],
                                            self.parameterAsSource(parameters, "POINTS", context),
                                            self.parameterAsString(parameters, "POINT_FIELD", context), formula, zf, crs,
                                            context.project(), axis)
        except (ValueError, SyntaxError) as exc:
            raise QgsProcessingException(self.tr("Bottom: %s") % exc)
        fr = parameters.get("FRICTION")
        opts = self._export_options(parameters, context, crs)
        opts.liq_path = self.parameterAsFileOutput(parameters, "OUTPUT_LIQ", context) or None
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        try:
            res = M.generate_and_export(layers, out, opts, crs, size, self.parameterAsDouble(parameters, "GROWTH", context),
                                        self.parameterAsInt(parameters, "SMOOTHING", context), sources, True,
                                        None if fr in (None, "") else float(fr), None, self._progress(feedback), context.project())
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
        if len(r) > 1 and not np.allclose(r[0], r[-1]):
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


MESH_ALGORITHMS = [CreateMeshProjectAlgorithm, ChannelTemplateAlgorithm, GenerateMeshAlgorithm, ExportMeshAlgorithm,
                   CutMeshAlgorithm, ChannelMeshAlgorithm, RegularGridAlgorithm, InterpolateAlgorithm,
                   BoundaryConditionsAlgorithm, MeshQualityAlgorithm]
