"""Terrain-driven meshing tools: element sizes from a DEM, breaklines from a DEM, and MMG optimisation or
adaptation of an existing mesh (to the terrain or to a TELEMAC result)."""

from __future__ import annotations

import os

import numpy as np
from qgis.core import (
    Qgis,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessing,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsSettings,
    QgsWkbTypes,
)

from .algorithms import SLF_FILTER, _Base, _field
from .mesh_algorithms import _DOUBLE, _INT, _advanced

MMG_KEY = "telemac_toolkit/mmg2d"


class _TerrainBase(_Base):
    def group(self):
        return self.tr("Mesh")

    def groupId(self):
        return "mesh"

    def _extent(self, parameters, context, dem_layer, crs, project_path):
        """(x0, y0, x1, y1) in ``crs``: the given extent, else the project's domain, else the DEM."""
        from ..qgis_bridge import meshing as M

        if parameters.get("EXTENT") not in (None, ""):
            r = self.parameterAsExtent(parameters, "EXTENT", context, crs)
            if not r.isEmpty():
                return (r.xMinimum(), r.yMinimum(), r.xMaximum(), r.yMaximum())
        if project_path:
            rings = M.domain_rings(M.open_mesh_project(project_path), crs, context.project())
            if rings:
                a = np.vstack(rings)
                return (float(a[:, 0].min()), float(a[:, 1].min()), float(a[:, 0].max()), float(a[:, 1].max()))
        r = dem_layer.extent()
        if dem_layer.crs() != crs and crs.isValid():
            from qgis.core import QgsCoordinateTransform

            r = QgsCoordinateTransform(dem_layer.crs(), crs, context.transformContext()).transformBoundingBox(r)
        return (r.xMinimum(), r.yMinimum(), r.xMaximum(), r.yMaximum())


# --------------------------------------------------------------------------
class TerrainSizesAlgorithm(_TerrainBase):
    def name(self):
        return "terrainsizes"

    def displayName(self):
        return self.tr("Element sizes from terrain (raster)")

    def shortHelpString(self):
        return self.tr(
            "The element sizes terrain refinement asks for, as a raster to look at before meshing: small where linear "
            "elements would miss the DEM by more than the tolerance (banks, levees, embankments), the largest size on "
            "flat ground, and never growing faster than the growth rate. TELEMAC interpolates the bed linearly in each "
            "element; over an element of size h that misses a curvature kappa by about h² kappa / 8, hence "
            "h = sqrt(8 tolerance / kappa). Generate mesh applies the same sizes with its terrain tolerance.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM")))
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (its domain sets the extent)"),
                                                     extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterExtent("EXTENT", self.tr("Extent (empty = the project's domain or the DEM)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("TOLERANCE", self.tr("Bed tolerance (m)"), _DOUBLE, 0.05, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("MIN_SIZE", self.tr("Smallest size (m)"), _DOUBLE, 2.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("MAX_SIZE", self.tr("Largest size (m)"), _DOUBLE, 30.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("GROWTH", self.tr("Growth rate"), _DOUBLE, 1.2, minValue=1.0, maxValue=3.0))
        self.addParameter(_advanced(QgsProcessingParameterNumber("CELL", self.tr("Cell size (m, 0 = half the smallest size)"),
                                                                 _DOUBLE, 0.0, minValue=0.0)))
        self.addParameter(QgsProcessingParameterRasterDestination("OUTPUT", self.tr("Sizes")))

    def processAlgorithm(self, parameters, context, feedback):
        from osgeo import gdal

        from ..mesher.demsize import terrain_sizes
        from ..qgis_bridge import meshing as M

        dem_layer = self.parameterAsRasterLayer(parameters, "DEM", context)
        crs = dem_layer.crs()
        project = self.parameterAsFile(parameters, "PROJECT", context)
        if project:
            crs = M.open_mesh_project(project)["domain"].crs()
        bbox = self._extent(parameters, context, dem_layer, crs, project)
        sampler = M.RasterSampler(dem_layer, 1, crs, project=context.project())
        try:
            gd = terrain_sizes(sampler.sample, bbox, self.parameterAsDouble(parameters, "TOLERANCE", context),
                               self.parameterAsDouble(parameters, "MIN_SIZE", context),
                               self.parameterAsDouble(parameters, "MAX_SIZE", context),
                               self.parameterAsDouble(parameters, "GROWTH", context),
                               self.parameterAsDouble(parameters, "CELL", context) or None)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        out = self.parameterAsOutputLayer(parameters, "OUTPUT", context)
        v = np.asarray(gd.values, dtype=np.float32)[::-1]  # north up
        ny, nx = v.shape
        ds = gdal.GetDriverByName("GTiff").Create(out, nx, ny, 1, gdal.GDT_Float32, options=["COMPRESS=DEFLATE"])
        ds.SetGeoTransform((gd.x0 - 0.5 * gd.dx, gd.dx, 0.0, gd.y0 + (ny - 0.5) * gd.dy, 0.0, -gd.dy))
        ds.SetProjection(crs.toWkt())
        ds.GetRasterBand(1).WriteArray(v)
        ds = None
        fine = float(np.mean(gd.values < 0.75 * float(np.max(gd.values))))
        feedback.pushInfo(self.tr("Sizes from %.3g m to %.3g m; %.0f%% of the area below three quarters of the largest size")
                          % (float(np.min(gd.values)), float(np.max(gd.values)), 100 * fine))
        return {"OUTPUT": out}


# --------------------------------------------------------------------------
class DemBreaklinesAlgorithm(_TerrainBase):
    KINDS = ["Crests and toes", "Crests (levee and bank tops, embankments)", "Toes (bank feet, channel edges)"]
    MODES = ["Guide lines (resampled)", "Breaklines (every vertex kept)"]

    def name(self):
        return "dembreaklines"

    def displayName(self):
        return self.tr("Breaklines from a DEM")

    def shortHelpString(self):
        return self.tr(
            "Lines along the sharp bends of the terrain, for the mesh to follow: crests (levee and bank tops, road "
            "embankments: the ground bends down on both sides) and toes (bank feet, channel edges: it bends up). A line "
            "is drawn where an element of the given size would miss the bend by more than the tolerance "
            "(|curvature| > 8 tolerance / size²). The cells beyond that are thinned to one-cell-wide lines, traced, "
            "simplified to the cell size, and lines shorter than the minimum length dropped. With a mesh project the "
            "lines are clipped to its domain and can be added to its constraint lines.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM")))
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (extent, and where to add the lines)"),
                                                     extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterExtent("EXTENT", self.tr("Extent (empty = the project's domain or the DEM)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("TOLERANCE", self.tr("Bed tolerance (m)"), _DOUBLE, 0.05, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("SIZE", self.tr("Element size the lines are for (m, 0 = the project's default size)"),
                                                       _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterEnum("KINDS", self.tr("Lines"), [self.tr(k) for k in self.KINDS], defaultValue=0))
        self.addParameter(_advanced(QgsProcessingParameterNumber("CELL", self.tr("Cell size (m, 0 = the DEM's pixel)"), _DOUBLE, 0.0,
                                                                 minValue=0.0)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("MIN_LENGTH", self.tr("Shortest line kept (m, 0 = three element sizes)"),
                                                                 _DOUBLE, 0.0, minValue=0.0)))
        self.addParameter(QgsProcessingParameterBoolean("ADD", self.tr("Add the lines to the project's constraint lines"), False))
        self.addParameter(QgsProcessingParameterEnum("MODE", self.tr("Added as"), [self.tr(m) for m in self.MODES], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSink("OUTPUT", self.tr("Breaklines"), QgsProcessing.SourceType.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.breaklines import breaklines_from_dem
        from ..qgis_bridge import meshing as M

        dem_layer = self.parameterAsRasterLayer(parameters, "DEM", context)
        project = self.parameterAsFile(parameters, "PROJECT", context)
        crs = dem_layer.crs()
        size = self.parameterAsDouble(parameters, "SIZE", context)
        layers = None
        if project:
            layers = M.open_mesh_project(project)
            crs = layers["domain"].crs()
            if not size:
                try:
                    size = float(M.read_settings(project).get("default_size") or 0)
                except ValueError:
                    size = 0.0
                size = size or M.auto_default_size(layers, context.project())
        if not size:
            raise QgsProcessingException(self.tr("Give the element size the lines are for."))
        bbox = self._extent(parameters, context, dem_layer, crs, project)
        cell = self.parameterAsDouble(parameters, "CELL", context) or max(dem_layer.rasterUnitsPerPixelX(), 1e-6)
        sampler = M.RasterSampler(dem_layer, 1, crs, project=context.project())
        res = breaklines_from_dem(sampler.sample, bbox, self.parameterAsDouble(parameters, "TOLERANCE", context), size, cell,
                                  self.parameterAsDouble(parameters, "MIN_LENGTH", context) or None)
        kind = self.parameterAsEnum(parameters, "KINDS", context)
        wanted = {0: ("crest", "toe"), 1: ("crest",), 2: ("toe",)}[kind]
        dom = None
        if layers is not None:
            geoms = [g for _f, g in M._features(layers["domain"], crs, context.project())]
            if geoms:
                # half an element inside the outline: a line along it would leave thin elements between them
                dom = QgsGeometry.unaryUnion(geoms) if len(geoms) > 1 else geoms[0]
                dom = dom.buffer(-0.5 * size, 8)
            # not inside structured channels or grids, nor within one of their edges of them
            subs = M.channel_submeshes(layers.get("channels"), crs, size, [], []) + M.grid_submeshes(layers.get("grid_zones"), crs, size, [])
            from ..telemac_io.mesh import boundary_loops

            for sm in subs:
                for lp in boundary_loops(sm.ikle, sm.x, sm.y):
                    ring = [QgsPointXY(float(sm.x[i]), float(sm.y[i])) for i in lp]
                    edge = float(np.median(np.hypot(np.diff(sm.x[lp]), np.diff(sm.y[lp])))) if len(lp) > 2 else size
                    hole = QgsGeometry.fromPolygonXY([ring + ring[:1]]).buffer(1.0 * edge, 8)
                    dom = dom.difference(hole) if dom is not None else None
            if subs:
                feedback.pushInfo(self.tr("Lines left out inside the %d structured part(s) and within one of their edges") % len(subs))
        fields = QgsFields()
        fields.append(_field("kind", "string"))
        fields.append(_field("length", "double"))
        sink, dest = self.parameterAsSink(parameters, "OUTPUT", context, fields, QgsWkbTypes.Type.LineString, crs)
        out_geoms = []
        for k in wanted:
            for c in res.get(k, []):
                g = QgsGeometry.fromPolylineXY([QgsPointXY(float(a), float(b)) for a, b in c])
                if dom is not None:
                    g = g.intersection(dom)
                    if g.isEmpty() or QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Line:
                        continue
                parts = g.asMultiPolyline() if g.isMultipart() else [g.asPolyline()]
                for p in parts:
                    if len(p) < 2:
                        continue
                    pg = QgsGeometry.fromPolylineXY(p)
                    ft = QgsFeature(fields)
                    ft.setGeometry(pg)
                    ft.setAttributes([k, pg.length()])
                    sink.addFeature(ft, QgsFeatureSink.Flag.FastInsert)
                    out_geoms.append((k, pg))
        feedback.pushInfo(self.tr("%d lines (%s), %.0f m in all") % (len(out_geoms), ", ".join(
            "%d %ss" % (sum(1 for kk, _g in out_geoms if kk == k), k) for k in wanted), sum(g.length() for _k, g in out_geoms)))
        if self.parameterAsBoolean(parameters, "ADD", context):
            if layers is None:
                raise QgsProcessingException(self.tr("Give the mesh project to add the lines to."))
            lyr = layers["constraint_lines"]
            mode = ("guide", "breakline")[self.parameterAsEnum(parameters, "MODE", context)]
            # each line's element size: no more than its distance to the nearest other line, so the elements between
            # two close lines (a bank top and its foot) are not squeezed
            sizes = []
            for i, (_k, g) in enumerate(out_geoms):
                near = [g.distance(h) for j, (_kk, h) in enumerate(out_geoms) if j != i and g.boundingBox().buffered(size).intersects(h.boundingBox())]
                sizes.append(max(min([size] + [d for d in near if d > 0]), 0.1 * size))
            lyr.startEditing()
            for (k, g), s in zip(out_geoms, sizes):
                f = QgsFeature(lyr.fields())
                f.setGeometry(g)
                for name, val in (("name", "DEM %s" % k), ("mode", mode), ("size", float(s) if s < size else None)):
                    if lyr.fields().indexOf(name) >= 0:
                        f.setAttribute(name, val)
                lyr.addFeature(f)
            if not lyr.commitChanges():
                raise QgsProcessingException("; ".join(lyr.commitErrors()))
            feedback.pushInfo(self.tr("%d lines added to the constraint lines as %s") % (len(out_geoms), mode))
            M.refresh_project_layers(project, context.project())
        return {"OUTPUT": dest}


# --------------------------------------------------------------------------
class MmgAdaptAlgorithm(_TerrainBase):
    MODES = ["Optimise (better shapes, the same sizes)", "Adapt to the terrain (DEM and bed tolerance)",
             "Adapt to a TELEMAC result (a variable and its tolerance)"]
    WHEN = ["Last frame", "Largest value over all frames"]

    def name(self):
        return "mmgadapt"

    def displayName(self):
        return self.tr("Optimise or adapt a mesh (MMG)")

    def shortHelpString(self):
        return self.tr(
            "Runs MMG's mmg2d on a TELEMAC mesh (MMG is not part of the plugin: install it from mmgtools.org or the "
            "conda-forge package mmgsuite, and give the path of mmg2d_O3 once; it is remembered).\n\n"
            "• Optimise: better element shapes with the same sizes.\n"
            "• Adapt to the terrain: sizes from the DEM, small where linear elements would miss it by more than the "
            "bed tolerance (see Element sizes from terrain).\n"
            "• Adapt to a result: sizes from the curvature of a variable of a TELEMAC result on this mesh (free "
            "surface, velocity, depth), with a tolerance in its unit, for a second run with elements where the flow "
            "needs them (solution-adaptive refinement).\n\n"
            "The boundary nodes do not move, so the boundary conditions carry over; elements inside the 'keep' polygons "
            "(structured channels, bridges) are left as they are, and sizes next to them are held to their edges. Every "
            "variable is carried over (unchanged nodes exactly, new ones interpolated); with a DEM the bottom is "
            "sampled from it. The .cli next to the mesh is rebuilt.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("INPUT", self.tr("TELEMAC mesh (SELAFIN, .cli next to it)"), fileFilter=SLF_FILTER))
        self.addParameter(QgsProcessingParameterEnum("MODE", self.tr("Mode"), [self.tr(m) for m in self.MODES], defaultValue=0))
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM (terrain mode; also resamples the bottom)"), optional=True))
        self.addParameter(QgsProcessingParameterNumber("TOLERANCE", self.tr("Tolerance (m for the terrain; the variable's unit for a result)"),
                                                       _DOUBLE, 0.05, minValue=1e-9))
        self.addParameter(QgsProcessingParameterFile("RESULT", self.tr("TELEMAC result on this mesh (result mode)"), fileFilter=SLF_FILTER,
                                                     optional=True))
        self.addParameter(QgsProcessingParameterString("VARIABLE", self.tr("Variable of the result"), "FREE SURFACE", optional=True))
        self.addParameter(QgsProcessingParameterEnum("WHEN", self.tr("Frames of the result"), [self.tr(w) for w in self.WHEN], defaultValue=1))
        self.addParameter(QgsProcessingParameterNumber("HMIN", self.tr("Smallest size (m)"), _DOUBLE, 1.0, minValue=1e-6))
        self.addParameter(QgsProcessingParameterNumber("HMAX", self.tr("Largest size (m, 0 = the mesh's longest edge)"), _DOUBLE, 0.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber("HGRAD", self.tr("Growth rate between neighbours"), _DOUBLE, 1.2, minValue=1.0, maxValue=3.0))
        self.addParameter(QgsProcessingParameterFeatureSource("KEEP", self.tr("Areas to keep as they are (polygons)"),
                                                              [Qgis.ProcessingSourceType.VectorPolygon], optional=True))
        self.addParameter(QgsProcessingParameterCrs("MESH_CRS", self.tr("CRS of the mesh (empty = project CRS)"), optional=True))
        self.addParameter(QgsProcessingParameterFile("MMG", self.tr("mmg2d program (empty = remembered, MMG2D variable or PATH)"),
                                                     optional=True))
        self.addParameter(QgsProcessingParameterFileDestination("OUTPUT", self.tr("New mesh"), SLF_FILTER))

    def processAlgorithm(self, parameters, context, feedback):
        from ..mesher.demsize import apply_floor, field_sizes, terrain_sizes
        from ..mesher.edit import match_nodes, region_elements
        from ..mesher.mmg import find_mmg2d, run_mmg2d
        from ..qgis_bridge import meshing as M
        from ..telemac_io import SelafinFile
        from ..telemac_io.mesh import boundary_loops
        from .edit_algorithms import _EditBase, _copy_prj, _read_slf

        path = self.parameterAsFile(parameters, "INPUT", context)
        src = _read_slf(path)
        x, y, ik = src["x"], src["y"], src["ikle"]
        exe = self.parameterAsFile(parameters, "MMG", context) or str(QgsSettings().value(MMG_KEY, "") or "") or None
        prog = find_mmg2d(exe)
        if prog is None:
            raise QgsProcessingException(self.tr("mmg2d was not found. Install MMG (mmgtools.org, or conda-forge 'mmgsuite') and give "
                                                 "the path of mmg2d_O3."))
        QgsSettings().setValue(MMG_KEY, prog)
        crs = self.parameterAsCrs(parameters, "MESH_CRS", context)
        if not crs.isValid():
            from ..qgis_bridge.crs import sidecar_crs

            crs = sidecar_crs(path) or (context.project().crs() if context.project() else crs)
        keep = np.zeros(len(ik), dtype=bool)
        ksrc = self.parameterAsSource(parameters, "KEEP", context)
        if ksrc is not None:
            for _f, g in M._features(ksrc, crs, context.project()):
                for part in M.polygon_parts(g):
                    keep |= region_elements(x, y, ik, part)
            feedback.pushInfo(self.tr("%d elements kept as they are") % int(keep.sum()))
        mode = self.parameterAsEnum(parameters, "MODE", context)
        hmin = self.parameterAsDouble(parameters, "HMIN", context)
        P = np.stack([x[ik], y[ik]], axis=2)
        longest = float(np.max(np.hypot(*(np.roll(P, -1, axis=1) - P).reshape(-1, 2).T)))
        hmax = self.parameterAsDouble(parameters, "HMAX", context) or longest
        hgrad = self.parameterAsDouble(parameters, "HGRAD", context)
        tol = self.parameterAsDouble(parameters, "TOLERANCE", context)
        dem_layer = self.parameterAsRasterLayer(parameters, "DEM", context)
        dem = M.RasterSampler(dem_layer, 1, crs, project=context.project()).sample if dem_layer is not None else None
        sizes = None
        if mode == 1:
            if dem is None:
                raise QgsProcessingException(self.tr("Terrain mode needs a DEM."))
            pad = 2 * hmax
            gd = terrain_sizes(dem, (x.min() - pad, y.min() - pad, x.max() + pad, y.max() + pad), tol, hmin, hmax, hgrad)
            if keep.any():
                lp = boundary_loops(ik[keep], x, y)
                a = np.vstack([np.stack([x[lo], y[lo]], axis=1) for lo in lp])
                b = np.vstack([np.roll(np.stack([x[lo], y[lo]], axis=1), -1, axis=0) for lo in lp])
                gd = apply_floor(gd, a, b, np.hypot(*(b - a).T), hgrad)
            sizes = gd(x, y)
        elif mode == 2:
            rpath = self.parameterAsFile(parameters, "RESULT", context)
            if not rpath:
                raise QgsProcessingException(self.tr("Result mode needs a TELEMAC result."))
            var = (self.parameterAsString(parameters, "VARIABLE", context) or "FREE SURFACE").strip().upper()
            with SelafinFile(rpath) as s:
                names = [n.strip().upper() for n in s.header.var_names]
                if var not in names:
                    raise QgsProcessingException(self.tr("The result has no %s (it has %s)") % (var, ", ".join(names)))
                if s.header.npoin != len(x):
                    raise QgsProcessingException(self.tr("The result is on another mesh (%d nodes, the mesh has %d)") % (s.header.npoin, len(x)))
                if self.parameterAsEnum(parameters, "WHEN", context) == 0:
                    f = s.read(var, s.nframes - 1)
                else:
                    f = np.max(np.array([s.read(var, k) for k in range(s.nframes)]), axis=0)
            sizes = field_sizes(x, y, ik, f, tol, hmin, hmax)
            feedback.pushInfo(self.tr("Sizes from the curvature of %s: %.3g m to %.3g m") % (var, float(sizes.min()), float(sizes.max())))
        try:
            nx, ny, nik, log = run_mmg2d(x, y, ik, sizes=sizes, optimise=(mode == 0), keep=keep if keep.any() else None,
                                         hmin=hmin, hmax=hmax, hgrad=hgrad, exe=prog)
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        for line in log[-4:]:
            feedback.pushInfo("mmg2d: " + line.strip())
        feedback.pushInfo(self.tr("Mesh: %d nodes, %d elements (before: %d, %d)") % (len(nx), len(nik), len(x), len(ik)))
        out = self.parameterAsFileOutput(parameters, "OUTPUT", context)
        match = match_nodes(x, y, nx, ny, max(1e-9 * max(np.ptp(x), np.ptp(y)), 1e-9))
        res = _EditBase._carry(self, feedback, src, nx, ny, nik, os.path.splitext(path)[0] + ".cli", out, match)
        if dem is not None:  # the bottom from the DEM at every node
            names = [n.upper() for n in src["names"]]
            for bname in ("BOTTOM", "FOND"):
                if bname in names:
                    _resample_variable(out, src["names"][names.index(bname)], np.asarray(dem(nx, ny), dtype=np.float64))
                    feedback.pushInfo(self.tr("%s sampled from the DEM at every node") % bname)
                    break
        _copy_prj(path, out)
        return res


def _resample_variable(path, name, values):
    """Overwrite one variable of every frame of a SELAFIN file (nodes where values are NaN keep theirs)."""
    from ..telemac_io import SelafinFile
    from ..telemac_io.selafin import SelafinWriter

    with SelafinFile(path) as s:
        h = s.header
        frames = [s.read_frame(k) for k in range(s.nframes)]
        times = list(s.times)
        names = [n.strip() for n in h.var_names]
        units = [u.strip() for u in h.var_units]
        x, y, ik, title, date = h.x.copy(), h.y.copy(), h.ikle.copy(), h.title, h.date
    i = names.index(name.strip())
    tmp = path + ".tmp"
    with SelafinWriter(tmp, names, units, x, y, ik, title=title.strip() if isinstance(title, str) else "", date=date) as w:
        for t, fr in zip(times, frames):
            fr = np.array(fr, dtype=np.float64)
            fr[i] = np.where(np.isfinite(values), values, fr[i])
            w.write_frame(float(t), fr)
    os.replace(tmp, path)


class SuggestSizesAlgorithm(_TerrainBase):
    def name(self):
        return "suggestsizes"

    def displayName(self):
        return self.tr("Suggest element sizes from a DEM")

    def shortHelpString(self):
        return self.tr(
            "Sizes for a mesh from what the DEM can show, the terrain, the rivers and a node budget:\n"
            "• smallest size: two pixels (finer elements only interpolate the same pixels);\n"
            "• terrain tolerance: about the vertical accuracy of a DEM with such pixels (1 m LiDAR 0.05 m, 5 m DEMs "
            "0.2 m, 15 m DEMs 0.5 m);\n"
            "• default size (floodplain): what the median terrain inside the domain asks for at that tolerance, at "
            "least four pixels, at most a twentieth of the domain, within half the node budget;\n"
            "• rivers (the project's structured channels): a spacing across of a eighth of the width (7 nodes or more), "
            "1.5 times that along;\n"
            "• the number of nodes with the terrain refinement (the tolerance is relaxed if it is over the budget) and "
            "the time step at Courant number 1.\n\n"
            "Without a project, or before the domain is drawn, the DEM's own extent is used. With 'Apply', the default "
            "size, terrain tolerance and smallest size go into the project's settings (step 1 of the panel).")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer("DEM", self.tr("DEM")))
        self.addParameter(QgsProcessingParameterFile("PROJECT", self.tr("Mesh project (domain and channels; optional)"),
                                                     extension="gpkg", optional=True))
        self.addParameter(QgsProcessingParameterNumber("NODE_BUDGET", self.tr("Node budget"), _INT, 300000, minValue=1000))
        self.addParameter(_advanced(QgsProcessingParameterNumber("DEPTH", self.tr("Typical depth for the time step (m)"), _DOUBLE, 2.0,
                                                                 minValue=0.01)))
        self.addParameter(_advanced(QgsProcessingParameterNumber("VELOCITY", self.tr("Typical velocity for the time step (m/s)"),
                                                                 _DOUBLE, 1.5, minValue=0.0)))
        self.addParameter(QgsProcessingParameterBoolean("APPLY", self.tr("Apply to the project's settings"), False))
        from qgis.core import QgsProcessingOutputNumber

        for key, label in (("DEFAULT_SIZE", "Default size"), ("SMALLEST", "Smallest size"), ("TOLERANCE", "Terrain tolerance"),
                           ("CHANNEL_ACROSS", "Spacing across rivers"), ("CHANNEL_ALONG", "Spacing along rivers"),
                           ("NODES", "Nodes"), ("TIME_STEP", "Time step")):
            self.addOutput(QgsProcessingOutputNumber(key, self.tr(label)))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_bridge import meshing as M

        dem = self.parameterAsRasterLayer(parameters, "DEM", context)
        if dem is None:
            raise QgsProcessingException(self.tr("Give a DEM."))
        path = self.parameterAsFile(parameters, "PROJECT", context)
        layers = M.open_mesh_project(path) if path else None
        crs = layers["domain"].crs() if layers and layers.get("domain") is not None else dem.crs()
        try:
            sug = M.suggest_sizes_for(layers, dem, crs, context.project(), self.parameterAsInt(parameters, "NODE_BUDGET", context),
                                      self.parameterAsDouble(parameters, "DEPTH", context),
                                      self.parameterAsDouble(parameters, "VELOCITY", context))
        except ValueError as exc:
            raise QgsProcessingException(str(exc))
        for line in sug.report():
            feedback.pushInfo(line)
        if self.parameterAsBoolean(parameters, "APPLY", context):
            if not path:
                raise QgsProcessingException(self.tr("Give the mesh project to apply the sizes to."))
            M.write_settings(path, {"default_size": "%g" % sug.default, "terrain_tolerance": "%g" % sug.tolerance,
                                    "terrain_min_size": "%g" % sug.smallest})
            feedback.pushInfo(self.tr("Default size, terrain tolerance and smallest size written to the project"))
        return {"DEFAULT_SIZE": sug.default, "SMALLEST": sug.smallest, "TOLERANCE": sug.tolerance,
                "CHANNEL_ACROSS": sug.channel_across or 0.0, "CHANNEL_ALONG": sug.channel_along or 0.0, "NODES": sug.nodes,
                "TIME_STEP": sug.time_step}


TERRAIN_ALGORITHMS = [TerrainSizesAlgorithm, DemBreaklinesAlgorithm, MmgAdaptAlgorithm, SuggestSizesAlgorithm]
