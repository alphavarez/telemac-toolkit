"""QGIS side of the mesher: drawing layers in, mesh files out.

A mesh project is a GeoPackage with one layer per kind of object:

- ``domain``: polygons to mesh (holes are islands), optional ``size``
- ``islands``: polygons left out of the mesh
- ``refinement_areas`` / ``refinement_points``: target element ``size``
- ``constraint_lines``: lines the mesh follows; ``mode`` = breakline (every
  vertex kept) or guide (resampled to the local size); optional ``size``
- ``fixed_points``: nodes that must exist
- ``channels``: left and right banks and an optional axis of a channel, meshed
  as a structured strip (``nodes_across``, ``spacing``)
- ``grid_zones``: rectangles meshed as regular grids (``dx``, ``dy``)
- ``roughness_zones``: ``material`` id and ``friction`` value
- ``boundaries``: boundary condition lines (``type``, ``name``, ``value``)
- ``settings``: a table of key/value settings (sizes, bottom formula, export)

Layers of 0.2.0 projects (outline, density, hard_lines, ...) are read too.
The functions read any vector layer or Processing source, build the NumPy
mesher input with GEOS, triangulate with GEOS and hand a neutral ``Mesh`` to
the writers (TELEMAC first, other programs from the same mesh).
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from qgis.core import (
    Qgis,
    QgsCategorizedSymbolRenderer,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsField,
    QgsFillSymbol,
    QgsGeometry,
    QgsLineSymbol,
    QgsMarkerSymbol,
    QgsMeshLayer,
    QgsPalLayerSettings,
    QgsPointXY,
    QgsProject,
    QgsRectangle,
    QgsRendererCategory,
    QgsVectorFileWriter,
    QgsVectorLayer,
    QgsVectorLayerSimpleLabeling,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QMetaType
from qgis.PyQt.QtGui import QColor

from ..mesher import MeshInput, MeshResult, Submesh, channel_mesh, generate
from ..mesher.backends import triangulate_geos
from ..mesher.boundary import BoundaryLine, build_boundary_table, liq_template, parse_type
from ..mesher.geometry import as_ring, snap_line_to_nodes, snap_rings_to_nodes
from ..mesher.interp import FormulaSource, PointSource, bilinear, fill_nearest, interpolate_nodes
from ..mesher.meshdata import Frame, Mesh, pairs_to_quads
from ..mesher.regular import regular_grid_mesh
from ..mesher.sizefield import LineDensity, PointDensity, PolygonDensity, SizeField
from ..mesher.writers import FORMATS, write_frame_note, write_mesh
from ..telemac_io import write_cli, write_selafin
from ..telemac_io.cli import BOUNDARY_TYPES
from ..telemac_io.mesh import boundary_loops
from .crs import LOCAL_CRS_WKT, is_local_crs, local_crs, write_prj  # noqa: F401 (re-exported)

__all__ = [
    "PROJECT_LAYERS", "LOCAL_CRS_WKT", "local_crs", "create_mesh_project", "open_mesh_project", "style_project_layer",
    "read_settings", "write_settings", "add_flume", "build_mesh_input", "run_mesher", "RasterSampler", "point_source",
    "interpolate_values", "boundary_lines_from_layer", "make_boundary_files", "write_geometry", "write_prj",
    "style_geometry_mesh", "export_mesh", "read_mesh_file", "zones_from_layer",
]

#: key, geometry, fields (name, kind), title in the layer tree, older names read as this layer
PROJECT_LAYERS = [
    ("domain", "Polygon", [("name", "string"), ("size", "double")], "Domain (holes are islands)", ["outline"]),
    ("islands", "Polygon", [("name", "string")], "Islands (not meshed)", []),
    ("refinement_areas", "Polygon", [("name", "string"), ("size", "double")], "Refinement areas", ["density"]),
    ("refinement_points", "Point", [("name", "string"), ("size", "double")], "Refinement points", ["density_points"]),
    ("constraint_lines", "LineString", [("name", "string"), ("mode", "string"), ("size", "double")],
     "Constraint lines (breakline or guide)", []),
    ("fixed_points", "Point", [("name", "string"), ("size", "double")], "Fixed points", ["hard_points"]),
    ("channels", "LineString", [("channel", "string"), ("role", "string"), ("nodes_across", "int"), ("spacing", "double")],
     "Structured channels (left, right, axis)", []),
    ("grid_zones", "Polygon", [("name", "string"), ("dx", "double"), ("dy", "double")], "Structured grid zones", ["grids"]),
    ("roughness_zones", "Polygon", [("name", "string"), ("material", "int"), ("friction", "double")], "Roughness zones", []),
    ("boundaries", "LineString", [("name", "string"), ("type", "string"), ("value", "double")], "Boundary conditions", []),
]
LAYER_KEYS = [p[0] for p in PROJECT_LAYERS]
_OLD_LAYERS = {"outline": "domain", "density": "refinement_areas", "density_points": "refinement_points",
               "hard_points": "fixed_points", "grids": "grid_zones"}
_FIELD_ALIASES = {"size": ["size", "edge_len"], "nodes_across": ["nodes_across", "cross_nodes"],
                  "spacing": ["spacing", "interval"]}
BOUNDARY_TYPE_NAMES = ["wall", "Q", "H", "QH", "UV", "UVH", "free", "waves"]
CONSTRAINT_MODES = ("breakline", "guide")

SETTINGS_DEFAULTS = {
    "default_size": "", "growth": "1.2", "smoothing": "8", "bottom_formula": "", "bottom_axis": "", "z_factor": "1",
    "format": "telemac", "quads": "0", "frame": "keep", "frame_x": "0", "frame_y": "0", "frame_angle": "0",
    "store_offset": "1", "mesh_file": "",
}


def _qfield(name, kind):
    t = {"int": QMetaType.Type.Int, "double": QMetaType.Type.Double, "string": QMetaType.Type.QString}[kind]
    return QgsField(name, t)


# ---------------------------------------------------------------------------
# mesh project (GeoPackage)
# ---------------------------------------------------------------------------
def _write_layer(mem, path, name, first, ctx):
    opts = QgsVectorFileWriter.SaveVectorOptions()
    opts.driverName = "GPKG"
    opts.layerName = name
    opts.actionOnExistingFile = (QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile if first
                                 else QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer)
    res = QgsVectorFileWriter.writeAsVectorFormatV3(mem, path, ctx, opts)
    if res[0] != QgsVectorFileWriter.WriterError.NoError:
        raise RuntimeError("cannot write layer %s: %s" % (name, res[1]))


def create_mesh_project(path: str, crs: QgsCoordinateReferenceSystem, project: Optional[QgsProject] = None,
                        group_name: Optional[str] = None, add_to_project: bool = True,
                        settings: Optional[Dict[str, str]] = None) -> Dict[str, QgsVectorLayer]:
    """Create the mesh project GeoPackage (empty layers + settings table), styled and added to the project."""
    project = project or QgsProject.instance()
    if os.path.exists(path):
        raise FileExistsError(path)
    ctx = project.transformContext()
    for i, (name, geom, fields, _title, _old) in enumerate(PROJECT_LAYERS):
        mem = QgsVectorLayer("%s?crs=%s" % (geom, crs.authid() or "EPSG:4326"), name, "memory")
        mem.setCrs(crs)
        mem.dataProvider().addAttributes([_qfield(n, k) for n, k in fields])
        mem.updateFields()
        _write_layer(mem, path, name, i == 0, ctx)
    st = QgsVectorLayer("None", "settings", "memory")
    st.dataProvider().addAttributes([_qfield("key", "string"), _qfield("value", "string")])
    st.updateFields()
    values = dict(SETTINGS_DEFAULTS)
    values.update(settings or {})
    feats = []
    for k, v in values.items():
        f = QgsFeature(st.fields())
        f.setAttributes([k, str(v)])
        feats.append(f)
    st.dataProvider().addFeatures(feats)
    _write_layer(st, path, "settings", False, ctx)
    layers = open_mesh_project(path)
    if add_to_project:
        add_project_to_map(layers, path, project, group_name)
    return layers


def open_mesh_project(path: str) -> Dict[str, QgsVectorLayer]:
    """The project's layers by key (0.2.0 layer names are mapped to the new keys), styled."""
    from qgis.core import QgsProviderRegistry

    if not os.path.exists(path):
        raise RuntimeError("cannot open %s" % path)
    names = [d.name() for d in QgsProviderRegistry.instance().querySublayers(path)]
    titles = {p[0]: p[3] for p in PROJECT_LAYERS}
    old_titles = {"hard_lines": "Breaklines (0.2.0 hard lines)", "soft_lines": "Guide lines (0.2.0 soft lines)"}
    layers = {}
    for name in names:
        if name == "settings":
            continue
        key = _OLD_LAYERS.get(name, name)
        if key not in titles and name not in old_titles:
            continue
        lyr = QgsVectorLayer("%s|layername=%s" % (path, name), titles.get(key, old_titles.get(name, name)), "ogr")
        if lyr.isValid():
            style_project_layer(lyr, key)
            layers[key if key in titles else name] = lyr
    return layers


def add_project_to_map(layers: Dict[str, QgsVectorLayer], path: str, project: Optional[QgsProject] = None,
                       group_name: Optional[str] = None):
    project = project or QgsProject.instance()
    root = project.layerTreeRoot()
    group = root.insertGroup(0, group_name or "Mesh project – %s" % os.path.splitext(os.path.basename(path))[0])
    group.setCustomProperty("telemac_toolkit/mesh_project", path)
    order = LAYER_KEYS + ["hard_lines", "soft_lines"]
    # boundaries and points on top, domain at the bottom
    for key in sorted(layers, key=lambda k: -order.index(k) if k in order else 0, reverse=True):
        project.addMapLayer(layers[key], False)
        group.addLayer(layers[key])
    return group


def read_settings(path: str) -> Dict[str, str]:
    out = dict(SETTINGS_DEFAULTS)
    lyr = QgsVectorLayer("%s|layername=settings" % path, "settings", "ogr")
    if lyr.isValid():
        for f in lyr.getFeatures():
            k, v = f["key"], f["value"]
            if k:
                out[str(k)] = "" if v is None or (hasattr(v, "isNull") and v.isNull()) else str(v)
    return out


def write_settings(path: str, values: Dict[str, str]):
    """Update (or add) settings in the project's settings table."""
    lyr = QgsVectorLayer("%s|layername=settings" % path, "settings", "ogr")
    if not lyr.isValid():
        st = QgsVectorLayer("None", "settings", "memory")
        st.dataProvider().addAttributes([_qfield("key", "string"), _qfield("value", "string")])
        st.updateFields()
        _write_layer(st, path, "settings", False, QgsProject.instance().transformContext())
        lyr = QgsVectorLayer("%s|layername=settings" % path, "settings", "ogr")
    existing = {f["key"]: f.id() for f in lyr.getFeatures()}
    lyr.startEditing()
    ki, vi = lyr.fields().indexOf("key"), lyr.fields().indexOf("value")
    for k, v in values.items():
        if k in existing:
            lyr.changeAttributeValue(existing[k], vi, str(v))
        else:
            f = QgsFeature(lyr.fields())
            f.setAttribute(ki, k)
            f.setAttribute(vi, str(v))
            lyr.addFeature(f)
    if not lyr.commitChanges():
        raise RuntimeError("; ".join(lyr.commitErrors()))


def _label(layer, expression):
    pal = QgsPalLayerSettings()
    pal.fieldName = expression
    pal.isExpression = True
    layer.setLabeling(QgsVectorLayerSimpleLabeling(pal))
    layer.setLabelsEnabled(True)


def style_project_layer(layer: QgsVectorLayer, kind: str):
    """Symbols for the mesh project layers."""
    size_field = "coalesce(\"size\", \"edge_len\")" if layer.fields().indexOf("edge_len") >= 0 and layer.fields().indexOf("size") >= 0 else (
        "\"edge_len\"" if layer.fields().indexOf("edge_len") >= 0 else "\"size\"")
    if kind == "domain":
        layer.renderer().setSymbol(QgsFillSymbol.createSimple({"color": "0,0,0,0", "outline_color": "#202020",
                                                               "outline_width": "0.6"}))
    elif kind == "islands":
        layer.renderer().setSymbol(QgsFillSymbol.createSimple({"color": "120,120,120,90", "outline_color": "#202020",
                                                               "outline_width": "0.5", "style": "b_diagonal"}))
    elif kind == "refinement_areas":
        layer.renderer().setSymbol(QgsFillSymbol.createSimple({"color": "255,170,0,45", "outline_color": "#e08a00",
                                                               "outline_width": "0.4", "outline_style": "dash"}))
        _label(layer, "'size ' || format_number(%s, 2)" % size_field)
    elif kind == "refinement_points":
        layer.renderer().setSymbol(QgsMarkerSymbol.createSimple({"name": "circle", "color": "#ffaa00", "size": "2.6",
                                                                 "outline_color": "#8a5a00"}))
        _label(layer, "format_number(%s, 2)" % size_field)
    elif kind == "constraint_lines":
        cats = [QgsRendererCategory("breakline", QgsLineSymbol.createSimple({"line_color": "#d62728", "line_width": "0.6"}),
                                    "Breakline (every vertex kept)"),
                QgsRendererCategory("guide", QgsLineSymbol.createSimple({"line_color": "#1f77b4", "line_width": "0.5",
                                                                         "line_style": "dash"}), "Guide line (resampled)")]
        # an empty mode is a breakline for the mesher, so it is drawn as one
        layer.setRenderer(QgsCategorizedSymbolRenderer(
            "CASE WHEN lower(\"mode\") IN ('guide', 'soft', 'resampled') THEN 'guide' ELSE 'breakline' END", cats))
    elif kind in ("hard_lines",):
        layer.renderer().setSymbol(QgsLineSymbol.createSimple({"line_color": "#d62728", "line_width": "0.6"}))
    elif kind in ("soft_lines",):
        layer.renderer().setSymbol(QgsLineSymbol.createSimple({"line_color": "#1f77b4", "line_width": "0.5",
                                                               "line_style": "dash"}))
    elif kind == "fixed_points":
        layer.renderer().setSymbol(QgsMarkerSymbol.createSimple({"name": "cross2", "color": "#d62728", "size": "3.2",
                                                                 "outline_color": "#d62728", "outline_width": "0.5"}))
    elif kind == "channels":
        cats = []
        for value, color, style, label in (("left", "#0b3d91", "solid", "Left bank"), ("right", "#7a1fa2", "solid", "Right bank"),
                                           ("axis", "#17a2c4", "dash", "Axis (thalweg)")):
            sym = QgsLineSymbol.createSimple({"line_color": color, "line_width": "0.7", "line_style": style})
            cats.append(QgsRendererCategory(value, sym, label))
        # every spelling the mesher accepts for the axis is drawn (and listed) once
        layer.setRenderer(QgsCategorizedSymbolRenderer(
            "CASE WHEN lower(\"role\") IN ('axis', 'thalweg', 'centre', 'center', 'centerline', 'centreline', 'a', 't') "
            "THEN 'axis' WHEN lower(\"role\") IN ('l', 'left') THEN 'left' WHEN lower(\"role\") IN ('r', 'right') "
            "THEN 'right' ELSE lower(\"role\") END", cats))
    elif kind == "grid_zones":
        layer.renderer().setSymbol(QgsFillSymbol.createSimple({"color": "148,103,189,40", "outline_color": "#6a3d9a",
                                                               "outline_width": "0.5"}))
        _label(layer, "format_number(\"dx\", 2) || ' x ' || format_number(coalesce(\"dy\", \"dx\"), 2)")
    elif kind == "roughness_zones":
        layer.renderer().setSymbol(QgsFillSymbol.createSimple({"color": "46,158,68,35", "outline_color": "#2e9e44",
                                                               "outline_width": "0.4", "outline_style": "dot"}))
        _label(layer, "coalesce(\"name\", '') || ' #' || coalesce(to_string(\"material\"), '') || "
                      "if(\"friction\" is null, '', ' n=' || format_number(\"friction\", 3))")
    elif kind == "boundaries":
        cats = []
        for name in BOUNDARY_TYPE_NAMES:
            codes = parse_type(name)
            label, color = BOUNDARY_TYPES.get(codes, (name, "#7f7f7f"))
            sym = QgsLineSymbol.createSimple({"line_color": color, "line_width": "1.4"})
            cats.append(QgsRendererCategory(name.lower(), sym, "%s (%d%d%d)" % (label, *codes)))
        layer.setRenderer(QgsCategorizedSymbolRenderer("lower(\"type\")", cats))
        _label(layer, "coalesce(\"name\", '') || ' [' || coalesce(\"type\", '') || ']'")
    layer.triggerRepaint()


# ---------------------------------------------------------------------------
# reading layers
# ---------------------------------------------------------------------------
def _crs_of(src) -> QgsCoordinateReferenceSystem:
    if hasattr(src, "sourceCrs"):
        return src.sourceCrs()
    return src.crs()


def _features(src, crs, project=None):
    """(feature, geometry in ``crs``, straight segments) for every non-empty feature."""
    if src is None:
        return
    project = project or QgsProject.instance()
    scrs = _crs_of(src)
    ct = None
    if scrs.isValid() and crs.isValid() and scrs != crs:
        ct = QgsCoordinateTransform(scrs, crs, project.transformContext())
    for f in src.getFeatures():
        g = f.geometry()
        if g is None or g.isNull() or g.isEmpty():
            continue
        g = QgsGeometry(g)
        if QgsWkbTypes.isCurvedType(g.wkbType()):
            g = QgsGeometry(g.constGet().segmentize())
        if ct is not None:
            g.transform(ct)
        yield f, g


def _value(f, name, default=None):
    names = _FIELD_ALIASES.get(name, [name])
    for n in names:
        try:
            idx = f.fields().indexOf(n)
        except Exception:
            idx = -1
        if idx < 0:
            continue
        v = f.attribute(idx)
        if v is None:
            continue
        if callable(getattr(v, "isNull", None)) and v.isNull():  # a NULL QVariant
            continue
        if isinstance(v, str) and not v.strip():
            continue
        return v
    return default


def _num(f, name, default=None):
    v = _value(f, name, None)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return default
    return v if math.isfinite(v) else default


def _xy(points) -> np.ndarray:
    return np.array([[p.x(), p.y()] for p in points], dtype=np.float64).reshape(-1, 2)


def polygon_parts(g: QgsGeometry) -> List[List[np.ndarray]]:
    """Polygon(s) as lists of rings [outer, hole, ...] (closing vertex removed)."""
    if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Polygon:
        parts = []
        if g.isMultipart() or QgsWkbTypes.flatType(g.wkbType()) == Qgis.WkbType.GeometryCollection:
            for p in g.asGeometryCollection():
                if QgsWkbTypes.geometryType(p.wkbType()) == Qgis.GeometryType.Polygon:
                    parts += polygon_parts(p)
        return parts
    polys = g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]
    out = []
    for poly in polys:
        rings = [as_ring(_xy(r)) for r in poly if len(r) >= 4]
        if rings:
            out.append(rings)
    return out


def line_parts(g: QgsGeometry) -> List[np.ndarray]:
    t = QgsWkbTypes.geometryType(g.wkbType())
    if t == Qgis.GeometryType.Line:
        lines = g.asMultiPolyline() if g.isMultipart() else [g.asPolyline()]
        return [_xy(line) for line in lines if len(line) >= 2]
    if t == Qgis.GeometryType.Polygon:  # polygons used as lines: their rings
        return [np.vstack([r, r[:1]]) for part in polygon_parts(g) for r in part]
    if g.isMultipart() or QgsWkbTypes.flatType(g.wkbType()) == Qgis.WkbType.GeometryCollection:
        out = []
        for p in g.asGeometryCollection():
            if QgsWkbTypes.geometryType(p.wkbType()) == Qgis.GeometryType.Line:
                out += line_parts(p)
        return out
    return []


def point_parts(g: QgsGeometry) -> np.ndarray:
    if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Point:
        return np.zeros((0, 2))
    pts = g.asMultiPoint() if g.isMultipart() else [g.asPoint()]
    return _xy(pts)


def _is_closed(c: np.ndarray) -> bool:
    return len(c) >= 4 and np.allclose(c[0], c[-1])


def _poly_geom(rings: Sequence[np.ndarray]) -> QgsGeometry:
    return QgsGeometry.fromPolygonXY([[QgsPointXY(float(a), float(b)) for a, b in np.vstack([r, r[:1]])] for r in rings])


def nice_number(v: float) -> float:
    """1, 2, 2.5 or 5 times a power of ten, the nearest below ``v``."""
    if not (v > 0 and math.isfinite(v)):
        return 0.0
    p = 10.0 ** math.floor(math.log10(v))
    for m in (5.0, 2.5, 2.0, 1.0):
        if m * p <= v * 1.0000001:
            return m * p
    return p


def auto_default_size(layers: Dict[str, object], project=None) -> float:
    """A default element size when none is set: about 1/40 of the square root of the domain area
    (a 20 m x 1 m flume gives 0.1 m, a 2 km² valley 25 m)."""
    dom = normalize_layers(layers).get("domain")
    if dom is None:
        return 0.0
    area = 0.0
    for _f, g in _features(dom, _crs_of(dom), project):
        area += abs(g.area())
    return nice_number(math.sqrt(area) / 40.0) if area > 0 else 0.0


def normalize_layers(layers: Dict[str, object]) -> Dict[str, object]:
    """Accept 0.2.0 keys (outline, density, hard_lines, ...) next to the current ones."""
    out = {}
    for k, v in (layers or {}).items():
        if v is None:
            continue
        out[_OLD_LAYERS.get(k, k)] = v if _OLD_LAYERS.get(k, k) not in out else out[_OLD_LAYERS.get(k, k)]
    return out


# ---------------------------------------------------------------------------
# building the mesher input
# ---------------------------------------------------------------------------
@dataclass
class MeshBuild:
    inp: MeshInput
    crs: QgsCoordinateReferenceSystem
    warnings: List[str] = field(default_factory=list)
    info: List[str] = field(default_factory=list)
    domain_geom: Optional[QgsGeometry] = None


def channel_submeshes(src, crs, default_edge: float, warnings: List[str], info: List[str],
                      cross_default: Optional[int] = None, interval_default: Optional[float] = None) -> List[Submesh]:
    """One structured channel strip per value of the 'channel' field (roles left, right, axis)."""
    groups: Dict[str, dict] = {}
    for f, g in _features(src, crs):
        parts = line_parts(g)
        if not parts:
            continue
        line = max(parts, key=lambda c: float(np.sum(np.hypot(*np.diff(c, axis=0).T))))
        cid = str(_value(f, "channel", "1"))
        role = str(_value(f, "role", "")).strip().lower()
        role = {"l": "left", "r": "right", "a": "axis", "t": "axis", "thalweg": "axis", "centre": "axis", "center": "axis",
                "centerline": "axis", "centreline": "axis"}.get(role, role)
        grp = groups.setdefault(cid, {})
        if role not in ("left", "right", "axis"):
            warnings.append("Channel %s: feature %d has role %r; use left, right or axis." % (cid, f.id(), role))
            continue
        grp[role] = line
        n = _num(f, "nodes_across")
        if n:
            grp["nodes_across"] = int(n)
        iv = _num(f, "spacing")
        if iv:
            grp["spacing"] = iv
    subs = []
    for cid, grp in groups.items():
        if "left" not in grp or "right" not in grp:
            warnings.append("Channel %s needs a left and a right bank line; skipped." % cid)
            continue
        n = int(grp.get("nodes_across", 0) or cross_default or 0)
        ax = grp.get("axis")
        if n < 2:
            n = 7 if ax is not None else 5
        if ax is not None and n % 2 == 0:
            n += 1
            warnings.append("Channel %s: nodes across raised to %d (odd with an axis line)." % (cid, n))
        iv = float(grp.get("spacing", 0.0) or 0.0) or interval_default or default_edge * 0.5
        x, y, ik = channel_mesh(grp["left"], grp["right"], ax, cross_nodes=n, interval=iv)
        subs.append(Submesh(x, y, ik, pairs_to_quads(x, y, ik)))
        info.append("Channel %s: %d nodes across, spacing %.3g, %d nodes, %d cells" % (cid, n, iv, x.size, len(ik) // 2))
    return subs


def grid_submeshes(src, crs, default_edge: float, info: List[str]) -> List[Submesh]:
    """Regular grids fitted to the minimum oriented rectangle of each grid zone."""
    subs = []
    for f, g in _features(src, crs):
        dx = _num(f, "dx") or default_edge * 0.5
        dy = _num(f, "dy") or dx
        bb = g.orientedMinimumBoundingBox()
        box = bb[0] if isinstance(bb, tuple) else bb
        ring = box.asPolygon()[0][:4]
        p = _xy(ring)
        k = int(np.argmin(p[:, 0] + p[:, 1]))  # lower-left corner first
        p = np.roll(p, -k, axis=0)
        u, v = p[1] - p[0], p[3] - p[0]
        if u[0] * v[1] - u[1] * v[0] < 0:
            u, v = v, u
        w, h = math.hypot(*u), math.hypot(*v)
        if w <= 0 or h <= 0:
            continue
        nx, ny = max(1, int(round(w / dx))), max(1, int(round(h / dy)))
        ang = math.degrees(math.atan2(u[1], u[0]))
        x, y, ik = regular_grid_mesh(p[0, 0], p[0, 1], w / nx, h / ny, nx, ny, angle_deg=ang)
        subs.append(Submesh(x, y, ik, pairs_to_quads(x, y, ik)))
        info.append("Grid %s: %d x %d cells of %.3g x %.3g, rotated %.1f deg" % (
            _value(f, "name", f.id()), nx, ny, w / nx, h / ny, ang))
    return subs


def _constraint_sources(layers: Dict[str, object]):
    """(source, forced mode or None) for every layer that holds constraint lines."""
    out = []
    if layers.get("constraint_lines") is not None:
        out.append((layers["constraint_lines"], None))
    if layers.get("hard_lines") is not None:
        out.append((layers["hard_lines"], "breakline"))
    if layers.get("soft_lines") is not None:
        out.append((layers["soft_lines"], "guide"))
    return out


def build_mesh_input(layers: Dict[str, object], crs: Optional[QgsCoordinateReferenceSystem], default_size: float,
                     growth: float = 1.2, smoothing: int = 8, project=None) -> MeshBuild:
    """Read mesh project layers (keys as in PROJECT_LAYERS, or 0.2.0 names); missing ones are skipped."""
    layers = normalize_layers(layers)
    domain_src = layers.get("domain")
    if domain_src is None:
        raise ValueError("a domain polygon layer is required")
    crs = crs if crs is not None and crs.isValid() else _crs_of(domain_src)
    warnings: List[str] = []
    info: List[str] = []
    geoms, polys_d = [], []
    for f, g in _features(domain_src, crs, project):
        if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Polygon:
            continue
        geoms.append(g)
        s = _num(f, "size")
        if s and s > 0:
            for part in polygon_parts(g):
                polys_d.append(PolygonDensity(part, s))
    if not geoms:
        raise ValueError("the domain layer has no polygon")
    out_geom = QgsGeometry.unaryUnion(geoms) if len(geoms) > 1 else geoms[0]
    if not out_geom.isGeosValid():
        out_geom = out_geom.makeValid()
        warnings.append("The domain was not a valid polygon and was repaired.")
    isl = [g for _f, g in _features(layers.get("islands"), crs, project)
           if QgsWkbTypes.geometryType(g.wkbType()) == Qgis.GeometryType.Polygon]
    if isl:
        out_geom = out_geom.difference(QgsGeometry.unaryUnion(isl) if len(isl) > 1 else isl[0])
        info.append("Islands: %d" % len(isl))

    # structured parts ----------------------------------------------------------
    subs = channel_submeshes(layers.get("channels"), crs, default_size, warnings, info)
    subs += grid_submeshes(layers.get("grid_zones"), crs, default_size, info)
    sub_nodes = np.zeros((0, 2))
    seg_a = seg_b = np.zeros((0, 2))
    domain_geom = out_geom
    snap_tol = on_tol = 0.0
    if subs:
        polys, nodes, sa, sb, emax = [], [], [], [], 0.0
        for sm in subs:
            loops = boundary_loops(sm.ikle, sm.x, sm.y)
            rings = [np.stack([sm.x[lp], sm.y[lp]], axis=1) for lp in loops]
            polys.append(_poly_geom(rings))
            for r in rings:
                nodes.append(r)
                sa.append(r)
                sb.append(np.roll(r, -1, axis=0))
                emax = max(emax, float(np.max(np.hypot(*(np.roll(r, -1, axis=0) - r).T))))
        sub_geom = QgsGeometry.unaryUnion(polys) if len(polys) > 1 else polys[0]
        if sum(p.area() for p in polys) > sub_geom.area() * (1 + 1e-9):
            warnings.append("Structured parts overlap each other; the result will not be a valid mesh.")
        if not out_geom.buffer(emax * 1e-3, 4).contains(sub_geom):
            warnings.append("A structured part extends beyond the domain; it is kept whole.")
        domain_geom = out_geom.difference(sub_geom)
        sub_nodes = np.vstack(nodes)
        seg_a, seg_b = np.vstack(sa), np.vstack(sb)
        snap_tol = 0.51 * emax
        on_tol = max(1e-7 * max(out_geom.boundingBox().width(), out_geom.boundingBox().height()), 1e-6)
    domain: List[List[np.ndarray]] = []
    for part in polygon_parts(domain_geom):
        if subs:
            outer = snap_rings_to_nodes(part[:1], sub_nodes, snap_tol, (seg_a, seg_b), on_tol)
            if not outer:
                continue
            part = outer + snap_rings_to_nodes(part[1:], sub_nodes, snap_tol, (seg_a, seg_b), on_tol)
        domain.append(part)
    if not domain and not subs:
        raise ValueError("the domain is empty")
    if domain:
        nholes = sum(len(p) - 1 for p in domain)
        info.append("Domain: %d polygon(s), %d hole(s) (islands and structured parts)" % (len(domain), nholes))
    else:
        info.append("Structured parts cover the whole domain")

    # refinement ------------------------------------------------------------------
    pts_d, lines_d = [], []
    for f, g in _features(layers.get("refinement_areas"), crs, project):
        s = _num(f, "size")
        if not s or s <= 0:
            warnings.append("Refinement area %d has no size; ignored." % f.id())
            continue
        for part in polygon_parts(g):
            polys_d.append(PolygonDensity(part, s))
    dxy, ds = [], []
    for f, g in _features(layers.get("refinement_points"), crs, project):
        s = _num(f, "size")
        if s and s > 0:
            for p in point_parts(g):
                dxy.append(p)
                ds.append(s)

    # constraints -------------------------------------------------------------------
    hard, soft = [], []
    for src, forced in _constraint_sources(layers):
        for f, g in _features(src, crs, project):
            mode = forced or str(_value(f, "mode", "breakline")).strip().lower()
            mode = {"hard": "breakline", "fixed": "breakline", "break": "breakline", "soft": "guide",
                    "resampled": "guide"}.get(mode, mode)
            if mode not in CONSTRAINT_MODES:
                warnings.append("Constraint line %d has mode %r; read as breakline." % (f.id(), mode))
                mode = "breakline"
            s = _num(f, "size")
            closed_src = any(_is_closed(c) for c in line_parts(g))
            gi = g.intersection(domain_geom)
            if gi.isEmpty():
                continue
            merged = gi.mergeLines() if QgsWkbTypes.geometryType(gi.wkbType()) == Qgis.GeometryType.Line else gi
            for c in line_parts(merged):
                if subs:
                    c = snap_line_to_nodes(c, sub_nodes, snap_tol, (seg_a, seg_b), on_tol)
                if len(c) < 2:
                    continue
                closed = _is_closed(c) and closed_src
                (hard if mode == "breakline" else soft).append((c[:-1] if closed else c, closed))
                if s and s > 0:
                    lines_d.append(LineDensity(c, s, closed))
    if hard:
        info.append("Breaklines: %d" % len(hard))
    if soft:
        info.append("Guide lines: %d" % len(soft))
    hp = []
    for f, g in _features(layers.get("fixed_points"), crs, project):
        s = _num(f, "size")
        for p in point_parts(g):
            if domain_geom.contains(QgsGeometry.fromPointXY(QgsPointXY(*p))):
                hp.append(p)
                if s and s > 0:
                    dxy.append(p)
                    ds.append(s)
            else:
                warnings.append("Fixed point %d is outside the area to mesh; ignored." % f.id())
    if dxy:
        pts_d.append(PointDensity(np.array(dxy), np.array(ds)))
    size = SizeField(default_size, growth, polys_d, pts_d, lines_d)
    inp = MeshInput(domain=domain, size=size, soft_lines=soft, hard_lines=hard,
                    hard_points=np.array(hp) if hp else None, submeshes=subs, smooth_iterations=int(smoothing))
    return MeshBuild(inp, crs, warnings, info, domain_geom)


def run_mesher(build: MeshBuild, progress: Optional[Callable[[float, str], bool]] = None) -> MeshResult:
    return generate(build.inp, triangulate_geos, progress)


def structured_quads(build: MeshBuild, res: MeshResult) -> np.ndarray:
    """Quads of the structured parts, as node ids of the generated mesh."""
    return res.quads


# ---------------------------------------------------------------------------
# flume / schematic channel template
# ---------------------------------------------------------------------------
def add_flume(layers: Dict[str, QgsVectorLayer], length: float, width: float, origin=(0.0, 0.0), angle_deg: float = 0.0,
              cell: Optional[float] = None, structured: bool = True, inflow: Optional[float] = None,
              outflow_level: Optional[float] = None, name: str = "flume") -> List[str]:
    """Draw a straight channel or flume: domain rectangle, inlet (Q) and outlet (H) lines, optional grid zone.

    The rectangle starts at ``origin`` and runs ``length`` along ``angle_deg``
    (counter-clockwise from the x axis); the inlet is the first short side.
    """
    a = math.radians(angle_deg)
    ux, uy = math.cos(a), math.sin(a)
    vx, vy = -uy, ux
    x0, y0 = origin
    c = [(x0, y0), (x0 + length * ux, y0 + length * uy), (x0 + length * ux + width * vx, y0 + length * uy + width * vy),
         (x0 + width * vx, y0 + width * vy)]
    ring = [QgsPointXY(*p) for p in c + c[:1]]
    report = []

    def add(key, geom, **attrs):
        lyr = layers[key]
        f = QgsFeature(lyr.fields())
        f.setGeometry(geom)
        for k, v in attrs.items():
            if lyr.fields().indexOf(k) >= 0:
                f.setAttribute(k, v)
        lyr.startEditing()
        lyr.addFeature(f)
        if not lyr.commitChanges():
            raise RuntimeError("; ".join(lyr.commitErrors()))

    add("domain", QgsGeometry.fromPolygonXY([ring]), name=name, size=cell)
    # boundary lines slightly inside the corners so the corner nodes stay walls
    m = 0.6 * (cell or width / 10.0)  # more than the snapping distance (half an edge)
    inlet = [QgsPointXY(x0 + m * vx, y0 + m * vy), QgsPointXY(x0 + (width - m) * vx, y0 + (width - m) * vy)]
    ex, ey = x0 + length * ux, y0 + length * uy
    outlet = [QgsPointXY(ex + m * vx, ey + m * vy), QgsPointXY(ex + (width - m) * vx, ey + (width - m) * vy)]
    add("boundaries", QgsGeometry.fromPolylineXY(inlet), name="inlet", type="Q", value=inflow)
    add("boundaries", QgsGeometry.fromPolylineXY(outlet), name="outlet", type="H", value=outflow_level)
    report.append("Domain %g m x %g m from (%g, %g), %g deg; inlet (Q) and outlet (H) lines" % (length, width, x0, y0, angle_deg))
    if structured and cell:
        add("grid_zones", QgsGeometry.fromPolygonXY([ring]), name=name, dx=cell, dy=cell)
        report.append("Structured grid zone, cells of %g m" % cell)
    return report


def add_planform(layers: Dict[str, QgsVectorLayer], pf, cell: float, structured: bool = True,
                 inflow: Optional[float] = None, outflow_level: Optional[float] = None, name: str = "flume") -> List[str]:
    """Draw a channel of any planform (bend, meanders): domain, inlet and outlet lines and, when
    ``structured``, its left bank, right bank and axis in the channels layer, vertex by vertex, so the
    strip's cross-sections are normal to the axis."""
    report = []

    def add(key, geom, **attrs):
        lyr = layers[key]
        f = QgsFeature(lyr.fields())
        f.setGeometry(geom)
        for k, v in attrs.items():
            if lyr.fields().indexOf(k) >= 0:
                f.setAttribute(k, v)
        lyr.startEditing()
        lyr.addFeature(f)
        if not lyr.commitChanges():
            raise RuntimeError("; ".join(lyr.commitErrors()))

    def pts(a):
        return [QgsPointXY(float(x), float(y)) for x, y in a]

    ring = pf.outline()
    add("domain", QgsGeometry.fromPolygonXY([pts(ring)]), name=name, size=cell)
    W = pf.width
    m = 0.6 * cell / W  # boundary lines stop short of the corners, which stay walls
    L, R = pf.left, pf.right
    for k, label, typ, val in ((0, "inlet", "Q", inflow), (-1, "outlet", "H", outflow_level)):
        a, b = R[k], L[k]
        line = [a + m * (b - a), b - m * (b - a)]
        add("boundaries", QgsGeometry.fromPolylineXY(pts(line)), name=label, type=typ, value=val)
    report.append("%s channel: axis length %.4g m, width %g m, %d cross-sections; inlet (Q) and outlet (H) lines"
                  % (pf.kind.capitalize(), pf.length, W, len(pf.axis)))
    if structured:
        across = max(3, int(round(W / cell)) + 1)
        if across % 2 == 0:
            across += 1
        for role, line in (("left", L), ("right", R), ("axis", pf.axis)):
            add("channels", QgsGeometry.fromPolylineXY(pts(line)), channel=name, role=role, nodes_across=across, spacing=cell)
        report.append("Structured strip: %d nodes across, cross-sections every %.3g m along the axis "
                      "(%.3g m on the inner bank, %.3g m on the outer one at the tightest bend)" % (
                          across, pf.length / (len(pf.axis) - 1), *_bank_spacing(pf)))
    return report


def _bank_spacing(pf):
    dl = np.hypot(*np.diff(pf.left, axis=0).T)
    dr = np.hypot(*np.diff(pf.right, axis=0).T)
    return float(min(dl.min(), dr.min())), float(max(dl.max(), dr.max()))


def channel_axis(layers: Dict[str, object], name: str = "", crs=None, project=None) -> Optional[np.ndarray]:
    """Vertices of the channel axis line used by a bed formula in s and n.

    The axis of channel ``name`` in the channels layer (role axis), or the only
    axis there when no name is given.
    """
    src = normalize_layers(layers).get("channels")
    if src is None:
        return None
    crs = crs if crs is not None and crs.isValid() else _crs_of(src)
    found = []
    for f, g in _features(src, crs, project):
        role = str(_value(f, "role", "")).strip().lower()
        if role not in ("axis", "thalweg", "centre", "center", "centerline", "centreline", "a", "t"):
            continue
        parts = line_parts(g)
        if parts:
            found.append((str(_value(f, "channel", "")), max(parts, key=len)))
    if name:
        found = [c for c in found if c[0] == name] or found
    if not found:
        return None
    if len(found) > 1 and not name:
        raise ValueError("the bed formula uses s or n but the project has %d channel axes: set the setting "
                         "bottom_axis to the channel name" % len(found))
    return np.asarray(found[0][1], dtype=np.float64)[:, :2]


# ---------------------------------------------------------------------------
# roughness zones
# ---------------------------------------------------------------------------
def zones_from_layer(src, crs, project=None):
    """[(rings, material, friction or None, name)] in feature order (later zones win)."""
    zones = []
    for k, (f, g) in enumerate(_features(src, crs, project)):
        mid = _num(f, "material")
        fr = _num(f, "friction")
        name = str(_value(f, "name", "") or "")
        mid = int(mid) if mid is not None else k + 2
        for part in polygon_parts(g):
            zones.append((part, mid, fr, name or "zone %d" % mid))
    return zones
# ---------------------------------------------------------------------------
# interpolation from QGIS layers
# ---------------------------------------------------------------------------
_DTYPES = {}
for _n, _d in (("Byte", "u1"), ("Int8", "i1"), ("UInt16", "u2"), ("Int16", "i2"), ("UInt32", "u4"), ("Int32", "i4"),
               ("Float32", "f4"), ("Float64", "f8")):
    _t = getattr(Qgis.DataType, _n, None)
    if _t is not None:
        _DTYPES[_t] = _d


class RasterSampler:
    """Bilinear sampling of a raster band at mesh nodes, read tile by tile.

    Points are transformed to the raster CRS when it differs from the mesh CRS.
    ``scale`` converts vertical units (for example 0.3048 for feet).
    """

    def __init__(self, layer, band: int = 1, mesh_crs: Optional[QgsCoordinateReferenceSystem] = None,
                 scale: float = 1.0, tile: int = 2048, project=None):
        self.layer = layer
        self.band = int(band)
        self.scale = float(scale)
        self.tile = int(tile)
        self.name = layer.name()
        self.ct = None
        if mesh_crs is not None and mesh_crs.isValid() and layer.crs().isValid() and layer.crs() != mesh_crs:
            self.ct = QgsCoordinateTransform(mesh_crs, layer.crs(), (project or QgsProject.instance()).transformContext())

    def _block(self, c0, r0, c1, r1):
        lyr = self.layer
        ext = lyr.extent()
        xres, yres = lyr.rasterUnitsPerPixelX(), lyr.rasterUnitsPerPixelY()
        rect = QgsRectangle(ext.xMinimum() + c0 * xres, ext.yMaximum() - r1 * yres,
                            ext.xMinimum() + c1 * xres, ext.yMaximum() - r0 * yres)
        prov = lyr.dataProvider()
        blk = prov.block(self.band, rect, c1 - c0, r1 - r0)
        dt = _DTYPES.get(blk.dataType())
        if dt is None:
            raise ValueError("unsupported raster data type in %s" % self.name)
        raw = bytes(blk.data())
        arr = np.frombuffer(raw, dtype=np.dtype(dt)).astype(np.float64)
        arr = arr[: (r1 - r0) * (c1 - c0)].reshape(r1 - r0, c1 - c0)
        if blk.hasNoDataValue():
            arr = np.where(arr == blk.noDataValue(), np.nan, arr)
        elif prov.sourceHasNoDataValue(self.band):
            arr = np.where(arr == prov.sourceNoDataValue(self.band), np.nan, arr)
        return arr, rect

    def sample(self, px, py) -> np.ndarray:
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        if self.ct is not None:
            rx = np.empty_like(px)
            ry = np.empty_like(py)
            for i in range(px.size):
                p = self.ct.transform(QgsPointXY(float(px[i]), float(py[i])))
                rx[i], ry[i] = p.x(), p.y()
        else:
            rx, ry = px, py
        lyr = self.layer
        ext = lyr.extent()
        xres, yres = lyr.rasterUnitsPerPixelX(), lyr.rasterUnitsPerPixelY()
        W, H = lyr.width(), lyr.height()
        col = (rx - ext.xMinimum()) / xres
        row = (ext.yMaximum() - ry) / yres
        out = np.full(px.shape, np.nan)
        inside = (col >= 0) & (col <= W) & (row >= 0) & (row <= H)
        if not inside.any():
            return out
        ic = np.clip(np.floor(col).astype(np.int64), 0, W - 1)
        ir = np.clip(np.floor(row).astype(np.int64), 0, H - 1)
        tc, tr = ic // self.tile, ir // self.tile
        keys = np.unique(np.stack([tr[inside], tc[inside]], axis=1), axis=0)
        for kr, kc in keys:
            sel = np.flatnonzero(inside & (tr == kr) & (tc == kc))
            c0 = max(int(kc) * self.tile - 1, 0)
            r0 = max(int(kr) * self.tile - 1, 0)
            c1 = min((int(kc) + 1) * self.tile + 1, W)
            r1 = min((int(kr) + 1) * self.tile + 1, H)
            arr, rect = self._block(c0, r0, c1, r1)
            out[sel] = bilinear(arr, rect.xMinimum(), rect.yMaximum(), xres, -yres, rx[sel], ry[sel])
        return out * self.scale


def point_source(layer, field_name: Optional[str], mesh_crs, scale: float = 1.0, project=None) -> PointSource:
    """Survey points (value from a field, or the geometry Z) triangulated into a TIN."""
    xy, z = [], []
    for f, g in _features(layer, mesh_crs, project):
        if field_name:
            v = _num(f, field_name)
            if v is None:
                continue
            for p in point_parts(g):
                xy.append(p)
                z.append(v)
        else:
            for part in (g.asGeometryCollection() if g.isMultipart() else [g]):
                pt = part.constGet()
                if hasattr(pt, "z") and QgsWkbTypes.hasZ(pt.wkbType()):
                    xy.append([pt.x(), pt.y()])
                    z.append(pt.z())
    if len(xy) < 3:
        raise ValueError("fewer than three survey points with a value")
    xy = np.array(xy, dtype=np.float64)
    ik = triangulate_geos(xy)
    name = layer.name() if hasattr(layer, "name") else (layer.sourceName() if hasattr(layer, "sourceName") else "points")
    return PointSource(xy, np.array(z, dtype=np.float64), ik, scale=scale, name=name)


def interpolate_values(x, y, sources, fill: bool = True):
    """Values at nodes from prioritised sources; returns (values, source index, report lines)."""
    vals, src = interpolate_nodes(x, y, sources)
    report = []
    for k, s in enumerate(sources):
        report.append("%s: %d nodes" % (getattr(s, "name", "source %d" % (k + 1)), int(np.sum(src == k))))
    miss = int(np.sum(src < 0))
    if miss:
        report.append("No source covers %d nodes%s" % (miss, ", filled from the nearest node" if fill else ""))
        if fill:
            vals = fill_nearest(x, y, vals)
    return vals, src, report


# ---------------------------------------------------------------------------
# boundary conditions
# ---------------------------------------------------------------------------
def boundary_lines_from_layer(src, crs, type_field: str = "type", name_field: str = "name",
                              value_field: Optional[str] = "value", project=None):
    lines, values, problems = [], {}, []
    for f, g in _features(src, crs, project):
        t = _value(f, type_field, None)
        name = str(_value(f, name_field, "") or "") or "boundary %d" % f.id()
        try:
            codes = parse_type(t)
        except ValueError as exc:
            problems.append("Feature %d: %s" % (f.id(), exc))
            continue
        for c in line_parts(g):
            lines.append(BoundaryLine(c, codes, name))
        if value_field:
            v = _num(f, value_field)
            if v is not None:
                values[name] = v
    return lines, values, problems


def make_boundary_files(x, y, ikle, lines, values: Dict[str, float], cli_path: str, liq_path: Optional[str] = None,
                        tolerance: Optional[float] = None):
    """Write the .cli (and optionally a .liq template); returns (table, report, cas_lines)."""
    from ..telemac_io.cli import liquid_boundaries

    table, report = build_boundary_table(x, y, ikle, lines, tolerance)
    write_cli(cli_path, table)
    lb = liquid_boundaries(table)
    q, h = [], []
    liq_values = {}
    for b in lb:
        name = table.comments[b["start"]]
        v = values.get(name)
        hcode, ucode, _ = b["codes"]
        q.append("%g" % v if (v is not None and ucode == 5) else "0.")
        h.append("%g" % v if (v is not None and hcode == 5 and ucode != 5) else "0.")
        if v is not None:
            if ucode == 5:
                liq_values["Q(%d)" % b["number"]] = v
            elif hcode == 5:
                liq_values["SL(%d)" % b["number"]] = v
    cas = []
    if lb:
        cas.append("PRESCRIBED FLOWRATES = " + ";".join(q))
        cas.append("PRESCRIBED ELEVATIONS = " + ";".join(h))
    if liq_path:
        with open(liq_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(liq_template(table, values=liq_values))
    return table, report, cas


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
def write_geometry(path: str, x, y, ikle, crs=None, variables: Optional[Dict[str, np.ndarray]] = None,
                   units: Optional[Dict[str, str]] = None, title: str = "TELEMAC Toolkit mesh") -> str:
    """Geometry SELAFIN (BOTTOM and any other variables) plus a .prj sidecar."""
    variables = dict(variables or {})
    if not variables:
        variables = {"BOTTOM": np.zeros(len(x))}
    units = dict(units or {})
    for k in variables:
        units.setdefault(k, "M" if k.upper() in ("BOTTOM", "FOND") else "")
    write_selafin(path, np.asarray(x, float), np.asarray(y, float), np.asarray(ikle, np.int64), variables,
                  units=units, title=title)
    write_prj(path, crs)
    return path


def style_geometry_mesh(layer: QgsMeshLayer, edges: bool = True):
    """Mesh edges visible, bottom elevation on the terrain ramp."""
    from . import styles
    from ..telemac_io.variables import STYLE_BOTTOM

    rs = layer.rendererSettings()
    ms = rs.nativeMeshSettings()
    ms.setEnabled(edges)
    ms.setColor(QColor("#3a3f4b"))
    ms.setLineWidth(0.15)
    rs.setNativeMeshSettings(ms)
    layer.setRendererSettings(rs)
    if layer.dataProvider() is not None and layer.dataProvider().datasetGroupCount() > 0:
        class _L:
            pass
        holder = _L()
        holder.layer = layer
        styles.style_scalar_group(holder, 0, STYLE_BOTTOM)
        rs = layer.rendererSettings()
        rs.setActiveScalarDatasetGroup(0)
        layer.setRendererSettings(rs)
    layer.triggerRepaint()


# ---------------------------------------------------------------------------
# meshes in, meshes out (any format)
# ---------------------------------------------------------------------------
@dataclass
class ExportOptions:
    fmt: str = "telemac"
    quads: bool = False
    frame: str = "keep"  # keep | lower_left | point
    frame_point: Tuple[float, float] = (0.0, 0.0)
    angle: float = 0.0
    store_offset: bool = True
    double: bool = False
    cli_path: Optional[str] = None
    liq_path: Optional[str] = None


FRAME_MODES = ["keep", "lower_left", "point"]


def mesh_from_result(res: MeshResult, keep_quads: bool, crs=None) -> Mesh:
    m = Mesh.from_result(res, res.quads, keep_quads)
    if crs is not None and crs.isValid():
        m.crs_wkt = crs.toWkt()
    return m


def read_mesh_file(path: str):
    """(Mesh, CRS or None) from a SELAFIN file or any mesh QGIS reads (2DM, UGRID, ...)."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".slf", ".res", ".geo", ".srf", ".sel", ".ser"):
        from ..telemac_io import SelafinFile
        from .crs import sidecar_crs

        with SelafinFile(path) as s:
            h = s.header
            if h.is_3d or h.ndp != 3:
                raise ValueError("a 2D triangular SELAFIN file is required")
            frame = s.read_frame(0) if s.nframes else np.zeros((0, h.npoin))
            names = [n.strip() for n in h.var_names]
            units = [u.strip() for u in h.var_units]
            m = Mesh(h.x.copy(), h.y.copy(), h.ikle.copy())
            for i, n in enumerate(names):
                if i >= len(frame):
                    break
                if n.upper() in ("BOTTOM", "FOND") and m.z is None:
                    m.z = frame[i].copy()
                else:
                    m.node_vars[n] = (units[i], frame[i].copy())
        crs = sidecar_crs(path)
        if crs is not None:
            m.crs_wkt = crs.toWkt()
        return m, crs
    from qgis.core import QgsMesh

    lyr = QgsMeshLayer(path, "mesh", "mdal")
    if not lyr.isValid():
        raise ValueError("QGIS cannot read %s as a mesh" % path)
    qm = QgsMesh()
    lyr.dataProvider().populateMesh(qm)
    n = qm.vertexCount()
    x = np.array([qm.vertex(i).x() for i in range(n)])
    y = np.array([qm.vertex(i).y() for i in range(n)])
    z = np.array([qm.vertex(i).z() for i in range(n)])
    tri, quads = [], []
    for k in range(qm.faceCount()):
        f = list(qm.face(k))
        if len(f) == 3:
            tri.append(f)
        elif len(f) == 4:
            quads.append(f)
        else:
            raise ValueError("element %d has %d nodes; only triangles and quadrilaterals are supported" % (k + 1, len(f)))
    m = Mesh(x, y, np.array(tri, dtype=np.int64).reshape(-1, 3), np.array(quads, dtype=np.int64).reshape(-1, 4),
             z=np.where(np.isfinite(z), z, 0.0))
    crs = lyr.crs() if lyr.crs().isValid() else None
    if crs is None:
        from .crs import sidecar_crs

        crs = sidecar_crs(path)  # the .prj written next to the mesh
    if crs is not None:
        m.crs_wkt = crs.toWkt()
    return m, crs


def resolve_frame(mesh: Mesh, opts: ExportOptions) -> Frame:
    if opts.frame == "keep" and not opts.angle:
        return Frame()
    if opts.frame == "lower_left":
        o = (float(np.min(mesh.x)), float(np.min(mesh.y)))
    elif opts.frame == "point":
        o = (float(opts.frame_point[0]), float(opts.frame_point[1]))
    else:
        o = (0.0, 0.0)
    if opts.fmt == "telemac" and opts.store_offset and not opts.angle:
        o = (float(math.floor(o[0])), float(math.floor(o[1])))  # SELAFIN stores the offset as integers
    return Frame(o, float(opts.angle or 0.0))


def _frame_geom(g: QgsGeometry, frame: Frame) -> QgsGeometry:
    if frame.is_identity:
        return g
    from qgis.PyQt.QtGui import QTransform

    a = math.radians(frame.angle_deg)
    c, s = math.cos(a), math.sin(a)
    x0, y0 = frame.origin
    g = QgsGeometry(g)
    g.transform(QTransform(c, -s, s, c, -(x0 * c + y0 * s), x0 * s - y0 * c))
    return g


def liquid_values(table, values: Dict[str, float]):
    """(.liq column values, steering-file lines) from boundary names and values."""
    from ..telemac_io.cli import liquid_boundaries

    lb = liquid_boundaries(table)
    q, h, liq = [], [], {}
    for b in lb:
        name = table.comments[b["start"]]
        v = values.get(name)
        hcode, ucode, _ = b["codes"]
        q.append("%g" % v if (v is not None and ucode == 5) else "0.")
        h.append("%g" % v if (v is not None and hcode == 5 and ucode != 5) else "0.")
        if v is not None:
            if ucode == 5:
                liq["Q(%d)" % b["number"]] = v
            elif hcode == 5:
                liq["SL(%d)" % b["number"]] = v
    cas = []
    if lb:
        cas = ["PRESCRIBED FLOWRATES = " + ";".join(q), "PRESCRIBED ELEVATIONS = " + ";".join(h)]
    return liq, cas


def write_ugrid(path: str, mesh: Mesh, crs) -> List[str]:
    """UGRID NetCDF through QGIS's MDAL (2DM written first, read back, saved as UGRID)."""
    import tempfile

    from qgis.core import QgsMesh, QgsProviderRegistry

    from ..mesher.writers import write_2dm

    with tempfile.TemporaryDirectory() as d:
        tmp = os.path.join(d, "mesh.2dm")
        write_2dm(tmp, mesh)
        lyr = QgsMeshLayer(tmp, "tmp", "mdal")
        qm = QgsMesh()
        lyr.dataProvider().populateMesh(qm)
        md = QgsProviderRegistry.instance().providerMetadata("mdal")
        ok = md.createMeshData(qm, path, "Ugrid", crs if crs is not None else QgsCoordinateReferenceSystem())
        del lyr
    if not ok or not os.path.exists(path):
        raise RuntimeError("MDAL could not write %s" % path)
    return [path]


def write_hecras(path: str, drawing: Dict[str, object], crs, frame: Frame, default_size: float, project=None) -> List[str]:
    """GeoPackage for RAS Mapper: 2D area perimeter, breaklines, refinement regions, boundary condition lines."""
    drawing = normalize_layers(drawing or {})
    if drawing.get("domain") is None:
        raise ValueError("the HEC-RAS export needs the drawing (domain, lines, refinement areas), not only a mesh")
    out_crs = crs if frame.is_identity else local_crs()
    ctx = (project or QgsProject.instance()).transformContext()
    if os.path.exists(path):
        os.remove(path)
    specs = []
    per = QgsVectorLayer("Polygon", "perimeter", "memory")
    per.setCrs(out_crs)
    per.dataProvider().addAttributes([_qfield("name", "string"), _qfield("cell_size", "double")])
    per.updateFields()
    isl = [g for _f, g in _features(drawing.get("islands"), crs, project)]
    for f, g in _features(drawing["domain"], crs, project):
        if isl:
            g = g.difference(QgsGeometry.unaryUnion(isl))
        nf = QgsFeature(per.fields())
        nf.setGeometry(_frame_geom(g, frame))
        nf.setAttributes([str(_value(f, "name", "2D area %d" % f.id())), float(_num(f, "size") or default_size)])
        per.dataProvider().addFeature(nf)
    specs.append(per)
    brk = QgsVectorLayer("LineString", "breaklines", "memory")
    brk.setCrs(out_crs)
    brk.dataProvider().addAttributes([_qfield("name", "string"), _qfield("kind", "string"), _qfield("cell_spacing", "double")])
    brk.updateFields()
    for src, forced in _constraint_sources(drawing):
        for f, g in _features(src, crs, project):
            for c in line_parts(g):
                nf = QgsFeature(brk.fields())
                nf.setGeometry(_frame_geom(QgsGeometry.fromPolylineXY([QgsPointXY(*p) for p in c]), frame))
                nf.setAttributes([str(_value(f, "name", "")), forced or str(_value(f, "mode", "breakline")),
                                  float(_num(f, "size") or default_size)])
                brk.dataProvider().addFeature(nf)
    for f, g in _features(drawing.get("channels"), crs, project):
        if str(_value(f, "role", "")).lower() in ("left", "right"):
            for c in line_parts(g):
                nf = QgsFeature(brk.fields())
                nf.setGeometry(_frame_geom(QgsGeometry.fromPolylineXY([QgsPointXY(*p) for p in c]), frame))
                nf.setAttributes(["%s %s bank" % (_value(f, "channel", ""), _value(f, "role", "")), "channel bank",
                                  float(_num(f, "spacing") or default_size)])
                brk.dataProvider().addFeature(nf)
    specs.append(brk)
    ref = QgsVectorLayer("Polygon", "refinement_regions", "memory")
    ref.setCrs(out_crs)
    ref.dataProvider().addAttributes([_qfield("name", "string"), _qfield("cell_size", "double")])
    ref.updateFields()
    for key in ("refinement_areas", "grid_zones"):
        for f, g in _features(drawing.get(key), crs, project):
            s = _num(f, "size") or _num(f, "dx")
            if not s:
                continue
            nf = QgsFeature(ref.fields())
            nf.setGeometry(_frame_geom(g, frame))
            nf.setAttributes([str(_value(f, "name", "")), float(s)])
            ref.dataProvider().addFeature(nf)
    specs.append(ref)
    bc = QgsVectorLayer("LineString", "bc_lines", "memory")
    bc.setCrs(out_crs)
    bc.dataProvider().addAttributes([_qfield("name", "string"), _qfield("type", "string"), _qfield("value", "double")])
    bc.updateFields()
    for f, g in _features(drawing.get("boundaries"), crs, project):
        for c in line_parts(g):
            nf = QgsFeature(bc.fields())
            nf.setGeometry(_frame_geom(QgsGeometry.fromPolylineXY([QgsPointXY(*p) for p in c]), frame))
            nf.setAttributes([str(_value(f, "name", "")), str(_value(f, "type", "")), _num(f, "value")])
            bc.dataProvider().addFeature(nf)
    specs.append(bc)
    for i, lyr in enumerate(specs):
        _write_layer(lyr, path, lyr.name(), i == 0, ctx)
    return [path]


def export_mesh(path: str, mesh: Mesh, crs, opts: ExportOptions, lines: Optional[List[BoundaryLine]] = None,
                values: Optional[Dict[str, float]] = None, tolerance: Optional[float] = None,
                drawing: Optional[Dict[str, object]] = None, default_size: float = 0.0, project=None):
    """Write ``mesh`` in ``opts.fmt``; returns (files, report lines, steering-file lines)."""
    if opts.fmt not in FORMATS:
        raise ValueError("unknown format %r" % opts.fmt)
    fmt = FORMATS[opts.fmt]
    base, ext = os.path.splitext(path)
    report, cas = [], []
    if fmt.suffix and not base.endswith(fmt.suffix):
        base += fmt.suffix
    if ext.lower() != fmt.ext or path != base + ext:
        path = base + fmt.ext
        report.append("File name set to %s for the %s format" % (os.path.basename(path), fmt.label))
    values = values or {}
    if not opts.quads or not fmt.quads:
        if mesh.has_quads and not fmt.quads:
            report.append("%s takes triangles only: quadrilaterals are split" % fmt.label)
        mesh = Mesh(mesh.x, mesh.y, mesh.all_triangles(), z=mesh.z, node_vars=mesh.node_vars,
                    material=None if mesh.material is None else np.concatenate(
                        [mesh.material[:len(mesh.triangles)], np.tile(mesh.material[len(mesh.triangles):], 2)]),
                    material_names=mesh.material_names, crs_wkt=mesh.crs_wkt, title=mesh.title) if mesh.has_quads else mesh
    frame = resolve_frame(mesh, opts)
    m = mesh if frame.is_identity else mesh.transformed(frame)
    if not frame.is_identity:
        report.append("Local coordinates: %s" % frame.describe())
    table = None
    if lines and opts.fmt != "hecras":
        lt = lines if frame.is_identity else [BoundaryLine(np.stack(frame.forward(bl.coords[:, 0], bl.coords[:, 1]), axis=1),
                                                          bl.codes, bl.name) for bl in lines]
        table, rep = build_boundary_table(m.x, m.y, m.all_triangles(), lt, tolerance)
        report += rep
        liq, cas = liquid_values(table, values)
    else:
        liq = {}
    stored = opts.fmt == "telemac" and opts.store_offset and not frame.is_identity and not frame.angle_deg
    if opts.fmt == "hecras":
        files = write_hecras(path, drawing, crs, frame, default_size, project)
    elif opts.fmt == "ugrid":
        files = write_ugrid(path, m, crs if (frame.is_identity or is_local_crs(crs)) else None)
    elif opts.fmt == "telemac":
        cli = opts.cli_path or (base + ".cli" if table is not None else None)
        files = write_mesh("telemac", path, m, table, liq, cli_path=cli, liq_path=opts.liq_path,
                           store_offset=opts.store_offset, double=opts.double)
    else:
        proj = ""
        if crs is not None and crs.isValid() and frame.is_identity and not is_local_crs(crs):
            try:
                proj = crs.toWkt(Qgis.CrsWktVariant.Wkt1Esri)
            except Exception:
                proj = crs.toWkt()
        files = write_mesh(opts.fmt, path, m, table, values, projection=proj)
    if opts.fmt != "hecras" and crs is not None and crs.isValid() and (frame.is_identity or stored or is_local_crs(crs)):
        prj = write_prj(path, crs)
        if prj:
            files.append(prj)
    if not frame.is_identity and not stored:
        note = write_frame_note(path, frame, crs.description() if crs is not None else "")
        if note:
            files.append(note)
    return files, report, cas


# ---------------------------------------------------------------------------
# the whole chain, shared by the Processing algorithm and the workflow panel
# ---------------------------------------------------------------------------
@dataclass
class ChainResult:
    files: List[str]
    report: List[str]
    warnings: List[str]
    cas: List[str]
    mesh: Mesh
    stats: dict
    crs: QgsCoordinateReferenceSystem
    path: str


def bottom_sources_from(rasters=(), points=None, point_field=None, formula: str = "", z_factor: float = 1.0,
                        crs=None, project=None, axis=None) -> list:
    sources = [RasterSampler(r, 1, crs, z_factor, project=project) for r in rasters or []]
    if points is not None:
        ps = point_source(points, point_field or None, crs, z_factor, project=project)
        ps.name = "survey points"
        sources.append(ps)
    if formula and formula.strip():
        sources.append(FormulaSource(formula.strip(), "formula %s" % formula.strip(), axis))
    return sources


def generate_and_export(layers: Dict[str, object], path: str, opts: ExportOptions, crs=None, default_size: float = 10.0,
                        growth: float = 1.2, smoothing: int = 8, bottom_sources=(), fill: bool = True,
                        default_friction: Optional[float] = None, tolerance: Optional[float] = None,
                        progress: Optional[Callable[[float, str], bool]] = None, project=None) -> ChainResult:
    """Drawing layers -> mesh -> bottom, roughness, boundaries -> file(s) in the chosen format."""
    layers = normalize_layers(layers)

    def sub(lo, hi):
        if progress is None:
            return None
        return lambda f, msg="": progress(lo + (hi - lo) * f, msg)

    build = build_mesh_input(layers, crs, default_size, growth, smoothing, project)
    res = run_mesher(build, sub(0.0, 0.8))
    mesh = mesh_from_result(res, opts.quads, build.crs)
    report = list(build.info) + list(res.log)
    if bottom_sources:
        if progress:
            progress(0.82, "Bottom elevation")
        z, _src, rep = interpolate_values(mesh.x, mesh.y, list(bottom_sources), fill)
        mesh.z = z
        report += ["Bottom: " + r for r in rep]
    zones = zones_from_layer(layers.get("roughness_zones"), build.crs, project) if layers.get("roughness_zones") is not None else []
    if zones or default_friction is not None:
        report += mesh.assign_zones(zones, default_friction=default_friction)
    lines, values, problems = [], {}, []
    if layers.get("boundaries") is not None:
        lines, values, problems = boundary_lines_from_layer(layers["boundaries"], build.crs, project=project)
    if progress:
        progress(0.9, "Writing %s" % opts.fmt)
    files, rep, cas = export_mesh(path, mesh, build.crs, opts, lines or None, values, tolerance, drawing=layers,
                                  default_size=default_size, project=project)
    report += rep
    st = dict(res.stats)
    if progress:
        progress(1.0, "Done")
    return ChainResult(files, report, list(build.warnings) + problems, cas, mesh, st, build.crs, files[0] if files else path)


def refresh_project_layers(path: str, project=None):
    """Reload the map layers that read from the project GeoPackage (after edits made elsewhere)."""
    project = project or QgsProject.instance()
    target = os.path.normcase(os.path.abspath(path))
    for lyr in project.mapLayers().values():
        if not isinstance(lyr, QgsVectorLayer):
            continue
        src = lyr.source().split("|")[0]
        if os.path.normcase(os.path.abspath(src)) == target:
            lyr.reload()
            lyr.triggerRepaint()


def project_path_of(layer) -> Optional[str]:
    """The mesh project GeoPackage a layer belongs to, if any."""
    if not isinstance(layer, QgsVectorLayer) or layer.providerType() != "ogr":
        return None
    src = layer.source().split("|")[0]
    return src if src.lower().endswith(".gpkg") else None
