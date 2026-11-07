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

from ..mesher import MeshInput, MeshResult, Submesh, generate
from ..mesher.channel import channel_strip
from ..mesher.backends import triangulate_geos
from ..mesher.boundary import BoundaryLine, build_boundary_table, liq_template, parse_codes, parse_tracer, parse_type
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
    "read_settings", "write_settings", "add_flume", "add_river", "add_river_network", "line_arrays", "crs_with_code",
    "suggest_sizes_for", "dem_pixel_size", "add_composition", "add_outline", "read_shape_file", "formula_axis",
    "layers_using", "free_output", "free_name",
    "assign_project_crs", "build_mesh_input", "run_mesher", "RasterSampler", "point_source",
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
    ("channels", "LineString", [("channel", "string"), ("role", "string"), ("nodes_across", "int"), ("spacing", "double"),
                                ("across_spacing", "double")],
     "Structured channels (left, right, axis)", []),
    ("grid_zones", "Polygon", [("name", "string"), ("dx", "double"), ("dy", "double")], "Structured grid zones", ["grids"]),
    ("roughness_zones", "Polygon", [("name", "string"), ("material", "int"), ("friction", "double")], "Roughness zones", []),
    ("boundaries", "LineString", [("name", "string"), ("type", "string"), ("value", "double"), ("tracer", "string"),
                                  ("codes", "string")], "Boundary conditions", []),
]
LAYER_KEYS = [p[0] for p in PROJECT_LAYERS]
_OLD_LAYERS = {"outline": "domain", "density": "refinement_areas", "density_points": "refinement_points",
               "hard_points": "fixed_points", "grids": "grid_zones"}
_FIELD_ALIASES = {"size": ["size", "edge_len"], "nodes_across": ["nodes_across", "cross_nodes"],
                  "spacing": ["spacing", "interval"]}
BOUNDARY_TYPE_NAMES = ["wall", "Q", "H", "QH", "UV", "UVH", "free", "waves"]
#: Blue Kenue's names for the boundary types (attribute forms, hover tips)
BOUNDARY_TYPE_LABELS = [("Closed boundary (wall)", "wall"), ("Open boundary with prescribed Q", "Q"),
                        ("Open boundary with prescribed H", "H"), ("Open boundary with prescribed Q and H", "QH"),
                        ("Open boundary with prescribed UV", "UV"), ("Open boundary with prescribed UV and H", "UVH"),
                        ("Open boundary with incident waves", "waves"), ("Open boundary, free (Thompson)", "free"),
                        ("Custom (codes field)", "custom")]
TRACER_LABELS = [("Closed boundary (wall)", "wall"), ("Open boundary with prescribed tracer", "prescribed"),
                 ("Open boundary with free tracer", "free")]
CONSTRAINT_MODES = ("breakline", "guide")

SETTINGS_DEFAULTS = {
    "default_size": "", "growth": "1.2", "smoothing": "8", "max_nodes": "2000000", "method": "balanced", "min_angle": "30",
    "outline_tolerance": "0",
    "terrain_tolerance": "0", "terrain_min_size": "0", "bottom_formula": "", "bottom_axis": "", "bottom_halfwidth": "",
    "z_factor": "1",
    "format": "telemac", "quads": "0", "frame": "keep", "frame_x": "0", "frame_y": "0", "frame_angle": "0",
    "store_offset": "1", "mesh_file": "",
}


def _qfield(name, kind):
    t = {"int": QMetaType.Type.Int, "double": QMetaType.Type.Double, "string": QMetaType.Type.QString}[kind]
    return QgsField(name, t)


# ---------------------------------------------------------------------------
# CRS without a code (QGIS lists such layers with an empty [ ])
# ---------------------------------------------------------------------------
def _proj_params(proj: str) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for tok in (proj or "").split():
        if tok.startswith("+"):
            k, _, v = tok[1:].partition("=")
            try:
                out[k] = round(float(v), 9)
            except ValueError:
                out[k] = v
    return out


def same_projection(a: QgsCoordinateReferenceSystem, b: QgsCoordinateReferenceSystem) -> bool:
    """Same map projection, parameters and ellipsoid (coordinates do not change between the two)."""
    pa, pb = _proj_params(a.toProj()), _proj_params(b.toProj())
    keys = ("proj", "zone", "south", "lat_0", "lon_0", "lat_1", "lat_2", "lat_ts", "k", "k_0", "x_0", "y_0", "units")

    def ell(p):
        e = p.get("ellps") or {"WGS84": "WGS84", "NAD83": "GRS80", "NAD27": "clrk66"}.get(str(p.get("datum", "")), p.get("datum"))
        return {"WGS84": "GRS80"}.get(e, e)  # WGS 84 and GRS 1980 differ by 0.1 mm in the semi-minor axis

    return all(pa.get(k) == pb.get(k) for k in keys) and ell(pa) == ell(pb)


def crs_with_code(crs: Optional[QgsCoordinateReferenceSystem]) -> Optional[QgsCoordinateReferenceSystem]:
    """An EPSG CRS with the same projection as ``crs`` when ``crs`` has no code (often the case of a CRS read from a
    DEM's own definition, e.g. INEGI's ITRF92 UTM), chosen by the datum named in its definition; None otherwise."""
    if crs is None or not crs.isValid() or crs.authid() or is_local_crs(crs):
        return None
    p = _proj_params(crs.toProj())
    wkt = crs.toWkt().upper().replace("_", "").replace(" ", "")
    zone = None
    if p.get("proj") == "utm":
        zone = int(p.get("zone", 0) or 0)
    elif p.get("proj") == "tmerc" and p.get("k", p.get("k_0")) == 0.9996 and p.get("x_0") == 500000:
        z = (float(p.get("lon_0", 0)) + 183.0) / 6.0
        zone = int(round(z)) if abs(z - round(z)) < 1e-9 else None
    south = "south" in p or p.get("y_0") == 10000000
    codes = []
    if zone:
        if "ITRF92" in wkt or "ITRF1992" in wkt:
            codes = [4473 + zone] if 11 <= zone <= 16 and not south else []
        elif "ITRF2008" in wkt:
            codes = [6355 + zone] if 11 <= zone <= 16 and not south else []
        elif "NAD83" in wkt or "NORTHAMERICAN1983" in wkt:
            codes = [26900 + zone] if not south else []
        elif "WGS84" in wkt or "WGS1984" in wkt:
            codes = [(32700 if south else 32600) + zone]
        elif not south:  # an unnamed datum on GRS 1980: Mexico's current frame, then the older one, then NAD83
            codes = ([6355 + zone, 4473 + zone] if 11 <= zone <= 16 else []) + [26900 + zone]
    elif p.get("proj") == "lcc" and ("ITRF2008" in wkt or "ITRF" not in wkt):
        codes = [6372]  # Mexico ITRF2008 / LCC (INEGI)
    for code in codes:
        cand = QgsCoordinateReferenceSystem("EPSG:%d" % code)
        if cand.isValid() and same_projection(crs, cand):
            return cand
    return None


def assign_project_crs(path: str, crs: QgsCoordinateReferenceSystem) -> int:
    """Give every layer of the mesh project GeoPackage ``crs`` (an EPSG CRS with the same projection, so the
    coordinates stay as they are). Returns the number of layers changed."""
    import sqlite3

    code = int(crs.authid().split(":")[1]) if crs.authid().upper().startswith("EPSG:") else None
    if code is None:
        raise ValueError("an EPSG CRS is needed")
    wkt = crs.toWkt(QgsCoordinateReferenceSystem.WktVariant.WKT1_GDAL)
    con = sqlite3.connect(path)
    try:
        cur = con.cursor()
        if not cur.execute("SELECT 1 FROM gpkg_spatial_ref_sys WHERE srs_id = ?", (code,)).fetchone():
            cur.execute("INSERT INTO gpkg_spatial_ref_sys (srs_name, srs_id, organization, organization_coordsys_id, definition, "
                        "description) VALUES (?, ?, 'EPSG', ?, ?, ?)", (crs.description(), code, code, wkt, crs.description()))
        n = cur.execute("UPDATE gpkg_geometry_columns SET srs_id = ?", (code,)).rowcount
        cur.execute("UPDATE gpkg_contents SET srs_id = ? WHERE data_type = 'features'", (code,))
        con.commit()
    finally:
        con.close()
    return n


def ensure_fields(lyr: QgsVectorLayer, key: str) -> None:
    """Add the fields of project layer ``key`` that a project made by an older version lacks."""
    fields = next((p[2] for p in PROJECT_LAYERS if p[0] == key), [])
    missing = [(n, k) for n, k in fields if lyr.fields().indexOf(n) < 0]
    if missing:
        lyr.dataProvider().addAttributes([_qfield(n, k) for n, k in missing])
        lyr.updateFields()


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
            if key in titles:
                try:
                    ensure_fields(lyr, key)  # fields added by later versions (empty for older features)
                except Exception:  # noqa: BLE001 - a read-only file stays as it is
                    pass
            style_project_layer(lyr, key)
            if not lyr.crs().authid() and is_local_crs(lyr.crs()):
                lyr.setCrs(local_crs())  # named in QGIS instead of "Unknown CRS"
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
        cats.append(QgsRendererCategory("custom", QgsLineSymbol.createSimple({"line_color": "#9467bd", "line_width": "1.4",
                                                                              "line_style": "dash"}), "Custom codes"))
        layer.setRenderer(QgsCategorizedSymbolRenderer("lower(\"type\")", cats))
        _label(layer, "coalesce(\"name\", '') || ' [' || coalesce(\"type\", '') || ']'")
    _forms(layer, kind)
    _map_tip(layer, kind)
    layer.triggerRepaint()


# ---------------------------------------------------------------------------
# attribute forms and hover tips
# ---------------------------------------------------------------------------
_NULL_VALUE = "{2839923C-8B7D-419E-B84B-CA2FE9B80EC7}"  # QgsValueMapFieldFormatter.NULL_VALUE

_ALIASES = {
    "name": "Name", "size": "Element size (m)", "edge_len": "Element size (m)", "mode": "Mode",
    "channel": "Channel name", "role": "Role", "nodes_across": "Nodes across (the fewest when they follow the width)",
    "spacing": "Cross-section spacing along (m)", "across_spacing": "Spacing across (m; empty = constant nodes across)",
    "dx": "Cell size along x (m)", "dy": "Cell size along y (m)", "material": "Material id",
    "friction": "Friction coefficient (law in the steering file; Manning n by default)", "type": "Boundary code",
    "value": "Value (m³/s for Q, m for H; optional)", "tracer": "Tracer code (empty = from the boundary code)",
    "codes": "Custom codes: LIHBOR LIUBOR LIVBOR [LITBOR], e.g. 5 4 4 4",
}


def _value_map(layer, field: str, items) -> None:
    from qgis.core import QgsEditorWidgetSetup

    i = layer.fields().indexOf(field)
    if i >= 0:
        layer.setEditorWidgetSetup(i, QgsEditorWidgetSetup("ValueMap", {"map": [{label: value} for label, value in items]}))


def _forms(layer, kind: str) -> None:
    """Friendly field names, and drop-down lists where a field takes a few values (as Blue Kenue's CONLIM editor)."""
    for name, alias in _ALIASES.items():
        i = layer.fields().indexOf(name)
        if i >= 0:
            layer.setFieldAlias(i, alias)
    if kind == "boundaries":
        _value_map(layer, "type", BOUNDARY_TYPE_LABELS)
        _value_map(layer, "tracer", [("Automatic (from the boundary code)", _NULL_VALUE)] + TRACER_LABELS)
    elif kind == "constraint_lines":
        _value_map(layer, "mode", [("Breakline (every vertex kept)", "breakline"), ("Guide line (resampled to the local size)", "guide")])
    elif kind == "channels":
        _value_map(layer, "role", [("Left bank", "left"), ("Right bank", "right"), ("Axis (thalweg)", "axis")])


def _fmt(expr: str, places: int = 1) -> str:
    return "format_number(%s, %d)" % (expr, places)


def _row(label: str, value_expr: str, when: str = "") -> str:
    body = "'%s: ' || %s || '<br>'" % (label, value_expr)
    return "[%% %s %%]" % ("CASE WHEN %s THEN %s ELSE '' END" % (when, body) if when else body)


def _map_tip(layer, kind: str) -> None:
    """Hover text: what each feature is and its measures (area, length, node estimate)."""
    has = layer.fields().indexOf
    size = "\"size\"" if has("size") >= 0 else ("\"edge_len\"" if has("edge_len") >= 0 else "NULL")
    area = "area($geometry)"
    area_txt = "%s || ' m²' || CASE WHEN %s >= 10000 THEN ' (' || %s || ' ha)' ELSE '' END" % (
        _fmt(area), area, _fmt(area + " / 10000", 2))
    nodes = _row("Nodes at this size", "'about ' || %s" % _fmt("1.155 * %s / (%s ^ 2)" % (area, size), 0), "%s > 0" % size)
    unset = "not set: this area is ignored" if kind in ("refinement_areas", "refinement_points") else "the default size of the project"
    size_row = _row("Element size", "CASE WHEN %s > 0 THEN %s || ' m' ELSE '%s' END" % (size, _fmt(size, 2), unset))
    length = _row("Length", "%s || ' m'" % _fmt("length($geometry)"))
    title_of = {"domain": "Domain", "islands": "Island (not meshed)", "refinement_areas": "Refinement area",
                "refinement_points": "Refinement point", "constraint_lines": "Constraint line", "fixed_points": "Fixed point",
                "grid_zones": "Structured grid zone", "roughness_zones": "Roughness zone", "boundaries": "Boundary line"}
    name = "coalesce(\"name\", '')" if has("name") >= 0 else "''"
    title = "[%% CASE WHEN %s <> '' THEN %s ELSE '%s' END %%]" % (name, name, title_of.get(kind, kind))
    rows = []
    if kind == "domain":
        rows = [size_row, _row("Area", area_txt), _row("Perimeter", "%s || ' m'" % _fmt("perimeter($geometry)")), nodes]
    elif kind == "islands":
        rows = [_row("Area", area_txt), _row("Perimeter", "%s || ' m'" % _fmt("perimeter($geometry)"))]
    elif kind == "refinement_areas":
        rows = [size_row, _row("Area", area_txt), nodes]
    elif kind in ("refinement_points", "fixed_points"):
        rows = [size_row, _row("X, Y", "%s || ', ' || %s" % (_fmt("$x", 2), _fmt("$y", 2)))]
    elif kind == "constraint_lines":
        rows = [_row("Mode", "CASE WHEN lower(\"mode\") IN ('guide', 'soft', 'resampled') THEN 'guide line (nodes at the local size)' "
                             "ELSE 'breakline (every vertex kept)' END"),
                _row("Element size along it", "%s || ' m'" % _fmt(size, 2), "%s > 0" % size), length,
                _row("Vertices", "num_points($geometry)")]
    elif kind == "channels":
        title = "[% 'Channel ' || coalesce(\"channel\", '') %]"
        across = ("CASE WHEN \"across_spacing\" > 0 THEN 'follow the width, about ' || %s || ' m apart (at least ' || "
                  "coalesce(to_string(\"nodes_across\"), '3') || ')' ELSE coalesce(to_string(\"nodes_across\"), 'automatic') END"
                  % _fmt("\"across_spacing\"", 2)) if has("across_spacing") >= 0 else "coalesce(to_string(\"nodes_across\"), 'automatic')"
        rows = [_row("Role", "CASE WHEN lower(\"role\") IN ('l', 'left') THEN 'left bank' WHEN lower(\"role\") IN ('r', 'right') "
                             "THEN 'right bank' ELSE 'axis (thalweg)' END"),
                _row("Nodes across", across), _row("Cross-sections every", "%s || ' m'" % _fmt("\"spacing\"", 2), "\"spacing\" > 0"),
                length,
                _row("Cross-sections", "'about ' || to_string(round(length($geometry) / \"spacing\") + 1)", "\"spacing\" > 0")]
    elif kind == "grid_zones":
        rows = [_row("Cells", "%s || ' x ' || %s || ' m'" % (_fmt("\"dx\"", 2), _fmt("coalesce(\"dy\", \"dx\")", 2))),
                _row("Area", area_txt),
                _row("Cells in the zone", "'about ' || %s" % _fmt("%s / (\"dx\" * coalesce(\"dy\", \"dx\"))" % area, 0),
                     "\"dx\" > 0")]
    elif kind == "roughness_zones":
        rows = [_row("Material id", "to_string(\"material\")", "\"material\" IS NOT NULL"),
                _row("Friction coefficient", _fmt("\"friction\"", 3), "\"friction\" IS NOT NULL"), _row("Area", area_txt)]
    elif kind == "boundaries":
        cases = " ".join("WHEN lower(\"type\") = '%s' THEN '%s'" % (code.lower(), label) for label, code in BOUNDARY_TYPE_LABELS)
        codes = " ".join("WHEN lower(\"type\") = '%s' THEN ' (%d %d %d)'" % (n.lower(), *parse_type(n)) for n in BOUNDARY_TYPE_NAMES)
        rows = [_row("Type", "CASE %s ELSE coalesce(\"type\", 'not set: wall') END || CASE %s ELSE '' END" % (cases, codes)),
                _row("Value", "%s || CASE WHEN lower(\"type\") IN ('q', 'qh') THEN ' m³/s' ELSE ' m' END" % _fmt("\"value\"", 3),
                     "\"value\" IS NOT NULL")]
        if has("tracer") >= 0:
            tcases = " ".join("WHEN lower(\"tracer\") = '%s' THEN '%s'" % (code, label) for label, code in TRACER_LABELS)
            rows.append(_row("Tracer", "CASE %s ELSE 'from the boundary code' END" % tcases))
        if has("codes") >= 0:
            rows.append(_row("Custom codes", "\"codes\"", "lower(\"type\") = 'custom'"))
        rows.append(length)
    bg, fg = tip_colors()
    layer.setMapTipTemplate('<div style="background-color:%s;color:%s;padding:3px 6px">'
                            "<b>%s</b><br>%s</div>" % (bg, fg, title, "".join(rows)))
    if has("name") >= 0:
        layer.setDisplayExpression("\"name\"")


def tip_colors():
    """Background and text colour for hover tips. QGIS draws the tip frame in the theme's colours but leaves
    the text black, unreadable on a dark theme: the tip sets both, from the theme, or white and near-black."""
    def lum(c):
        return 0.2126 * c.redF() + 0.7152 * c.greenF() + 0.0722 * c.blueF()

    try:
        from qgis.PyQt.QtGui import QPalette
        from qgis.PyQt.QtWidgets import QApplication

        pal = QApplication.palette()
        bg, fg = pal.color(QPalette.ColorRole.Base), pal.color(QPalette.ColorRole.Text)
        if abs(lum(bg) - lum(fg)) >= 0.5:
            return bg.name(), fg.name()
    except Exception:  # noqa: BLE001 - no application palette (headless): the light pair
        pass
    return "#ffffff", "#1a1a1a"


def refresh_map_tips(layers: Dict[str, object]) -> None:
    """Rewrite the hover tips of project layers already on the map (saved with an older template or theme)."""
    for key, lyr in normalize_layers(layers).items():
        try:
            if key in LAYER_KEYS and isinstance(lyr, QgsVectorLayer):
                _map_tip(lyr, key)
        except Exception:  # noqa: BLE001 - a tip never stops a project from opening
            pass


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
    from ..mesher.geometry import same_point

    return len(c) >= 4 and same_point(c[0], c[-1])


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
        a = _num(f, "across_spacing")
        if a and a > 0:
            grp["across_spacing"] = a
    subs = []
    for cid, grp in groups.items():
        if "left" not in grp and "right" not in grp and "axis" in grp:
            continue  # an axis alone: the reference line of a bed formula in channel coordinates
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
        # with a spacing across, nodes_across is the fewest nodes across
        st = channel_strip(grp["left"], grp["right"], ax, cross_nodes=n, interval=iv, across_spacing=grp.get("across_spacing"),
                           min_nodes=int(grp.get("nodes_across", 3)))
        subs.append(Submesh(st.x, st.y, st.ikle, st.quads))
        lo, hi = int(st.counts.min()), int(st.counts.max())
        across = "%d nodes across" % lo if lo == hi else "%d to %d nodes across (about %.3g m apart)" % (lo, hi, grp["across_spacing"])
        q = st.quality()
        info.append("Channel %s: %s, spacing %.3g, %d nodes, %d cells; cells along/across %.2g to %.2g, cross-sections "
                    "within %.0f deg of normal to the banks" % (cid, across, iv, st.x.size, len(st.quads),
                                                                 q["ratio_min"], q["ratio_max"], q["skew_max"]))
        if q["ratio_min"] < 0.25:
            warnings.append("Channel %s: cells on a bank near (%.1f, %.1f) are %.2g times as long as they are wide (an inner "
                            "bank in a tight bend): turn the cross-sections there (river helper) or use fewer nodes across."
                            % (cid, *q["ratio_min_at"], q["ratio_min"]))
        if q["ratio_max"] > 6.0:
            warnings.append("Channel %s: cells near (%.1f, %.1f) are %.2g times longer along the channel than across: use a "
                            "smaller spacing along or fewer nodes across." % (cid, *q["ratio_max_at"], q["ratio_max"]))
        if q["skew_max"] > 35.0:
            warnings.append("Channel %s: a cross-section near (%.1f, %.1f) is %.0f deg from normal to a bank: smooth the banks "
                            "there or use the river helper." % (cid, *q["skew_at"], q["skew_max"]))
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


def _check_lines_on_outline(layers, out_geom: QgsGeometry, crs, default_size: float, warnings: List[str], project=None):
    """Warn about boundary lines off the domain outline and channel ends short of it.

    Only boundary nodes within half an edge of a boundary line take its type, so a line drawn beside the
    outline can leave a boundary of walls; snapping (panel step 2) puts the lines on it.
    """
    if out_geom is None or out_geom.isEmpty():
        return
    outline = QgsGeometry(out_geom.constGet().boundary())
    tol = max(0.25 * float(default_size or 0.0), 1e-3)
    blines = []
    for f, g in _features(layers.get("boundaries"), crs, project):
        pts = [QgsGeometry.fromPointXY(QgsPointXY(v.x(), v.y())) for v in g.vertices()]
        if not pts:
            continue
        blines.append(g)
        dmax = max(outline.distance(p) for p in pts)
        if dmax > tol:
            warnings.append("Boundary line %r lies up to %.3g m from the domain outline (more than a quarter of the default "
                            "size): snap it to the outline, or its nodes stay walls." % (str(_value(f, "name", f.id())), dmax))
    banks: Dict[str, Dict[str, tuple]] = {}
    across: Dict[str, int] = {}
    for f, g in _features(layers.get("channels"), crs, project):
        role = str(_value(f, "role", "")).strip().lower()
        v = list(g.vertices())
        if role in ("left", "right") and len(v) >= 2:
            name = str(_value(f, "channel", ""))
            banks.setdefault(name, {})[role] = (QgsPointXY(v[0].x(), v[0].y()), QgsPointXY(v[-1].x(), v[-1].y()))
            n = _num(f, "nodes_across")
            if n and n > 2:
                across[name] = int(n)
    for name, b in banks.items():
        if "left" not in b or "right" not in b:
            continue
        (l0, l1), (r0, r1) = b["left"], b["right"]
        if l0.distance(r0) + l1.distance(r1) > l0.distance(r1) + l1.distance(r0):
            r0, r1 = r1, r0  # banks drawn in opposite directions
        for which, a, c in (("first", l0, r0), ("last", l1, r1)):
            da = outline.distance(QgsGeometry.fromPointXY(a))
            dc = outline.distance(QgsGeometry.fromPointXY(c))
            if min(da, dc) <= tol < max(da, dc):
                warnings.append("Channel %r: at one end a bank reaches the domain outline and the other stops %.3g m from it; "
                                "snap both ends to the outline (or both inside, for a channel that ends in the domain)."
                                % (name, max(da, dc)))
            elif max(da, dc) <= tol:
                # an end on the outline: a boundary line should cover its inner nodes
                n = across.get(name, 5)

                def covered(frac):
                    p = QgsGeometry.fromPointXY(QgsPointXY(a.x() + frac * (c.x() - a.x()), a.y() + frac * (c.y() - a.y())))
                    return any(bl.distance(p) <= tol for bl in blines)

                if not covered(0.5):
                    warnings.append("Channel %r: its %s cross-section lies on the domain outline but no boundary line covers "
                                    "it, so it is a closed end (a wall); draw a boundary line there if water flows through."
                                    % (name, which))
                elif not (covered(1.0 / (n - 1)) and covered(1.0 - 1.0 / (n - 1))):
                    warnings.append("Channel %r: the boundary line at its %s cross-section covers only part of it; its "
                                    "nodes next to the banks stay walls. Snap the line to the outline across the whole "
                                    "channel (it may stop short of the two bank nodes)." % (name, which))


def _check_transitions(subs, seg_a: np.ndarray, seg_b: np.ndarray, size, warnings: List[str], ratio: float = 2.5) -> None:
    """Warn where an edge of a structured part is much longer than the elements the size field asks for next to it:
    the triangles there would be stretched (the strip's nodes are fixed)."""
    if not subs or not len(seg_a):
        return
    L = np.hypot(*(seg_b - seg_a).T)
    mid = 0.5 * (seg_a + seg_b)
    h = np.asarray(size(mid[:, 0], mid[:, 1]), dtype=np.float64)
    r = L / np.maximum(h, 1e-300)
    bad = r > ratio
    if bad.any():
        k = int(np.argmax(r))
        warnings.append("%d edges of structured parts are more than %g times the element size next to them (up to %.3g m "
                        "against %.3g m, near (%.1f, %.1f)): the triangles there will be stretched. Use a smaller spacing in "
                        "the strip, or a larger size around it." % (int(bad.sum()), ratio, L[k], h[k], *mid[k]))


def build_mesh_input(layers: Dict[str, object], crs: Optional[QgsCoordinateReferenceSystem], default_size: float,
                     growth: float = 1.2, smoothing: int = 8, project=None, max_nodes: int = 2_000_000,
                     min_angle: Optional[float] = None, terrain=None, outline_tolerance: Optional[float] = None) -> MeshBuild:
    """Read mesh project layers (keys as in PROJECT_LAYERS, or 0.2.0 names); missing ones are skipped.

    ``terrain``: (dem callable, tolerance, smallest size) for terrain-driven sizes (mesher.demsize), or None.
    ``outline_tolerance``: how far (m) the mesh boundary may stray from the drawn outline and lines between
    their kept vertices; None or 0 = 5% of the local element size.
    """
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

    _check_lines_on_outline(layers, out_geom, crs, default_size, warnings, project)

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
    grids = []
    if terrain is not None and terrain[0] is not None and terrain[1] and terrain[1] > 0:
        from ..mesher.demsize import apply_floor, terrain_sizes

        dem, tol, hmin = terrain
        hmin = float(hmin) if hmin and hmin > 0 else default_size / 15.0
        bb = out_geom.boundingBox()
        pad = 2 * default_size
        gd = terrain_sizes(dem, (bb.xMinimum() - pad, bb.yMinimum() - pad, bb.xMaximum() + pad, bb.yMaximum() + pad),
                           float(tol), min(hmin, default_size), default_size, growth)
        if subs:  # elements next to structured parts keep their edge length
            sa, sb, sl = [], [], []
            for sm in subs:
                for lp in boundary_loops(sm.ikle, sm.x, sm.y):
                    a = np.stack([sm.x[lp], sm.y[lp]], axis=1)
                    b = np.roll(a, -1, axis=0)
                    sa.append(a)
                    sb.append(b)
                    sl.append(np.hypot(*(b - a).T))
            gd = apply_floor(gd, np.vstack(sa), np.vstack(sb), np.concatenate(sl), growth)
        grids.append(gd)
        n_fine = int(np.sum(gd.values < 0.75 * default_size))
        info.append("Terrain sizes: tolerance %g m, smallest %.3g m; grid of %d x %d cells of %.3g m, %.0f%% of it refined"
                    % (tol, float(np.min(gd.values)), gd.values.shape[1], gd.values.shape[0], gd.dx,
                       100.0 * n_fine / gd.values.size))
    size = SizeField(default_size, growth, polys_d, pts_d, lines_d, grids)
    _check_transitions(subs, seg_a, seg_b, size, warnings)
    inp = MeshInput(domain=domain, size=size, soft_lines=soft, hard_lines=hard,
                    hard_points=np.array(hp) if hp else None, submeshes=subs, smooth_iterations=int(smoothing),
                    max_nodes=int(max_nodes) if max_nodes else None, min_angle=min_angle or None,
                    outline_tolerance=outline_tolerance or None)
    return MeshBuild(inp, crs, warnings, info, domain_geom)


#: meshing methods: built-in ones, then Gmsh's (optional package)
METHOD_KEYS = ["balanced", "quality", "gmsh_delaunay", "gmsh_frontal", "gmsh_meshadapt", "gmsh_quads"]


def gmsh_missing_message() -> str:
    return ("The Gmsh methods need the gmsh Python package in QGIS's Python. Install it with "
            "'python -m pip install gmsh' (on Windows from the OSGeo4W Shell, on macOS with QGIS's own python3), "
            "then restart QGIS.")


def run_mesher(build: MeshBuild, progress: Optional[Callable[[float, str], bool]] = None, method: str = "balanced") -> MeshResult:
    """Mesh with a method of METHOD_KEYS (the quality-guaranteed angle comes with build.inp.min_angle)."""
    if method.startswith("gmsh_"):
        from ..mesher import gmsh_engine

        if not gmsh_engine.available():
            raise ValueError(gmsh_missing_message())
        return generate(build.inp, triangulate_geos, progress, interior=gmsh_engine.gmsh_interior(method[5:]))
    if method not in ("balanced", "quality"):
        raise ValueError("unknown meshing method %r" % method)
    return generate(build.inp, triangulate_geos, progress)


def domain_rings(layers: Dict[str, object], crs, project=None) -> List[np.ndarray]:
    """Rings (outer and holes) of the domain minus the islands, as drawn: for measuring how well a mesh follows it."""
    lay = normalize_layers(layers)
    geoms = [g for _f, g in _features(lay.get("domain"), crs, project)]
    if not geoms:
        return []
    dom = QgsGeometry.unaryUnion(geoms) if len(geoms) > 1 else geoms[0]
    isl = [g for _f, g in _features(lay.get("islands"), crs, project)]
    if isl:
        dom = dom.difference(QgsGeometry.unaryUnion(isl) if len(isl) > 1 else isl[0])
    out = []
    for part in polygon_parts(dom):
        out += [np.asarray(r, dtype=np.float64) for r in part]
    return out


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


def profile_text(s, v) -> str:
    """'s:v;s:v' for a profile along a channel axis (the project's settings)."""
    return ";".join("%.10g:%.10g" % (a, b) for a, b in zip(np.asarray(s, float), np.asarray(v, float)))


def parse_profile(text: str):
    """(s, v) arrays from 's:v;s:v', or None."""
    pts = [p.split(":") for p in (text or "").split(";") if ":" in p]
    if len(pts) < 2:
        return None
    a = np.array([[float(u), float(w)] for u, w in pts])
    return a[:, 0], a[:, 1]


def formula_axis(layers, settings: Dict[str, str], crs=None, project=None):
    """(axis, half-width profile) for a bed formula in channel coordinates (s, n, b)."""
    axis = channel_axis(layers, settings.get("bottom_axis", ""), crs, project)
    return axis, parse_profile(settings.get("bottom_halfwidth", ""))


def _joint_extent(comp, i: int, end: int, cell: float):
    """(extension of piece i's outline beyond that end, length of its strip to leave out there) for the joint at its
    start (end = 0) or end (end = -1); (0, 0) at the channel's own ends."""
    j = i - 1 if end == 0 else i
    if j < 0 or j >= len(comp.joints):
        return 0.0, 0.0
    joint = comp.joints[j]
    p = comp.pieces[i]
    other = comp.pieces[i - 1] if end == 0 else comp.pieces[i + 1]
    w_here = float(p.widths[0 if end == 0 else -1])
    w_there = float(other.widths[-1 if end == 0 else 0])
    if joint["kind"] == "elbow":
        w = max(w_here, w_there)
        e = min(0.5 * w * math.tan(math.radians(min(abs(joint["angle"]), 170.0)) / 2.0), 3.0 * w)
        return e, e + cell
    # a sudden width change: the narrower piece reaches into the wider one; the wider one's strip stops short
    if w_here < w_there:
        return 0.25 * w_here, 0.0
    return 0.0, max(2.0 * cell, 0.25 * (w_here - w_there))


def add_composition(layers: Dict[str, QgsVectorLayer], comp, cell: float, structured: bool = True,
                    inflow: Optional[float] = None, outflow_level: Optional[float] = None, name: str = "flume") -> List[str]:
    """Write a composed or drawn channel (mesher.composer.Composition) into a mesh project: the domain (the pieces
    joined, with mitred sharp corners and steps), one structured strip per piece (stopped short of corners and
    expansions, which get triangles), the whole axis for the bed formula, and inlet (Q) and outlet (H) lines."""
    report = list(comp.log)
    adder = _Adder(layers)
    polys = []
    trims = []
    for i, p in enumerate(comp.pieces):
        e0, g0 = _joint_extent(comp, i, 0, cell)
        e1, g1 = _joint_extent(comp, i, -1, cell)
        L, R = p.left, p.right
        t0 = np.array([math.cos(p.heading[0]), math.sin(p.heading[0])])
        t1 = np.array([math.cos(p.heading[-1]), math.sin(p.heading[-1])])
        if e0:
            L = np.vstack([L[:1] - e0 * t0, L])
            R = np.vstack([R[:1] - e0 * t0, R])
        if e1:
            L = np.vstack([L, L[-1:] + e1 * t1])
            R = np.vstack([R, R[-1:] + e1 * t1])
        ring = np.vstack([R, L[::-1], R[:1]])
        polys.append(QgsGeometry.fromPolygonXY([_pts(ring)]).buffer(0.0, 1))
        trims.append((g0, g1))
    dom = QgsGeometry.unaryUnion(polys) if len(polys) > 1 else polys[0]
    if dom.isMultipart():
        parts = dom.asMultiPolygon()
        dom = QgsGeometry.fromPolygonXY(max(parts, key=lambda pp: QgsGeometry.fromPolygonXY(pp).area()))
    adder.add("domain", dom, name=name, size=cell)
    report.append("Domain: %.4g m² in %d piece(s)" % (dom.area(), len(comp.pieces)))
    if structured:
        ensure_fields(layers["channels"], "channels")
        for i, (p, (g0, g1)) in enumerate(zip(comp.pieces, trims)):
            step = float(np.median(np.diff(p.s))) if len(p.s) > 1 else cell
            a = int(math.ceil(g0 / step - 1e-9)) if g0 else 0
            b = int(math.ceil(g1 / step - 1e-9)) if g1 else 0
            if len(p.s) - a - b < 3:
                report.append("Warning: piece %d (%.3g m) is too short for a structured strip next to its corners or steps; "
                              "it is meshed with triangles" % (i + 1, p.length))
                continue
            st = p.trimmed(a, b)
            cname = "%s_%d" % (name, i + 1) if len(comp.pieces) > 1 else name
            for role, line in (("left", st.left), ("right", st.right), ("axis", st.axis)):
                adder.add("channels", QgsGeometry.fromPolylineXY(_pts(line)), channel=cname, role=role,
                          nodes_across=MIN_NODES_ACROSS, spacing=cell, across_spacing=cell)
            w = st.widths
            report.append("Structured strip %r: %.4g m, %s wide, nodes across follow the width (about %g m apart)" % (
                cname, st.length, ("%.4g m" % w[0]) if np.ptp(w) < 1e-9 else "%.4g to %.4g m" % (w.min(), w.max()), cell))
    axis_name = "%s_axis" % name
    adder.add("channels", QgsGeometry.fromPolylineXY(_pts(comp.axis)), channel=axis_name, role="axis")
    for k, label, typ, val in ((0, "inlet", "Q", inflow), (-1, "outlet", "H", outflow_level)):
        p = comp.pieces[k]
        a, b = p.right[k], p.left[k]
        W = float(p.widths[k])
        m = min(0.6 * cell, 0.25 * W) / W
        adder.add("boundaries", QgsGeometry.fromPolylineXY(_pts([a + m * (b - a), b - m * (b - a)])), name=label, type=typ, value=val)
    report.append("Inlet (Q) and outlet (H) lines across the end sections")
    return report


def add_outline(layers: Dict[str, QgsVectorLayer], rings: List[np.ndarray], lines: List[np.ndarray], cell: float,
                name: str = "basin") -> List[str]:
    """An imported outline: the largest closed line is the domain, closed lines inside it islands, open lines
    breaklines (walls, steps, edges of features drawn in the file)."""
    if not rings:
        raise ValueError("the file has no closed line for the outline (a polyline that ends where it starts)")
    adder = _Adder(layers)
    geoms = sorted((QgsGeometry.fromPolygonXY([_pts(r)]).makeValid() for r in rings), key=lambda g: -g.area())
    dom = geoms[0]
    adder.add("domain", dom, name=name, size=cell)
    report = ["Domain: the largest closed line, %.4g m²" % dom.area()]
    isl = 0
    for g in geoms[1:]:
        if dom.contains(g.pointOnSurface()):
            adder.add("islands", g, name="island %d" % (isl + 1))
            isl += 1
    if isl:
        report.append("Islands: %d closed lines inside the outline" % isl)
    if len(geoms) - 1 > isl:
        report.append("Warning: %d closed lines outside the outline were left out" % (len(geoms) - 1 - isl))
    for k, c in enumerate(lines):
        adder.add("constraint_lines", QgsGeometry.fromPolylineXY(_pts(c)), name="line %d" % (k + 1), mode="breakline")
    if lines:
        report.append("Breaklines: %d open lines" % len(lines))
    report.append("Draw the inflow and outflow lines on the outline in the boundaries layer (tracing, T, follows it)")
    return report


def read_shape_file(path: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """(closed lines, open lines) from a TXT/CSV (x y [width] per line; a line without numbers ends a polyline)
    or a DXF file (lines, polylines and closed polylines; Z is dropped)."""
    from ..mesher.composer import read_xy_file

    if path.lower().endswith(".dxf"):
        lyr = QgsVectorLayer(path, "dxf", "ogr")
        if not lyr.isValid():
            raise ValueError("%s could not be read as a DXF file" % path)
        polys = []
        for f in lyr.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            gt = QgsWkbTypes.geometryType(g.wkbType())
            if gt == Qgis.GeometryType.Line:
                polys += [c for c in line_parts(g)]
            elif gt == Qgis.GeometryType.Polygon:
                for part in polygon_parts(g):
                    polys += [np.asarray(r, dtype=np.float64) for r in part]
    else:
        polys = read_xy_file(path)
    closed, open_ = [], []
    for c in polys:
        c = np.asarray(c, dtype=np.float64)
        xy = c[:, :2]
        if len(xy) >= 4 and np.hypot(*(xy[0] - xy[-1])) <= 1e-9 * max(float(np.ptp(xy)), 1.0):
            closed.append(c)
        else:
            open_.append(c)
    return closed, open_


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


def line_arrays(src, crs, project=None) -> List[np.ndarray]:
    """Every line of a layer or feature source as an (n, 2) array in ``crs``; the parts of a
    multi-line feature that touch end to end are merged first."""
    out = []
    for _f, g in _features(src, crs, project):
        if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Line:
            continue
        if g.isMultipart():
            m = g.mergeLines()
            parts = m.asMultiPolyline() if m.isMultipart() else [m.asPolyline()]
        else:
            parts = [g.asPolyline()]
        for part in parts:
            if len(part) >= 2:
                out.append(np.array([[p.x(), p.y()] for p in part], dtype=np.float64))
    return out


def line_arrays_with(src, crs, project=None, field: Optional[str] = None) -> List[Tuple[np.ndarray, Optional[float]]]:
    """Like line_arrays, with the number in ``field`` of the feature each line comes from (None without one)."""
    out = []
    for f, g in _features(src, crs, project):
        if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Line:
            continue
        v = _num(f, field) if field else None
        if g.isMultipart():
            m = g.mergeLines()
            parts = m.asMultiPolyline() if m.isMultipart() else [m.asPolyline()]
        else:
            parts = [g.asPolyline()]
        for part in parts:
            if len(part) >= 2:
                out.append((np.array([[p.x(), p.y()] for p in part], dtype=np.float64), v))
    return out


def add_river(layers: Dict[str, QgsVectorLayer], reach, cell: float, nodes_across: int = 0, floodplain: float = 200.0,
              floodplain_size: Optional[float] = None, use_domain: bool = False, boundaries: bool = True,
              inflow_span: str = "channel", outflow_span: str = "valley", inflow: Optional[float] = None,
              outflow_level: Optional[float] = None, name: str = "river", project=None, vary_across: bool = False) -> List[str]:
    """Write a river reach (mesher.river.Reach) into a mesh project.

    ``vary_across`` (with ``nodes_across`` 0): the nodes across follow the width, about ``cell`` apart.

    The left bank, right bank and axis go to the channels layer with matched vertices, so the
    structured strip's cross-sections are those of the reach. Without ``use_domain`` the domain is
    the reach widened by ``floodplain`` on each side, with straight ends normal to the axis (the
    end cross-sections lie on them); with it the project's domain is kept and the reach is cut to
    the cross-sections inside it. Inflow (Q, upstream) and outflow (H, downstream) lines go where
    the reach meets the domain outline.
    """
    from ..mesher.river import end_boundary_lines

    report = list(reach.log)
    warnings = list(reach.warnings)

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

    taken = {str(f["channel"]) for f in layers["channels"].getFeatures()} if layers["channels"].fields().indexOf("channel") >= 0 else set()
    base, k = name or "river", 2
    while name in taken or not name:
        name = "%s_%d" % (base, k)
        k += 1
    crs = layers["domain"].crs()
    if use_domain:
        geoms = [g for _f, g in _features(layers["domain"], crs, project)]
        if not geoms:
            raise ValueError("the mesh project has no domain: draw one, or let the river helper create it")
        dom = QgsGeometry.unaryUnion(geoms) if len(geoms) > 1 else geoms[0]
        eng = QgsGeometry.createGeometryEngine(dom.constGet())
        eng.prepareGeometry()

        def inside(p):
            g = QgsGeometry.fromPointXY(QgsPointXY(float(p[0]), float(p[1])))  # keep it alive while GEOS reads it
            return eng.contains(g.constGet())

        ok = np.array([inside(a) and inside(b) for a, b in zip(reach.left, reach.right)], dtype=bool)
        if not ok.any():
            raise ValueError("no cross-section of the reach lies inside the domain")
        # longest run of sections inside, one cross-section of margin from the outline
        runs, start = [], None
        for i, v in enumerate(np.r_[ok, False]):
            if v and start is None:
                start = i
            elif not v and start is not None:
                runs.append((start, i))
                start = None
        a, b = max(runs, key=lambda r: r[1] - r[0])
        if b - a < 3:
            raise ValueError("fewer than three cross-sections of the reach lie inside the domain")
        if a > 0 or b < len(ok):
            report.append("Reach cut to the %d cross-sections inside the domain (of %d)" % (b - a, len(ok)))
            reach = reach.trimmed(np.arange(a, b))
        polys = dom.asMultiPolygon() if dom.isMultipart() else [dom.asPolygon()]
        mid = QgsGeometry.fromPointXY(QgsPointXY(*reach.axis[len(reach.axis) // 2]))
        poly = next((pp for pp in polys if QgsGeometry.fromPolygonXY(pp).contains(mid)), polys[0])
        ring = np.array([[p.x(), p.y()] for p in poly[0]])
    else:
        half = float(np.max(np.maximum(reach.left_w, reach.right_w)))
        axis_g = QgsGeometry.fromPolylineXY(pts(reach.axis))
        dom = axis_g.buffer(half + float(floodplain), 16, Qgis.EndCapStyle.Flat, Qgis.JoinStyle.Round, 2.0)
        if dom.isMultipart():
            parts = dom.asMultiPolygon()
            dom = QgsGeometry.fromPolygonXY(max(parts, key=lambda pp: QgsGeometry.fromPolygonXY(pp).area()))
        add("domain", dom, name=name, size=floodplain_size if floodplain_size else None)
        ring = np.array([[p.x(), p.y()] for p in dom.asPolygon()[0]])
        report.append("Domain: the reach and %g m of floodplain on each side, %.4g km²" % (floodplain, dom.area() / 1e6))
    across, a_sp = _nodes_across(reach, cell, nodes_across, vary_across)
    if a_sp:
        ensure_fields(layers["channels"], "channels")
    for role, line in (("left", reach.left), ("right", reach.right), ("axis", reach.axis)):
        add("channels", QgsGeometry.fromPolylineXY(pts(line)), channel=name, role=role, nodes_across=across, spacing=cell,
            across_spacing=a_sp)
    report.append(_strip_report(name, reach, across, a_sp, cell))
    if boundaries:
        lines = end_boundary_lines(reach, ring, inflow_span, outflow_span, margin=0.6 * cell)
        for lname, typ, coords, dist in lines:
            add("boundaries", QgsGeometry.fromPolylineXY(pts(coords)), name="%s %s" % (name, lname), type=typ,
                value=inflow if typ == "Q" else outflow_level)
            length = float(np.hypot(*np.diff(coords, axis=0).T).sum())
            report.append("Boundary %s (%s): %.4g m of the outline%s" % (lname, typ, length,
                          "" if dist < 0.5 * cell else ", %.3g m from the end cross-section" % dist))
        if use_domain:
            report.append("The boundary lines follow the domain outline at the ends of the reach; check them on the map.")
    return report + ["Warning: " + w for w in warnings]


#: fewest nodes across a river strip whose nodes follow the width: an inflow across it then has 5 liquid nodes
MIN_NODES_ACROSS = 7


def _nodes_across(reach, cell: float, nodes_across: int, vary: bool):
    """(nodes across, spacing across or None): fixed when given; otherwise from the median width, and following
    the width along the reach when ``vary`` (nodes across is then the fewest, MIN_NODES_ACROSS)."""
    W = float(np.median(reach.width))
    if vary and not (nodes_across and nodes_across > 0):
        return MIN_NODES_ACROSS, float(cell)
    across = int(nodes_across) if nodes_across and nodes_across > 0 else max(3, int(round(W / cell)) + 1)
    if across % 2 == 0:
        across += 1
    return across, None


def _strip_report(name, reach, across, a_sp, cell) -> str:
    W = float(np.median(reach.width))
    if a_sp:
        from ..mesher.channel import _counts

        c = np.maximum(2 * _counts(reach.width, 2, a_sp, 2) + 1, across)
        return ("Structured channel %r: nodes across follow the width (%d to %d, about %.3g m apart), cross-sections every %.3g m"
                % (name, int(c.min()), int(c.max()), a_sp, cell))
    return ("Structured channel %r: %d nodes across (%.3g m apart at the median width %.3g m), cross-sections every %.3g m"
            % (name, across, W / (across - 1), W, cell))


def _junction_size(net, J, cell: float, nodes_across: int, vary: bool) -> float:
    """Element size in a junction: the finest strip spacing (along or across) of the branches there."""
    from ..mesher.channel import _counts

    h = float(cell)
    for r, e in J.members:
        st = net.reaches[r]
        across, a_sp = _nodes_across(st, cell, nodes_across, vary)
        k = 0 if e == 0 else -1
        if a_sp:
            across = max(int(2 * _counts(st.width, 2, a_sp, 2)[k] + 1), across)
        h = min(h, float(st.width[k]) / max(across - 1, 1))
    return h


def add_river_network(layers: Dict[str, QgsVectorLayer], net, cell: float, nodes_across: int = 0, floodplain: float = 200.0,
                      floodplain_size: Optional[float] = None, boundaries: bool = True, inflow_span: str = "channel",
                      outflow_span: str = "valley", inflow: Optional[float] = None, outflow_level: Optional[float] = None,
                      name: str = "river", vary_across: bool = True) -> List[str]:
    """Write a river network (mesher.network.Network) into a mesh project.

    One structured channel per reach (cut short of its junctions); per junction a refinement area and its banks as
    guide lines, so the junction is meshed with triangles of the strips' size; the domain is every reach widened
    by ``floodplain``, with straight ends normal to the axis at the free ends; inflow (Q) lines at the upstream
    ends and outflow (H) lines at the downstream ones.
    """
    from ..mesher.river import end_line

    report = list(net.log)
    warnings = list(net.warnings)
    adder = _Adder(layers)
    taken = {str(f["channel"]) for f in layers["channels"].getFeatures()} if layers["channels"].fields().indexOf("channel") >= 0 else set()
    base = name or "river"
    names = []
    k = 1
    for _r in net.reaches:
        while "%s_%d" % (base, k) in taken:
            k += 1
        names.append("%s_%d" % (base, k))
        k += 1
    # domain: each reach widened, round at the junctions, cut straight at the free ends
    parts = []
    free = {}
    for r, e, kind in net.ends:
        free.setdefault(r, []).append((e, kind))
    for r, full in enumerate(net.full):
        half = float(np.max(np.maximum(full.left_w, full.right_w)))
        R = half + float(floodplain)
        g = QgsGeometry.fromPolylineXY(_pts(full.axis)).buffer(R, 16, Qgis.EndCapStyle.Round, Qgis.JoinStyle.Round, 2.0)
        for e, _kind in free.get(r, []):
            p = full.axis[e]
            t = full.axis[0] - full.axis[1] if e == 0 else full.axis[-1] - full.axis[-2]
            t = t / max(float(np.hypot(*t)), 1e-300)
            nrm = np.array([-t[1], t[0]])
            h = R * 1.5 + 1.0
            box = [p + nrm * h, p + nrm * h + t * h, p - nrm * h + t * h, p - nrm * h, p + nrm * h]
            g = g.difference(QgsGeometry.fromPolygonXY([_pts(box)]))
        parts.append(g)
    for J in net.junctions:
        parts.append(QgsGeometry.fromPolygonXY([_pts(J.ring)]).buffer(0.0, 1))
    dom = QgsGeometry.unaryUnion(parts)
    if dom.isMultipart():
        polys = dom.asMultiPolygon()
        if len(polys) > 1:
            warnings.append("The widened reaches do not form one area; only the largest is kept")
        dom = QgsGeometry.fromPolygonXY(max(polys, key=lambda pp: QgsGeometry.fromPolygonXY(pp).area()))
    adder.add("domain", dom, name=base, size=floodplain_size if floodplain_size else None)
    report.append("Domain: the network and %g m of floodplain on each side, %.4g km²" % (floodplain, dom.area() / 1e6))
    # strips
    a_sp = None
    for r, st in enumerate(net.reaches):
        across, a_sp = _nodes_across(st, cell, nodes_across, vary_across)
        if a_sp:
            ensure_fields(layers["channels"], "channels")
        for role, line in (("left", st.left), ("right", st.right), ("axis", st.axis)):
            adder.add("channels", QgsGeometry.fromPolylineXY(_pts(line)), channel=names[r], role=role, nodes_across=across,
                      spacing=cell, across_spacing=a_sp)
        report.append(_strip_report(names[r], st, across, a_sp, cell))
    # junctions
    for k, J in enumerate(net.junctions):
        h = _junction_size(net, J, cell, nodes_across, vary_across)
        poly = QgsGeometry.fromPolygonXY([_pts(J.ring)]).buffer(0.5 * h, 4)
        jname = "%s junction %d" % (base, k + 1)
        adder.add("refinement_areas", poly, name=jname, size=h)
        for b in J.banks:
            adder.add("constraint_lines", QgsGeometry.fromPolylineXY(_pts(b)), name=jname + " bank", mode="guide", size=h)
        report.append("Junction %d at (%.1f, %.1f): %s; channel area %.4g m² meshed with triangles of %.3g m, %d bank lines"
                      % (k + 1, J.point[0], J.point[1], ", ".join(names[r] for r, _e in J.members), abs(_ring_area(J.ring)), h,
                         len(J.banks)))
    if boundaries:
        ring = np.array([[p.x(), p.y()] for p in dom.asPolygon()[0]])
        n_in = sum(1 for _r, _e, kind in net.ends if kind == "inflow")
        n_out = len(net.ends) - n_in
        for r, e, kind in net.ends:
            span = inflow_span if kind == "inflow" else outflow_span
            coords, dist = end_line(net.reaches[r], e, ring, span, margin=0.6 * cell)
            typ = "Q" if kind == "inflow" else "H"
            if kind == "inflow":
                val = inflow if n_in == 1 else None
            else:
                val = outflow_level
            adder.add("boundaries", QgsGeometry.fromPolylineXY(_pts(coords)), name="%s %s" % (names[r], kind), type=typ, value=val)
            length = float(np.hypot(*np.diff(coords, axis=0).T).sum())
            report.append("Boundary %s %s (%s): %.4g m of the outline%s" % (names[r], kind, typ, length,
                          "" if dist < 0.5 * cell else ", %.3g m from the end cross-section" % dist))
        if n_in > 1 and inflow is not None:
            warnings.append("%d inflows: set the discharge of each one in the boundary layer (the %g m³/s given was not "
                            "split)" % (n_in, inflow))
        if n_out > 1 and outflow_level is not None:
            warnings.append("%d outflows got the same level %g m: check each one" % (n_out, outflow_level))
    return report + ["Warning: " + w for w in warnings]


def _ring_area(r: np.ndarray) -> float:
    x, y = r[:, 0], r[:, 1]
    return 0.5 * float(np.dot(x[:-1], y[1:]) - np.dot(x[1:], y[:-1]))


def _pts(a):
    return [QgsPointXY(float(x), float(y)) for x, y in a]


class _Adder:
    """Adds features with attributes (those the layer has) to project layers, one commit each."""

    def __init__(self, layers):
        self.layers = layers

    def add(self, key, geom, **attrs):
        lyr = self.layers[key]
        f = QgsFeature(lyr.fields())
        f.setGeometry(geom)
        for k, v in attrs.items():
            if lyr.fields().indexOf(k) >= 0:
                f.setAttribute(k, v)
        lyr.startEditing()
        lyr.addFeature(f)
        if not lyr.commitChanges():
            raise RuntimeError("; ".join(lyr.commitErrors()))


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


def dem_pixel_size(layer, crs: Optional[QgsCoordinateReferenceSystem] = None, project=None) -> float:
    """The DEM's pixel size in the units of ``crs`` (its own CRS when none): the mean of x and y."""
    p = 0.5 * (layer.rasterUnitsPerPixelX() + layer.rasterUnitsPerPixelY())
    if crs is None or not crs.isValid() or not layer.crs().isValid() or layer.crs() == crs:
        return float(p)
    ct = QgsCoordinateTransform(layer.crs(), crs, (project or QgsProject.instance()).transformContext())
    ext = ct.transformBoundingBox(layer.extent())
    return float(0.5 * (ext.width() / max(layer.width(), 1) + ext.height() / max(layer.height(), 1)))


def channel_widths(layers: Dict[str, object], crs, project=None) -> List[float]:
    """Median width of each structured channel of the project (distance between its banks)."""
    banks: Dict[str, Dict[str, np.ndarray]] = {}
    for f, g in _features(normalize_layers(layers).get("channels"), crs, project):
        role = str(_value(f, "role", "")).strip().lower()
        parts = line_parts(g)
        if role in ("left", "right") and parts:
            banks.setdefault(str(_value(f, "channel", "")), {})[role] = max(parts, key=len)
    out = []
    for b in banks.values():
        if "left" in b and "right" in b:
            from ..mesher.river import resample

            L = resample(b["left"], 1.0)
            R = resample(b["right"], 1.0)
            f = np.linspace(0, 1, 50)
            li = L[np.minimum((f * (len(L) - 1)).astype(int), len(L) - 1)]
            ri = R[np.minimum((f * (len(R) - 1)).astype(int), len(R) - 1)]
            if np.hypot(*(L[0] - R[0])) > np.hypot(*(L[0] - R[-1])):
                ri = ri[::-1]
            out.append(float(np.median(np.hypot(*(li - ri).T))))
    return out


def suggest_sizes_for(layers: Optional[Dict[str, object]], dem_layer, crs: Optional[QgsCoordinateReferenceSystem] = None,
                      project=None, node_budget: int = 300_000, depth: float = 2.0, velocity: float = 1.5):
    """Sizes suggested from a DEM for a mesh project (its domain and channels), or for the DEM's extent without one.
    Returns a mesher.suggest.SizeSuggestion."""
    from ..mesher.geometry import points_in_rings
    from ..mesher.suggest import suggest_sizes

    crs = crs if crs is not None and crs.isValid() else dem_layer.crs()
    rings = domain_rings(layers, crs, project) if layers else []
    sampler = RasterSampler(dem_layer, 1, crs, project=project)
    if rings:
        allp = np.vstack(rings)
        bbox = (float(allp[:, 0].min()), float(allp[:, 1].min()), float(allp[:, 0].max()), float(allp[:, 1].max()))

        def inside(x, y):
            return points_in_rings(x, y, rings)
    else:
        ext = dem_layer.extent()
        if dem_layer.crs().isValid() and crs.isValid() and dem_layer.crs() != crs:
            ext = QgsCoordinateTransform(dem_layer.crs(), crs, (project or QgsProject.instance()).transformContext()).transformBoundingBox(ext)
        bbox = (ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum())
        inside = None
    widths = channel_widths(layers, crs, project) if layers else []
    sug = suggest_sizes(sampler.sample, bbox, dem_pixel_size(dem_layer, crs, project), inside, widths, node_budget, depth, velocity)
    if not rings:
        sug.notes.append("No domain drawn yet: the area is the DEM's (where it has data).")
    return sug


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
            custom = _value(f, "codes", None)
            tracer = None
            if str(t or "").strip().lower() == "custom" or (not t and custom):
                codes, tracer = parse_codes(custom)
            else:
                codes = parse_type(t)
            tr_text = _value(f, "tracer", None)
            if tr_text not in (None, ""):
                tracer = parse_tracer(tr_text)
        except ValueError as exc:
            problems.append("Feature %d: %s" % (f.id(), exc))
            continue
        for c in line_parts(g):
            lines.append(BoundaryLine(c, codes, name, tracer))
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
        cas.append("VELOCITY PROFILES = " + ";".join("4" if tuple(b["codes"]) == (4, 5, 5) else "1" for b in lb))
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
    if ext == ".t3s":  # Blue Kenue T3 mesh
        from ..telemac_io.bluekenue import read_t3s
        from .crs import sidecar_crs

        x, y, z, ik = read_t3s(path)
        m = Mesh(x, y, ik)
        m.z = z
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
        cas = ["PRESCRIBED FLOWRATES = " + ";".join(q), "PRESCRIBED ELEVATIONS = " + ";".join(h),
               "VELOCITY PROFILES = " + ";".join("4" if tuple(b["codes"]) == (4, 5, 5) else "1" for b in lb)]
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
                                                          bl.codes, bl.name, bl.tracer) for bl in lines]
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
                        crs=None, project=None, axis=None, half_width=None) -> list:
    sources = [RasterSampler(r, 1, crs, z_factor, project=project) for r in rasters or []]
    if points is not None:
        ps = point_source(points, point_field or None, crs, z_factor, project=project)
        ps.name = "survey points"
        sources.append(ps)
    if formula and formula.strip():
        sources.append(FormulaSource(formula.strip(), "formula %s" % formula.strip(), axis, half_width))
    return sources


def generate_and_export(layers: Dict[str, object], path: str, opts: ExportOptions, crs=None, default_size: float = 10.0,
                        growth: float = 1.2, smoothing: int = 8, bottom_sources=(), fill: bool = True,
                        default_friction: Optional[float] = None, tolerance: Optional[float] = None,
                        progress: Optional[Callable[[float, str], bool]] = None, project=None,
                        max_nodes: int = 2_000_000, min_angle: Optional[float] = None, method: str = "balanced",
                        terrain=None, outline_tolerance: Optional[float] = None) -> ChainResult:
    """Drawing layers -> mesh -> bottom, roughness, boundaries -> file(s) in the chosen format."""
    layers = normalize_layers(layers)

    def sub(lo, hi):
        if progress is None:
            return None
        return lambda f, msg="": progress(lo + (hi - lo) * f, msg)

    build = build_mesh_input(layers, crs, default_size, growth, smoothing, project, max_nodes, min_angle, terrain,
                             outline_tolerance)
    res = run_mesher(build, sub(0.0, 0.8), method)
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


def _file_of(layer) -> str:
    """The file behind a layer (mesh, raster or vector), normalised for comparing; '' when none."""
    src = layer.source() or ""
    if src.upper().startswith("SELAFIN:"):
        src = src.split(":", 1)[1]
    src = src.split("|")[0].strip().strip('"')
    if src.count('":') == 1:  # 'DRIVER:"path":mesh'
        src = src.split('":')[0].strip('"')
    return os.path.normcase(os.path.abspath(src)) if src and os.path.splitext(src)[1] else ""


def layers_using(paths, project=None) -> list:
    """Layers of the project that read one of ``paths``."""
    project = project or QgsProject.instance()
    want = {os.path.normcase(os.path.abspath(p)) for p in paths if p}
    return [lyr for lyr in project.mapLayers().values() if _file_of(lyr) in want]


def free_output(paths, project=None) -> List[str]:
    """Remove the layers that read ``paths`` (main thread only): a file loaded in QGIS cannot be overwritten on Windows.
    Returns their names."""
    project = project or QgsProject.instance()
    used = layers_using(paths, project)
    names = [lyr.name() for lyr in used]
    if used:
        project.removeMapLayers([lyr.id() for lyr in used])
    return names


def free_name(path: str) -> str:
    """``path``, or ``stem_2.ext``, ``stem_3.ext``... the first that is not loaded in QGIS."""
    stem, ext = os.path.splitext(path)
    k = 2
    cand = path
    while layers_using([cand]):
        cand = "%s_%d%s" % (stem, k, ext)
        k += 1
    return cand


def project_path_of(layer) -> Optional[str]:
    """The mesh project GeoPackage a layer belongs to, if any."""
    if not isinstance(layer, QgsVectorLayer) or layer.providerType() != "ogr":
        return None
    src = layer.source().split("|")[0]
    return src if src.lower().endswith(".gpkg") else None
