"""Readable names and default styles for TELEMAC dataset groups."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from qgis.core import (
    QgsColorRampShader,
    QgsGradientColorRamp,
    QgsGradientStop,
    QgsMeshDatasetIndex,
    QgsStyle,
)
from qgis.PyQt.QtGui import QColor

from ..telemac_io.variables import (
    STYLE_BOTTOM,
    STYLE_DEPTH,
    STYLE_ELEVATION,
    STYLE_FROUDE,
    STYLE_TIME,
    STYLE_VELOCITY,
    display_name,
    lookup,
)


def _terrain_ramp() -> QgsGradientColorRamp:
    r = QgsGradientColorRamp(QColor("#1d5c3a"), QColor("#f4f1ea"))
    r.setStops([
        QgsGradientStop(0.30, QColor("#7fa95a")),
        QgsGradientStop(0.60, QColor("#d9c88f")),
        QgsGradientStop(0.85, QColor("#a5774f")),
    ])
    return r


def _named_ramp(name: str, fallback: Tuple[str, str]):
    try:
        r = QgsStyle.defaultStyle().colorRamp(name)
    except Exception:
        r = None
    return r if r is not None else QgsGradientColorRamp(QColor(fallback[0]), QColor(fallback[1]))


def ramp_for(style: str):
    if style == STYLE_DEPTH:
        return _named_ramp("Blues", ("#deebf7", "#08306b"))
    if style == STYLE_VELOCITY:
        return _named_ramp("Viridis", ("#440154", "#fde725"))
    if style == STYLE_ELEVATION:
        return _named_ramp("YlGnBu", ("#ffffd9", "#081d58"))
    if style == STYLE_BOTTOM:
        return _terrain_ramp()
    if style == STYLE_TIME:
        return _named_ramp("Magma", ("#000004", "#fcfdbf"))
    return _named_ramp("Spectral", ("#d7191c", "#2b83ba"))


def _items_from_ramp(ramp, vmin: float, vmax: float, n: int = 11):
    items = []
    if vmax <= vmin:
        vmax = vmin + 1e-6
    for k in range(n):
        t = k / (n - 1)
        v = vmin + t * (vmax - vmin)
        items.append(QgsColorRampShader.ColorRampItem(v, ramp.color(t), "%.3g" % v))
    return items


def _froude_items(vmax: float):
    top = max(1.5, vmax)
    return [
        QgsColorRampShader.ColorRampItem(0.0, QColor("#2c7bb6"), "0"),
        QgsColorRampShader.ColorRampItem(0.6, QColor("#abd9e9"), "0.6"),
        QgsColorRampShader.ColorRampItem(1.0, QColor("#f7f7f7"), "1 (critical)"),
        QgsColorRampShader.ColorRampItem(1.0 + 0.5 * (top - 1.0), QColor("#fdae61"), "%.2g" % (1.0 + 0.5 * (top - 1.0))),
        QgsColorRampShader.ColorRampItem(top, QColor("#d7191c"), "%.2g" % top),
    ]


def group_mapping(rs) -> Tuple[Dict[int, int], Dict[Tuple[int, int], int]]:
    """Map SELAFIN variable indices to MDAL dataset group indices.

    MDAL names groups "<name lower-case> <unit>" and merges U/V pairs into one
    vector group named after the common prefix ("velocity").
    """
    dp = rs.layer.dataProvider()
    names = [n.strip().lower() for n in rs.names]
    scalar: Dict[int, int] = {}
    vector: Dict[Tuple[int, int], int] = {}
    for g in range(dp.datasetGroupCount()):
        meta = dp.datasetGroupMetadata(g)
        gname = meta.name().lower().strip()
        if meta.isVector():
            for i, j, _label in rs.pairs:
                a, b = names[i], names[j]
                k = 0
                while k < min(len(a), len(b)) and a[k] == b[k]:
                    k += 1
                prefix = a[:k].strip()
                if prefix and gname.startswith(prefix) and (i, j) not in vector:
                    vector[(i, j)] = g
                    break
        else:
            best, blen = None, -1
            for i, n in enumerate(names):
                if gname.startswith(n) and len(n) > blen and i not in scalar:
                    best, blen = i, len(n)
            if best is not None:
                scalar[best] = g
    return scalar, vector


def rename_groups(rs):
    """Give dataset groups readable names: 'Water depth (m)', 'Velocity (m/s)'."""
    scalar, vector = group_mapping(rs)
    root = rs.layer.datasetGroupTreeRootItem()
    labels = {}
    for i, g in scalar.items():
        labels[g] = display_name(rs.names[i], rs.units[i])
    for (i, j), g in vector.items():
        label = next((lab for a, b, lab in rs.pairs if (a, b) == (i, j)), "Vector")
        u = rs.units[i].strip().lower()
        labels[g] = "%s (%s)" % (label, u) if u else label
    if root is not None:
        for k in range(root.childCount()):
            item = root.child(k)
            g = item.datasetGroupIndex()
            if g in labels:
                item.setName(labels[g])
    rs.group_of_var = scalar
    rs.group_of_pair = vector
    return scalar, vector


def style_scalar_group(rs, g: int, style: str, dry_threshold: Optional[float] = None):
    layer = rs.layer
    meta = layer.dataProvider().datasetGroupMetadata(g)
    vmin, vmax = float(meta.minimum()), float(meta.maximum())
    settings = layer.rendererSettings()
    ss = settings.scalarSettings(g)
    if style == STYLE_FROUDE:
        items = _froude_items(vmax)
        lo, hi = 0.0, items[-1].value
        # a source ramp is required: QGIS reclassifies from it (a missing one crashes QGIS)
        ramp = QgsGradientColorRamp(items[0].color, items[-1].color)
        ramp.setStops([QgsGradientStop((it.value - lo) / (hi - lo), it.color) for it in items[1:-1]])
    else:
        lo = vmin
        if style == STYLE_DEPTH and dry_threshold is not None:
            lo = max(dry_threshold, 0.0)
            vmax = max(vmax, lo + 1e-3)
        hi = vmax if vmax > lo else lo + 1e-6
        ramp = ramp_for(style)
        items = _items_from_ramp(ramp, lo, hi)
    # min/max first: it may reclassify the current shader; our items are set last
    ss.setClassificationMinimumMaximum(lo, hi)
    shader = QgsColorRampShader(lo, hi)
    shader.setSourceColorRamp(ramp.clone())
    shader.setColorRampItemList(items)
    if style == STYLE_DEPTH:
        shader.setClip(True)  # dry nodes (below threshold) are not drawn
    ss.setColorRampShader(shader)
    settings.setScalarSettings(g, ss)
    layer.setRendererSettings(settings)


#: QGIS draws the user grid in output pixels and fixed arrows in millimetres; at the
#: 96 dpi of a screen one millimetre is this many pixels
PX_PER_MM = 96.0 / 25.4


@dataclass
class ArrowStyle:
    """How vector arrows are drawn (velocity, unit discharge, wind ...).

    placement  "grid": one arrow per cell of a regular screen grid, interpolated
               from the mesh, so the density is the same in fine and coarse parts
               of the mesh and at every zoom; "nodes": one arrow per mesh node.
    spacing    grid spacing in screen pixels.
    length     "speed": proportional to the magnitude, the fastest value of the
               whole result gets ``longest`` mm (0 = 90 % of the spacing, so arrows
               never overlap); "scale": proportional with a fixed scale of
               ``scale`` mm per unit (m/s), the same for every result, to compare
               runs; "fixed": every arrow ``longest`` mm long, direction only.
    min_value  arrows below this magnitude are not drawn (dry or still water).
    """

    placement: str = "grid"
    spacing: int = 40
    length: str = "speed"
    longest: float = 0.0
    scale: float = 10.0
    min_value: float = 0.01
    color: str = "#1b1b1b"
    by_value: bool = False
    width: float = 0.3

    KEYS = ("placement", "spacing", "length", "longest", "scale", "min_value", "color", "by_value", "width")

    def longest_mm(self) -> float:
        if self.longest > 0:
            return float(self.longest)
        if self.placement == "grid":
            return 0.9 * self.spacing / PX_PER_MM
        return 5.0

    def length_of(self, value: float, vmax: float) -> float:
        """Arrow length in mm for a magnitude (what the panel shows as reference)."""
        if self.length == "fixed":
            return self.longest_mm()
        if self.length == "scale":
            return self.scale * value
        return self.longest_mm() * value / vmax if vmax > 0 else 0.0

    @classmethod
    def load(cls) -> "ArrowStyle":
        from qgis.core import QgsSettings

        st = cls()
        s = QgsSettings()
        for k in cls.KEYS:
            v = s.value("telemac_toolkit/arrows/" + k, None)
            if v is None:
                continue
            cur = getattr(st, k)
            try:
                if isinstance(cur, bool):
                    v = v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes")
                else:
                    v = type(cur)(v)
            except (TypeError, ValueError):
                continue
            setattr(st, k, v)
        if st.placement not in ("grid", "nodes"):
            st.placement = "grid"
        if st.length not in ("speed", "scale", "fixed"):
            st.length = "speed"
        st.spacing = int(min(max(st.spacing, 8), 400))
        return st

    def save(self):
        from qgis.core import QgsSettings

        s = QgsSettings()
        for k in self.KEYS:
            s.setValue("telemac_toolkit/arrows/" + k, getattr(self, k))


def group_max_magnitude(layer, g: int) -> float:
    try:
        return float(layer.dataProvider().datasetGroupMetadata(g).maximum())
    except Exception:
        return 0.0


def style_vector_group(rs, g: int, style: Optional[ArrowStyle] = None):
    """Arrows of vector group ``g``: placement, length, filter and colour from ``style``.

    Lengths use QGIS units: the grid and the "Scaled" method are in pixels, the
    "Fixed" method in millimetres; ``ArrowStyle`` is in millimetres throughout.
    """
    from qgis.core import QgsMeshRendererVectorArrowSettings

    style = style or ArrowStyle.load()
    layer = rs.layer
    settings = layer.rendererSettings()
    vs = settings.vectorSettings(g)
    vmax = group_max_magnitude(layer, g)
    try:
        vs.setSymbology(vs.Symbology.Arrows)
    except Exception:
        pass
    vs.setOnUserDefinedGrid(style.placement == "grid")
    vs.setUserGridCellWidth(int(style.spacing))
    vs.setUserGridCellHeight(int(style.spacing))
    vs.setFilterMin(float(style.min_value) if style.min_value > 0 else -1.0)
    vs.setFilterMax(-1.0)
    vs.setLineWidth(float(style.width))
    vs.setColor(QColor(style.color))
    try:
        from qgis.core import QgsInterpolatedLineColor

        if style.by_value and vmax > 0:
            ramp = ramp_for(STYLE_VELOCITY)
            shader = QgsColorRampShader(0.0, vmax)
            shader.setSourceColorRamp(ramp.clone())
            shader.setColorRampItemList(_items_from_ramp(ramp, 0.0, vmax))
            vs.setColorRampShader(shader)
            vs.setColoringMethod(QgsInterpolatedLineColor.ColoringMethod.ColorRamp)
        else:
            vs.setColoringMethod(QgsInterpolatedLineColor.ColoringMethod.SingleColor)
    except Exception:
        pass
    arrows = vs.arrowSettings()
    method = QgsMeshRendererVectorArrowSettings.ArrowScalingMethod
    if style.length == "fixed":
        arrows.setShaftLengthMethod(method.Fixed)
        arrows.setFixedShaftLength(style.longest_mm())
    else:
        # "Scaled" multiplies the magnitude into screen pixels (QGIS's MinMax would rescale every
        # frame to that frame's own range, so equal lengths would mean different speeds over time)
        mm_per_unit = style.scale if style.length == "scale" else (style.longest_mm() / vmax if vmax > 0 else 0.0)
        arrows.setShaftLengthMethod(method.Scaled)
        arrows.setScaleFactor(mm_per_unit * PX_PER_MM)
    arrows.setArrowHeadWidthRatio(0.25)
    arrows.setArrowHeadLengthRatio(0.35)
    vs.setArrowsSettings(arrows)
    settings.setVectorSettings(g, vs)
    layer.setRendererSettings(settings)


def apply_arrow_style(rs, style: Optional[ArrowStyle] = None):
    """Restyle the arrows of every vector group of a result."""
    style = style or ArrowStyle.load()
    for g in getattr(rs, "group_of_pair", {}).values():
        style_vector_group(rs, g, style)
    rs.layer.triggerRepaint()


def apply_default_style(rs):
    """Readable names, a style per variable kind, depth shown with velocity arrows."""
    scalar, vector = rename_groups(rs)
    for i, g in scalar.items():
        info = lookup(rs.names[i])
        style_scalar_group(rs, g, info.style if info else "generic", rs.dry_threshold)
    arrows = ArrowStyle.load()
    for (i, j), g in vector.items():
        style_scalar_group(rs, g, STYLE_VELOCITY)
        style_vector_group(rs, g, arrows)
    settings = rs.layer.rendererSettings()
    dv = rs.default_variable()
    if dv in scalar:
        settings.setActiveScalarDatasetGroup(scalar[dv])
    vel = next((g for (i, j), g in vector.items() if rs.names[i].strip().upper() in ("VELOCITY U", "VITESSE U")), None)
    if vel is not None:
        settings.setActiveVectorDatasetGroup(vel)
    rs.layer.setRendererSettings(settings)


def set_active_variable(rs, var):
    """Show a variable (index or (i, j) vector pair) as the coloured surface."""
    settings = rs.layer.rendererSettings()
    g = None
    if isinstance(var, tuple):
        g = getattr(rs, "group_of_pair", {}).get(var)
    else:
        g = getattr(rs, "group_of_var", {}).get(var)
    if g is None:
        return False
    settings.setActiveScalarDatasetGroup(g)
    rs.layer.setRendererSettings(settings)
    rs.layer.triggerRepaint()
    return True


def set_arrows_visible(rs, visible: bool):
    settings = rs.layer.rendererSettings()
    if visible:
        vel = next(iter(getattr(rs, "group_of_pair", {}).values()), None)
        if vel is None:
            return False
        settings.setActiveVectorDatasetGroup(vel)
    else:
        settings.setActiveVectorDatasetGroup(-1)
    rs.layer.setRendererSettings(settings)
    rs.layer.triggerRepaint()
    return True
