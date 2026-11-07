"""Check the drawing of a mesh project before meshing, show the problems on the map and fix the safe ones.

The measures come from ``mesher.drawcheck`` (NumPy), against the element size wanted at each place. QGIS adds
what only it knows: invalid geometries (self-intersections), islands crossing the outline, domain polygons
overlapping each other. The safe fixes (near-duplicate vertices, small spikes, invalid polygons, line ends
snapped onto the outline) are made in the layers' edit buffers, so they can be reviewed and undone.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from qgis.core import (
    Qgis,
    QgsCategorizedSymbolRenderer,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsGeometry,
    QgsMarkerSymbol,
    QgsPointXY,
    QgsProject,
    QgsRendererCategory,
    QgsVectorLayer,
    QgsWkbTypes,
)

from ..mesher.drawcheck import Problem, clean_vertices, narrow_gaps, near_ends, vertex_problems
from . import meshing as M

__all__ = ["DrawingIssue", "DrawingCheck", "check_drawing", "fix_drawing", "issues_layer", "CHECK_LAYER_NAME"]

CHECK_LAYER_NAME = "Drawing check"
#: at most this many problems of one kind on one feature are listed one by one
PER_FEATURE = 20

_NAMES = {"domain": "domain", "islands": "island", "grid_zones": "grid zone", "constraint_lines": "line",
          "channels": "channel line"}


@dataclass
class DrawingIssue:
    problem: Problem
    layer: str          # project layer key
    fid: int = -1
    label: str = ""     # the object, e.g. "island 'pier'"

    @property
    def text(self) -> str:
        m = self.problem.message
        body = m if (self.label and self.label in m) or not self.label else "%s: %s" % (self.label, m)
        return "%s%s (at %.6g, %.6g)" % (body[:1].upper(), body[1:], self.problem.x, self.problem.y)


@dataclass
class DrawingCheck:
    issues: List[DrawingIssue] = field(default_factory=list)
    hidden: Dict[str, int] = field(default_factory=dict)  # problems not listed one by one, per kind
    crs: Optional[QgsCoordinateReferenceSystem] = None

    @property
    def fixable(self) -> List[DrawingIssue]:
        return [i for i in self.issues if i.problem.fixable]

    def summary(self) -> List[str]:
        if not self.issues:
            return ["Drawing check: no problem found."]
        counts: Dict[str, int] = {}
        for i in self.issues:
            counts[i.problem.code] = counts.get(i.problem.code, 0) + 1
        for k, n in self.hidden.items():
            counts[k] = counts.get(k, 0) + n
        words = {"duplicate": "near-duplicate vertices", "spike": "spikes", "sharp": "very sharp corners",
                 "corners": "corners closer than the elements", "breakline": "dense breakline vertices",
                 "neck": "narrow parts of the outline", "gap": "narrow gaps between outlines",
                 "near": "line ends or vertices just off the outline", "invalid": "invalid geometries",
                 "crosses": "islands crossing the outline", "overlap": "overlapping domain polygons"}
        parts = ["%d %s" % (n, words.get(k, k)) for k, n in sorted(counts.items(), key=lambda kv: -kv[1])]
        nfix = len(self.fixable)
        return ["Drawing check: " + ", ".join(parts) + "." + (" %d can be fixed safely (Fix the safe problems)." % nfix
                                                             if nfix else "")]


def _name(f, key: str) -> str:
    try:
        n = f["name"] if f.fields().indexOf("name") >= 0 else None
    except KeyError:
        n = None
    if key == "channels":
        try:
            n = "%s %s" % (f["channel"] or "", f["role"] or "")
        except KeyError:
            pass
    what = _NAMES.get(key, key)
    if key == "constraint_lines":
        try:
            what = "guide line" if str(f["mode"] or "").lower() in ("guide", "soft", "resampled") else "breakline"
        except KeyError:
            what = "breakline"
    return "%s %r" % (what, str(n).strip()) if n not in (None, "") and str(n).strip() and str(n) != "NULL" else \
        "%s %d" % (what, f.id())


def _invalid(g: QgsGeometry):
    """Where a geometry is invalid (self-intersection, ring crossing): [(x, y, what)]."""
    out = []
    try:
        errs = g.validateGeometry(Qgis.GeometryValidationEngine.Geos)
    except (AttributeError, TypeError):
        errs = g.validateGeometry()
    for e in errs:
        p = e.where() if e.hasWhere() else g.centroid().asPoint()
        out.append((p.x(), p.y(), e.what()))
    return out


def check_drawing(layers: Dict[str, object], size: Callable, crs=None, project=None) -> DrawingCheck:
    """Problems of the drawing, measured against ``size`` (the project's size field, f(x, y) -> m)."""
    lay = M.normalize_layers(layers)
    dom_src = lay.get("domain")
    crs = crs if crs is not None and crs.isValid() else (M._crs_of(dom_src) if dom_src is not None else QgsCoordinateReferenceSystem())
    res = DrawingCheck(crs=crs)
    counts: Dict[Tuple[str, int, str], int] = {}

    def add(p: Problem, key: str, fid: int, label: str):
        k = (key, fid, p.code)
        counts[k] = counts.get(k, 0) + 1
        if counts[k] > PER_FEATURE:
            res.hidden[p.code] = res.hidden.get(p.code, 0) + 1
            return
        res.issues.append(DrawingIssue(p, key, fid, label))

    # ---- polygons: invalid geometry, vertices --------------------------------------------------
    dom_geoms = []
    for key in ("domain", "islands", "grid_zones", "refinement_areas", "roughness_zones"):
        for f, g in M._features(lay.get(key), crs, project):
            if QgsWkbTypes.geometryType(g.wkbType()) != Qgis.GeometryType.Polygon:
                continue
            label = _name(f, key) if key in _NAMES else "%s %d" % (key.replace("_", " ")[:-1], f.id())
            for x, y, what in _invalid(g):
                add(Problem("invalid", x, y, "not a valid polygon (%s)" % what, fixable=True), key, f.id(), label)
            if key == "domain":
                dom_geoms.append((f, g))
            if key in ("domain", "islands", "grid_zones"):
                for part in M.polygon_parts(g):
                    for ring in part:
                        for p in vertex_problems(ring, True, size, kind="outline"):
                            add(p, key, f.id(), label)
    # ---- lines --------------------------------------------------------------------------------
    ends, end_labels, end_keys = [], [], []
    for key in ("constraint_lines", "channels"):
        for f, g in M._features(lay.get(key), crs, project):
            label = _name(f, key)
            kind = "breakline" if label.startswith("breakline") else "line"
            for line in M.line_parts(g):
                closed = len(line) > 3 and np.allclose(line[0], line[-1], rtol=0.0, atol=1e-12 * max(1.0, float(np.abs(line).max())))
                c = line[:-1] if closed else line
                for p in vertex_problems(c, closed, size, kind=kind):
                    add(p, key, f.id(), label)
                if not closed:
                    for which, pt in (("start", c[0]), ("end", c[-1])):
                        ends.append(pt)
                        end_labels.append("%s (%s)" % (label, which))
                        end_keys.append((key, f.id()))
    for f, g in M._features(lay.get("grid_zones"), crs, project):
        for part in M.polygon_parts(g):
            for r in part:
                for pt in r:
                    ends.append(pt)
                    end_labels.append("%s (corner)" % _name(f, "grid_zones"))
                    end_keys.append(("grid_zones", f.id()))
    # ---- the area meshed: domain minus islands --------------------------------------------------
    if not dom_geoms:
        return res
    for i in range(len(dom_geoms)):
        for j in range(i + 1, len(dom_geoms)):
            inter = dom_geoms[i][1].intersection(dom_geoms[j][1])
            if inter is not None and not inter.isEmpty() and inter.area() > 0:
                c = inter.centroid().asPoint()
                add(Problem("overlap", c.x(), c.y(), "two domain polygons overlap (%.3g m²): they are merged; draw one "
                                                     "polygon, or let them share an edge" % inter.area()),
                    "domain", dom_geoms[j][0].id(), _name(dom_geoms[j][0], "domain"))
    dom = QgsGeometry.unaryUnion([g for _f, g in dom_geoms]) if len(dom_geoms) > 1 else QgsGeometry(dom_geoms[0][1])
    if not dom.isGeosValid():
        dom = dom.makeValid()
    islands = [(f, g) for f, g in M._features(lay.get("islands"), crs, project)
               if QgsWkbTypes.geometryType(g.wkbType()) == Qgis.GeometryType.Polygon]
    bnd = QgsGeometry(dom.constGet().boundary()) if dom.constGet() is not None else None
    for f, g in islands:
        if bnd is not None and g.intersects(bnd) and not g.touches(bnd):
            inter = g.intersection(bnd)
            c = (inter if inter is not None and not inter.isEmpty() else g).centroid().asPoint()
            add(Problem("crosses", c.x(), c.y(), "crosses the domain outline: only its part inside is cut out"),
                "islands", f.id(), _name(f, "islands"))
    area = dom
    if islands:
        isl = QgsGeometry.unaryUnion([g for _f, g in islands])
        area = dom.difference(isl)
    rings, names = [], []
    for part in M.polygon_parts(area):
        for k, ring in enumerate(part):
            if k == 0:
                rings.append(ring)
                names.append("the domain outline")
                continue
            # name a hole after the island it comes from
            pt = QgsGeometry.fromPointXY(QgsPointXY(*ring[0]))
            nm = "a hole of the domain"
            for f, g in islands:
                if g.buffer(1e-9 * max(1.0, float(np.abs(ring).max())), 2).intersects(pt):
                    nm = _name(f, "islands")
                    break
            rings.append(ring)
            names.append(nm)
    for p in narrow_gaps(rings, names, size):
        key, label = ("islands", p.other) if p.other.startswith("island") else ("domain", "the domain outline")
        add(p, key, -1, label)
    if ends and rings:
        for p in near_ends(np.array(ends), end_labels, rings, size):
            key, fid = end_keys[p.index]
            add(p, key, fid, end_labels[p.index])
    return res


# ---------------------------------------------------------------------------
# fixes
# ---------------------------------------------------------------------------
def _to_layer(lyr, crs, project):
    lcrs = lyr.crs()
    if crs is None or not crs.isValid() or not lcrs.isValid() or lcrs == crs:
        return None, None
    ctx = (project or QgsProject.instance()).transformContext()
    return QgsCoordinateTransform(crs, lcrs, ctx), QgsCoordinateTransform(lcrs, crs, ctx)


def _cleaned(g: QgsGeometry, size: Callable, fwd, back) -> Optional[QgsGeometry]:
    """The geometry without near-duplicate vertices and small spikes (None when nothing changes)."""
    g = QgsGeometry(g)
    if QgsWkbTypes.isCurvedType(g.wkbType()):
        g = QgsGeometry(g.constGet().segmentize())
    if back is not None:
        g.transform(back)  # into the check CRS, where sizes are measured
    changed = False
    gt = QgsWkbTypes.geometryType(g.wkbType())
    if gt == Qgis.GeometryType.Polygon:
        polys = []
        for part in M.polygon_parts(g):
            rings = []
            for r in part:
                c = clean_vertices(r, True, size)
                changed |= len(c) != len(r)
                if len(c) >= 3:
                    rings.append([QgsPointXY(float(a), float(b)) for a, b in np.vstack([c, c[:1]])])
            if rings:
                polys.append(rings)
        if not changed or not polys:
            return None
        out = QgsGeometry.fromMultiPolygonXY(polys) if g.isMultipart() else QgsGeometry.fromPolygonXY(polys[0])
    elif gt == Qgis.GeometryType.Line:
        lines = []
        for line in M.line_parts(g):
            c = clean_vertices(line, False, size)
            changed |= len(c) != len(line)
            if len(c) >= 2:
                lines.append([QgsPointXY(float(a), float(b)) for a, b in c])
        if not changed or not lines:
            return None
        out = QgsGeometry.fromMultiPolylineXY(lines) if g.isMultipart() else QgsGeometry.fromPolylineXY(lines[0])
    else:
        return None
    if fwd is not None:
        out.transform(fwd)
    return out


def fix_drawing(layers: Dict[str, object], check: DrawingCheck, size: Callable, project=None, commit: bool = False) -> List[str]:
    """Make the safe fixes of ``check`` in the layers' edit buffers (``commit``: save them). Returns a report."""
    lay = M.normalize_layers(layers)
    report: List[str] = []
    by_feature: Dict[Tuple[str, int], List[DrawingIssue]] = {}
    for i in check.fixable:
        by_feature.setdefault((i.layer, i.fid), []).append(i)
    touched = set()
    n_clean = n_valid = n_snap = 0
    for (key, fid), items in by_feature.items():
        lyr = lay.get(key)
        if not isinstance(lyr, QgsVectorLayer) or fid < 0:
            continue
        f = lyr.getFeature(fid)
        if not f.isValid():
            continue
        if not lyr.isEditable():
            lyr.startEditing()
        fwd, back = _to_layer(lyr, check.crs, project)
        g = QgsGeometry(f.geometry())
        codes = {i.problem.code for i in items}
        if "invalid" in codes:
            v = g.makeValid()
            if v is not None and not v.isEmpty() and QgsWkbTypes.geometryType(v.wkbType()) == QgsWkbTypes.geometryType(g.wkbType()):
                g = v
                n_valid += 1
        if codes & {"duplicate", "spike"}:
            c = _cleaned(g, size, fwd, back)
            if c is not None:
                g = c
                n_clean += 1
        for i in items:
            if i.problem.code != "near" or i.problem.target is None:
                continue
            px, py = i.problem.x, i.problem.y
            tx, ty = i.problem.target
            if fwd is not None:
                a = fwd.transform(QgsPointXY(px, py))
                b = fwd.transform(QgsPointXY(tx, ty))
                px, py, tx, ty = a.x(), a.y(), b.x(), b.y()
            _pt, at, _b, _a, d2 = g.closestVertex(QgsPointXY(px, py))
            if at >= 0 and d2 <= (1e-7 * max(1.0, abs(px), abs(py))) ** 2:
                if g.moveVertex(tx, ty, at):
                    n_snap += 1
        extra = []
        if g.isMultipart() and not QgsWkbTypes.isMultiType(lyr.wkbType()):
            # a repaired bow-tie is two polygons: the layer holds one per feature, so the others become new
            # features with the same attributes
            parts = [QgsGeometry(p.clone()) for p in g.constParts()]
            g = parts[0]
            for extra_g in parts[1:]:
                nf = QgsFeature(f)
                nf.setGeometry(extra_g)
                extra.append(nf)
        if lyr.changeGeometry(fid, g):
            touched.add(lyr)
        for nf in extra:
            nf.setId(-1)
            if "fid" in [fl.name() for fl in lyr.fields()]:
                nf.setAttribute("fid", None)
            lyr.addFeature(nf)
    if n_clean:
        report.append("Near-duplicate vertices and small spikes removed from %d objects." % n_clean)
    if n_valid:
        report.append("%d invalid polygons repaired." % n_valid)
    if n_snap:
        report.append("%d line ends or vertices snapped onto the outline." % n_snap)
    for lyr in touched:
        if commit:
            lyr.commitChanges(False)
        lyr.triggerRepaint()
    if touched and not commit:
        report.append("The changes are in the layers' edit buffers: review them and save the layers (or undo).")
    if not report:
        report.append("Nothing to fix safely.")
    return report


# ---------------------------------------------------------------------------
# the problems on the map
# ---------------------------------------------------------------------------
def issues_layer(check: DrawingCheck, name: str = CHECK_LAYER_NAME) -> QgsVectorLayer:
    """A memory point layer with one point per problem (fixable ones orange, the others red)."""
    crs = check.crs.authid() if check.crs is not None and check.crs.authid() else ""
    lyr = QgsVectorLayer("Point" + ("?crs=%s" % crs if crs else ""), name, "memory")
    if not crs and check.crs is not None and check.crs.isValid():
        lyr.setCrs(check.crs)
    pr = lyr.dataProvider()
    pr.addAttributes([M._qfield("problem", "string"), M._qfield("object", "string"), M._qfield("message", "string"),
                      M._qfield("distance", "double"), M._qfield("size", "double"), M._qfield("fixable", "int")])
    lyr.updateFields()
    feats = []
    for i in check.issues:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(i.problem.x, i.problem.y)))
        d = i.problem.distance
        f.setAttributes([i.problem.code, i.label, i.problem.message, None if d != d else float(d),
                         None if i.problem.size != i.problem.size else float(i.problem.size), 1 if i.problem.fixable else 0])
        feats.append(f)
    pr.addFeatures(feats)
    lyr.updateExtents()
    cats = []
    for value, color, label in ((1, "#f28e2b", "can be fixed safely"), (0, "#d62728", "to look at")):
        sym = QgsMarkerSymbol.createSimple({"name": "circle", "color": color, "size": "3.2", "outline_color": "#ffffff",
                                            "outline_width": "0.4"})
        cats.append(QgsRendererCategory(value, sym, label))
    lyr.setRenderer(QgsCategorizedSymbolRenderer("fixable", cats))
    bg, fg = M.tip_colors()
    lyr.setMapTipTemplate('<div style="background-color:%s;color:%s;padding:3px 6px"><b>[%% "object" %%]</b><br>'
                          '[%% "message" %%]</div>' % (bg, fg))
    lyr.setCustomProperty("telemac_toolkit/drawing_check", "1")
    return lyr
