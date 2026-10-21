"""Coordinate systems for TELEMAC files, which store none.

SELAFIN keeps plain x, y numbers. The plugin decides what they mean, in this
order: a .prj file next to the mesh or result; the project CRS when the
numbers make sense in it; otherwise the user chooses (map CRS or the local
laboratory grid) and may save the choice as a .prj, so the file opens directly
next time. "Make sense" is checked against the CRS's area of use, which is
what tells a 20 m flume (x from 0 to 20) apart from UTM coordinates.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
from qgis.core import (
    Qgis,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsPointXY,
    QgsProject,
    QgsRectangle,
)

#: a local metric grid for laboratory flumes and schematic channels: a transverse Mercator plane with
#: scale 1 and no false origin, so (0, 0) is the origin and lengths are exact metres near it. It is a
#: normal projected CRS for QGIS (measuring, labels, every format), not a place on Earth.
LOCAL_CRS_NAME = "Local metric grid (laboratory)"
LOCAL_CRS_WKT = (
    'PROJCRS["Local metric grid (laboratory)",BASEGEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
    'ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,'
    'ANGLEUNIT["degree",0.0174532925199433]]],CONVERSION["Local plane",METHOD["Transverse Mercator",ID["EPSG",9807]],'
    'PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433]],'
    'PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433]],'
    'PARAMETER["Scale factor at natural origin",1,SCALEUNIT["unity",1]],'
    'PARAMETER["False easting",0,LENGTHUNIT["metre",1]],PARAMETER["False northing",0,LENGTHUNIT["metre",1]]],'
    'CS[Cartesian,2],AXIS["easting (X)",east,ORDER[1],LENGTHUNIT["metre",1]],'
    'AXIS["northing (Y)",north,ORDER[2],LENGTHUNIT["metre",1]]]')

#: models whose coordinates all lie within this distance of (0, 0) can be local (laboratory, schematic)
LOCAL_LIMIT = 100_000.0


def local_crs() -> QgsCoordinateReferenceSystem:
    return QgsCoordinateReferenceSystem.fromWkt(LOCAL_CRS_WKT)


def is_local_crs(crs) -> bool:
    if crs is None or not crs.isValid():
        return False
    if crs.authid() and not crs.authid().upper().startswith("USER:"):
        return False
    return "local metric grid" in crs.description().replace("_", " ").lower()


def crs_label(crs) -> str:
    """Short text for the panels: 'local coordinates', 'EPSG:32614', 'no CRS'."""
    if crs is None or not crs.isValid():
        return "no CRS"
    if is_local_crs(crs):
        return "local coordinates (laboratory grid)"
    return crs.authid() or crs.description() or "custom CRS"


@dataclass
class Extent:
    xmin: float
    ymin: float
    xmax: float
    ymax: float

    @classmethod
    def of(cls, x, y) -> "Extent":
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if x.size == 0:
            return cls(0.0, 0.0, 0.0, 0.0)
        return cls(float(np.nanmin(x)), float(np.nanmin(y)), float(np.nanmax(x)), float(np.nanmax(y)))

    def rect(self) -> QgsRectangle:
        return QgsRectangle(self.xmin, self.ymin, self.xmax, self.ymax)

    @property
    def fits_lonlat(self) -> bool:
        return -180.0 <= self.xmin and self.xmax <= 360.0 and -90.0 <= self.ymin and self.ymax <= 90.0

    @property
    def small(self) -> bool:
        """All coordinates near (0, 0): a laboratory or schematic model is possible."""
        return max(abs(self.xmin), abs(self.xmax), abs(self.ymin), abs(self.ymax)) < LOCAL_LIMIT


def plausible(crs, ext: Extent, project: Optional[QgsProject] = None) -> bool:
    """True when coordinates in ``ext`` make sense in ``crs`` (inside its area of use)."""
    if crs is None or not crs.isValid():
        return False
    if is_local_crs(crs):
        return ext.small
    if crs.isGeographic():
        return ext.fits_lonlat
    try:
        bounds = crs.bounds()  # area of use, longitude/latitude
    except Exception:
        bounds = QgsRectangle()
    if bounds.isNull() or bounds.isEmpty() or bounds.width() <= 0:
        return True  # a custom CRS without area of use: nothing to check
    geo = QgsCoordinateReferenceSystem("EPSG:4326")
    tr = QgsCoordinateTransform(crs, geo, project or QgsProject.instance())
    tr.setBallparkTransformsAreAppropriate(True)
    pad_x, pad_y = 0.1 * bounds.width() + 0.1, 0.1 * bounds.height() + 0.1  # models may overhang a little
    inside = 0
    corners = [(ext.xmin, ext.ymin), (ext.xmax, ext.ymax), (0.5 * (ext.xmin + ext.xmax), 0.5 * (ext.ymin + ext.ymax))]
    for cx, cy in corners:
        try:
            p = tr.transform(QgsPointXY(cx, cy))
        except Exception:
            return False
        if not (math.isfinite(p.x()) and math.isfinite(p.y())):
            return False
        if (bounds.xMinimum() - pad_x <= p.x() <= bounds.xMaximum() + pad_x
                and bounds.yMinimum() - pad_y <= p.y() <= bounds.yMaximum() + pad_y):
            inside += 1
    return inside == len(corners)


def sidecar_crs(path: str) -> Optional[QgsCoordinateReferenceSystem]:
    """CRS from a .prj file next to the file (``name.prj`` or ``name.slf.prj``)."""
    for cand in (os.path.splitext(path)[0] + ".prj", path + ".prj"):
        if os.path.exists(cand):
            try:
                with open(cand, "r", encoding="utf-8", errors="replace") as f:
                    txt = f.read().strip()
            except OSError:
                continue
            crs = QgsCoordinateReferenceSystem.fromWkt(txt)
            if not crs.isValid():
                crs = QgsCoordinateReferenceSystem(txt)  # "EPSG:32614" also accepted
            if crs.isValid():
                return crs
    return None


def write_prj(path: str, crs: QgsCoordinateReferenceSystem) -> Optional[str]:
    """``name.prj`` next to ``name.slf`` (or any mesh file)."""
    if crs is None or not crs.isValid():
        return None
    prj = os.path.splitext(path)[0] + ".prj"
    try:
        # the local grid keeps its name only in WKT2; map CRSs get the widely read ESRI WKT1
        wkt = crs.toWkt(Qgis.CrsWktVariant.Preferred if is_local_crs(crs) else Qgis.CrsWktVariant.Wkt1Esri)
    except Exception:
        wkt = crs.toWkt()
    if not wkt:
        wkt = crs.toWkt()
    with open(prj, "w", encoding="utf-8") as f:
        f.write(wkt)
    return prj


def _finite_rect(r) -> bool:
    return (r is not None and not r.isNull() and all(math.isfinite(v) for v in
            (r.xMinimum(), r.yMinimum(), r.xMaximum(), r.yMaximum())) and r.width() < 1e9 and r.height() < 1e9)


def show_layer(layer, canvas, project: Optional[QgsProject] = None, first: bool = False, force_zoom: bool = False) -> str:
    """Make sure a freshly loaded layer is on screen; returns what was done ('' if nothing).

    * first layer of the project, or a project with no CRS: the project takes
      the layer's CRS and the map zooms to it;
    * a local laboratory layer in a map project (or the reverse): the two cannot
      be drawn together, so the project switches to the layer's CRS;
    * otherwise the map zooms to the layer only when none of it is in view
      (QGIS itself only zooms for the first layer of a project, which is why a
      result opened after a mesh project could stay off screen).
    """
    project = project or QgsProject.instance()
    if canvas is None or layer is None:
        return ""
    lcrs = layer.crs()
    pcrs = project.crs()
    done = ""
    if lcrs.isValid() and (first or not pcrs.isValid() or is_local_crs(lcrs) != is_local_crs(pcrs)):
        if pcrs != lcrs:
            project.setCrs(lcrs)
            done = "project CRS set to %s" % crs_label(lcrs)
        canvas.setDestinationCrs(lcrs)
        force_zoom = True
    try:
        ext = canvas.mapSettings().layerExtentToOutputExtent(layer, layer.extent())
    except Exception:
        ext = None
    if not _finite_rect(ext) or ext.isEmpty():
        # the layer cannot be transformed to the map CRS: show it in its own CRS
        if lcrs.isValid() and pcrs != lcrs:
            project.setCrs(lcrs)
            canvas.setDestinationCrs(lcrs)
            done = "project CRS set to %s" % crs_label(lcrs)
        ext = layer.extent()
        force_zoom = True
    if not _finite_rect(ext):
        return done
    if force_zoom or not canvas.extent().intersects(ext):
        r = QgsRectangle(ext)
        if r.width() == 0 or r.height() == 0:
            r.grow(1.0)
        r.scale(1.05)
        canvas.setExtent(r)
        canvas.refresh()
        done = (done + "; " if done else "") + "zoomed to the layer"
    return done


@dataclass
class CrsChoice:
    crs: Optional[QgsCoordinateReferenceSystem]
    reason: str
    ask: bool  # no safe automatic choice: ask the user (``crs`` is the suggestion)


def suggest_crs(path: str, ext: Extent, project: Optional[QgsProject] = None,
                last: Optional[QgsCoordinateReferenceSystem] = None, ignore=()) -> CrsChoice:
    """What the coordinates of a TELEMAC file most likely are.

    ``ignore`` holds layer ids not to count as "the project already has a map"
    (the layer being loaded).
    """
    project = project or QgsProject.instance()
    side = sidecar_crs(path)
    if side is not None:
        return CrsChoice(side, "from %s" % os.path.basename(os.path.splitext(path)[0] + ".prj"), False)
    pcrs = project.crs()
    has_map = any(lid not in ignore for lid in project.mapLayers())
    # the project CRS, when the numbers fit it; but not the lon/lat default of an empty project for
    # numbers near (0, 0), which are far more often a flume than a place in the Gulf of Guinea
    if plausible(pcrs, ext, project) and (has_map or not (pcrs.isGeographic() and ext.small)):
        return CrsChoice(pcrs, "project CRS", False)
    # nothing certain: a suggestion for the question
    if last is not None and last.isValid() and plausible(last, ext, project):
        return CrsChoice(last, "last choice", True)
    if ext.small:
        # near (0, 0): a laboratory or schematic model is the likely case (an empty project's CRS is only a default)
        return CrsChoice(local_crs(), "coordinates near (0, 0)", True)
    if pcrs.isValid() and plausible(pcrs, ext, project):
        return CrsChoice(pcrs, "project CRS", True)
    return CrsChoice(None, "unknown", True)
