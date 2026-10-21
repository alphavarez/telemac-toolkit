"""A TELEMAC result loaded in QGIS: the SELAFIN reader + the QGIS mesh layer.

QGIS (MDAL) draws the mesh and handles the temporal controller. The plugin
keeps its own lazy reader of the same file for the live cursor, probes and
profiles, so values come straight from the file at full precision.
"""

from __future__ import annotations

import math
import os
from typing import List, Optional, Tuple

import numpy as np
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsMeshLayer,
    QgsProject,
)
from qgis.PyQt.QtCore import QDateTime, QObject, Qt, pyqtSignal

from ..telemac_io import SelafinFile, TriangleLocator, extract_2d, sample_polyline
from ..telemac_io.selafin import precision_warning
from ..telemac_io.variables import canonical, display_name, lookup, vector_pairs
from . import styles
from .crs import Extent, crs_label, sidecar_crs, suggest_crs  # noqa: F401 (re-exported)

def _qdatetime(d):
    from qgis.PyQt.QtCore import QDate, QTime

    return QDateTime(QDate(d.year, d.month, d.day), QTime(d.hour, d.minute, d.second), Qt.TimeSpec.UTC)


def format_duration(seconds: float) -> str:
    """'2 h 46 min', '35 min 10 s', '12.5 s', '3 d 4 h'."""
    s = float(seconds)
    sign = "-" if s < 0 else ""
    s = abs(s)
    if s < 60:
        return "%s%g s" % (sign, round(s, 2))
    d, rem = divmod(int(round(s)), 86400)
    h, rem = divmod(rem, 3600)
    m, sec = divmod(rem, 60)
    if d:
        return "%s%d d %d h" % (sign, d, h) if not m else "%s%d d %d h %d min" % (sign, d, h, m)
    if h:
        return "%s%d h %02d min" % (sign, h, m) if m else "%s%d h" % (sign, h)
    return "%s%d min %02d s" % (sign, m, sec) if sec else "%s%d min" % (sign, m)


class ProfileSampler:
    """Points along a line with their triangles and weights, reused every frame."""

    def __init__(self, line_xy, locator: TriangleLocator, step: float):
        self.line_xy = np.asarray(line_xy, dtype=float)
        self.px, self.py, self.s = sample_polyline(self.line_xy, step=step)
        self.tri, self.nodes, self.w = locator.locate_many(self.px, self.py)
        self.inside = self.tri >= 0

        # only these nodes are read from the file at each frame
        self.unique_nodes, inv = np.unique(self.nodes.ravel(), return_inverse=True)
        self._inv = inv.reshape(self.nodes.shape)

    def apply(self, node_values: np.ndarray, local: bool = False) -> np.ndarray:
        """node_values for all mesh nodes, or (local=True) for ``unique_nodes`` only."""
        nv = np.asarray(node_values)
        idx = self._inv if local else self.nodes
        v = np.einsum("ij,ij->i", nv[idx], self.w)
        v[~self.inside] = np.nan
        return v


class ResultSet(QObject):
    """One loaded TELEMAC result."""

    reloaded = pyqtSignal()

    def __init__(self, source_path: str, display_path: str, layer: QgsMeshLayer, derived_from_3d: str = ""):
        super().__init__()
        self.source_path = os.path.abspath(source_path)
        self.path = os.path.abspath(display_path)
        self.layer = layer
        self.layer_id = layer.id()
        self.derived_from_3d = derived_from_3d  # "", "average", "plane 0", ...
        self.slf = SelafinFile(self.path)
        self._locator: Optional[TriangleLocator] = None
        self._node_area: Optional[np.ndarray] = None
        self.dry_threshold = 0.01

    # ---- metadata -----------------------------------------------------------
    @property
    def header(self):
        return self.slf.header

    @property
    def names(self) -> List[str]:
        return self.slf.var_names

    @property
    def units(self) -> List[str]:
        return self.slf.var_units

    def label(self, i: int) -> str:
        return display_name(self.names[i], self.units[i])

    @property
    def pairs(self):
        return vector_pairs(self.names)

    @property
    def times(self) -> np.ndarray:
        return self.slf.times

    @property
    def nframes(self) -> int:
        return self.slf.nframes

    @property
    def has_calendar_date(self) -> bool:
        return self.header.has_calendar_date

    def reference_datetime(self) -> QDateTime:
        d = self.header.start_datetime
        if d is None:
            from qgis.PyQt.QtCore import QDate, QTime

            return QDateTime(QDate(1900, 1, 1), QTime(0, 0, 0), Qt.TimeSpec.UTC)
        return _qdatetime(d)

    def time_label(self, i: int) -> str:
        t = float(self.times[i])
        txt = format_duration(t)
        if abs(t) >= 60:
            txt += "  (%s s)" % format(int(round(t)), ",").replace(",", " ")
        if self.has_calendar_date:
            dts = self.slf.datetimes
            txt = dts[i].strftime("%Y-%m-%d %H:%M:%S") + "   " + txt
        return txt

    def var_index(self, cid_or_name: str) -> Optional[int]:
        for i, n in enumerate(self.names):
            if canonical(n) == cid_or_name or n.strip().upper() == str(cid_or_name).strip().upper():
                return i
        return None

    def default_variable(self) -> int:
        for cid in ("water_depth", "free_surface", "scalar_velocity", "bottom"):
            i = self.var_index(cid)
            if i is not None:
                return i
        return 0

    def warnings(self) -> List[str]:
        out = []
        w = precision_warning(self.header)
        if w:
            out.append(w)
        if self.slf.partial_bytes:
            out.append("The last frame is incomplete; the simulation may still be running. Use Refresh to read new frames.")
        if self.header.is_partitioned:
            out.append("This is a sub-domain file of a parallel run; merge the partitions (gretel) to see the whole model.")
        return out

    # ---- geometry -----------------------------------------------------------
    @property
    def locator(self) -> TriangleLocator:
        if self._locator is None:
            h = self.header
            self._locator = TriangleLocator(h.x, h.y, h.ikle)
        return self._locator

    @property
    def node_area(self) -> np.ndarray:
        """Area associated with each node (one third of adjacent triangles)."""
        if self._node_area is None:
            h = self.header
            x, y, ik = h.x, h.y, h.ikle
            a = 0.5 * np.abs((x[ik[:, 1]] - x[ik[:, 0]]) * (y[ik[:, 2]] - y[ik[:, 0]])
                             - (x[ik[:, 2]] - x[ik[:, 0]]) * (y[ik[:, 1]] - y[ik[:, 0]]))
            na = np.zeros(h.npoin)
            for k in range(3):
                np.add.at(na, ik[:, k], a / 3.0)
            self._node_area = na
        return self._node_area

    def mean_edge(self) -> float:
        h = self.header
        ik = h.ikle
        e = np.hypot(h.x[ik[:, 1]] - h.x[ik[:, 0]], h.y[ik[:, 1]] - h.y[ik[:, 0]])
        return float(np.mean(e)) if e.size else 1.0

    # ---- values -------------------------------------------------------------
    def locate(self, px: float, py: float):
        return self.locator.locate(px, py)

    def values_at(self, px: float, py: float, frame: int) -> Optional[List[Tuple[str, float, str]]]:
        """[(label, value, unit)] for every variable at a point, plus vector magnitudes."""
        r = self.locate(px, py)
        if r is None:
            return None
        _, nodes, w = r
        fr = np.asarray(self.slf.values[self.slf._frame_index(frame)][:, nodes], dtype=np.float64) @ w
        rows = []
        for i, n in enumerate(self.names):
            info = lookup(n)
            rows.append((info.label if info else n.strip().capitalize(), float(fr[i]), self.units[i].lower()))
        for i, j, label in self.pairs:
            mag = math.hypot(fr[i], fr[j])
            ang = (math.degrees(math.atan2(fr[i], fr[j])) + 360.0) % 360.0  # bearing, 0 = north
            if not self.has_file_magnitude((i, j)):
                rows.append(("%s magnitude" % label, mag, self.units[i].lower()))
            rows.append(("%s direction (bearing)" % label, ang, "deg"))
        return rows

    def has_file_magnitude(self, pair) -> bool:
        """True when TELEMAC already wrote the magnitude of this vector (SCALAR VELOCITY)."""
        i, _j = pair
        return canonical(self.names[i]) == "velocity_u" and self.var_index("scalar_velocity") is not None

    def node_series(self, var: int, nodes, w) -> np.ndarray:
        return self.slf.series_weighted(var, nodes, w)

    def series_at(self, px: float, py: float, var) -> Optional[np.ndarray]:
        r = self.locate(px, py)
        if r is None:
            return None
        if isinstance(var, tuple):  # vector magnitude (i, j)
            a = self.slf.series_weighted(var[0], r[1], r[2])
            b = self.slf.series_weighted(var[1], r[1], r[2])
            return np.hypot(a, b)
        return self.slf.series_weighted(var, r[1], r[2])

    def values_at_nodes(self, var, frame: int, nodes: np.ndarray) -> np.ndarray:
        """Values of a variable (or vector magnitude) at selected nodes only."""
        f = self.slf._frame_index(frame)
        vals = self.slf.values
        if isinstance(var, tuple):
            a = np.asarray(vals[f, var[0]][nodes], dtype=np.float64)
            b = np.asarray(vals[f, var[1]][nodes], dtype=np.float64)
            return np.hypot(a, b)
        return np.asarray(vals[f, self.header.var_index(var)][nodes], dtype=np.float64)

    def profile(self, sampler: "ProfileSampler", var, frame: int) -> np.ndarray:
        return sampler.apply(self.values_at_nodes(var, frame, sampler.unique_nodes), local=True)

    def frame_values(self, var, frame: int) -> np.ndarray:
        if isinstance(var, tuple):
            return np.hypot(self.slf.read(var[0], frame), self.slf.read(var[1], frame))
        return self.slf.read(var, frame)

    def sampler(self, line_xy, step: Optional[float] = None) -> ProfileSampler:
        if step is None:
            step = self.mean_edge() / 4.0
            total = float(np.sum(np.hypot(*np.diff(np.asarray(line_xy, float), axis=0).T)))
            step = min(max(step, total / 4000.0), total / 50.0) if total > 0 else 1.0
        return ProfileSampler(line_xy, self.locator, step)

    def domain_mean_series(self, var) -> np.ndarray:
        """Area-weighted mean over the mesh for every frame."""
        na = self.node_area
        tot = na.sum()
        out = np.empty(self.nframes)
        for f in range(self.nframes):
            out[f] = float(self.frame_values(var, f) @ na) / tot
        return out

    def peak_frame(self, var, point=None) -> int:
        if point is not None:
            s = self.series_at(point[0], point[1], var)
            if s is not None and np.isfinite(s).any():
                return int(np.nanargmax(s))
        return int(np.argmax(self.domain_mean_series(var)))

    # ---- life cycle -----------------------------------------------------------
    def refresh(self) -> int:
        """Re-read the file (new frames from a running simulation)."""
        n = self.slf.refresh()
        try:
            self.layer.dataProvider().reloadData()
            self.layer.reload()
            self.layer.triggerRepaint()
        except Exception:
            pass
        self.reloaded.emit()
        return n

    def close(self):
        self.slf.close()


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def coordinates_look_geographic(header) -> bool:
    if not header.npoin:
        return False
    return Extent.of(header.x, header.y).fits_lonlat


def file_extent(path: str) -> Extent:
    with SelafinFile(path) as probe:
        h = probe.header
        return Extent.of(h.x2 if h.is_3d else h.x, h.y2 if h.is_3d else h.y)


def derived_2d_path(source: str, mode: str) -> str:
    base, _ = os.path.splitext(source)
    tag = mode.replace(" ", "")
    return "%s_2d_%s.slf" % (base, tag)


def load_result(
    path: str,
    plane: str = "average",
    crs: Optional[QgsCoordinateReferenceSystem] = None,
    name: Optional[str] = None,
    project: Optional[QgsProject] = None,
    add_to_project: bool = True,
    out_dir: Optional[str] = None,
    progress=None,
) -> ResultSet:
    """Open a SELAFIN result as a styled QGIS mesh layer.

    3D files are converted to a 2D file first (QGIS cannot read 3D SELAFIN):
    ``plane`` is "average", "bottom", "surface" or a plane index. The 2D file is
    written next to the source (or in ``out_dir``) and reused if up to date.
    """
    path = os.path.abspath(path)
    project = project or QgsProject.instance()
    with SelafinFile(path) as probe:
        is3d = probe.header.is_3d
        ndp = probe.header.ndp
    if ndp != 3 and not is3d:
        raise ValueError(
            "This file has %d node(s) per element, not triangles (a point history or a "
            "quadratic file). It cannot be shown as a mesh." % ndp
        )
    display = path
    mode = ""
    if is3d:
        mode = plane if isinstance(plane, str) else "plane %d" % int(plane)
        target = derived_2d_path(path, mode)
        if out_dir:
            target = os.path.join(out_dir, os.path.basename(target))
        if not os.path.isdir(os.path.dirname(target)) or not os.access(os.path.dirname(target), os.W_OK):
            import tempfile

            target = os.path.join(tempfile.gettempdir(), os.path.basename(target))
        if not (os.path.exists(target) and os.path.getmtime(target) >= os.path.getmtime(path)):
            extract_2d(path, target, plane if isinstance(plane, str) else int(plane), progress=progress)
        display = target

    lname = name or os.path.splitext(os.path.basename(path))[0] + (" (%s)" % mode if mode else "")
    layer = QgsMeshLayer(display, lname, "mdal")
    if not layer.isValid():
        raise RuntimeError("QGIS could not open %s as a mesh layer:\n%s" % (display, layer.error().summary()))
    if crs is None or not crs.isValid():
        # SELAFIN stores no CRS: a .prj next to the file, else the project CRS if the numbers fit it
        choice = suggest_crs(path, file_extent(path), project)
        crs = choice.crs  # the best guess when nobody can be asked (scripts)
    if crs is not None and crs.isValid():
        layer.setCrs(crs)

    rs = ResultSet(path, display, layer, derived_from_3d=mode)
    # temporal reference: file date, or TELEMAC's 1900-01-01 default
    tp = layer.temporalProperties()
    ref = tp.referenceTime()
    if not ref.isValid():
        try:
            tp.setReferenceTime(rs.reference_datetime(), layer.dataProvider().temporalCapabilities())
        except Exception:
            pass
    styles.apply_default_style(rs)
    if add_to_project:
        project.addMapLayer(layer)
    return rs
