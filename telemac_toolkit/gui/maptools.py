"""Map tools: live cursor, probe points, profile line."""

from __future__ import annotations

import contextlib
from qgis.core import (
    Qgis,
    QgsCoordinateTransform,
    QgsGeometry,
    QgsPointXY,
    QgsProject,
)
from qgis.gui import QgsMapTool, QgsRubberBand, QgsVertexMarker
from qgis.PyQt.QtCore import Qt, pyqtSignal
from qgis.PyQt.QtGui import QColor, QCursor


class _ResultTool(QgsMapTool):
    """Base: converts canvas points to the result layer's CRS."""

    def __init__(self, canvas, get_result):
        super().__init__(canvas)
        self.get_result = get_result  # callable -> ResultSet or None
        self.setCursor(QCursor(Qt.CursorShape.CrossCursor))

    def to_layer(self, map_point):
        rs = self.get_result()
        if rs is None:
            return None, None
        canvas_crs = self.canvas().mapSettings().destinationCrs()
        lyr_crs = rs.layer.crs()
        pt = QgsPointXY(map_point)
        if canvas_crs.isValid() and lyr_crs.isValid() and canvas_crs != lyr_crs:
            tr = QgsCoordinateTransform(canvas_crs, lyr_crs, QgsProject.instance())
            try:
                pt = tr.transform(pt)
            except Exception:
                return rs, None
        return rs, pt

    def from_layer(self, rs, pt):
        canvas_crs = self.canvas().mapSettings().destinationCrs()
        lyr_crs = rs.layer.crs()
        if canvas_crs.isValid() and lyr_crs.isValid() and canvas_crs != lyr_crs:
            tr = QgsCoordinateTransform(lyr_crs, canvas_crs, QgsProject.instance())
            return tr.transform(QgsPointXY(pt))
        return QgsPointXY(pt)


class LiveCursorTool(_ResultTool):
    """Reports every variable under the mouse at the current time."""

    moved = pyqtSignal(float, float)  # layer coordinates
    left_mesh = pyqtSignal()

    def canvasMoveEvent(self, e):
        rs, pt = self.to_layer(e.mapPoint())
        if pt is None:
            self.left_mesh.emit()
            return
        self.moved.emit(pt.x(), pt.y())


class ProbeTool(_ResultTool):
    """Click to add probe points (time series); right-click clears them."""

    added = pyqtSignal(float, float)  # layer coordinates
    cleared = pyqtSignal()

    def __init__(self, canvas, get_result, colors):
        super().__init__(canvas, get_result)
        self.colors = colors
        self.markers = []

    def canvasReleaseEvent(self, e):
        if e.button() == Qt.MouseButton.RightButton:
            self.clear()
            self.cleared.emit()
            return
        rs, pt = self.to_layer(e.mapPoint())
        if pt is None or rs.locate(pt.x(), pt.y()) is None:
            return
        m = QgsVertexMarker(self.canvas())
        m.setCenter(e.mapPoint())
        m.setIconType(QgsVertexMarker.IconType.ICON_CIRCLE)
        m.setIconSize(11)
        m.setPenWidth(3)
        m.setColor(QColor(self.colors[len(self.markers) % len(self.colors)]))
        self.markers.append(m)
        self.added.emit(pt.x(), pt.y())

    def clear(self):
        for m in self.markers:
            self.canvas().scene().removeItem(m)
        self.markers = []


class ProfileTool(_ResultTool):
    """Left-click vertices, right-click (or double-click) to finish the line."""

    finished = pyqtSignal(list)  # [(x, y), ...] in layer coordinates

    def __init__(self, canvas, get_result):
        super().__init__(canvas, get_result)
        self.band = QgsRubberBand(canvas, Qgis.GeometryType.Line)
        self.band.setColor(QColor("#D55E00"))
        self.band.setWidth(2)
        self.points = []  # map coordinates
        self.done = False

    def reset(self):
        self.band.reset(Qgis.GeometryType.Line)
        self.points = []
        self.done = False

    def canvasMoveEvent(self, e):
        if self.points and not self.done:
            self.band.setToGeometry(QgsGeometry.fromPolylineXY(self.points + [QgsPointXY(e.mapPoint())]), None)

    def canvasReleaseEvent(self, e):
        if self.done:
            self.reset()
        if e.button() == Qt.MouseButton.LeftButton:
            self.points.append(QgsPointXY(e.mapPoint()))
            if len(self.points) > 1:
                self.band.setToGeometry(QgsGeometry.fromPolylineXY(self.points), None)
        elif e.button() == Qt.MouseButton.RightButton:
            self._finish()

    def canvasDoubleClickEvent(self, e):
        self._finish()

    def _finish(self):
        if len(self.points) < 2:
            self.reset()
            return
        self.done = True
        self.band.setToGeometry(QgsGeometry.fromPolylineXY(self.points), None)
        rs = self.get_result()
        if rs is None:
            return
        out = []
        for p in self.points:
            _, q = self.to_layer(p)
            if q is not None:
                out.append((q.x(), q.y()))
        if len(out) >= 2:
            self.finished.emit(out)

    def keyPressEvent(self, e):
        if e.key() == Qt.Key.Key_Escape:
            self.reset()
        elif e.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self._finish()

    def deactivate(self):
        super().deactivate()

    def remove(self):
        self.reset()
        with contextlib.suppress(RuntimeError):  # the canvas may already be gone
            self.canvas().scene().removeItem(self.band)
