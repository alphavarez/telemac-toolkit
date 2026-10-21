"""A small line-plot widget drawn with QPainter (no matplotlib needed).

Supports several series, NaN gaps, a vertical time cursor, hover read-out,
a filled band between two series (bed and water surface in profiles), and
export to PNG / CSV / clipboard from the right-click menu.
"""

from __future__ import annotations

import math
from typing import List, Optional, Sequence

import numpy as np
from qgis.PyQt.QtCore import QPointF, QRectF, Qt, pyqtSignal
from qgis.PyQt.QtGui import QBrush, QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPalette, QPen, QPolygonF
from qgis.PyQt.QtWidgets import QApplication, QFileDialog, QMenu, QSizePolicy, QWidget

from ..compat import menu_exec

#: Colour-blind-safe categorical palette (Okabe-Ito).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#7a6a00", "#555555"]


def nice_ticks(lo: float, hi: float, n: int = 6) -> List[float]:
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return []
    if hi == lo:
        hi = lo + 1.0
    span = hi - lo
    raw = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if span / step <= n:
            break
    start = math.ceil(lo / step - 1e-9) * step
    ticks = []
    v = start
    while v <= hi + step * 1e-6:
        ticks.append(0.0 if abs(v) < step * 1e-9 else v)
        v += step
    return ticks


def fmt_tick(v: float, step: float) -> str:
    if abs(v) >= 1e5 or (abs(v) < 1e-3 and v != 0):
        return "%.3g" % v
    if not step:
        return "%g" % v
    e = math.floor(math.log10(abs(step)))
    mant = abs(step) / 10 ** e
    dec = max(0, -e + (1 if math.isclose(mant, 2.5) else 0))
    return "%.*f" % (dec, v)


class Series:
    def __init__(self, label: str, x, y, color: str, width: float = 1.6, dashed: bool = False):
        self.label = label
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.color = color
        self.width = width
        self.dashed = dashed


class PlotWidget(QWidget):
    """Line plot with a vertical cursor and hover values."""

    clickedX = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.series: List[Series] = []
        self.band = None  # (x, y_low, y_high, color)
        self.title = ""
        self.xlabel = ""
        self.ylabel = ""
        self.cursor_x: Optional[float] = None
        self.x_formatter = None
        self.empty_text = "No data"
        self._hover: Optional[QPointF] = None
        self._plot_rect = QRectF()
        self._xr = (0.0, 1.0)
        self._yr = (0.0, 1.0)
        self.setMouseTracking(True)
        self.setMinimumHeight(170)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.DefaultContextMenu)

    # ---- data -------------------------------------------------------------
    def clear(self, empty_text: Optional[str] = None):
        self.series = []
        self.band = None
        if empty_text is not None:
            self.empty_text = empty_text
        self.update()

    def set_series(self, series: Sequence[Series], title: str = "", xlabel: str = "", ylabel: str = "", band=None):
        self.series = list(series)
        self.band = band
        self.title, self.xlabel, self.ylabel = title, xlabel, ylabel
        self.update()

    def set_cursor(self, x: Optional[float]):
        self.cursor_x = x
        self.update()

    # ---- geometry ---------------------------------------------------------
    def _ranges(self):
        xs = [s.x[np.isfinite(s.x)] for s in self.series if s.x.size]
        ys = [s.y[np.isfinite(s.y)] for s in self.series if s.y.size]
        if self.band is not None:
            ys += [np.asarray(self.band[1])[np.isfinite(self.band[1])], np.asarray(self.band[2])[np.isfinite(self.band[2])]]
        xs = np.concatenate(xs) if xs else np.array([])
        ys = np.concatenate(ys) if ys else np.array([])
        if xs.size == 0 or ys.size == 0:
            return None
        x0, x1 = float(xs.min()), float(xs.max())
        y0, y1 = float(ys.min()), float(ys.max())
        if x1 == x0:
            x1 = x0 + 1.0
        if y1 == y0:
            pad = abs(y0) * 0.05 or 1.0
            y0, y1 = y0 - pad, y1 + pad
        else:
            pad = (y1 - y0) * 0.06
            y0, y1 = y0 - pad, y1 + pad
        return (x0, x1), (y0, y1)

    def _to_px(self, x, y):
        r = self._plot_rect
        (x0, x1), (y0, y1) = self._xr, self._yr
        px = r.left() + (np.asarray(x) - x0) / (x1 - x0) * r.width()
        py = r.bottom() - (np.asarray(y) - y0) / (y1 - y0) * r.height()
        return px, py

    def _x_from_px(self, px: float) -> float:
        r = self._plot_rect
        x0, x1 = self._xr
        return x0 + (px - r.left()) / max(r.width(), 1) * (x1 - x0)

    # ---- painting -----------------------------------------------------------
    def paintEvent(self, _event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.render_to(p, self.rect())
        p.end()

    def render_to(self, p: QPainter, rect):
        pal = self.palette()
        bg = pal.color(QPalette.ColorRole.Base)
        fg = pal.color(QPalette.ColorRole.Text)
        quiet = QColor(fg)
        quiet.setAlpha(150)
        grid = QColor(fg)
        grid.setAlpha(28)
        p.fillRect(rect, bg)
        font = QFont(self.font())
        fm = QFontMetrics(font)
        rng = self._ranges()
        if rng is None:
            p.setPen(quiet)
            p.drawText(QRectF(rect), int(Qt.AlignmentFlag.AlignCenter), self.empty_text)
            return
        self._xr, self._yr = rng
        yt = nice_ticks(*self._yr, n=6)
        xt = nice_ticks(*self._xr, n=7)
        ystep = (yt[1] - yt[0]) if len(yt) > 1 else 1.0
        xstep = (xt[1] - xt[0]) if len(xt) > 1 else 1.0
        ylabels = [fmt_tick(v, ystep) for v in yt]
        left = 14 + max([fm.horizontalAdvance(s) for s in ylabels] + [20]) + (fm.height() + 4 if self.ylabel else 0)
        top = 10 + (fm.height() + 6 if self.title else 0)
        bottom = 10 + fm.height() * 2 + (fm.height() if self.xlabel else 0)
        right = 16
        r = QRectF(rect.left() + left, rect.top() + top, rect.width() - left - right, rect.height() - top - bottom)
        if r.width() < 40 or r.height() < 30:
            return
        self._plot_rect = r

        # title, axis labels
        p.setPen(fg)
        if self.title:
            f2 = QFont(font)
            f2.setBold(True)
            p.setFont(f2)
            p.drawText(QRectF(rect.left() + left, rect.top() + 6, r.width(), fm.height() + 2),
                       int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter), self.title)
            p.setFont(font)
        p.setPen(quiet)
        if self.xlabel:
            p.drawText(QRectF(r.left(), rect.bottom() - fm.height() - 4, r.width(), fm.height()),
                       int(Qt.AlignmentFlag.AlignCenter), self.xlabel)
        if self.ylabel:
            p.save()
            p.translate(rect.left() + 4 + fm.height() / 2, r.center().y())
            p.rotate(-90)
            p.drawText(QRectF(-r.height() / 2, -fm.height() / 2, r.height(), fm.height()),
                       int(Qt.AlignmentFlag.AlignCenter), self.ylabel)
            p.restore()

        # grid and ticks
        for v, s in zip(yt, ylabels):
            _, py = self._to_px(self._xr[0], v)
            if r.top() - 1 <= py <= r.bottom() + 1:
                p.setPen(QPen(grid, 1))
                p.drawLine(QPointF(r.left(), py), QPointF(r.right(), py))
                p.setPen(quiet)
                p.drawText(QRectF(r.left() - 8 - fm.horizontalAdvance(s), py - fm.height() / 2, fm.horizontalAdvance(s), fm.height()),
                           int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter), s)
        for v in xt:
            px, _ = self._to_px(v, self._yr[0])
            if r.left() - 1 <= px <= r.right() + 1:
                p.setPen(QPen(grid, 1))
                p.drawLine(QPointF(px, r.top()), QPointF(px, r.bottom()))
                s = self.x_formatter(v, xstep) if self.x_formatter else fmt_tick(v, xstep)
                w = fm.horizontalAdvance(s)
                p.setPen(quiet)
                p.drawText(QRectF(px - w / 2, r.bottom() + 4, w, fm.height()), int(Qt.AlignmentFlag.AlignCenter), s)
        p.setPen(QPen(quiet, 1))
        p.drawLine(QPointF(r.left(), r.bottom()), QPointF(r.right(), r.bottom()))
        p.drawLine(QPointF(r.left(), r.top()), QPointF(r.left(), r.bottom()))

        p.save()
        p.setClipRect(r.adjusted(-1, -1, 1, 1))
        # band (e.g. water between bed and surface)
        if self.band is not None:
            bx, blo, bhi, bcol = self.band
            bx, blo, bhi = (np.asarray(a, float) for a in (bx, blo, bhi))
            ok = np.isfinite(bx) & np.isfinite(blo) & np.isfinite(bhi) & (bhi > blo)
            col = QColor(bcol)
            col.setAlpha(70)
            for seg in _runs(ok):
                xs, ylo, yhi = bx[seg], blo[seg], bhi[seg]
                pxs, plo = self._to_px(xs, ylo)
                _, phi = self._to_px(xs, yhi)
                poly = QPolygonF([QPointF(a, b) for a, b in zip(pxs, phi)] + [QPointF(a, b) for a, b in zip(pxs[::-1], plo[::-1])])
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QBrush(col))
                p.drawPolygon(poly)
        # series
        for s in self.series:
            pen = QPen(QColor(s.color), s.width)
            if s.dashed:
                pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            ok = np.isfinite(s.x) & np.isfinite(s.y)
            px, py = self._to_px(s.x, s.y)
            for seg in _runs(ok):
                path = QPainterPath(QPointF(px[seg][0], py[seg][0]))
                for a, b in zip(px[seg][1:], py[seg][1:]):
                    path.lineTo(QPointF(a, b))
                if len(px[seg]) == 1:
                    p.drawEllipse(QPointF(px[seg][0], py[seg][0]), 2, 2)
                p.drawPath(path)
        # time cursor
        if self.cursor_x is not None and self._xr[0] <= self.cursor_x <= self._xr[1]:
            cx, _ = self._to_px(self.cursor_x, 0)
            pen = QPen(QColor("#D55E00"), 1.2)
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.drawLine(QPointF(cx, r.top()), QPointF(cx, r.bottom()))
        p.restore()

        # legend
        labels = [s for s in self.series if s.label]
        if labels:
            lw = max(fm.horizontalAdvance(s.label) for s in labels) + 30
            lh = fm.height() * len(labels) + 8
            lx = r.right() - lw - 6
            ly = r.top() + 6
            back = QColor(bg)
            back.setAlpha(215)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(back)
            p.drawRoundedRect(QRectF(lx, ly, lw, lh), 4, 4)
            for k, s in enumerate(labels):
                yy = ly + 4 + fm.height() * k + fm.height() / 2
                pen = QPen(QColor(s.color), 2)
                if s.dashed:
                    pen.setStyle(Qt.PenStyle.DashLine)
                p.setPen(pen)
                p.drawLine(QPointF(lx + 6, yy), QPointF(lx + 22, yy))
                p.setPen(fg)
                p.drawText(QRectF(lx + 26, yy - fm.height() / 2, lw - 28, fm.height()),
                           int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter), s.label)

        # hover read-out
        if self._hover is not None and r.contains(self._hover):
            xh = self._x_from_px(self._hover.x())
            lines = [(self.x_formatter(xh, 0) if self.x_formatter else "x = %.4g" % xh, fg)]
            for s in self.series:
                if s.x.size == 0:
                    continue
                k = int(np.nanargmin(np.abs(s.x - xh)))
                if np.isfinite(s.y[k]):
                    lines.append(("%s: %.4g" % (s.label or "value", s.y[k]), QColor(s.color)))
            hx = self._hover.x() + 12
            w = max(fm.horizontalAdvance(t) for t, _ in lines) + 12
            h = fm.height() * len(lines) + 8
            if hx + w > r.right():
                hx = self._hover.x() - w - 12
            hy = max(r.top(), min(self._hover.y() - h / 2, r.bottom() - h))
            back = QColor(bg)
            back.setAlpha(235)
            p.setPen(QPen(grid, 1))
            p.setBrush(back)
            p.drawRoundedRect(QRectF(hx, hy, w, h), 4, 4)
            for k, (t, c) in enumerate(lines):
                p.setPen(c)
                p.drawText(QRectF(hx + 6, hy + 4 + k * fm.height(), w - 8, fm.height()),
                           int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter), t)
            p.setPen(QPen(grid, 1))
            p.drawLine(QPointF(self._hover.x(), r.top()), QPointF(self._hover.x(), r.bottom()))

    # ---- interaction --------------------------------------------------------
    def mouseMoveEvent(self, e):
        pos = e.position() if hasattr(e, "position") else QPointF(e.pos())
        self._hover = QPointF(pos)
        self.update()

    def leaveEvent(self, _e):
        self._hover = None
        self.update()

    def mouseReleaseEvent(self, e):
        if e.button() == Qt.MouseButton.LeftButton and self._plot_rect.contains(self._hover or QPointF(-1, -1)):
            self.clickedX.emit(self._x_from_px(self._hover.x()))

    def contextMenuEvent(self, e):
        m = QMenu(self)
        a_copy = m.addAction(self.tr("Copy data (CSV)"))
        a_csv = m.addAction(self.tr("Save data as CSV…"))
        a_png = m.addAction(self.tr("Save image as PNG…"))
        act = menu_exec(m, e.globalPos())
        if act == a_copy:
            QApplication.clipboard().setText(self.to_csv())
        elif act == a_csv:
            fn, _ = QFileDialog.getSaveFileName(self, self.tr("Save CSV"), "", "CSV (*.csv)")
            if fn:
                with open(fn, "w", encoding="utf-8") as f:
                    f.write(self.to_csv())
        elif act == a_png:
            fn, _ = QFileDialog.getSaveFileName(self, self.tr("Save PNG"), "", "PNG (*.png)")
            if fn:
                self.grab().save(fn, "PNG")

    def to_csv(self) -> str:
        """All series side by side; series with a different x get their own columns."""
        if not self.series:
            return ""
        same_x = all(s.x.shape == self.series[0].x.shape and np.allclose(s.x, self.series[0].x, equal_nan=True)
                     for s in self.series)
        out = []
        if same_x:
            head = [self.xlabel or "x"] + [s.label or "y%d" % k for k, s in enumerate(self.series)]
            out.append(",".join('"%s"' % h for h in head))
            for k in range(self.series[0].x.size):
                row = [self.series[0].x[k]] + [s.y[k] for s in self.series]
                out.append(",".join("" if not np.isfinite(v) else repr(float(v)) for v in row))
        else:
            for s in self.series:
                out.append('"%s","%s"' % (self.xlabel or "x", s.label or "y"))
                for a, b in zip(s.x, s.y):
                    out.append("%s,%s" % (repr(float(a)), "" if not np.isfinite(b) else repr(float(b))))
                out.append("")
        return "\n".join(out) + "\n"


def _runs(mask: np.ndarray):
    """Index arrays of consecutive True runs."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) > 1)
    return np.split(idx, breaks + 1)
