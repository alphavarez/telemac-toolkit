"""The TELEMAC Results dock: time controls, cursor values, time series, profiles."""

from __future__ import annotations

import os
from typing import List, Optional, Tuple

import numpy as np
from qgis.gui import QgsDockWidget
from qgis.PyQt.QtCore import QCoreApplication, QSize, Qt, QTimer
from qgis.PyQt.QtGui import QColor, QIcon
from qgis.PyQt.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QSlider,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..qgis_bridge import styles
from ..qgis_bridge.results import format_duration
from ..telemac_io.variables import canonical
from .plot_widget import PALETTE, PlotWidget, Series

ICON_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "icons")


def icon(name: str) -> QIcon:
    return QIcon(os.path.join(ICON_DIR, name))


def time_axis(times: np.ndarray):
    """(factor, unit) so the time axis reads in s, min, h or days."""
    span = float(times[-1] - times[0]) if times.size > 1 else 0.0
    if span > 4 * 86400:
        return 86400.0, "d"
    if span > 3 * 3600:
        return 3600.0, "h"
    if span > 180:
        return 60.0, "min"
    return 1.0, "s"


class ResultsDock(QgsDockWidget):
    WATER_PROFILE = "__water_profile__"

    def __init__(self, plugin, parent=None):
        super().__init__(QCoreApplication.translate("TelemacToolkit", "TELEMAC Results"), parent)
        self.setObjectName("TelemacToolkitResultsDock")
        self.plugin = plugin
        self.rs = None
        self.probes: List[Tuple[float, float]] = []
        self.cursor_pt: Optional[Tuple[float, float]] = None
        self.profile_line = None
        self.sampler = None
        self._play = QTimer(self)
        self._play.timeout.connect(self._advance)
        self._building = False

        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)

        # --- result picker ------------------------------------------------
        row = QHBoxLayout()
        self.result_combo = QComboBox()
        self.result_combo.setToolTip(self.tr("TELEMAC results loaded in this project"))
        self.result_combo.currentIndexChanged.connect(self._on_result_changed)
        row.addWidget(self.result_combo, 1)
        self.btn_open = self._tool("open.svg", self.tr("Open a TELEMAC result (SELAFIN)…"), plugin.open_result_dialog)
        self.btn_refresh = self._tool("refresh.svg", self.tr("Re-read the file (picks up new frames of a running simulation)"), self.refresh)
        self.btn_zoom = self._tool("zoom.svg", self.tr("Zoom to this result (also re-centres the map after opening "
                                                       "another project)"), lambda: plugin.zoom_to_result(self.rs))
        self.btn_crs = self._tool("crs.svg", self.tr("Coordinate system of this result (local laboratory coordinates "
                                                     "or a map CRS); can save it next to the file"),
                                  lambda: plugin.change_result_crs(self.rs))
        self.btn_calc = self._tool("calc.svg", self.tr("Calculator: new variables from formulas, with gradients, "
                                                       "divergence, vorticity and time maxima"),
                                   lambda: plugin.open_calculator(self.rs))
        for b in (self.btn_open, self.btn_refresh, self.btn_zoom, self.btn_crs, self.btn_calc):
            row.addWidget(b)
        lay.addLayout(row)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.info)

        # --- time controls -----------------------------------------------
        trow = QHBoxLayout()
        self.btn_first = self._tool("first.svg", self.tr("First frame"), lambda: self.set_frame(0))
        self.btn_prev = self._tool("prev.svg", self.tr("Previous frame"), lambda: self.set_frame(self.frame - 1))
        self.btn_play = self._tool("play.svg", self.tr("Play / pause"), self.toggle_play)
        self.btn_play.setCheckable(True)
        self.btn_next = self._tool("next.svg", self.tr("Next frame"), lambda: self.set_frame(self.frame + 1))
        self.btn_last = self._tool("last.svg", self.tr("Last frame"), lambda: self.set_frame(self.nframes - 1))
        for b in (self.btn_first, self.btn_prev, self.btn_play, self.btn_next, self.btn_last):
            trow.addWidget(b)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.valueChanged.connect(self._on_slider)
        trow.addWidget(self.slider, 1)
        self.btn_peak = self._tool("peak.svg", self.tr("Jump to the peak of the shown variable: at the first probe if there is one, "
                                                       "otherwise of its area-weighted mean over the mesh"), self.jump_to_peak)
        trow.addWidget(self.btn_peak)
        lay.addLayout(trow)
        srow = QHBoxLayout()
        self.time_label = QLabel("—")
        f = self.time_label.font()
        f.setBold(True)
        self.time_label.setFont(f)
        srow.addWidget(self.time_label, 1)
        srow.addWidget(QLabel(self.tr("Speed")))
        self.speed = QSpinBox()
        self.speed.setRange(1, 30)
        self.speed.setValue(4)
        self.speed.setSuffix(self.tr(" fps"))
        self.speed.valueChanged.connect(lambda v: self._play.setInterval(int(1000 / v)))
        srow.addWidget(self.speed)
        self.loop = QCheckBox(self.tr("Loop"))
        srow.addWidget(self.loop)
        lay.addLayout(srow)

        # --- map variable ---------------------------------------------------
        form = QFormLayout()
        mrow = QHBoxLayout()
        self.map_var = QComboBox()
        self.map_var.currentIndexChanged.connect(self._on_map_var)
        mrow.addWidget(self.map_var, 1)
        self.arrows = QCheckBox(self.tr("Arrows"))
        self.arrows.setChecked(True)
        self.arrows.toggled.connect(self._on_arrows)
        mrow.addWidget(self.arrows)
        form.addRow(self.tr("Map shows"), mrow)
        self.dry = QDoubleSpinBox()
        self.dry.setDecimals(3)
        self.dry.setRange(0.0, 10.0)
        self.dry.setSingleStep(0.005)
        self.dry.setValue(0.01)
        self.dry.setSuffix(" m")
        self.dry.setToolTip(self.tr("Water depth below this value is drawn as dry (transparent)"))
        self.dry.editingFinished.connect(self._on_dry)
        form.addRow(self.tr("Dry below"), self.dry)
        lay.addLayout(form)
        lay.addWidget(self._build_arrow_box())

        # --- tabs -----------------------------------------------------------
        self.tabs = QTabWidget()
        lay.addWidget(self.tabs, 1)

        # cursor tab
        ct = QWidget()
        cl = QVBoxLayout(ct)
        crow = QHBoxLayout()
        self.btn_cursor = self._tool("cursor.svg", self.tr("Live cursor: values under the mouse"), plugin.toggle_live_cursor)
        self.btn_cursor.setCheckable(True)
        crow.addWidget(self.btn_cursor)
        self.cursor_xy = QLabel(self.tr("Move the mouse over the mesh"))
        crow.addWidget(self.cursor_xy, 1)
        cl.addLayout(crow)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels([self.tr("Variable"), self.tr("Value"), self.tr("Unit")])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        cl.addWidget(self.table)
        self.tabs.addTab(ct, icon("cursor.svg"), self.tr("Cursor"))

        # time series tab
        tt = QWidget()
        tl = QVBoxLayout(tt)
        trow2 = QHBoxLayout()
        self.btn_probe = self._tool("probe.svg", self.tr("Probe: click points on the map (right-click clears)"), plugin.toggle_probe)
        self.btn_probe.setCheckable(True)
        trow2.addWidget(self.btn_probe)
        self.ts_var = QComboBox()
        self.ts_var.currentIndexChanged.connect(self.update_timeseries)
        trow2.addWidget(self.ts_var, 1)
        self.btn_clear_probes = self._tool("clear.svg", self.tr("Clear probe points"), self.clear_probes)
        trow2.addWidget(self.btn_clear_probes)
        tl.addLayout(trow2)
        self.ts_plot = PlotWidget()
        self.ts_plot.empty_text = self.tr("Turn on the probe tool and click the mesh")
        self.ts_plot.clickedX.connect(self._on_ts_click)
        tl.addWidget(self.ts_plot, 1)
        self.tabs.addTab(tt, icon("probe.svg"), self.tr("Time series"))

        # profile tab
        pt = QWidget()
        pl = QVBoxLayout(pt)
        prow = QHBoxLayout()
        self.btn_profile = self._tool("profile.svg", self.tr("Profile: click vertices, right-click to finish"), plugin.toggle_profile)
        self.btn_profile.setCheckable(True)
        prow.addWidget(self.btn_profile)
        self.pr_var = QComboBox()
        self.pr_var.currentIndexChanged.connect(self.update_profile)
        prow.addWidget(self.pr_var, 1)
        pl.addLayout(prow)
        self.pr_plot = PlotWidget()
        self.pr_plot.empty_text = self.tr("Turn on the profile tool and draw a line")
        pl.addWidget(self.pr_plot, 1)
        self.tabs.addTab(pt, icon("profile.svg"), self.tr("Profile"))

        self.setWidget(w)
        self._set_enabled(False)

    # ---- arrow settings -----------------------------------------------------------
    def _build_arrow_box(self):
        from qgis.gui import QgsCollapsibleGroupBox, QgsColorButton

        st = styles.ArrowStyle.load()
        box = QgsCollapsibleGroupBox(self.tr("Arrow settings"))
        box.setObjectName("TelemacToolkitArrowBox")
        box.setCollapsed(True)
        box.setSaveCollapsedState(False)
        f = QFormLayout(box)
        self.ar_place = QComboBox()
        self.ar_place.addItem(self.tr("Regular grid (same density everywhere)"), "grid")
        self.ar_place.addItem(self.tr("Every mesh node"), "nodes")
        self.ar_place.setToolTip(self.tr("Grid: one arrow per grid cell on the screen, interpolated from the mesh, so fine "
                                         "parts of the mesh are not cluttered and the density stays the same when you zoom."))
        f.addRow(self.tr("Placement"), self.ar_place)
        self.ar_spacing = QSpinBox()
        self.ar_spacing.setRange(8, 400)
        self.ar_spacing.setSuffix(" px")
        self.ar_spacing.setToolTip(self.tr("Distance between arrows on the screen; smaller shows more arrows "
                                           "(use 15 to 25 px for narrow channels)"))
        f.addRow(self.tr("Spacing"), self.ar_spacing)
        self.ar_length = QComboBox()
        self.ar_length.addItem(self.tr("Proportional (fastest = longest)"), "speed")
        self.ar_length.addItem(self.tr("Proportional, fixed scale"), "scale")
        self.ar_length.addItem(self.tr("Same length (direction only)"), "fixed")
        f.addRow(self.tr("Length"), self.ar_length)
        self.ar_longest = QDoubleSpinBox()
        self.ar_longest.setRange(0.0, 100.0)
        self.ar_longest.setDecimals(1)
        self.ar_longest.setSuffix(" mm")
        self.ar_longest.setSpecialValueText(self.tr("auto (fits the spacing)"))
        self.ar_longest.setToolTip(self.tr("Length of the fastest arrow of the whole result (or of every arrow with "
                                           "'Same length'). Auto: 90 % of the spacing, so arrows never overlap."))
        f.addRow(self.tr("Longest arrow"), self.ar_longest)
        self.ar_scale = QDoubleSpinBox()
        self.ar_scale.setRange(0.01, 10000.0)
        self.ar_scale.setDecimals(2)
        self.ar_scale.setSuffix(self.tr(" mm per unit"))
        self.ar_scale.setToolTip(self.tr("Fixed scale: 1 m/s (or 1 unit of the vector) is drawn this long. "
                                         "Use the same value to compare results."))
        f.addRow(self.tr("Scale"), self.ar_scale)
        self.ar_min = QDoubleSpinBox()
        self.ar_min.setRange(0.0, 1000.0)
        self.ar_min.setDecimals(3)
        self.ar_min.setSingleStep(0.01)
        self.ar_min.setToolTip(self.tr("Arrows below this magnitude are not drawn (dry or still water)"))
        f.addRow(self.tr("Hide below"), self.ar_min)
        crow = QHBoxLayout()
        self.ar_color = QgsColorButton()
        self.ar_color.setAllowOpacity(True)
        crow.addWidget(self.ar_color)
        self.ar_byvalue = QCheckBox(self.tr("Coloured by magnitude"))
        crow.addWidget(self.ar_byvalue, 1)
        f.addRow(self.tr("Colour"), crow)
        self.ar_width = QDoubleSpinBox()
        self.ar_width.setRange(0.05, 3.0)
        self.ar_width.setDecimals(2)
        self.ar_width.setSingleStep(0.05)
        self.ar_width.setSuffix(" mm")
        f.addRow(self.tr("Line width"), self.ar_width)
        self.ar_ref = QLabel()
        self.ar_ref.setWordWrap(True)
        f.addRow(self.ar_ref)
        self._arrow_timer = QTimer(self)
        self._arrow_timer.setSingleShot(True)
        self._arrow_timer.setInterval(150)
        self._arrow_timer.timeout.connect(self._apply_arrows)
        self._set_arrow_widgets(st)
        for sig in (self.ar_place.currentIndexChanged, self.ar_spacing.valueChanged, self.ar_length.currentIndexChanged,
                    self.ar_longest.valueChanged, self.ar_scale.valueChanged, self.ar_min.valueChanged,
                    self.ar_color.colorChanged, self.ar_byvalue.toggled, self.ar_width.valueChanged):
            sig.connect(lambda *_: self._arrow_timer.start())
        self.arrow_box = box
        return box

    def _set_arrow_widgets(self, st):
        for w in (self.ar_place, self.ar_spacing, self.ar_length, self.ar_longest, self.ar_scale, self.ar_min,
                  self.ar_color, self.ar_byvalue, self.ar_width):
            w.blockSignals(True)
        self.ar_place.setCurrentIndex(max(0, self.ar_place.findData(st.placement)))
        self.ar_spacing.setValue(int(st.spacing))
        self.ar_length.setCurrentIndex(max(0, self.ar_length.findData(st.length)))
        self.ar_longest.setValue(float(st.longest))
        self.ar_scale.setValue(float(st.scale))
        self.ar_min.setValue(float(st.min_value))
        self.ar_color.setColor(QColor(st.color))
        self.ar_byvalue.setChecked(bool(st.by_value))
        self.ar_width.setValue(float(st.width))
        for w in (self.ar_place, self.ar_spacing, self.ar_length, self.ar_longest, self.ar_scale, self.ar_min,
                  self.ar_color, self.ar_byvalue, self.ar_width):
            w.blockSignals(False)
        self._update_arrow_enabling()

    def arrow_style(self):
        c = self.ar_color.color()
        return styles.ArrowStyle(placement=self.ar_place.currentData(), spacing=self.ar_spacing.value(),
                                 length=self.ar_length.currentData(), longest=self.ar_longest.value(),
                                 scale=self.ar_scale.value(), min_value=self.ar_min.value(),
                                 color=c.name(QColor.NameFormat.HexArgb) if c.alpha() < 255 else c.name(),
                                 by_value=self.ar_byvalue.isChecked(), width=self.ar_width.value())

    def _update_arrow_enabling(self):
        grid = self.ar_place.currentData() == "grid"
        mode = self.ar_length.currentData()
        self.ar_spacing.setEnabled(grid)
        self.ar_longest.setEnabled(mode in ("speed", "fixed"))
        self.ar_scale.setEnabled(mode == "scale")
        self.ar_color.setEnabled(not self.ar_byvalue.isChecked())
        self.ar_longest.setSpecialValueText(self.tr("auto (fits the spacing)") if grid else self.tr("auto (5 mm)"))

    def _apply_arrows(self):
        self._update_arrow_enabling()
        st = self.arrow_style()
        st.save()
        rs = self.rs
        if rs is None:
            self.ar_ref.setText("")
            return
        styles.apply_arrow_style(rs, st)
        self._update_arrow_reference(st)

    def _update_arrow_reference(self, st=None):
        rs = self.rs
        st = st or self.arrow_style()
        pairs = getattr(rs, "group_of_pair", {}) if rs is not None else {}
        if not pairs:
            self.ar_ref.setText(self.tr("This result has no vector variable."))
            return
        (i, _j), g = next(iter(pairs.items()))
        vmax = styles.group_max_magnitude(rs.layer, g)
        unit = rs.units[i].strip().lower() or self.tr("units")
        if st.length == "fixed":
            txt = self.tr("Every arrow %.1f mm; the colour or the map shows the magnitude.") % st.longest_mm()
        else:
            ref = vmax if vmax > 0 else 1.0
            txt = self.tr("Reference: %.3g %s = %.1f mm (fastest in the result: %.3g %s)") % (
                ref, unit, st.length_of(ref, vmax), vmax, unit)
        self.ar_ref.setText("<span style='color:gray'>%s</span>" % txt)

    # ---- helpers --------------------------------------------------------------
    def _tool(self, ico, tip, slot):
        b = QToolButton()
        b.setIcon(icon(ico))
        b.setIconSize(QSize(18, 18))
        b.setToolTip(tip)
        b.setAutoRaise(True)
        b.clicked.connect(lambda *_: slot())
        return b

    def _set_enabled(self, on: bool):
        for wdg in (self.btn_refresh, self.btn_first, self.btn_prev, self.btn_play, self.btn_next, self.btn_last,
                    self.slider, self.btn_peak, self.map_var, self.arrows, self.dry, self.tabs, self.speed, self.loop,
                    self.arrow_box, self.btn_zoom, self.btn_crs, self.btn_calc):
            wdg.setEnabled(on)

    @property
    def frame(self) -> int:
        return self.plugin.timesync.frame

    @property
    def nframes(self) -> int:
        return self.rs.nframes if self.rs else 0

    # ---- result list ------------------------------------------------------------
    def refresh_result_list(self, select=None):
        self._building = True
        cur = select or (self.rs.layer_id if self.rs else None)
        self.result_combo.clear()
        for lid, rs in self.plugin.results.items():
            self.result_combo.addItem(icon("telemac.svg"), rs.layer.name(), lid)
        self._building = False
        if self.result_combo.count() == 0:
            self.set_result(None)
            return
        idx = max(0, self.result_combo.findData(cur)) if cur else 0
        if self.result_combo.currentIndex() != idx:
            self.result_combo.setCurrentIndex(idx)
        else:
            self._on_result_changed(idx)

    def _on_result_changed(self, idx):
        if self._building:
            return
        lid = self.result_combo.itemData(idx)
        self.set_result(self.plugin.results.get(lid))

    def set_result(self, rs):
        if rs is self.rs and rs is not None:
            return
        self._play.stop()
        self.btn_play.setChecked(False)
        self.rs = rs
        self.probes = []
        self.plugin.clear_probe_markers()
        self.sampler = None
        self.profile_line = None
        self.table.setRowCount(0)
        self.ts_plot.clear()
        self.pr_plot.clear()
        if rs is None:
            self.info.setText(self.tr("Open a TELEMAC result (.slf) to start."))
            self.time_label.setText("—")
            self._set_enabled(False)
            self.plugin.timesync.bind(None)
            return
        self._set_enabled(True)
        self._fill_variable_combos()
        self.slider.blockSignals(True)
        self.slider.setRange(0, max(rs.nframes - 1, 0))
        self.slider.blockSignals(False)
        self._update_info()
        self.plugin.timesync.bind(rs)
        self.on_frame(self.plugin.timesync.frame)

    def _update_info(self):
        rs = self.rs
        h = rs.header
        t = rs.times
        dt = float(np.median(np.diff(t))) if t.size > 1 else 0.0
        parts = [
            "%s nodes · %s elements" % (format(h.npoin, ",").replace(",", " "), format(h.nelem, ",").replace(",", " ")),
            self.tr("%d frames every %s") % (rs.nframes, format_duration(dt)) if dt else self.tr("%d frame(s)") % rs.nframes,
        ]
        if rs.derived_from_3d:
            parts.append(self.tr("3D file shown as: %s") % rs.derived_from_3d)
        if rs.has_calendar_date:
            parts.append(self.tr("start %s") % h.start_datetime.strftime("%Y-%m-%d %H:%M"))
        from ..qgis_bridge.results import crs_label

        parts.append(crs_label(rs.layer.crs()))
        html = " · ".join(parts)
        for wmsg in rs.warnings():
            html += "<br><span style='color:#b35900'>⚠ %s</span>" % wmsg
        self.info.setText(html)
        self.info.setToolTip(rs.source_path)

    def _variable_items(self):
        rs = self.rs
        items = []
        for i in range(len(rs.names)):
            items.append((rs.label(i), i))
        for i, j, label in rs.pairs:
            if rs.has_file_magnitude((i, j)):
                continue  # TELEMAC already wrote it (SCALAR VELOCITY)
            u = rs.units[i].strip().lower()
            items.append(("%s magnitude%s" % (label, " (%s)" % u if u else ""), (i, j)))
        return items

    def _fill_variable_combos(self):
        rs = self.rs
        items = self._variable_items()
        for combo in (self.map_var, self.ts_var, self.pr_var):
            combo.blockSignals(True)
            combo.clear()
        # map: only variables QGIS knows as a dataset group
        for label, v in items:
            if (isinstance(v, tuple) and v in getattr(rs, "group_of_pair", {})) or \
               (not isinstance(v, tuple) and v in getattr(rs, "group_of_var", {})):
                self.map_var.addItem(label, v)
        for label, v in items:
            self.ts_var.addItem(label, v)
        if self._water_profile_vars() is not None:
            self.pr_var.addItem(self.tr("Water profile (bed and free surface)"), self.WATER_PROFILE)
        for label, v in items:
            self.pr_var.addItem(label, v)
        dv = rs.default_variable()
        k = self.map_var.findData(dv)
        self.map_var.setCurrentIndex(max(k, 0))
        k = self.ts_var.findData(dv)
        self.ts_var.setCurrentIndex(max(k, 0))
        self.pr_var.setCurrentIndex(0)
        for combo in (self.map_var, self.ts_var, self.pr_var):
            combo.blockSignals(False)
        self.arrows.blockSignals(True)
        self.arrows.setEnabled(bool(getattr(rs, "group_of_pair", {})))
        self.arrows.setChecked(rs.layer.rendererSettings().activeVectorDatasetGroup() >= 0)
        self.arrows.blockSignals(False)
        self.dry.setValue(rs.dry_threshold)
        self._update_arrow_reference()

    def _water_profile_vars(self):
        """(bottom, surface) indices, surface may be computed from depth."""
        rs = self.rs
        b = rs.var_index("bottom")
        s = rs.var_index("free_surface")
        h = rs.var_index("water_depth")
        if b is not None and (s is not None or h is not None):
            return b, s, h
        return None

    # ---- time -------------------------------------------------------------------
    def set_frame(self, i: int):
        if not self.rs:
            return
        self.plugin.timesync.set_frame(int(np.clip(i, 0, self.nframes - 1)))

    def _on_slider(self, v):
        self.set_frame(v)

    def toggle_play(self):
        if self._play.isActive():
            self._play.stop()
            self.btn_play.setChecked(False)
            self.btn_play.setIcon(icon("play.svg"))
        else:
            if self.frame >= self.nframes - 1:
                self.set_frame(0)
            self._play.setInterval(int(1000 / self.speed.value()))
            self._play.start()
            self.btn_play.setChecked(True)
            self.btn_play.setIcon(icon("pause.svg"))

    def _advance(self):
        if self.frame >= self.nframes - 1:
            if self.loop.isChecked():
                self.set_frame(0)
            else:
                self.toggle_play()
            return
        self.set_frame(self.frame + 1)

    def on_frame(self, i: int):
        """Called by TimeSync whenever the clock moves."""
        if not self.rs:
            return
        self.slider.blockSignals(True)
        self.slider.setValue(i)
        self.slider.blockSignals(False)
        self.time_label.setText("%s   ·   %d / %d" % (self.rs.time_label(i), i + 1, self.nframes))
        f, _u = time_axis(self.rs.times)
        self.ts_plot.set_cursor(float(self.rs.times[i]) / f)
        if self.cursor_pt is not None:
            self.update_cursor(*self.cursor_pt)
        self.update_profile()

    def jump_to_peak(self):
        if not self.rs:
            return
        v = self.map_var.currentData()
        if v is None:
            v = self.rs.default_variable()
        self.set_frame(self.rs.peak_frame(v, self.probes[0] if self.probes else None))

    def _on_ts_click(self, x):
        if not self.rs:
            return
        f, _u = time_axis(self.rs.times)
        self.set_frame(int(np.argmin(np.abs(self.rs.times - x * f))))

    # ---- map style ------------------------------------------------------------------
    def _on_map_var(self, _idx):
        if self.rs is None:
            return
        v = self.map_var.currentData()
        if v is not None:
            styles.set_active_variable(self.rs, v)

    def _on_arrows(self, on):
        if self.rs is not None:
            styles.set_arrows_visible(self.rs, on)

    def _on_dry(self):
        rs = self.rs
        if rs is None:
            return
        rs.dry_threshold = float(self.dry.value())
        h = rs.var_index("water_depth")
        g = getattr(rs, "group_of_var", {}).get(h) if h is not None else None
        if g is not None:
            styles.style_scalar_group(rs, g, "depth", rs.dry_threshold)
            rs.layer.triggerRepaint()
        self.update_profile()

    # ---- cursor -----------------------------------------------------------------------
    def update_cursor(self, x: float, y: float):
        rs = self.rs
        if rs is None:
            return
        self.cursor_pt = (x, y)
        rows = rs.values_at(x, y, self.frame)
        self.cursor_xy.setText("x = %.2f   y = %.2f" % (x, y))
        if rows is None:
            self.table.setRowCount(0)
            self.cursor_xy.setText(self.cursor_xy.text() + self.tr("   (outside the mesh)"))
            return
        self.table.setRowCount(len(rows))
        for r, (label, value, unit) in enumerate(rows):
            self.table.setItem(r, 0, QTableWidgetItem(label))
            it = QTableWidgetItem("%.4g" % value if abs(value) < 1e6 else "%.6g" % value)
            it.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
            self.table.setItem(r, 1, it)
            self.table.setItem(r, 2, QTableWidgetItem(unit))

    def cursor_left(self):
        self.cursor_xy.setText(self.tr("Move the mouse over the mesh"))

    # ---- time series ---------------------------------------------------------------------
    def add_probe(self, x: float, y: float):
        if self.rs is None:
            return
        self.probes.append((x, y))
        self.tabs.setCurrentIndex(1)
        self.update_timeseries()

    def clear_probes(self):
        self.probes = []
        self.plugin.clear_probe_markers()
        self.update_timeseries()

    def update_timeseries(self, *_):
        rs = self.rs
        if rs is None or not self.probes:
            self.ts_plot.clear()
            return
        v = self.ts_var.currentData()
        if v is None:
            return
        f, unit = time_axis(rs.times)
        t = rs.times / f
        series = []
        for k, (x, y) in enumerate(self.probes):
            s = rs.series_at(x, y, v)
            if s is None:
                continue
            series.append(Series("P%d" % (k + 1), t, s, PALETTE[k % len(PALETTE)]))
        label = self.ts_var.currentText()
        xlabel = self.tr("Time (%s)") % unit
        if rs.has_calendar_date:
            xlabel += self.tr(" since %s") % rs.header.start_datetime.strftime("%Y-%m-%d %H:%M")
        self.ts_plot.set_series(series, title=label, xlabel=xlabel, ylabel=label)
        self.ts_plot.set_cursor(float(rs.times[self.frame]) / f)

    # ---- profile -----------------------------------------------------------------------------
    def set_profile_line(self, line_xy):
        if self.rs is None:
            return
        self.profile_line = line_xy
        self.sampler = self.rs.sampler(line_xy)
        self.tabs.setCurrentIndex(2)
        self.update_profile()

    def update_profile(self, *_):
        rs = self.rs
        if rs is None or self.sampler is None:
            return
        v = self.pr_var.currentData()
        sp = self.sampler
        fr = self.frame
        if v == self.WATER_PROFILE:
            b, s, h = self._water_profile_vars()
            bed = rs.profile(sp, b, fr)
            if s is not None:
                surf = rs.profile(sp, s, fr)
            else:
                surf = bed + rs.profile(sp, h, fr)
            if h is not None:
                depth = rs.profile(sp, h, fr)
                dry = ~(depth > rs.dry_threshold)
            else:
                dry = ~(surf - bed > rs.dry_threshold)
            wet_surf = surf.copy()
            wet_surf[dry] = np.nan
            series = [Series(self.tr("Bed"), sp.s, bed, "#7a5230", 1.8),
                      Series(self.tr("Free surface"), sp.s, wet_surf, "#0072B2", 1.8)]
            band = (sp.s, bed, np.where(dry, np.nan, surf), "#56B4E9")
            self.pr_plot.set_series(series, title=self.tr("Water profile"), xlabel=self.tr("Distance along line (m)"),
                                    ylabel=self.tr("Elevation (m)"), band=band)
        else:
            vals = rs.profile(sp, v, fr)
            label = self.pr_var.currentText()
            self.pr_plot.set_series([Series(label, sp.s, vals, PALETTE[0])], title=label,
                                    xlabel=self.tr("Distance along line (m)"), ylabel=label)

    # ---- refresh ---------------------------------------------------------------------------------
    def refresh(self):
        if self.rs is None:
            return
        n_before = self.rs.nframes
        self.rs.refresh()
        self.slider.setRange(0, max(self.rs.nframes - 1, 0))
        self._update_info()
        self.plugin.timesync.bind(self.rs)
        if self.rs.nframes > n_before:
            self.set_frame(self.rs.nframes - 1)
        self.update_timeseries()
