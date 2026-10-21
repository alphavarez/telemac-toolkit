"""The Mesh workflow panel: project, drawing, bottom and roughness, mesh and export, checks.

Everything the panel does goes through the Processing algorithms of the
provider, so a run here and a run from the toolbox, a model or a script give
the same result. Settings are kept in the project GeoPackage.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from qgis.core import (
    Qgis,
    QgsApplication,
    QgsCoordinateTransform,
    QgsPointXY,
    QgsProcessingAlgRunnerTask,
    QgsProcessingContext,
    QgsProcessingFeedback,
    QgsProject,
    QgsRasterLayer,
    QgsRectangle,
    QgsVectorLayer,
)
from qgis.gui import QgsDockWidget, QgsFieldComboBox, QgsFileWidget, QgsMapLayerComboBox, QgsVertexMarker
from qgis.PyQt.QtCore import QCoreApplication, Qt, pyqtSignal
from qgis.PyQt.QtGui import QColor, QIcon
from qgis.PyQt.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QToolBox,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..mesher.writers import FORMATS
from ..qgis_bridge import meshing as M

ICON_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "icons")
FORMAT_KEYS = list(FORMATS)

#: what each layer is for, shown in the Draw step
HINTS = {
    "domain": "Area to mesh. Holes are islands. Optional size per polygon.",
    "islands": "Areas left out of the mesh (buildings, piers, dry islands).",
    "refinement_areas": "Finer (or coarser) elements inside: size in metres.",
    "refinement_points": "Element size at a point (structures, gauges).",
    "constraint_lines": "Lines the mesh follows. mode = breakline (every vertex kept) or guide (resampled).",
    "fixed_points": "Nodes that must exist.",
    "channels": "Structured strip: one left and one right bank line (role), optional axis; same 'channel' value.",
    "grid_zones": "Rectangles meshed as regular grids (dx, dy).",
    "roughness_zones": "Material id and friction value per zone.",
    "boundaries": "Boundary lines: type wall, Q, H, QH, UV, UVH, free, waves or a code like 455; name; value.",
}


def tr(text):
    return QCoreApplication.translate("TelemacToolkit", text)


class _Feedback(QgsProcessingFeedback):
    """Processing feedback that forwards messages through a queued Qt signal (safe from the task thread)."""

    message = pyqtSignal(str, int)  # text, level (0 info, 1 warning, 2 error)

    def pushInfo(self, info):
        super().pushInfo(info)
        self.message.emit(info, 0)

    def pushWarning(self, warning):
        super().pushWarning(warning)
        self.message.emit(warning, 1)

    def reportError(self, error, fatalError=False):
        super().reportError(error, fatalError)
        self.message.emit(error, 2)


class MeshPanel(QgsDockWidget):
    projectChanged = pyqtSignal(str)

    def __init__(self, iface, parent=None):
        super().__init__(tr("Mesh workflow"), parent)
        self.setObjectName("TelemacToolkitMeshPanel")
        self.iface = iface
        self.path: Optional[str] = None
        self.layers: Dict[str, QgsVectorLayer] = {}
        self.task = None
        self.feedback = None
        self.marker = None
        self.last_output: Optional[str] = None
        self._mesh_cache = None
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)
        top = QHBoxLayout()
        self.project_label = QLabel(tr("No mesh project"))
        self.project_label.setWordWrap(True)
        top.addWidget(self.project_label, 1)
        self.btn_new = QToolButton()
        self.btn_new.setText(tr("New…"))
        self.btn_open = QToolButton()
        self.btn_open.setText(tr("Open…"))
        top.addWidget(self.btn_new)
        top.addWidget(self.btn_open)
        v.addLayout(top)
        self.toolbox = QToolBox()
        v.addWidget(self.toolbox, 1)
        self._build_project_page()
        self._build_draw_page()
        self._build_bottom_page()
        self._build_export_page()
        self._build_check_page()
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setMinimumHeight(90)
        v.addWidget(self.progress)
        v.addWidget(self.log)
        self.setWidget(w)
        self.btn_new.clicked.connect(self.new_project)
        self.btn_open.clicked.connect(self.open_project_dialog)
        self._set_enabled(False)

    # ------------------------------------------------------------------ pages
    def _build_project_page(self):
        p = QWidget()
        f = QFormLayout(p)
        self.crs_label = QLabel("—")
        f.addRow(tr("Coordinates"), self.crs_label)
        self.size = QDoubleSpinBox()
        self.size.setRange(0.0, 1e6)
        self.size.setDecimals(3)
        self.size.setSuffix(" m")
        self.size.setSpecialValueText(tr("auto (from the domain area)"))
        self.size.setToolTip(tr("Element size where no refinement applies. Auto: about 1/40 of the square root of "
                                "the domain area (0.1 m for a 20 m x 1 m flume), shown in the log when you mesh."))
        f.addRow(tr("Default size"), self.size)
        self.growth = QDoubleSpinBox()
        self.growth.setRange(1.0, 3.0)
        self.growth.setSingleStep(0.05)
        self.growth.setValue(1.2)
        self.growth.setToolTip(tr("How fast the element size grows away from refinement objects (1.2 = 20% per element)"))
        f.addRow(tr("Growth rate"), self.growth)
        self.smoothing = QSpinBox()
        self.smoothing.setRange(0, 60)
        self.smoothing.setValue(8)
        f.addRow(tr("Smoothing"), self.smoothing)
        self.btn_flume = QPushButton(tr("Channel or flume template (straight, bend, meanders)…"))
        self.btn_flume.setToolTip(tr("Draw a laboratory flume or schematic channel (straight, circular bend or sine-generated "
                                     "meanders) from (0, 0), with its bed formula; for natural rivers draw the domain on the map"))
        f.addRow(self.btn_flume)
        self.btn_flume.clicked.connect(self.flume_dialog)
        self.toolbox.addItem(p, tr("1 · Project"))

    def _build_draw_page(self):
        p = QWidget()
        v = QVBoxLayout(p)
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels([tr("Layer"), tr("Objects")])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.cellDoubleClicked.connect(lambda r, c: self.edit_layer(r))
        v.addWidget(self.table)
        self.hint = QLabel()
        self.hint.setWordWrap(True)
        v.addWidget(self.hint)
        self.table.currentCellChanged.connect(lambda r, *a: self._show_hint(r))
        row = QHBoxLayout()
        self.btn_edit = QPushButton(tr("Draw in this layer"))
        self.btn_check = QPushButton(tr("Check the drawing"))
        row.addWidget(self.btn_edit)
        row.addWidget(self.btn_check)
        v.addLayout(row)
        self.btn_edit.clicked.connect(lambda: self.edit_layer(self.table.currentRow()))
        self.btn_check.clicked.connect(self.check_drawing)
        self.toolbox.addItem(p, tr("2 · Draw"))

    def _build_bottom_page(self):
        p = QWidget()
        f = QFormLayout(p)
        self.rasters = QListWidget()
        self.rasters.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.rasters.setToolTip(tr("Ticked rasters are used; the first one that covers a node wins. Drag to reorder."))
        self.rasters.setMaximumHeight(110)
        f.addRow(tr("Rasters"), self.rasters)
        self.points = QgsMapLayerComboBox()
        self.points.setFilters(Qgis.LayerFilter.PointLayer)
        self.points.setAllowEmptyLayer(True)
        self.points.setLayer(None)
        f.addRow(tr("Survey points"), self.points)
        self.point_field = QgsFieldComboBox()
        self.point_field.setAllowEmptyFieldName(True)
        self.points.layerChanged.connect(self.point_field.setLayer)
        f.addRow(tr("Value field"), self.point_field)
        self.formula = QLineEdit()
        self.formula.setPlaceholderText("10 - 0.001*x")
        self.formula.setToolTip(tr("Bottom as a formula of x and y (or of s and n along a channel axis), used where no raster or point covers a node"))
        f.addRow(tr("Formula"), self.formula)
        self.friction = QLineEdit()
        self.friction.setPlaceholderText(tr("none"))
        self.friction.setToolTip(tr("Friction value outside the roughness zones (empty = none)"))
        f.addRow(tr("Default friction"), self.friction)
        self.toolbox.addItem(p, tr("3 · Bottom and roughness"))

    def _build_export_page(self):
        p = QWidget()
        f = QFormLayout(p)
        self.fmt = QComboBox()
        for k in FORMAT_KEYS:
            self.fmt.addItem(FORMATS[k].label, k)
        f.addRow(tr("Format"), self.fmt)
        self.fmt_note = QLabel()
        self.fmt_note.setWordWrap(True)
        f.addRow(self.fmt_note)
        self.quads = QCheckBox(tr("Keep structured cells as quadrilaterals"))
        f.addRow(self.quads)
        self.frame = QComboBox()
        self.frame.addItems([tr("Map coordinates as drawn"), tr("Lower-left corner at (0, 0)"), tr("A chosen point at (0, 0)")])
        f.addRow(tr("Coordinates"), self.frame)
        pt = QHBoxLayout()
        self.fx = QDoubleSpinBox()
        self.fy = QDoubleSpinBox()
        for s in (self.fx, self.fy):
            s.setRange(-1e9, 1e9)
            s.setDecimals(3)
            pt.addWidget(s)
        f.addRow(tr("Point (x, y)"), pt)
        self.angle = QDoubleSpinBox()
        self.angle.setRange(-360, 360)
        self.angle.setDecimals(4)
        self.angle.setSuffix(" °")
        f.addRow(tr("Rotation"), self.angle)
        self.store = QCheckBox(tr("TELEMAC: keep the shift in the file header"))
        self.store.setChecked(True)
        f.addRow(self.store)
        self.output = QgsFileWidget()
        self.output.setStorageMode(QgsFileWidget.StorageMode.SaveFile)
        f.addRow(tr("Mesh file"), self.output)
        self.btn_generate = QPushButton(QIcon(os.path.join(ICON_DIR, "mesh.svg")), tr("Generate mesh"))
        f.addRow(self.btn_generate)
        self.btn_generate.clicked.connect(lambda: self.generate())
        self.fmt.currentIndexChanged.connect(self._format_changed)
        self.frame.currentIndexChanged.connect(self._frame_changed)
        self._format_changed()
        self._frame_changed()
        self.toolbox.addItem(p, tr("4 · Mesh and export"))

    def _build_check_page(self):
        p = QWidget()
        f = QFormLayout(p)
        self.btn_quality = QPushButton(tr("Mesh quality layer"))
        f.addRow(self.btn_quality)
        self.btn_quality.clicked.connect(self.quality)
        row = QHBoxLayout()
        self.find_kind = QComboBox()
        self.find_kind.addItems([tr("Node"), tr("Element")])
        self.find_num = QSpinBox()
        self.find_num.setRange(1, 2 ** 31 - 1)
        self.btn_find = QPushButton(tr("Find"))
        row.addWidget(self.find_kind)
        row.addWidget(self.find_num, 1)
        row.addWidget(self.btn_find)
        f.addRow(tr("Number"), row)
        self.found = QLabel()
        self.found.setWordWrap(True)
        f.addRow(self.found)
        self.btn_find.clicked.connect(lambda: self.find(self.find_kind.currentIndex() == 1, self.find_num.value()))
        self.toolbox.addItem(p, tr("5 · Check"))

    # ------------------------------------------------------------------ state
    def _set_enabled(self, on: bool):
        for i in range(self.toolbox.count()):
            self.toolbox.widget(i).setEnabled(on if i < 4 else True)

    def _format_changed(self):
        k = self.fmt.currentData()
        fmt = FORMATS[k]
        self.fmt_note.setText(fmt.note)
        self.quads.setEnabled(fmt.quads)
        self.store.setEnabled(k == "telemac")
        filt = "%s (*%s)" % (fmt.label, fmt.ext)
        self.output.setFilter(filt)
        cur = self.output.filePath()
        if cur:
            base = os.path.splitext(cur)[0]
            for f in FORMATS.values():
                if f.suffix and base.endswith(f.suffix):
                    base = base[: -len(f.suffix)]
            self.output.setFilePath(base + fmt.suffix + fmt.ext)

    def _frame_changed(self):
        point = self.frame.currentIndex() == 2
        self.fx.setEnabled(point)
        self.fy.setEnabled(point)
        self.angle.setEnabled(self.frame.currentIndex() > 0)

    def _show_hint(self, row):
        if 0 <= row < self.table.rowCount():
            key = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
            self.hint.setText(HINTS.get(key, ""))

    def say(self, text, level=0):
        prefix = {1: "⚠ ", 2: "✖ "}.get(level, "")
        self.log.appendPlainText(prefix + text)

    # ------------------------------------------------------------------ project
    def new_project(self):
        from qgis.gui import QgsProjectionSelectionDialog

        path, _ = QFileDialog.getSaveFileName(self, tr("New mesh project"), "", "GeoPackage (*.gpkg)")
        if not path:
            return
        if not path.lower().endswith(".gpkg"):
            path += ".gpkg"
        from qgis.PyQt.QtWidgets import QMessageBox

        local = QMessageBox.question(self, tr("Coordinates"), tr("Use local coordinates (laboratory flume or schematic "
                                                                  "channel, drawn from (0, 0))?\n\nNo = choose a map CRS."))
        if local == QMessageBox.StandardButton.Yes:
            crs = M.local_crs()
        else:
            dlg = QgsProjectionSelectionDialog(self)
            dlg.setCrs(QgsProject.instance().crs())
            from ..compat import exec_dialog

            if not exec_dialog(dlg):
                return
            crs = dlg.crs()
        self.create_project(path, crs)

    def create_project(self, path, crs):
        if os.path.exists(path):
            os.remove(path)
        M.create_mesh_project(path, crs)
        if not QgsProject.instance().crs().isValid() or M.is_local_crs(crs):
            QgsProject.instance().setCrs(crs)
        self.set_project(path, add_to_map=False)

    def open_project_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, tr("Open mesh project"), "", "GeoPackage (*.gpkg)")
        if path:
            self.set_project(path, add_to_map=True)

    def set_project(self, path: str, add_to_map: bool = False):
        """Make ``path`` the current project; its layers on the map are reused, else added."""
        self.path = path
        on_map = self._map_layers(path)
        if on_map:
            self.layers = on_map
        else:
            self.layers = M.open_mesh_project(path)
            if add_to_map:
                M.add_project_to_map(self.layers, path)
        dom = self.layers.get("domain")
        crs = dom.crs() if dom is not None else None
        self.project_label.setText("<b>%s</b>" % os.path.basename(path))
        self.crs_label.setText((tr("local (laboratory)") if M.is_local_crs(crs) else crs.authid() or crs.description()) if crs else "—")
        st = M.read_settings(path)
        self._load_settings(st)
        self._fill_layers()
        self._fill_rasters()
        if not self.output.filePath():
            fmt = FORMATS[self.fmt.currentData()]
            self.output.setFilePath(os.path.splitext(path)[0] + "_mesh" + fmt.suffix + fmt.ext)
        self._set_enabled(True)
        self.projectChanged.emit(path)

    def _map_layers(self, path) -> Dict[str, QgsVectorLayer]:
        out = {}
        target = os.path.normcase(os.path.abspath(path))
        for lyr in QgsProject.instance().mapLayers().values():
            pp = M.project_path_of(lyr)
            if pp and os.path.normcase(os.path.abspath(pp)) == target:
                name = lyr.source().split("layername=")[-1].split("|")[0] if "layername=" in lyr.source() else ""
                key = M._OLD_LAYERS.get(name, name)
                if key:
                    out[key] = lyr
        return out

    def _load_settings(self, st: Dict[str, str]):
        def num(k, d):
            try:
                return float(st.get(k) or d)
            except ValueError:
                return d

        self.size.setValue(num("default_size", 0.0))
        self.growth.setValue(num("growth", 1.2))
        self.smoothing.setValue(int(num("smoothing", 8)))
        self.formula.setText(st.get("bottom_formula", ""))
        k = st.get("format", "telemac")
        if k in FORMAT_KEYS:
            self.fmt.setCurrentIndex(FORMAT_KEYS.index(k))
        self.quads.setChecked(st.get("quads") == "1")
        self.frame.setCurrentIndex(M.FRAME_MODES.index(st.get("frame", "keep")) if st.get("frame", "keep") in M.FRAME_MODES else 0)
        self.fx.setValue(num("frame_x", 0.0))
        self.fy.setValue(num("frame_y", 0.0))
        self.angle.setValue(num("frame_angle", 0.0))
        self.store.setChecked(st.get("store_offset", "1") != "0")
        if st.get("mesh_file"):
            self.output.setFilePath(st["mesh_file"])

    def settings(self) -> Dict[str, str]:
        return {"default_size": "%g" % self.size.value(), "growth": "%g" % self.growth.value(),
                "smoothing": str(self.smoothing.value()), "bottom_formula": self.formula.text().strip(),
                "format": self.fmt.currentData(), "quads": "1" if self.quads.isChecked() else "0",
                "frame": M.FRAME_MODES[self.frame.currentIndex()], "frame_x": "%.6f" % self.fx.value(),
                "frame_y": "%.6f" % self.fy.value(), "frame_angle": "%g" % self.angle.value(),
                "store_offset": "1" if self.store.isChecked() else "0", "mesh_file": self.output.filePath()}

    def _fill_layers(self):
        keys = [k for k in M.LAYER_KEYS if k in self.layers] + [k for k in ("hard_lines", "soft_lines") if k in self.layers]
        self.table.setRowCount(len(keys))
        for r, k in enumerate(keys):
            lyr = self.layers[k]
            it = QTableWidgetItem(lyr.name())
            it.setData(Qt.ItemDataRole.UserRole, k)
            it.setToolTip(HINTS.get(k, ""))
            self.table.setItem(r, 0, it)
            self.table.setItem(r, 1, QTableWidgetItem(str(lyr.featureCount())))
        if keys:
            self.table.setCurrentCell(0, 0)

    def _fill_rasters(self):
        self.rasters.clear()
        for lyr in QgsProject.instance().mapLayers().values():
            if isinstance(lyr, QgsRasterLayer):
                it = QListWidgetItem(lyr.name())
                it.setData(Qt.ItemDataRole.UserRole, lyr.id())
                it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                it.setCheckState(Qt.CheckState.Unchecked)
                self.rasters.addItem(it)

    def refresh(self):
        if self.path:
            for lyr in self.layers.values():
                lyr.reload()
            self._fill_layers()
            self._fill_rasters()

    def edit_layer(self, row):
        if row < 0 or row >= self.table.rowCount():
            return
        key = self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        lyr = self.layers.get(key)
        if lyr is None:
            return
        try:
            self.iface.setActiveLayer(lyr)
            if not lyr.isEditable():
                lyr.startEditing()
            act = self.iface.actionAddFeature()
            if act is not None:
                act.trigger()
        except AttributeError:
            pass
        self.say(tr("Drawing in %s. %s") % (lyr.name(), HINTS.get(key, "")))

    # ------------------------------------------------------------------ actions
    def _commit_edits(self):
        for lyr in self.layers.values():
            if lyr.isEditable() and lyr.isModified():
                lyr.commitChanges(False)

    def check_drawing(self):
        self._commit_edits()
        try:
            b = M.build_mesh_input(self.layers, None, self.size.value() or M.auto_default_size(self.layers) or 1.0, self.growth.value(), self.smoothing.value())
        except Exception as exc:
            self.say(str(exc), 2)
            return False
        for line in b.info:
            self.say(line)
        for w in b.warnings:
            self.say(w, 1)
        if not self.size.value():
            auto = M.auto_default_size(self.layers)
            self.say(tr("No default size set: %g will be used (step 1 to change it).") % auto if auto
                     else tr("Set a default size in step 1."), 1 if not auto else 0)
        self.say(tr("The drawing is ready to mesh.") if not b.warnings else tr("The drawing can be meshed; check the warnings."))
        return True

    def flume_dialog(self):
        if not self.path:
            return
        try:
            import processing

            processing.execAlgorithmDialog("telemactoolkit:channeltemplate", {"PROJECT": self.path})
        except Exception as exc:
            self.say(str(exc), 2)
        self.set_project(self.path)

    def _params(self) -> dict:
        checked = [self.rasters.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.rasters.count())
                   if self.rasters.item(i).checkState() == Qt.CheckState.Checked]
        params = {
            "PROJECT": self.path, "DEFAULT_SIZE": self.size.value(), "GROWTH": self.growth.value(),
            "SMOOTHING": self.smoothing.value(), "RASTERS": [QgsProject.instance().mapLayer(i) for i in checked],
            "FORMULA": self.formula.text().strip(), "FORMAT": self.fmt.currentIndex(), "QUADS": self.quads.isChecked(),
            "FRAME": self.frame.currentIndex(), "FRAME_ANGLE": self.angle.value(), "STORE_OFFSET": self.store.isChecked(),
            "OUTPUT": self.output.filePath(), "LOAD": True,
        }
        if self.fmt.currentData() == "telemac" and self.output.filePath():
            # a .liq template next to the mesh, with the columns TELEMAC expects for these boundaries
            params["OUTPUT_LIQ"] = os.path.splitext(self.output.filePath())[0] + ".liq"
        dom = self.layers.get("domain")
        if self.frame.currentIndex() == 2 and dom is not None:
            params["FRAME_POINT"] = "%.6f,%.6f" % (self.fx.value(), self.fy.value())
        if self.points.currentLayer() is not None:
            params["POINTS"] = self.points.currentLayer()
            params["POINT_FIELD"] = self.point_field.currentField()
        if self.friction.text().strip():
            try:
                params["FRICTION"] = float(self.friction.text())
            except ValueError:
                self.say(tr("Default friction is not a number; ignored."), 1)
        return params

    def generate(self, wait: bool = False):
        """Run Generate mesh in the background (or synchronously with ``wait``)."""
        if not self.path:
            self.say(tr("Open or create a mesh project first."), 2)
            return None
        if not self.output.filePath():
            self.say(tr("Choose a mesh file in step 4."), 2)
            return None
        self._commit_edits()
        try:
            M.write_settings(self.path, self.settings())
        except Exception as exc:
            self.say(tr("Settings not saved: %s") % exc, 1)
        params = self._params()
        alg = QgsApplication.processingRegistry().createAlgorithmById("telemactoolkit:generatemesh")
        if alg is None:
            self.say(tr("The TELEMAC Toolkit Processing provider is not loaded."), 2)
            return None
        self.feedback = _Feedback()
        self.feedback.message.connect(self.say)
        self.feedback.progressChanged.connect(lambda p: self.progress.setValue(int(p)))
        self.log.clear()
        self.say(tr("Generating the mesh…"))
        self.btn_generate.setEnabled(False)
        if wait:
            import processing

            try:
                res = processing.run(alg, params, feedback=self.feedback)
            except Exception as exc:
                self.say(str(exc), 2)
                res = None
            self._done(res is not None, res or {})
            return res
        context = QgsProcessingContext()
        context.setProject(QgsProject.instance())
        context.setFeedback(self.feedback)
        self._context = context  # keep it alive while the task runs
        self.task = QgsProcessingAlgRunnerTask(alg, params, context, self.feedback)
        self.task.executed.connect(self._done)
        QgsApplication.taskManager().addTask(self.task)
        return None

    def _done(self, ok, results):
        self.btn_generate.setEnabled(True)
        self.progress.setValue(100 if ok else 0)
        if ok and results.get("OUTPUT"):
            self.last_output = results["OUTPUT"]
            self._mesh_cache = None
            self.say(tr("Done: %s") % self.last_output)
        elif not ok:
            self.say(tr("The mesh was not generated (see the messages above)."), 2)
        self.task = None

    def quality(self):
        if not self.last_output or os.path.splitext(self.last_output)[1].lower() not in (".slf", ".2dm", ".nc", ".mesh"):
            self.say(tr("Generate a mesh in a format QGIS reads first (SELAFIN, 2DM, UGRID, MIKE)."), 1)
            return None
        import processing

        res = processing.run("telemactoolkit:meshquality", {"INPUT": self.last_output, "OUTPUT": "memory:"})
        lyr = res["OUTPUT"]
        if isinstance(lyr, QgsVectorLayer):
            lyr.setName(tr("Mesh quality"))
            QgsProject.instance().addMapLayer(lyr)
        return lyr

    def find(self, element: bool, number: int):
        """Centre the map on node or element ``number`` (1-based, as in TELEMAC listings)."""
        if not self.last_output:
            self.say(tr("Generate or open a mesh first."), 1)
            return None
        if self._mesh_cache is None:
            self._mesh_cache = M.read_mesh_file(self.last_output)
        mesh, crs = self._mesh_cache
        if element:
            if not 1 <= number <= mesh.nelem:
                self.found.setText(tr("The mesh has %d elements.") % mesh.nelem)
                return None
            cx, cy = mesh.element_centroids()
            x, y = float(cx[number - 1]), float(cy[number - 1])
        else:
            if not 1 <= number <= mesh.npoin:
                self.found.setText(tr("The mesh has %d nodes.") % mesh.npoin)
                return None
            x, y = float(mesh.x[number - 1]), float(mesh.y[number - 1])
        z = "" if element or mesh.z is None else tr(", bottom %.3f") % mesh.z[number - 1]
        self.found.setText("%s %d: x %.3f, y %.3f%s" % (tr("Element") if element else tr("Node"), number, x, y, z))
        canvas = self.iface.mapCanvas()
        pt = QgsPointXY(x, y)
        dest = canvas.mapSettings().destinationCrs()
        if crs is not None and crs.isValid() and dest.isValid() and crs != dest:
            pt = QgsCoordinateTransform(crs, dest, QgsProject.instance()).transform(pt)
        if self.marker is None:
            self.marker = QgsVertexMarker(canvas)
            self.marker.setIconType(QgsVertexMarker.IconType.ICON_CIRCLE)
            self.marker.setColor(QColor("#d62728"))
            self.marker.setIconSize(16)
            self.marker.setPenWidth(3)
        self.marker.setCenter(pt)
        canvas.setCenter(pt)
        canvas.refresh()
        return pt

    def autodetect(self):
        """Pick up a mesh project already on the map (its group remembers the GeoPackage)."""
        if self.path:
            return self.path
        root = QgsProject.instance().layerTreeRoot()
        for g in root.findGroups(True) if hasattr(root, "findGroups") else root.children():
            path = g.customProperty("telemac_toolkit/mesh_project") if hasattr(g, "customProperty") else None
            if path and os.path.exists(str(path)):
                self.set_project(str(path))
                return self.path
        return None

    def cleanup(self):
        if self.marker is not None:
            self.iface.mapCanvas().scene().removeItem(self.marker)
            self.marker = None
