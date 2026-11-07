"""The Mesh workflow panel: project, drawing, bottom and roughness, mesh and export, checks.

Everything the panel does goes through the Processing algorithms of the
provider, so a run here and a run from the toolbox, a model or a script give
the same result. Settings are kept in the project GeoPackage.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

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
    QgsVectorLayer,
)
from qgis.gui import QgsDockWidget, QgsFieldComboBox, QgsFileWidget, QgsMapLayerComboBox, QgsVertexMarker
from qgis.PyQt.QtCore import QCoreApplication, Qt, QTimer, pyqtSignal
from qgis.PyQt.QtGui import QColor, QIcon
from qgis.PyQt.QtWidgets import (
    QSizePolicy,
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


def _setting(key: str, default: str = "") -> str:
    from qgis.core import QgsSettings

    return str(QgsSettings().value("telemac_toolkit/" + key, default))


def _set_setting(key: str, value: str) -> None:
    from qgis.core import QgsSettings

    QgsSettings().setValue("telemac_toolkit/" + key, value)

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


#: mesh types in step 4: what the method (step 1) and the quadrilateral option (step 4) add up to
MESH_TYPES = [("triangles", "Triangles"), ("rectangles", "Rectangles where drawn"),
              ("quads", "Quadrilaterals everywhere (Gmsh)")]


def _narrow(combo, chars: int = 14):
    """A combo box that does not widen the panel to its longest item (the list still shows them whole)."""
    try:
        combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        combo.setMinimumContentsLength(chars)
    except AttributeError:
        pass
    return combo


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


class _LogFeedback(QgsProcessingFeedback):
    """Processing messages to the panel log."""

    def __init__(self, say):
        super().__init__()
        self._say = say

    def pushInfo(self, info):
        self._say(info)

    def pushWarning(self, warning):
        self._say(warning, 1)

    def reportError(self, error, fatalError=False):
        self._say(error, 2)


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
        self._build_run_page()
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
        self.outline_tol = QDoubleSpinBox()
        self.outline_tol.setRange(0.0, 1e4)
        self.outline_tol.setDecimals(3)
        self.outline_tol.setSingleStep(0.1)
        self.outline_tol.setSuffix(" m")
        self.outline_tol.setSpecialValueText(tr("auto (5% of the element size)"))
        self.outline_tol.setToolTip(tr(
            "How far the mesh boundary may stray from the drawn outline, islands and guide lines. The outline is "
            "resampled to the element size: its corners stay nodes, and other vertices stay too where leaving them "
            "out would move the boundary more than this. Smaller = closer fit, with more (and smaller) elements "
            "along wiggly outlines."))
        f.addRow(tr("Follow the outline within"), self.outline_tol)
        self.method = QComboBox()
        from ..mesher import gmsh_engine

        has_gmsh = gmsh_engine.available()
        for key, label in (("balanced", tr("Balanced")), ("quality", tr("Quality-guaranteed")),
                           ("gmsh_delaunay", tr("Gmsh Delaunay (fastest)")), ("gmsh_frontal", tr("Gmsh Frontal-Delaunay (best shapes)")),
                           ("gmsh_meshadapt", tr("Gmsh MeshAdapt (most robust)")), ("gmsh_quads", tr("Gmsh quadrilaterals"))):
            self.method.addItem(label, key)
            if key.startswith("gmsh_") and not has_gmsh:
                item = self.method.model().item(self.method.count() - 1)
                item.setEnabled(False)
                item.setToolTip(M.gmsh_missing_message())
        self.method.setToolTip(tr("Balanced: nodes on regular patterns, then smoothing (fast, smooth size changes). "
                                  "Quality-guaranteed: also refines until every triangle angle reaches the target, "
                                  "except next to sharper corners of the drawing; no final smoothing. "
                                  "Gmsh: Gmsh fills the inside (boundaries, lines and structured parts stay the plugin's); "
                                  "needs the gmsh Python package. Compare methods (step 4) measures them on your project."))
        f.addRow(tr("Method"), _narrow(self.method))
        self.min_angle = QDoubleSpinBox()
        self.min_angle.setRange(20.0, 30.0)
        self.min_angle.setSingleStep(1.0)
        self.min_angle.setDecimals(0)
        self.min_angle.setSuffix(" °")
        self.min_angle.setValue(30.0)
        self.min_angle.setToolTip(tr("Smallest triangle angle for the quality-guaranteed method (above 30° the refinement "
                                     "does not converge)"))
        f.addRow(tr("Smallest angle"), self.min_angle)
        self.method.currentIndexChanged.connect(lambda _i: self.min_angle.setEnabled(self.method.currentData() == "quality"))
        self.min_angle.setEnabled(False)
        self.terrain_tol = QDoubleSpinBox()
        self.terrain_tol.setRange(0.0, 100.0)
        self.terrain_tol.setDecimals(3)
        self.terrain_tol.setSingleStep(0.01)
        self.terrain_tol.setSpecialValueText(tr("off"))
        self.terrain_tol.setSuffix(" m")
        self.terrain_tol.setToolTip(tr("Terrain refinement: smaller elements where linear elements would miss the DEM (step 3, "
                                       "first ticked raster) by more than this; flat ground keeps the default size"))
        f.addRow(tr("Terrain tolerance"), self.terrain_tol)
        self.terrain_min = QDoubleSpinBox()
        self.terrain_min.setRange(0.0, 1e6)
        self.terrain_min.setDecimals(3)
        self.terrain_min.setSpecialValueText(tr("auto (default size / 15)"))
        self.terrain_min.setToolTip(tr("Smallest element the terrain refinement may ask for"))
        f.addRow(tr("Terrain smallest size"), self.terrain_min)
        self.max_nodes = QSpinBox()
        self.max_nodes.setRange(0, 2000000000)
        self.max_nodes.setSingleStep(500000)
        self.max_nodes.setValue(2000000)
        self.max_nodes.setSpecialValueText(tr("no limit"))
        self.max_nodes.setToolTip(tr("Generate mesh stops at once, before meshing, when the estimated number of nodes is above this. "
                                     "Raise it for a large model; 0 removes the check."))
        f.addRow(tr("Node limit"), self.max_nodes)
        self.btn_suggest = QPushButton(tr("Suggest sizes from the DEM…"))
        self.btn_suggest.setToolTip(tr("From the DEM ticked in step 3 (or the only DEM on the map): smallest size of two pixels, "
                                       "terrain tolerance for such a DEM, a floodplain size from the terrain and the domain, "
                                       "spacing across and along the rivers, the node count and the time step; Apply fills the "
                                       "sizes above"))
        f.addRow(self.btn_suggest)
        self.btn_suggest.clicked.connect(lambda: self.suggest_sizes())
        self.btn_flume = QPushButton(tr("Channel or flume template…"))
        self.btn_flume.setToolTip(tr("Draw a laboratory flume or schematic channel (straight, circular bend or sine-generated "
                                     "meanders) from (0, 0), with its bed formula; for natural rivers use the river reach helper"))
        f.addRow(self.btn_flume)
        self.btn_flume.clicked.connect(self.flume_dialog)
        self.btn_river = QPushButton(tr("River reach…"))
        self.btn_river.setToolTip(tr("A natural river from a centreline (constant width, or banks found in a DEM) or from two "
                                     "bank lines: structured channel, domain with floodplain, inflow and outflow lines. "
                                     "Several centrelines that join make a river network with its confluences"))
        f.addRow(self.btn_river)
        self.btn_river.clicked.connect(self.river_dialog)
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
        self.btn_fix = QPushButton(tr("Fix the safe problems"))
        self.btn_fix.setToolTip(tr("Remove near-duplicate vertices and small spikes, repair invalid polygons and snap line "
                                   "ends onto the outline, in the layers' edit buffers (review, then save or undo)"))
        self.btn_fix.setEnabled(False)
        v.addWidget(self.btn_fix)
        self.btn_snap = QPushButton(tr("Snapping for this project"))
        self.btn_snap.setToolTip(tr("Turn on QGIS snapping to the vertices and edges of the domain, islands, channels, "
                                    "lines and boundary layers, with snapping at intersections and topological editing, "
                                    "so lines end exactly on the outline and on each other; also tracing (T) to follow the "
                                    "outline when drawing boundary lines, the Advanced Digitizing panel and hover tips"))
        v.addWidget(self.btn_snap)
        self.auto_snap = QCheckBox(tr("Snapping and hover tips on opening"))
        self.auto_snap.setToolTip(tr("Set snapping, tracing, the Advanced Digitizing panel and hover tips when a project opens"))
        self.auto_snap.setChecked(_setting("auto_snapping", "1") == "1")
        self.auto_snap.toggled.connect(lambda on: _set_setting("auto_snapping", "1" if on else "0"))
        v.addWidget(self.auto_snap)
        self.btn_breaklines = QPushButton(tr("Breaklines from the DEM…"))
        self.btn_breaklines.setToolTip(tr("Crest and toe lines where the terrain bends sharply, added to the constraint lines"))
        v.addWidget(self.btn_breaklines)
        self.btn_breaklines.clicked.connect(self.breaklines_dialog)
        self.btn_edit.clicked.connect(lambda: self.edit_layer(self.table.currentRow()))
        self.btn_check.clicked.connect(lambda: self.check_drawing())
        self.btn_fix.clicked.connect(self.fix_drawing)
        self.btn_snap.clicked.connect(self.snapping_preset)
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
        from ..mesher import gmsh_engine

        self.mesh_type = QComboBox()
        for key, label in MESH_TYPES:
            self.mesh_type.addItem(tr(label), key)
        if not gmsh_engine.available():
            item = self.mesh_type.model().item(self.mesh_type.count() - 1)
            item.setEnabled(False)
            item.setToolTip(M.gmsh_missing_message())
        self.mesh_type.setToolTip(tr(
            "Triangles: unstructured triangles; structured channels and grid zones drawn in step 2 are still built "
            "cell by cell, and written as triangles. Rectangles where drawn: those structured parts keep their "
            "rectangular cells in the formats that accept them. Quadrilaterals everywhere: Gmsh fills the inside "
            "with quadrilaterals (needs the gmsh package). TELEMAC-2D computes on triangles only: in a SELAFIN "
            "file every rectangle is split in two along its shorter diagonal."))
        f.addRow(tr("Mesh type"), _narrow(self.mesh_type))
        self.mesh_type_note = QLabel()
        self.mesh_type_note.setWordWrap(True)
        f.addRow(self.mesh_type_note)
        self.fmt = QComboBox()
        for k in FORMAT_KEYS:
            self.fmt.addItem(FORMATS[k].label, k)
        f.addRow(tr("Format"), _narrow(self.fmt))
        self.fmt_note = QLabel()
        self.fmt_note.setWordWrap(True)
        f.addRow(self.fmt_note)
        self.quads = QCheckBox(tr("Keep structured cells as quadrilaterals"))
        f.addRow(self.quads)
        self.frame = QComboBox()
        self.frame.addItems([tr("Map coordinates as drawn"), tr("Lower-left corner at (0, 0)"), tr("A chosen point at (0, 0)")])
        f.addRow(tr("Coordinates"), _narrow(self.frame))
        pt = QHBoxLayout()
        self.fx = QDoubleSpinBox()
        self.fy = QDoubleSpinBox()
        for s in (self.fx, self.fy):
            s.setRange(-1e9, 1e9)
            s.setDecimals(3)
            s.setMinimumWidth(70)
            s.setSizePolicy(QSizePolicy.Policy.Ignored, s.sizePolicy().verticalPolicy())  # never widen the panel
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
        self.double = QCheckBox(tr("TELEMAC: double precision (SERAFIND)"))
        self.double.setToolTip(tr(
            "Coordinates and values in 8-byte floats. Not needed for accuracy: single precision with the origin in the "
            "header already resolves 0.01 mm. Files are twice as large, and the steering file must say "
            "GEOMETRY FILE FORMAT = 'SERAFIND' (the steering file tool in step 6 writes it)."))
        f.addRow(self.double)
        self.output = QgsFileWidget()
        self.output.setStorageMode(QgsFileWidget.StorageMode.SaveFile)
        f.addRow(tr("Mesh file"), self.output)
        self.btn_generate = QPushButton(QIcon(os.path.join(ICON_DIR, "mesh.svg")), tr("Generate mesh"))
        f.addRow(self.btn_generate)
        self.btn_generate.clicked.connect(lambda: self.generate())
        self.btn_compare = QPushButton(tr("Compare methods…"))
        self.btn_compare.setToolTip(tr("Mesh this project with several methods (and add meshes made elsewhere) and compare "
                                       "element shape, size changes, outline, bed error against a DEM, time step and nodes"))
        f.addRow(self.btn_compare)
        self.btn_compare.clicked.connect(self.compare_dialog)
        self.fmt.currentIndexChanged.connect(self._format_changed)
        self.frame.currentIndexChanged.connect(self._frame_changed)
        self._format_changed()
        self._frame_changed()
        self.mesh_type.currentIndexChanged.connect(self._mesh_type_changed)
        self.method.currentIndexChanged.connect(lambda _i: self._sync_mesh_type())
        self.quads.toggled.connect(lambda _on: self._sync_mesh_type())
        self._sync_mesh_type()
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
        ed = QLabel(tr("Edit the generated mesh:"))
        f.addRow(ed)
        self.btn_edit_mesh = QPushButton(tr("Edit in QGIS (Mesh Digitizing)"))
        self.btn_edit_mesh.setToolTip(tr("Open a copy of the mesh in edit mode: add, remove or move vertices, flip edges, split "
                                         "faces with QGIS's Mesh Digitizing toolbar; then Rebuild after editing"))
        self.btn_rebuild = QPushButton(tr("Rebuild after editing"))
        self.btn_rebuild.setToolTip(tr("Save the edits and write the mesh back with its bottom, friction and .cli"))
        self.btn_remesh = QPushButton(tr("Re-mesh inside polygons…"))
        self.btn_remesh.setToolTip(tr("Replace the elements inside polygons by a new mesh, joined node for node to the rest"))
        self.btn_mmg = QPushButton(tr("Optimise or adapt (MMG)…"))
        self.btn_mmg.setToolTip(tr("Better shapes, or sizes from the terrain or from a TELEMAC result, with MMG; the boundary stays"))
        self.btn_cli_lines = QPushButton(tr("Boundary lines from the .cli"))
        self.btn_cli_lines.setToolTip(tr("The liquid boundaries of the mesh as lines to edit, then Boundary conditions from lines"))
        for b in (self.btn_edit_mesh, self.btn_rebuild, self.btn_remesh, self.btn_mmg, self.btn_cli_lines):
            f.addRow(b)
        self.btn_edit_mesh.clicked.connect(self.edit_mesh_start)
        self.btn_rebuild.clicked.connect(self.edit_mesh_finish)
        self.btn_remesh.clicked.connect(lambda: self._open_tool("remeshregion", "INPUT", "_remeshed.slf"))
        self.btn_cli_lines.clicked.connect(self.cli_lines)
        self.btn_mmg.clicked.connect(lambda: self._open_tool("mmgadapt", "INPUT", "_mmg.slf"))
        self._editing = None
        self.toolbox.addItem(p, tr("5 · Check"))

    def _build_run_page(self):
        p = QWidget()
        f = QFormLayout(p)
        note = QLabel(tr("For the TELEMAC mesh generated in step 4 (or any SELAFIN mesh with its .cli):"))
        note.setWordWrap(True)
        f.addRow(note)
        self.run_buttons = {}
        tools = (
            ("stagedischarge", tr("Rating curve (outflow)…"),
             tr("Normal-depth stage-discharge curve of a prescribed elevation boundary, in TELEMAC's curves file format")),
            ("liquidboundariesfile", tr("Hydrographs (.liq) from a table…"),
             tr("Liquid boundaries file from a CSV or Excel table, columns in TELEMAC's boundary numbering")),
            ("initialconditions", tr("Initial state…"),
             tr("Initial conditions file: constant level, level along a channel axis, or depth above the bed")),
            ("steeringfile", tr("Steering file (.cas)…"),
             tr("TELEMAC-2D steering file for the mesh, its boundaries and the files above")),
            ("checkboundaries", tr("Check boundary conditions…"),
             tr("Liquid boundaries labelled with their TELEMAC numbers, checked against the .liq, curves and steering file")),
        )
        for key, label, tip in tools:
            b = QPushButton(label)
            b.setToolTip(tip)
            b.clicked.connect(lambda _c=False, k=key: self.run_tool(k))
            f.addRow(b)
            self.run_buttons[key] = b
        self.toolbox.addItem(p, tr("6 · Run files"))

    # ------------------------------------------------------------ snapping and editing
    def snapping_preset(self, quiet: bool = False):
        """QGIS snapping to the vertices and segments of the project's layers, at intersections, with topological editing;
        tracing, the Advanced Digitizing panel and hover tips on."""
        from qgis.core import Qgis, QgsSnappingConfig, QgsTolerance

        proj = QgsProject.instance()
        cfg = proj.snappingConfig()
        cfg.setEnabled(True)
        cfg.setMode(Qgis.SnappingMode.AdvancedConfiguration)
        cfg.setIntersectionSnapping(True)
        if hasattr(cfg, "setSelfSnapping"):
            cfg.setSelfSnapping(True)  # close a boundary line or an island on its own first vertex
        types = Qgis.SnappingTypes(Qgis.SnappingType.Vertex | Qgis.SnappingType.Segment)
        n = 0
        for key, lyr in self.layers.items():
            if lyr is None or key in ("refinement_points", "fixed_points"):
                continue
            s = QgsSnappingConfig.IndividualLayerSettings(True, types, 12.0, QgsTolerance.UnitType.Pixels, 0.0, 0.0)
            cfg.setIndividualLayerSettings(lyr, s)
            n += 1
        proj.setSnappingConfig(cfg)
        proj.setTopologicalEditing(True)
        extra = self._digitizing_aids()
        self.say(tr("Snapping on for %d layers: vertices and edges within 12 pixels, intersections, topological editing%s.")
                 % (n, (", " + ", ".join(extra)) if extra else ""))
        return n

    def _digitizing_aids(self):
        """Tracing (follow existing lines with T), the Advanced Digitizing panel and hover tips; what was turned on."""
        from qgis.PyQt.QtWidgets import QAction, QDockWidget

        done = []
        getter = getattr(self.iface, "mainWindow", None)
        win = getter() if callable(getter) else None
        if win is None:
            return done
        act = win.findChild(QAction, "mActionEnableTracing")
        if act is not None and act.isCheckable():
            if not act.isChecked():
                act.trigger()
            done.append(tr("tracing (T)"))
        dock = win.findChild(QDockWidget, "AdvancedDigitizingTools")
        if dock is not None:
            dock.show()
            done.append(tr("Advanced Digitizing panel"))
        tips = getattr(self.iface, "actionMapTips", None)
        act = tips() if callable(tips) else win.findChild(QAction, "mActionMapTips")
        if act is not None and act.isCheckable():
            if not act.isChecked():
                act.trigger()
            done.append(tr("hover tips"))
        return done

    def _free(self, path: str):
        """Remove the layers that read ``path`` before it is written again (Windows keeps loaded files locked)."""
        try:
            names = M.free_output([path])
        except Exception:  # noqa: BLE001 - writing will report it
            return
        if names:
            self.say(tr("Removed the layer %s from the map so its file can be written again; the new one is loaded.")
                     % ", ".join(names))

    def _current_mesh(self) -> Optional[str]:
        for cand in (self.last_output, self.output.filePath()):
            if cand and cand.lower().endswith(".slf") and os.path.exists(cand):
                return cand
        self.say(tr("Generate a TELEMAC mesh in step 4 first."), 1)
        return None

    def _open_tool(self, alg: str, key: str, suffix: str):
        mesh = self._current_mesh()
        params = {key: mesh, "OUTPUT": os.path.splitext(mesh)[0] + suffix} if mesh else {}
        try:
            import processing

            processing.execAlgorithmDialog("telemactoolkit:" + alg, params)
        except Exception as exc:
            self.say(str(exc), 2)

    def edit_mesh_start(self):
        """A copy of the mesh in edit mode for QGIS's Mesh Digitizing tools."""
        import shutil

        from qgis.core import QgsMeshLayer

        mesh = self._current_mesh()
        if not mesh:
            return None
        copy = os.path.splitext(mesh)[0] + "_editing.slf"
        shutil.copyfile(mesh, copy)
        lyr = QgsMeshLayer(copy, os.path.basename(copy), "mdal")
        if not lyr.isValid():
            self.say(tr("QGIS could not open the copy of the mesh."), 2)
            return None
        from ..qgis_bridge.crs import sidecar_crs

        crs = sidecar_crs(mesh)
        if crs is not None:
            lyr.setCrs(crs)
        QgsProject.instance().addMapLayer(lyr)
        canvas = self.iface.mapCanvas()
        ct = QgsCoordinateTransform(lyr.crs(), canvas.mapSettings().destinationCrs(), QgsProject.instance())
        if not lyr.startFrameEditing(ct):
            self.say(tr("QGIS cannot edit this mesh layer."), 2)
            return None
        try:
            self.iface.setActiveLayer(lyr)
        except AttributeError:
            pass
        self._editing = (lyr, copy, mesh)
        self.say(tr("Editing a copy of the mesh: use the Mesh Digitizing toolbar (View → Toolbars), then Rebuild after editing. "
                    "The bottom and other values are restored then; do not save the copy with another tool."))
        return lyr

    def edit_mesh_finish(self):
        """Save the edits, then rebuild the mesh with its values and .cli."""
        if not self._editing:
            self.say(tr("Start with Edit in QGIS."), 1)
            return None
        lyr, copy, mesh = self._editing
        try:
            if lyr.isEditable():
                ct = QgsCoordinateTransform(lyr.crs(), self.iface.mapCanvas().mapSettings().destinationCrs(), QgsProject.instance())
                if not lyr.commitFrameEditing(ct, False):
                    self.say(tr("QGIS could not save the edits."), 2)
                    return None
        except RuntimeError:  # the layer was removed
            pass
        import processing

        out = os.path.splitext(mesh)[0] + "_edited.slf"
        self._free(out)
        try:
            res = processing.run("telemactoolkit:rebuildedited", {"EDITED": copy, "ORIGINAL": mesh, "OUTPUT": out},
                                 feedback=_LogFeedback(self.say))
        except Exception as exc:
            self.say(str(exc), 2)
            return None
        try:
            QgsProject.instance().removeMapLayer(lyr.id())
        except RuntimeError:
            pass
        self._editing = None
        self.last_output = res["OUTPUT"]
        self._mesh_cache = None
        self.say(tr("Rebuilt: %s") % res["OUTPUT"])
        self._load_result(res["OUTPUT"])
        return res

    def cli_lines(self):
        mesh = self._current_mesh()
        if not mesh:
            return None
        import processing

        try:
            res = processing.run("telemactoolkit:clilines", {"MESH": mesh, "OUTPUT": "memory:"})
        except Exception as exc:
            self.say(str(exc), 2)
            return None
        lyr = res["OUTPUT"]
        lyr.setName(tr("Boundary lines of %s") % os.path.basename(mesh))
        QgsProject.instance().addMapLayer(lyr)
        self.say(tr("Edit the lines (snap to the mesh boundary), then run Boundary conditions from lines on the mesh."))
        return lyr

    def _load_result(self, path):
        try:
            from qgis.core import QgsMeshLayer

            from ..qgis_bridge import meshing as M
            from ..qgis_bridge.crs import sidecar_crs

            lyr = QgsMeshLayer(path, os.path.basename(path), "mdal")
            if lyr.isValid():
                crs = sidecar_crs(path)
                if crs is not None:
                    lyr.setCrs(crs)
                M.style_geometry_mesh(lyr)
                QgsProject.instance().addMapLayer(lyr)
        except Exception as exc:  # loading is a convenience
            self.say(str(exc), 1)

    def run_tool(self, alg: str):
        """Open a boundary or run-file tool with the last generated TELEMAC mesh filled in."""
        mesh = self.last_output if self.last_output and self.last_output.lower().endswith(".slf") else None
        if mesh is None and self.output.filePath().lower().endswith(".slf") and os.path.exists(self.output.filePath()):
            mesh = self.output.filePath()
        params = {"MESH": mesh} if mesh else {}
        if mesh:
            base = os.path.splitext(mesh)[0]
            params.update({
                "stagedischarge": {"OUTPUT": base + "_rating.txt"},
                "liquidboundariesfile": {"OUTPUT": base + ".liq"},
                "initialconditions": {"OUTPUT": base + "_init.slf"},
                "steeringfile": {"OUTPUT": os.path.join(os.path.dirname(mesh), "t2d_" + os.path.basename(base) + ".cas")},
            }.get(alg, {}))
            known = {"steeringfile": (("LIQ", ".liq"), ("STAGE", "_rating.txt"), ("INITIAL", "_init.slf")),
                     "checkboundaries": (("LIQ", ".liq"), ("STAGE", "_rating.txt"))}
            for key, suffix in known.get(alg, ()):
                if os.path.exists(base + suffix):
                    params[key] = base + suffix
        try:
            import processing

            processing.execAlgorithmDialog("telemactoolkit:" + alg, params)
        except Exception as exc:
            self.say(str(exc), 2)

    # ------------------------------------------------------------------ state
    def _set_enabled(self, on: bool):
        for i in range(self.toolbox.count()):
            self.toolbox.widget(i).setEnabled(on if i < 4 else True)

    def _structured_count(self) -> int:
        n = 0
        for key in ("channels", "grid_zones"):
            lyr = (self.layers or {}).get(key)
            try:
                n += int(lyr.featureCount()) if lyr is not None else 0
            except RuntimeError:
                pass
        return n

    def _mesh_type_changed(self, _i=None):
        """Mesh type sets the method (step 1) and the quadrilateral option (step 4)."""
        if getattr(self, "_syncing", False):
            return
        key = self.mesh_type.currentData()
        self._syncing = True
        try:
            if key == "quads":
                k = self.method.findData("gmsh_quads")
                if k >= 0:
                    self.method.setCurrentIndex(k)
                self.quads.setChecked(True)
            else:
                if self.method.currentData() == "gmsh_quads":
                    self.method.setCurrentIndex(self.method.findData("balanced"))
                self.quads.setChecked(key == "rectangles")
        finally:
            self._syncing = False
        self._mesh_type_note()

    def _sync_mesh_type(self):
        """The mesh type that the method and the quadrilateral option add up to."""
        if getattr(self, "_syncing", False):
            return
        key = ("quads" if self.method.currentData() == "gmsh_quads"
               else "rectangles" if self.quads.isChecked() else "triangles")
        self._syncing = True
        try:
            self.mesh_type.setCurrentIndex(max(0, self.mesh_type.findData(key)))
        finally:
            self._syncing = False
        self._mesh_type_note()

    def _mesh_type_note(self):
        key = self.mesh_type.currentData()
        fmt = FORMATS.get(self.fmt.currentData())
        notes = []
        if key in ("rectangles", "quads") and fmt is not None and not fmt.quads:
            notes.append(tr("%s takes triangles only: each rectangle is written as two triangles (shorter diagonal). "
                            "The rectangles are kept in 2DM, Gmsh, GiD and UGRID.") % fmt.label.split(" – ")[0].split(" - ")[0])
        if key == "rectangles" and self.layers and not self._structured_count():
            notes.append(tr("Nothing structured is drawn yet: draw Structured grid zones (any four-sided area, cells dx by "
                            "dy) or Structured channels (two bank lines) in step 2, or use the channel template."))
        if key == "triangles" and self.layers and self._structured_count():
            notes.append(tr("The structured channels and grid zones drawn are still built cell by cell (rectangles split "
                            "in two). Choose Rectangles where drawn to keep their rectangles in formats that accept them."))
        self.mesh_type_note.setText(" ".join(notes))

    def _format_changed(self):
        k = self.fmt.currentData()
        fmt = FORMATS[k]
        self.fmt_note.setText(fmt.note)
        self.quads.setEnabled(fmt.quads)
        self.store.setEnabled(k == "telemac")
        if hasattr(self, "double"):
            self.double.setEnabled(k == "telemac")
        if hasattr(self, "mesh_type_note"):
            self._mesh_type_note()
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
            cand = M.crs_with_code(crs)
            if cand is not None:
                ans = QMessageBox.question(self, tr("Coordinates"), tr(
                    "The CRS '%s' has no EPSG code (often the case of a CRS read from a DEM), so QGIS lists the project's "
                    "layers with an empty [ ]. It has the same projection as %s (%s); coordinates do not change.\n\n"
                    "Use %s?") % (crs.description() or tr("unnamed"), cand.authid(), cand.description(), cand.authid()))
                if ans == QMessageBox.StandardButton.Yes:
                    crs = cand
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
            M.refresh_map_tips(on_map)  # a project saved with an older tip, or with another theme
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
        if self.auto_snap.isChecked():
            try:
                self.snapping_preset()
            except Exception as exc:  # noqa: BLE001 - a convenience, never a reason to fail
                self.say(tr("Snapping not set: %s") % exc, 1)
        self._crs_check(crs)
        self.projectChanged.emit(path)

    def _crs_check(self, crs):
        """A project CRS without an EPSG code: offer the EPSG CRS with the same projection (once per project)."""
        from qgis.PyQt.QtWidgets import QMessageBox

        if crs is None or not crs.isValid() or crs.authid() or M.is_local_crs(crs):
            return
        self.crs_label.setText(tr("%s (no EPSG code)") % (crs.description() or tr("custom CRS")))
        cand = M.crs_with_code(crs)
        if cand is None:
            self.say(tr("The project's CRS has no EPSG code, so QGIS lists its layers with an empty [ ]. The coordinates are "
                        "right; to get a code, create the project with the EPSG CRS of your data (Project > Properties > CRS)."), 1)
            return
        key = "crs_asked/" + os.path.normcase(os.path.abspath(self.path))
        if _setting(key, "") == "1":
            self.say(tr("The project's CRS has no EPSG code; it has the same projection as %s (%s).") % (cand.authid(), cand.description()), 1)
            return
        _set_setting(key, "1")
        ans = QMessageBox.question(self, tr("Coordinates"), tr(
            "The project's CRS '%s' has no EPSG code, so QGIS lists its layers with an empty [ ].\n\nIt has the same projection "
            "as %s (%s). Give the project that code? The coordinates do not change.") % (
                crs.description() or tr("unnamed"), cand.authid(), cand.description()))
        if ans != QMessageBox.StandardButton.Yes:
            return
        try:
            M.assign_project_crs(self.path, cand)
        except Exception as exc:  # noqa: BLE001 - reported, the project keeps its CRS
            self.say(tr("The CRS could not be changed: %s") % exc, 2)
            return
        for lyr in self.layers.values():
            if lyr is not None:
                lyr.setCrs(cand)
        self.crs_label.setText(cand.authid())
        self.say(tr("Project CRS set to %s (%s).") % (cand.authid(), cand.description()))

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
        self.outline_tol.setValue(num("outline_tolerance", 0.0))
        self.max_nodes.setValue(int(num("max_nodes", 2000000)))
        self.terrain_tol.setValue(num("terrain_tolerance", 0.0))
        self.terrain_min.setValue(num("terrain_min_size", 0.0))
        k = self.method.findData(st.get("method", "balanced"))
        self.method.setCurrentIndex(k if k >= 0 and self.method.model().item(k).isEnabled() else 0)
        self.min_angle.setValue(min(30.0, max(20.0, num("min_angle", 30.0))))
        self.min_angle.setEnabled(self.method.currentData() == "quality")
        self.formula.setText(st.get("bottom_formula", ""))
        k = st.get("format", "telemac")
        if k in FORMAT_KEYS:
            self.fmt.setCurrentIndex(FORMAT_KEYS.index(k))
        self.quads.setChecked(st.get("quads") == "1")
        self.double.setChecked(st.get("double") == "1")
        self._sync_mesh_type()
        self.frame.setCurrentIndex(M.FRAME_MODES.index(st.get("frame", "keep")) if st.get("frame", "keep") in M.FRAME_MODES else 0)
        self.fx.setValue(num("frame_x", 0.0))
        self.fy.setValue(num("frame_y", 0.0))
        self.angle.setValue(num("frame_angle", 0.0))
        self.store.setChecked(st.get("store_offset", "1") != "0")
        if st.get("mesh_file"):
            self.output.setFilePath(st["mesh_file"])

    def settings(self) -> Dict[str, str]:
        return {"default_size": "%g" % self.size.value(), "growth": "%g" % self.growth.value(),
                "smoothing": str(self.smoothing.value()), "max_nodes": str(self.max_nodes.value()),
                "outline_tolerance": "%g" % self.outline_tol.value(),
                "terrain_tolerance": "%g" % self.terrain_tol.value(), "terrain_min_size": "%g" % self.terrain_min.value(),
                "method": self.method.currentData(), "min_angle": "%g" % self.min_angle.value(), "bottom_formula": self.formula.text().strip(),
                "format": self.fmt.currentData(), "quads": "1" if self.quads.isChecked() else "0",
                "double": "1" if self.double.isChecked() else "0", "mesh_type": self.mesh_type.currentData(),
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
        self._watch_layers()

    def _watch_layers(self):
        """Keep the object counts live: they follow every feature added or deleted, saved or not."""
        if getattr(self, "_count_timer", None) is None:
            self._count_timer = QTimer(self)
            self._count_timer.setSingleShot(True)
            self._count_timer.setInterval(150)
            self._count_timer.timeout.connect(self.update_counts)
            self._watched = {}
        for lyr in self.layers.values():
            try:
                lid = lyr.id()
            except RuntimeError:  # the layer was deleted
                continue
            if self._watched.get(lid) is lyr:
                continue
            self._watched[lid] = lyr
            for sig in (lyr.featureAdded, lyr.featuresDeleted, lyr.committedFeaturesAdded, lyr.committedFeaturesRemoved,
                        lyr.editingStopped, lyr.afterRollBack, lyr.dataChanged):
                try:
                    sig.connect(self._count_timer.start)
                except (AttributeError, TypeError):
                    pass

    def update_counts(self):
        """Object counts of the project layers, edits not yet saved included."""
        for r in range(self.table.rowCount()):
            it = self.table.item(r, 0)
            lyr = self.layers.get(it.data(Qt.ItemDataRole.UserRole)) if it is not None else None
            if lyr is None:
                continue
            try:
                n = lyr.featureCount()
            except RuntimeError:  # the layer was removed from the project
                continue
            cell = self.table.item(r, 1)
            if cell is None:
                self.table.setItem(r, 1, QTableWidgetItem(str(n)))
            elif cell.text() != str(n):
                cell.setText(str(n))
        if hasattr(self, "mesh_type_note"):
            self._mesh_type_note()

    def _dem_layer(self, layer_id: Optional[str] = None):
        """The DEM to suggest sizes from: ``layer_id``, the first raster ticked in step 3, the only raster on the map,
        or the active layer when it is a raster."""
        proj = QgsProject.instance()
        if layer_id:
            return proj.mapLayer(layer_id)
        for i in range(self.rasters.count()):
            it = self.rasters.item(i)
            if it.checkState() == Qt.CheckState.Checked:
                lyr = proj.mapLayer(it.data(Qt.ItemDataRole.UserRole))
                if lyr is not None:
                    return lyr
        rasters = [lyr for lyr in proj.mapLayers().values() if isinstance(lyr, QgsRasterLayer)]
        if len(rasters) == 1:
            return rasters[0]
        getter = getattr(self.iface, "activeLayer", None)
        act = getter() if callable(getter) else None
        return act if isinstance(act, QgsRasterLayer) else None

    def suggest_sizes(self, layer_id: Optional[str] = None, quiet: bool = False):
        """Sizes suggested from a DEM (mesher.suggest); with ``quiet`` only written to the log."""
        from qgis.PyQt.QtWidgets import QMessageBox

        dem = self._dem_layer(layer_id)
        if dem is None:
            self.say(tr("Add the DEM to the map and tick it in step 3 first."), 1)
            return None
        dom = self.layers.get("domain") if self.path else None
        crs = dom.crs() if dom is not None else dem.crs()
        try:
            sug = M.suggest_sizes_for(self.layers if self.path else None, dem, crs, QgsProject.instance(),
                                      min(300000, int(self.max_nodes.value()) or 300000))
        except Exception as exc:  # noqa: BLE001 - reported in the log
            self.say(tr("No suggestion from %s: %s") % (dem.name(), exc), 1)
            return None
        self.say(tr("Sizes suggested from the DEM %s:") % dem.name())
        for line in sug.report():
            self.say("  " + line)
        if quiet:
            self.say(tr("Suggest sizes from the DEM… (step 1) applies them."))
            return sug
        ans = QMessageBox.question(self, tr("Sizes from the DEM"), "\n".join(sug.report()) + "\n\n" + tr(
            "Apply the default size, terrain tolerance and smallest size to the project?"))
        if ans == QMessageBox.StandardButton.Yes:
            self.size.setValue(sug.default)
            self.terrain_tol.setValue(sug.tolerance)
            self.terrain_min.setValue(sug.smallest)
            if self.path:
                try:
                    M.write_settings(self.path, self.settings())
                except Exception as exc:  # noqa: BLE001
                    self.say(tr("Settings not saved: %s") % exc, 1)
            self.say(tr("Applied: default size %g m, terrain tolerance %g m, smallest size %g m.")
                     % (sug.default, sug.tolerance, sug.smallest))
        return sug

    def _raster_ticked(self, item):
        """A DEM ticked in step 3 while the project has no default size: its suggestion goes to the log (once)."""
        if item.checkState() != Qt.CheckState.Checked or self.size.value() > 0:
            return
        lid = item.data(Qt.ItemDataRole.UserRole)
        if lid in self._suggested:
            return
        self._suggested.add(lid)
        self.suggest_sizes(lid, quiet=True)

    def _fill_rasters(self):
        if not hasattr(self, "_suggested"):
            self._suggested = set()
            self.rasters.itemChanged.connect(self._raster_ticked)
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

    def check_drawing(self, commit: bool = True):
        if commit:
            self._commit_edits()
        try:
            b = M.build_mesh_input(self.layers, None, self.size.value() or M.auto_default_size(self.layers) or 1.0, self.growth.value(), self.smoothing.value(), None, self.max_nodes.value(),
                                   self.min_angle.value() if self.method.currentData() == "quality" else None,
                                   outline_tolerance=self.outline_tol.value() or None)
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
        problems = self._drawing_problems(b)
        self.say(tr("The drawing is ready to mesh.") if not (b.warnings or problems)
                 else tr("The drawing can be meshed; check the warnings."))
        return True

    def _drawing_problems(self, build) -> int:
        """Geometry checks against the element size; the problems as points on the map (layer Drawing check)."""
        from ..qgis_bridge import drawing as D

        self._check = None
        self.btn_fix.setEnabled(False)
        try:
            chk = D.check_drawing(self.layers, build.inp.size, build.crs)
        except Exception as exc:  # noqa: BLE001 - the checks never stop the meshing
            self.say(tr("Drawing check skipped: %s") % exc, 1)
            return 0
        proj = QgsProject.instance()
        for lyr in list(proj.mapLayers().values()):
            if lyr.customProperty("telemac_toolkit/drawing_check") == "1":
                proj.removeMapLayer(lyr.id())
        for line in chk.summary():
            self.say(line, 1 if chk.issues else 0)
        for i in chk.issues[:15]:
            self.say(i.text, 1)
        if len(chk.issues) > 15:
            self.say(tr("… and %d more: see the layer %s (one point per problem, with its message).")
                     % (len(chk.issues) - 15, D.CHECK_LAYER_NAME))
        if chk.issues:
            lyr = D.issues_layer(chk)
            proj.addMapLayer(lyr, False)
            proj.layerTreeRoot().insertLayer(0, lyr)
            self._check = (chk, build.inp.size)
            self.btn_fix.setEnabled(bool(chk.fixable))
        return len(chk.issues)

    def fix_drawing(self):
        from ..qgis_bridge import drawing as D

        if not getattr(self, "_check", None):
            return
        chk, size = self._check
        for line in D.fix_drawing(self.layers, chk, size):
            self.say(line)
        self.update_counts()
        self.say(tr("Checking again…"))
        self.check_drawing(commit=False)

    def flume_dialog(self):
        if not self.path:
            return
        try:
            import processing

            processing.execAlgorithmDialog("telemactoolkit:channeltemplate", {"PROJECT": self.path})
        except Exception as exc:
            self.say(str(exc), 2)
        self.set_project(self.path)

    def breaklines_dialog(self):
        if not self.path:
            return
        try:
            import processing

            params = {"PROJECT": self.path, "ADD": True}
            checked = [self.rasters.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.rasters.count())
                       if self.rasters.item(i).checkState() == Qt.CheckState.Checked]
            if checked:
                params["DEM"] = QgsProject.instance().mapLayer(checked[0])
            processing.execAlgorithmDialog("telemactoolkit:dembreaklines", params)
        except Exception as exc:
            self.say(str(exc), 2)
        self.set_project(self.path)

    def compare_dialog(self):
        if not self.path:
            return
        try:
            import processing

            params = {"PROJECT": self.path, "REPORT": os.path.splitext(self.path)[0] + "_methods.html"}
            checked = [self.rasters.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.rasters.count())
                       if self.rasters.item(i).checkState() == Qt.CheckState.Checked]
            if checked:
                params["DEM"] = QgsProject.instance().mapLayer(checked[0])
            processing.execAlgorithmDialog("telemactoolkit:comparemeshing", params)
        except Exception as exc:
            self.say(str(exc), 2)

    def river_dialog(self):
        if not self.path:
            return
        try:
            import processing

            processing.execAlgorithmDialog("telemactoolkit:riverreach", {"PROJECT": self.path})
        except Exception as exc:
            self.say(str(exc), 2)
        self.set_project(self.path)

    def _params(self) -> dict:
        checked = [self.rasters.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.rasters.count())
                   if self.rasters.item(i).checkState() == Qt.CheckState.Checked]
        params = {
            "PROJECT": self.path, "DEFAULT_SIZE": self.size.value(), "GROWTH": self.growth.value(),
            "SMOOTHING": self.smoothing.value(), "MAX_NODES": self.max_nodes.value(), "OUTLINE_TOLERANCE": self.outline_tol.value(),
            "METHOD": M.METHOD_KEYS.index(self.method.currentData()), "MIN_ANGLE": self.min_angle.value(),
            "TERRAIN_TOLERANCE": self.terrain_tol.value(), "TERRAIN_MIN_SIZE": self.terrain_min.value(), "RASTERS": [QgsProject.instance().mapLayer(i) for i in checked],
            "FORMULA": self.formula.text().strip(), "FORMAT": self.fmt.currentIndex(), "QUADS": self.quads.isChecked(), "DOUBLE": self.double.isChecked(),
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
        self._free(self.output.filePath())
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
