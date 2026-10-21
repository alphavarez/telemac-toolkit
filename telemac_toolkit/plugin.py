"""TELEMAC Toolkit for QGIS: plugin entry point."""

from __future__ import annotations

import os
import contextlib
import traceback

from qgis.core import Qgis, QgsApplication, QgsCoordinateReferenceSystem, QgsMeshLayer, QgsProject, QgsSettings
from qgis.gui import QgsCustomDropHandler
from qgis.PyQt.QtCore import QCoreApplication, Qt
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QFileDialog, QInputDialog, QMenu

from .compat import QAction, exec_dialog
from .log import log_exception
from .gui.maptools import LiveCursorTool, ProbeTool, ProfileTool
from .gui.plot_widget import PALETTE
from .gui.mesh_panel import MeshPanel
from .gui.results_dock import ResultsDock
from .tt_processing.provider import TelemacProvider
from .qgis_bridge.results import ResultSet, file_extent, load_result
from .qgis_bridge.crs import crs_label, show_layer, suggest_crs, write_prj
from .qgis_bridge import styles
from .qgis_bridge.timesync import TimeSync
from .telemac_io import SelafinFile

PLUGIN_DIR = os.path.dirname(__file__)
SLF_EXT = (".slf", ".res", ".srf")
MENU = "&TELEMAC Toolkit"
MESH_MENU = "Mesh"


def tr(text):
    return QCoreApplication.translate("TelemacToolkit", text)


def icon(name):
    return QIcon(os.path.join(PLUGIN_DIR, "icons", name))


class SelafinDropHandler(QgsCustomDropHandler):
    """Dropping a .slf on QGIS opens it through the plugin (handles 3D files too)."""

    def __init__(self, plugin):
        super().__init__()
        self.plugin = plugin

    def handleFileDrop(self, file):
        if file.lower().endswith(SLF_EXT):
            self.plugin.open_path(file)
            return True
        return False


class TelemacToolkitPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.results = {}  # layer id -> ResultSet
        self.provider = None
        self.dock = None
        self.timesync = None
        self.actions = []
        self.toolbar = None
        self.drop_handler = None
        self.mesh_menu = None
        self.interactive = True  # False: never ask (scripts, tests); the best CRS guess is used

    # ------------------------------------------------------------------ setup
    def initProcessing(self):
        provider = TelemacProvider()
        # addProvider refuses (and deletes) a provider whose id is already registered
        self.provider = provider if QgsApplication.processingRegistry().addProvider(provider) else None

    def initGui(self):
        self.initProcessing()
        canvas = self.iface.mapCanvas()
        self.timesync = TimeSync(canvas)
        if self.timesync.navigation() is None:
            # some builds return the base QgsTemporalController: use the panel's navigation object
            try:
                from qgis.gui import QgsTemporalControllerWidget

                w = self.iface.mainWindow().findChild(QgsTemporalControllerWidget)
                if w is not None:
                    self.timesync.nav_override = w.temporalController()
                    self.timesync.nav_override.updateTemporalRange.connect(self.timesync._on_range)
            except (ImportError, AttributeError, RuntimeError) as exc:
                log_exception("Temporal Controller panel not found; the clock follows the canvas only", exc)
        self.dock = ResultsDock(self, self.iface.mainWindow())
        self.iface.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.dock)
        self.dock.hide()
        self.timesync.frameChanged.connect(self.dock.on_frame)

        get = lambda: self.dock.rs  # noqa: E731
        self.cursor_tool = LiveCursorTool(canvas, get)
        self.cursor_tool.moved.connect(self.dock.update_cursor)
        self.cursor_tool.left_mesh.connect(self.dock.cursor_left)
        self.probe_tool = ProbeTool(canvas, get, PALETTE)
        self.probe_tool.added.connect(self.dock.add_probe)
        self.probe_tool.cleared.connect(self.dock.clear_probes)
        self.profile_tool = ProfileTool(canvas, get)
        self.profile_tool.finished.connect(self.dock.set_profile_line)
        canvas.mapToolSet.connect(self._on_map_tool_set)

        self.toolbar = self.iface.addToolBar("TELEMAC Toolkit")
        self.toolbar.setObjectName("TelemacToolkitToolbar")
        self.a_open = self._action("open.svg", tr("Open TELEMAC result…"), self.open_result_dialog)
        self.a_dock = self._action("results.svg", tr("TELEMAC Results panel"), self.toggle_dock, checkable=True)
        self.a_cursor = self._action("cursor.svg", tr("Live cursor"), self.toggle_live_cursor, checkable=True)
        self.a_probe = self._action("probe.svg", tr("Probe time series"), self.toggle_probe, checkable=True)
        self.a_profile = self._action("profile.svg", tr("Draw profile"), self.toggle_profile, checkable=True)
        self.a_extract = self._action("extract.svg", tr("Extract 2D from a 3D result…"),
                                      lambda: self._run_alg("telemactoolkit:extract2d"), toolbar=False)
        self.a_info = self._action("telemac.svg", tr("SELAFIN file information…"),
                                   lambda: self._run_alg("telemactoolkit:selafininfo"), toolbar=False)
        self.a_calc = self._action("calc.svg", tr("Calculator (gradients, vorticity, time maxima)…"),
                                   lambda: self.open_calculator(), toolbar=False)
        self.a_integrals = self._action("integral.svg", tr("Integrals and fluxes (volume, discharge through lines)…"),
                                        lambda: self._run_alg("telemactoolkit:integrals", self._current_input()), toolbar=False)
        # mesh: workflow panel + every tool in the Mesh submenu
        self.toolbar.addSeparator()
        self.mesh_panel = MeshPanel(self.iface, self.iface.mainWindow())
        self.iface.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.mesh_panel)
        self.mesh_panel.hide()
        self.a_meshpanel = self._action("mesh.svg", tr("Mesh workflow panel"), self.toggle_mesh_panel, checkable=True)
        self.mesh_panel.visibilityChanged.connect(self.a_meshpanel.setChecked)
        for ico, text, alg, tb in (
                ("meshproject.svg", tr("New mesh project…"), "createmeshproject", True),
                ("meshproject.svg", tr("Channel or flume template (straight, bend, meanders)…"), "channeltemplate", False),
                ("mesh.svg", tr("Generate mesh…"), "generatemesh", False),
                ("extract.svg", tr("Export or convert mesh…"), "exportmesh", False),
                ("extract.svg", tr("Extract part of a mesh…"), "cutmesh", False),
                ("mesh.svg", tr("Structured channel mesh…"), "channelmesh", False),
                ("mesh.svg", tr("Regular grid mesh…"), "regulargrid", False),
                ("extract.svg", tr("Edit node values of a mesh (bottom, friction)…"), "interpolatebottom", False),
                ("boundary.svg", tr("Boundary conditions from lines…"), "boundaryconditions", False),
                ("telemac.svg", tr("Mesh quality…"), "meshquality", False)):
            self._action(ico, text, (lambda a=alg: self._run_alg("telemactoolkit:" + a)), toolbar=tb, submenu=MESH_MENU)
        self.dock.visibilityChanged.connect(self.a_dock.setChecked)

        self.drop_handler = SelafinDropHandler(self)
        self.iface.registerCustomDropHandler(self.drop_handler)

        prj = QgsProject.instance()
        prj.layersWillBeRemoved.connect(self._on_layers_removed)
        prj.readProject.connect(self._adopt_project_layers)
        prj.layersAdded.connect(self._adopt_project_layers)
        prj.cleared.connect(self._on_project_cleared)
        self._adopt_project_layers()

    def _action(self, ico, text, slot, checkable=False, toolbar=True, submenu=None):
        a = QAction(icon(ico), text, self.iface.mainWindow())
        a.setCheckable(checkable)
        a.triggered.connect(lambda *_: slot())
        if submenu is None:
            self.iface.addPluginToMenu(MENU, a)
        else:
            if self.mesh_menu is None:
                self.mesh_menu = QMenu(tr(submenu), self.iface.mainWindow())
                self.mesh_menu.setIcon(icon("mesh.svg"))
                self.iface.addPluginToMenu(MENU, self.mesh_menu.menuAction())
            self.mesh_menu.addAction(a)
        if toolbar:
            self.toolbar.addAction(a)
        self.actions.append(a)
        return a

    def unload(self):
        canvas = self.iface.mapCanvas()
        for tool in (self.cursor_tool, self.probe_tool, self.profile_tool):
            if canvas.mapTool() is tool:
                canvas.unsetMapTool(tool)
        self.probe_tool.clear()
        self.profile_tool.remove()
        with contextlib.suppress(TypeError, RuntimeError):  # not connected
            canvas.mapToolSet.disconnect(self._on_map_tool_set)
        prj = QgsProject.instance()
        for sig, slot in ((prj.layersWillBeRemoved, self._on_layers_removed), (prj.readProject, self._adopt_project_layers),
                          (prj.layersAdded, self._adopt_project_layers),
                          (prj.cleared, self._on_project_cleared)):
            with contextlib.suppress(TypeError, RuntimeError):  # not connected
                sig.disconnect(slot)
        if self.timesync:
            self.timesync.disconnect_all()
        if self.drop_handler:
            self.iface.unregisterCustomDropHandler(self.drop_handler)
        for a in self.actions:
            self.iface.removePluginMenu(MENU, a)
        if self.mesh_menu is not None:
            self.iface.removePluginMenu(MENU, self.mesh_menu.menuAction())
            self.mesh_menu.deleteLater()
            self.mesh_menu = None
        self.actions = []
        if self.toolbar:
            self.toolbar.deleteLater()
        if self.dock:
            self.iface.removeDockWidget(self.dock)
            self.dock.deleteLater()
        if getattr(self, "mesh_panel", None) is not None:
            self.mesh_panel.cleanup()
            self.iface.removeDockWidget(self.mesh_panel)
            self.mesh_panel.deleteLater()
            self.mesh_panel = None
        if self.provider:
            try:
                QgsApplication.processingRegistry().removeProvider(self.provider)
            except RuntimeError:
                pass
            self.provider = None
        for rs in self.results.values():
            rs.close()
        self.results = {}

    # ------------------------------------------------------------------ opening
    def open_result_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self.iface.mainWindow(), tr("Open TELEMAC result"), "",
            tr("TELEMAC SELAFIN (*.slf *.res *.srf *.SLF);;All files (*)"))
        if path:
            self.open_path(path)

    def open_path(self, path, plane=None):
        bar = self.iface.messageBar()
        try:
            with SelafinFile(path) as probe:
                is3d = probe.header.is_3d
                nplan = probe.header.nplan
        except Exception as exc:
            bar.pushMessage(tr("TELEMAC Toolkit"), tr("Cannot read %s: %s") % (os.path.basename(path), exc), Qgis.MessageLevel.Critical, 10)
            return None
        if is3d and plane is None:
            choices = [tr("Depth average"), tr("Bottom plane"), tr("Surface plane")] + [tr("Plane %d") % k for k in range(nplan)]
            item, ok = QInputDialog.getItem(
                self.iface.mainWindow(), tr("3D result"),
                tr("QGIS cannot display 3D SELAFIN files directly.\nA 2D file will be written next to the result.\nShow:"),
                choices, 0, False)
            if not ok:
                return None
            k = choices.index(item)
            plane = ["average", "bottom", "surface"][k] if k < 3 else k - 3
        prj = QgsProject.instance()
        first_layer = len(prj.mapLayers()) == 0
        crs = self.choose_crs(path)
        if crs is False:
            return None
        try:
            rs = load_result(path, plane=plane if plane is not None else "average", crs=crs)
        except Exception as exc:
            bar.pushMessage(tr("TELEMAC Toolkit"), str(exc), Qgis.MessageLevel.Critical, 15)
            QgsApplication.messageLog().logMessage(traceback.format_exc(), "TELEMAC Toolkit", Qgis.MessageLevel.Critical)
            return None
        self.results[rs.layer_id] = rs
        self.dock.show()
        self.dock.raise_()
        self.dock.refresh_result_list(select=rs.layer_id)
        done = show_layer(rs.layer, self.iface.mapCanvas(), prj, first=first_layer)
        msgs = rs.warnings()
        text = tr("Opened %s (%d frames) in %s.") % (os.path.basename(path), rs.nframes, crs_label(rs.layer.crs()))
        if done and not first_layer:
            text += " " + tr("The map was adjusted: %s.") % done
        bar.pushMessage(tr("TELEMAC Toolkit"), text + ("  " + " ".join(msgs) if msgs else ""),
                        Qgis.MessageLevel.Warning if msgs else Qgis.MessageLevel.Info, 8)
        return rs

    def choose_crs(self, path, force_ask=False, current=None):
        """CRS for a SELAFIN file (they carry none).

        A .prj next to the file wins; then the project CRS, when the project
        already has a map and the coordinates lie in that CRS's area of use.
        Otherwise the user chooses between local laboratory coordinates and a
        map CRS (the suggestion follows the coordinates), and can save the
        choice as a .prj. Returns a CRS, None (leave unset) or False (cancelled).
        """
        from .gui.crs_dialog import CrsDialog

        try:
            ext = file_extent(path)
        except Exception:
            return None
        settings = QgsSettings()
        last = QgsCoordinateReferenceSystem(settings.value("telemac_toolkit/last_crs", "", type=str))
        choice = suggest_crs(path, ext, QgsProject.instance(), last)
        if (not choice.ask and not force_ask) or not self.interactive:
            return choice.crs
        dlg = CrsDialog(self.iface.mainWindow(), path, ext, current if current is not None else choice.crs, choice.reason)
        if not exec_dialog(dlg):
            return False
        crs = dlg.crs()
        if not crs.isValid():
            return None
        settings.setValue("telemac_toolkit/last_crs", crs.authid() or crs.toWkt())
        if dlg.save.isChecked():
            try:
                write_prj(path, crs)
            except OSError as exc:
                self.iface.messageBar().pushMessage(tr("TELEMAC Toolkit"), tr("Could not write the .prj file: %s") % exc,
                                                    Qgis.MessageLevel.Warning, 6)
        return crs

    def zoom_to_result(self, rs=None, quiet=False):
        rs = rs or (self.dock.rs if self.dock else None)
        if rs is None:
            return
        done = show_layer(rs.layer, self.iface.mapCanvas(), QgsProject.instance(), force_zoom=True)
        if done and not quiet and "CRS" in done:
            self.iface.messageBar().pushMessage(tr("TELEMAC Toolkit"), tr("Map adjusted: %s.") % done, Qgis.MessageLevel.Info, 5)

    def change_result_crs(self, rs=None):
        rs = rs or (self.dock.rs if self.dock else None)
        if rs is None:
            return
        crs = self.choose_crs(rs.source_path, force_ask=True, current=rs.layer.crs())
        if crs in (False, None):
            return
        rs.layer.setCrs(crs)
        others = [lid for lid in QgsProject.instance().mapLayers() if lid != rs.layer_id]
        show_layer(rs.layer, self.iface.mapCanvas(), QgsProject.instance(), first=not others, force_zoom=True)
        rs.layer.triggerRepaint()
        self.dock._update_info()

    def open_calculator(self, rs=None):
        """The Calculator dialog for the current result, with its variables listed in the expression box."""
        rs = rs or (self.dock.rs if self.dock else None)
        params = {}
        if rs is not None:
            from .tt_processing.calc_algorithms import template_for

            base = os.path.splitext(rs.path)[0]  # the 2D file for a 3D result
            params = {"INPUT": rs.path, "EXPRESSIONS": template_for(rs.path), "OUTPUT": base + "_calc.slf"}
            if rs.layer.crs().isValid():
                params["MESH_CRS"] = rs.layer.crs()
        self._run_alg("telemactoolkit:calculator", params)

    def _adopt_project_layers(self, *_):
        """Re-attach SELAFIN mesh layers of a reopened project to the plugin."""
        for lyr in QgsProject.instance().mapLayers().values():
            if not isinstance(lyr, QgsMeshLayer) or lyr.id() in self.results:
                continue
            src = lyr.source().strip('"')
            if lyr.providerType() != "mdal" or not src.lower().endswith(SLF_EXT) or not os.path.exists(src):
                continue
            try:
                rs = ResultSet(src, src, lyr)
                rs.group_of_var, rs.group_of_pair = styles.group_mapping(rs)
                self.results[lyr.id()] = rs
            except Exception as exc:  # a broken layer must not stop the others
                log_exception("Could not attach %s to the Results panel" % src, exc)
                continue
        if self.dock:
            self.dock.refresh_result_list()

    def _on_layers_removed(self, ids):
        changed = False
        for lid in ids:
            rs = self.results.pop(lid, None)
            if rs is not None:
                if self.dock.rs is rs:
                    self.dock.set_result(None)
                rs.close()
                changed = True
        if changed:
            self.dock.refresh_result_list()

    def _on_project_cleared(self):
        for rs in self.results.values():
            rs.close()
        self.results = {}
        if self.dock:
            self.dock.set_result(None)
            self.dock.refresh_result_list()

    # ------------------------------------------------------------------ tools
    def toggle_dock(self):
        self.dock.setVisible(not self.dock.isVisible())

    def toggle_mesh_panel(self):
        show = not self.mesh_panel.isVisible()
        self.mesh_panel.setVisible(show)
        if show:
            self.mesh_panel.autodetect()

    def _toggle_tool(self, tool, require_result=True):
        canvas = self.iface.mapCanvas()
        if canvas.mapTool() is tool:
            canvas.unsetMapTool(tool)
            return
        if require_result and self.dock.rs is None:
            self.iface.messageBar().pushMessage(tr("TELEMAC Toolkit"), tr("Open a TELEMAC result first."), Qgis.MessageLevel.Info, 4)
            self._on_map_tool_set(canvas.mapTool(), None)
            return
        self.dock.show()
        canvas.setMapTool(tool)

    def toggle_live_cursor(self):
        self._toggle_tool(self.cursor_tool)
        self.dock.tabs.setCurrentIndex(0)

    def toggle_probe(self):
        self._toggle_tool(self.probe_tool)
        self.dock.tabs.setCurrentIndex(1)

    def toggle_profile(self):
        self._toggle_tool(self.profile_tool)
        self.dock.tabs.setCurrentIndex(2)

    def clear_probe_markers(self):
        if getattr(self, "probe_tool", None):
            self.probe_tool.clear()

    def _on_map_tool_set(self, new, old=None):
        pairs = ((self.cursor_tool, self.a_cursor, self.dock.btn_cursor),
                 (self.probe_tool, self.a_probe, self.dock.btn_probe),
                 (self.profile_tool, self.a_profile, self.dock.btn_profile))
        for tool, action, button in pairs:
            on = new is tool
            action.setChecked(on)
            button.setChecked(on)

    def _current_input(self):
        rs = self.dock.rs if self.dock else None
        return {"INPUT": rs.path} if rs is not None else {}

    def _run_alg(self, alg_id, params=None):
        try:
            import processing

            processing.execAlgorithmDialog(alg_id, params or {})
        except Exception as exc:
            self.iface.messageBar().pushMessage(tr("TELEMAC Toolkit"), str(exc), Qgis.MessageLevel.Critical, 8)
