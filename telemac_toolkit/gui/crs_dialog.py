"""Ask what the coordinates of a TELEMAC file are (SELAFIN stores no CRS)."""

from __future__ import annotations

import os

from qgis.core import QgsCoordinateReferenceSystem
from qgis.gui import QgsProjectionSelectionWidget
from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QRadioButton,
    QVBoxLayout,
)

from ..qgis_bridge.crs import Extent, is_local_crs, local_crs, plausible


def tr(text):
    return QCoreApplication.translate("TelemacToolkit", text)


def _num(v: float) -> str:
    return ("%.3f" % v).rstrip("0").rstrip(".") if abs(v) < 1e6 else "%.0f" % v


class CrsDialog(QDialog):
    """Local laboratory grid or a map CRS, with a check against the coordinates."""

    def __init__(self, parent, path: str, ext: Extent, suggestion=None, reason: str = ""):
        super().__init__(parent)
        self.ext = ext
        self.setWindowTitle(tr("Coordinates of %s") % os.path.basename(path))
        lay = QVBoxLayout(self)
        intro = QLabel(tr("TELEMAC files do not store a coordinate system. In this file x runs from %s to %s and "
                          "y from %s to %s. What are these coordinates?") % (_num(ext.xmin), _num(ext.xmax),
                                                                             _num(ext.ymin), _num(ext.ymax)))
        intro.setWordWrap(True)
        lay.addWidget(intro)
        self.r_local = QRadioButton(tr("Local coordinates: laboratory flume or schematic channel, not a place on Earth"))
        self.r_map = QRadioButton(tr("Map coordinates in this CRS:"))
        grp = QButtonGroup(self)
        grp.addButton(self.r_local)
        grp.addButton(self.r_map)
        lay.addWidget(self.r_local)
        lay.addWidget(self.r_map)
        self.sel = QgsProjectionSelectionWidget(self)
        lay.addWidget(self.sel)
        self.check = QLabel()
        self.check.setWordWrap(True)
        lay.addWidget(self.check)
        self.save = QCheckBox(tr("Save the choice next to the file (%s), so it opens directly next time")
                              % (os.path.splitext(os.path.basename(path))[0] + ".prj"))
        self.save.setChecked(os.access(os.path.dirname(os.path.abspath(path)) or ".", os.W_OK))
        lay.addWidget(self.save)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

        if suggestion is not None and suggestion.isValid() and not is_local_crs(suggestion):
            self.sel.setCrs(suggestion)
            self.r_map.setChecked(True)
        else:
            from qgis.core import QgsProject

            pcrs = QgsProject.instance().crs()
            self.sel.setCrs(pcrs if pcrs.isValid() and not is_local_crs(pcrs) else QgsCoordinateReferenceSystem("EPSG:4326"))
            self.r_local.setChecked(suggestion is not None or ext.small)
            self.r_map.setChecked(not self.r_local.isChecked())
        self.r_local.toggled.connect(self._update)
        self.sel.crsChanged.connect(self._update)
        self._update()
        self.resize(520, self.sizeHint().height())

    def crs(self) -> QgsCoordinateReferenceSystem:
        return local_crs() if self.r_local.isChecked() else self.sel.crs()

    def _update(self, *_):
        self.sel.setEnabled(self.r_map.isChecked())
        crs = self.crs()
        if not crs.isValid():
            txt, color = tr("Choose a CRS."), "#b35900"
        elif plausible(crs, self.ext):
            txt, color = tr("✓ The coordinates fit this choice."), "#2e7d32"
        elif is_local_crs(crs):
            txt, color = tr("⚠ The coordinates are far from (0, 0): they look like map coordinates."), "#b35900"
        else:
            txt, color = tr("⚠ The coordinates fall outside the area of use of this CRS: the mesh would be drawn in "
                            "the wrong place (or not at all)."), "#b35900"
        self.check.setText("<span style='color:%s'>%s</span>" % (color, txt))
