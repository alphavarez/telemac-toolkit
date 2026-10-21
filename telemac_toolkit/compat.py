"""Small helpers so one codebase runs on QGIS 3.44 (Qt5) and QGIS 4.x (Qt6)."""

from qgis.PyQt import QtGui, QtWidgets

try:  # Qt6 moved QAction to QtGui; QGIS's shim may keep both
    QAction = QtWidgets.QAction
except AttributeError:  # pragma: no cover
    QAction = QtGui.QAction


def exec_dialog(widget):
    """QDialog/QMenu.exec() on Qt6, exec_() on old Qt5 builds."""
    fn = getattr(widget, "exec", None) or getattr(widget, "exec_")
    return fn()


def menu_exec(menu, pos):
    fn = getattr(menu, "exec", None) or getattr(menu, "exec_")
    return fn(pos)
