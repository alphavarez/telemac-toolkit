"""Messages for problems the plugin can work around: written to Python logging and, inside QGIS, to the
"TELEMAC Toolkit" tab of the Log Messages panel, so nothing is swallowed silently."""

from __future__ import annotations

import importlib.util
import logging
import traceback

LOG = logging.getLogger("telemac_toolkit")
_HAS_QGIS = importlib.util.find_spec("qgis") is not None


def log_exception(where: str, exc: BaseException, level: int = logging.WARNING):
    """Record an exception that was handled (the operation went on without it)."""
    text = "%s: %s: %s" % (where, type(exc).__name__, exc)
    LOG.log(level, text)
    LOG.debug("".join(traceback.format_exception(type(exc), exc, exc.__traceback__)))
    if _HAS_QGIS:
        from qgis.core import Qgis, QgsApplication

        if QgsApplication.instance() is not None:
            QgsApplication.messageLog().logMessage(text, "TELEMAC Toolkit", Qgis.MessageLevel.Info)
