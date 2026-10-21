"""One clock: keep the plugin's frame in step with QGIS's Temporal Controller."""

from __future__ import annotations

import numpy as np
from qgis.core import Qgis, QgsDateTimeRange, QgsInterval
from qgis.PyQt.QtCore import QObject, pyqtSignal


def _matching_closest_from_start():
    """Enum value for 'closest dataset from start of range' on QGIS 3.44 and 4.x."""
    for owner, name in (
        ("QgsMeshDataProviderTemporalCapabilities", "FindClosestDatasetFromStartRangeTime"),
    ):
        try:
            import qgis.core as qc

            cls = getattr(qc, owner)
            enum = getattr(cls, "MatchingTemporalDatasetMethod", cls)
            return getattr(enum, name)
        except AttributeError:
            continue
    try:
        return Qgis.MeshTemporalDatasetMatchingMethod.FindClosestDatasetFromStartRangeTime  # pragma: no cover
    except AttributeError:  # pragma: no cover
        return None


class TimeSync(QObject):
    """Bridges a ResultSet's file frames and the canvas temporal controller.

    ``set_frame(i)`` moves the QGIS clock; when the clock moves (play button,
    slider, another plugin), ``frameChanged(i)`` fires with the nearest file
    frame, so every view follows the same time.
    """

    frameChanged = pyqtSignal(int)

    def __init__(self, canvas=None):
        super().__init__()
        self.canvas = canvas
        self.rs = None
        self.frame = 0
        self._t0_ms = 0
        self._dt = 1.0
        self._ref = None
        ctrl = self.controller()
        if ctrl is not None:
            ctrl.updateTemporalRange.connect(self._on_range)

    nav_override = None  # set by the plugin when the canvas returns the base class

    def controller(self):
        if self.nav_override is not None:
            return self.nav_override
        if self.canvas is None:
            return None
        try:
            return self.canvas.temporalController()
        except Exception:  # pragma: no cover
            return None

    def disconnect_all(self):
        ctrl = self.controller()
        if ctrl is not None:
            try:
                ctrl.updateTemporalRange.disconnect(self._on_range)
            except Exception:
                pass
        self.rs = None

    def navigation(self):
        ctrl = self.controller()
        return ctrl if ctrl is not None and hasattr(ctrl, "setCurrentFrameNumber") else None

    # ------------------------------------------------------------------
    def bind(self, rs, configure: bool = True):
        self.rs = rs
        self.frame = min(self.frame, max(rs.nframes - 1, 0)) if rs else 0
        if rs is None:
            return
        tp = rs.layer.temporalProperties()
        m = _matching_closest_from_start()
        if m is not None:
            try:
                tp.setMatchingMethod(m)
            except Exception:  # pragma: no cover
                pass
        ref = tp.referenceTime()
        self._ref = ref if ref.isValid() else rs.reference_datetime()
        t = rs.times
        steps = np.diff(t)
        steps = steps[steps > 0]
        self._dt = float(np.median(steps)) if steps.size else 1.0
        self._t0_ms = int(round(float(t[0]) * 1000)) if t.size else 0
        if configure:
            self.configure_controller()
        self.set_frame(self.frame)

    def configure_controller(self):
        """Set the QGIS clock to this result's time span and output step."""
        nav = self.navigation()
        rs = self.rs
        if nav is None or rs is None or rs.nframes == 0:
            return
        t = rs.times
        begin = self._ref.addMSecs(int(round(float(t[0]) * 1000)))
        # extend by one step so the last output time is a frame of its own
        end = self._ref.addMSecs(int(round((float(t[-1]) + self._dt) * 1000)))
        try:
            nav.setNavigationMode(Qgis.TemporalNavigationMode.Animated)
        except AttributeError:  # pragma: no cover - very old API
            nav.setNavigationMode(nav.Animated)
        nav.setTemporalExtents(QgsDateTimeRange(begin, end))
        nav.setFrameDuration(QgsInterval(self._dt))

    # ------------------------------------------------------------------
    def set_frame(self, i: int):
        rs = self.rs
        if rs is None or rs.nframes == 0:
            return
        i = int(np.clip(i, 0, rs.nframes - 1))
        nav = self.navigation()
        if nav is not None and nav.navigationMode() != Qgis.TemporalNavigationMode.Animated:
            # the user switched the Temporal Controller off or to a fixed range:
            # turn animation back on so the map follows the panel
            self.configure_controller()
        if nav is not None and nav.navigationMode() == Qgis.TemporalNavigationMode.Animated:
            n = int(round((float(rs.times[i]) * 1000 - self._t0_ms) / (self._dt * 1000)))
            if n != nav.currentFrameNumber():
                nav.setCurrentFrameNumber(n)  # -> updateTemporalRange -> _on_range
                if self.frame == i:
                    return
            if self.frame != i:
                self.frame = i
                self.frameChanged.emit(i)
            return
        if i != self.frame:
            self.frame = i
        self.frameChanged.emit(i)

    def frame_for_range(self, rng) -> int:
        rs = self.rs
        secs = self._ref.msecsTo(rng.begin()) / 1000.0
        t = rs.times
        return int(np.argmin(np.abs(t - secs)))

    def _on_range(self, rng):
        if self.rs is None or self.rs.nframes == 0 or self._ref is None:
            return
        try:
            if not rng.begin().isValid():
                return
        except Exception:
            return
        i = self.frame_for_range(rng)
        if i != self.frame:
            self.frame = i
            self.frameChanged.emit(i)
