"""TELEMAC-2D steering file (.cas) template.

The numerical settings are those of the demo cases of this plugin that ran in TELEMAC-2D v9.1
(wetting and drying with positive depths: TREATMENT OF NEGATIVE DEPTHS = 2 and the settings it
requires; velocities advected with the NERD scheme, 14, which kept a river inflow stable on fine
meshes where the method of characteristics let it dry out). Everything specific to a case (files, boundaries, friction, times) comes from the
arguments; the rest is a sound starting point, to be reviewed for each model.
"""

from __future__ import annotations

import math
import os
from typing import List, Optional, Sequence

import numpy as np

__all__ = ["courant_time_step", "fitting_time_step", "steering_template", "FRICTION_LAWS"]

#: LAW OF BOTTOM FRICTION codes
FRICTION_LAWS = {"chezy": 2, "strickler": 3, "manning": 4, "nikuradse": 5}


def courant_time_step(x, y, ikle, depth: float, velocity: float = 1.0, courant: float = 1.0) -> float:
    """Time step for a given Courant number: smallest (element height) / (|u| + sqrt(g h)).

    The element height (twice the area over the longest edge) is the distance a wave crosses an element.
    """
    x, y, ikle = (np.asarray(v) for v in (x, y, ikle))
    P = np.stack([x[ikle], y[ikle]], axis=2)
    e = [np.hypot(*(P[:, (i + 1) % 3] - P[:, i]).T) for i in range(3)]
    lmax = np.maximum(np.maximum(e[0], e[1]), e[2])
    area = 0.5 * np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
                        - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
    hmin = float(np.min(2.0 * area / np.maximum(lmax, 1e-300)))
    c = abs(velocity) + math.sqrt(9.81 * max(depth, 1e-6))
    return courant * hmin / c


def fitting_time_step(dt_max: float, duration: float, decimals: int = 4) -> float:
    """The largest time step not above ``dt_max`` that divides ``duration`` exactly, with at most ``decimals`` decimals.

    TELEMAC runs round(DURATION / TIME STEP) steps; when the step does not divide the duration the last step
    can end past it, and a liquid boundaries file ending at the duration then stops the run with
    "OUT OF RANGE OF THE FILE OF LIQUID BOUNDARIES" (its tolerance is 1 ms). A step with few decimals that
    divides the duration also keeps the output periods whole numbers of steps.
    """
    if dt_max <= 0 or duration <= 0:
        return dt_max
    scale = 10 ** int(decimals)
    whole = round(duration * scale)
    k0 = max(1, int(math.ceil(duration / dt_max - 1e-9)))
    if abs(whole - duration * scale) < 1e-6:  # the duration itself has at most ``decimals`` decimals
        for k in range(k0, 4 * k0 + 1):
            if whole % k == 0:
                return (whole // k) / scale
    return duration / k0


def _rel(path: Optional[str], base: str) -> Optional[str]:
    if not path:
        return None
    try:
        r = os.path.relpath(path, base)
    except ValueError:  # another drive on Windows
        return path.replace("\\", "/")
    return r.replace("\\", "/")


def _list(values: Sequence[float]) -> str:
    return ";".join("%g" % v for v in values)


def _is_double(path: Optional[str]) -> bool:
    """True for a SERAFIND (double precision) file; False when it is single precision, missing or unreadable."""
    if not path or not os.path.exists(path):
        return False
    try:
        from .selafin import read_header

        return bool(read_header(path).is_double)
    except Exception:  # noqa: BLE001 - not a SELAFIN file: no format keyword
        return False


def steering_template(cas_path: str, geometry: str, boundary: str, n_liquid: int, results: Optional[str] = None,
                      liquid_file: Optional[str] = None, stage_file: Optional[str] = None,
                      stage_boundaries: Sequence[int] = (), prescribed_flowrates: Optional[Sequence[float]] = None,
                      prescribed_elevations: Optional[Sequence[float]] = None, initial_file: Optional[str] = None,
                      initial_level: Optional[float] = None, duration: float = 3600.0, time_step: float = 1.0,
                      output_period: float = 600.0, friction_law: str = "manning", friction: Optional[float] = 0.03,
                      friction_in_geometry: bool = False, title: str = "TELEMAC-2D model",
                      start: Optional[Sequence[int]] = None, notes: Sequence[str] = (),
                      kinds: Sequence[str] = ()) -> str:
    """The text of a TELEMAC-2D steering file; file names are written relative to the steering file.

    ``kinds``: the type of each liquid boundary in TELEMAC order (Q, H, QH, ...); discharge boundaries
    get a velocity profile in sqrt(h) (VELOCITY PROFILES = 4), which keeps the inflow stable where the
    depth varies across the section.
    """
    base = os.path.dirname(os.path.abspath(cas_path))
    law = FRICTION_LAWS[friction_law.lower()]
    out: List[str] = [
        "/ TELEMAC-2D steering file written by TELEMAC Toolkit for QGIS",
        "/ Review every value before running; the numerical settings below are a tested starting point.",
    ]
    out += ["/ " + n for n in notes]
    out += ["/", "/ ---- files", "TITLE = '%s'" % title.replace("'", " "),
            "GEOMETRY FILE = '%s'" % _rel(geometry, base),
            "BOUNDARY CONDITIONS FILE = '%s'" % _rel(boundary, base),
            "RESULTS FILE = '%s'" % _rel(results or os.path.splitext(cas_path)[0] + "_res.slf", base)]
    if _is_double(geometry):
        # a double-precision geometry (SERAFIND) must be declared; the results follow it
        out += ["GEOMETRY FILE FORMAT = 'SERAFIND'", "RESULTS FILE FORMAT = 'SERAFIND'"]
    if liquid_file:
        out.append("LIQUID BOUNDARIES FILE = '%s'" % _rel(liquid_file, base))
    if stage_file:
        out.append("STAGE-DISCHARGE CURVES FILE = '%s'" % _rel(stage_file, base))
    out += ["/", "/ ---- time"]
    if start:
        out += ["ORIGINAL DATE OF TIME = %d;%d;%d" % tuple(start[:3]), "ORIGINAL HOUR OF TIME = %d;%d;%d" % tuple(start[3:6])]
    out += ["TIME STEP = %.10g" % time_step, "DURATION = %.10g" % duration,
            "GRAPHIC PRINTOUT PERIOD = %d" % max(1, int(round(output_period / time_step))),
            "LISTING PRINTOUT PERIOD = %d" % max(1, int(round(output_period / time_step))),
            "VARIABLES FOR GRAPHIC PRINTOUTS = 'U,V,H,S,B,F,Q,M'", "MASS-BALANCE = YES"]
    out += ["/", "/ ---- initial state"]
    if initial_file:
        if _is_double(initial_file):
            out.append("PREVIOUS COMPUTATION FILE FORMAT = 'SERAFIND'")
        out += ["PREVIOUS COMPUTATION FILE = '%s'" % _rel(initial_file, base),
                "INITIAL TIME SET TO ZERO = YES"]
    elif initial_level is not None:
        out += ["INITIAL CONDITIONS = 'CONSTANT ELEVATION'", "INITIAL ELEVATION = %g" % initial_level]
    else:
        out += ["INITIAL CONDITIONS = 'ZERO DEPTH'"]
    if n_liquid:
        out += ["/", "/ ---- boundary conditions (%d liquid boundaries, numbered as TELEMAC numbers them)" % n_liquid]
        if prescribed_flowrates is not None:
            out.append("PRESCRIBED FLOWRATES = " + _list(prescribed_flowrates))
        if prescribed_elevations is not None:
            out.append("PRESCRIBED ELEVATIONS = " + _list(prescribed_elevations))
        if kinds:
            # discharge boundaries: velocity in sqrt(h), so shallow nodes at the banks do not run dry
            out.append("VELOCITY PROFILES = " + ";".join("4" if k == "Q" else "1" for k in kinds))
        if stage_file and stage_boundaries:
            sb = set(stage_boundaries)
            out.append("STAGE-DISCHARGE CURVES = " + ";".join("1" if i + 1 in sb else "0" for i in range(n_liquid)))
    out += ["/", "/ ---- friction and turbulence", "LAW OF BOTTOM FRICTION = %d" % law]
    if friction_in_geometry:
        out.append("/ FRICTION COEFFICIENT: read per node from the BOTTOM FRICTION variable of the geometry file")
        out.append("FRICTION COEFFICIENT = %g" % (friction if friction is not None else 0.03))
    else:
        out.append("FRICTION COEFFICIENT = %g" % (friction if friction is not None else 0.03))
    out += ["TURBULENCE MODEL = 1", "VELOCITY DIFFUSIVITY = 0.5",
            "/", "/ ---- numerics: wetting and drying with positive depths",
            "TIDAL FLATS = YES", "OPTION FOR THE TREATMENT OF TIDAL FLATS = 1", "TREATMENT OF NEGATIVE DEPTHS = 2",
            "MASS-LUMPING ON H = 1.", "CONTINUITY CORRECTION = YES", "SUPG OPTION = 0;0",
            "MAXIMUM NUMBER OF ITERATIONS FOR ADVECTION SCHEMES = 200",
            "/ velocities advected with the NERD scheme: with the method of characteristics (1) a river inflow",
            "/ could turn supercritical at a shallow bank node and dry out on fine meshes (small time steps)",
            "SCHEME FOR ADVECTION OF VELOCITIES = 14", "FREE SURFACE GRADIENT COMPATIBILITY = 0.9",
            "SOLVER = 1", "SOLVER ACCURACY = 1.E-6", "MAXIMUM NUMBER OF ITERATIONS FOR SOLVER = 200",
            "INFORMATION ABOUT SOLVER = NO"]
    return "\n".join(out) + "\n"
