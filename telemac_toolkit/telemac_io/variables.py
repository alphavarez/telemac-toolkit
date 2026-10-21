"""TELEMAC variable catalogue: English and French names, units, vector pairs.

TELEMAC writes variable names in the language chosen in the steering file
(``LANGUAGE = 1`` French, ``2`` English). The catalogue maps both to one
canonical id so the plugin can style, pair and compute regardless of language.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

__all__ = ["VarInfo", "CATALOGUE", "canonical", "lookup", "vector_pairs", "display_name", "STYLE_DEPTH",
           "STYLE_VELOCITY", "STYLE_ELEVATION", "STYLE_BOTTOM", "STYLE_GENERIC", "STYLE_FROUDE", "STYLE_TIME"]

STYLE_DEPTH = "depth"
STYLE_VELOCITY = "velocity"
STYLE_ELEVATION = "elevation"
STYLE_BOTTOM = "bottom"
STYLE_FROUDE = "froude"
STYLE_TIME = "time"
STYLE_GENERIC = "generic"


@dataclass(frozen=True)
class VarInfo:
    cid: str            # canonical id
    en: Tuple[str, ...]  # English names as written by TELEMAC (upper case)
    fr: Tuple[str, ...]  # French names
    unit: str
    label: str           # readable label
    style: str = STYLE_GENERIC
    symbol: str = ""     # short symbol for the calculator (H, U, V, S, B ...)


_V = VarInfo
CATALOGUE: Tuple[VarInfo, ...] = (
    # --- TELEMAC-2D ---------------------------------------------------------
    _V("velocity_u", ("VELOCITY U",), ("VITESSE U",), "M/S", "Velocity U", STYLE_VELOCITY, "U"),
    _V("velocity_v", ("VELOCITY V",), ("VITESSE V",), "M/S", "Velocity V", STYLE_VELOCITY, "V"),
    _V("velocity_w", ("VELOCITY W",), ("VITESSE W",), "M/S", "Velocity W", STYLE_VELOCITY, "W"),
    _V("celerity", ("CELERITY",), ("CELERITE",), "M/S", "Wave celerity", STYLE_GENERIC, "C"),
    _V("water_depth", ("WATER DEPTH",), ("HAUTEUR D'EAU", "HAUTEUR D EAU"), "M", "Water depth", STYLE_DEPTH, "H"),
    _V("free_surface", ("FREE SURFACE",), ("SURFACE LIBRE",), "M", "Free surface", STYLE_ELEVATION, "S"),
    _V("bottom", ("BOTTOM",), ("FOND",), "M", "Bottom", STYLE_BOTTOM, "B"),
    _V("froude", ("FROUDE NUMBER",), ("FROUDE",), "", "Froude number", STYLE_FROUDE, "F"),
    _V("scalar_flowrate", ("SCALAR FLOWRATE",), ("DEBIT SCALAIRE",), "M2/S", "Unit discharge", STYLE_VELOCITY, "Q"),
    _V("flowrate_x", ("FLOWRATE ALONG X",), ("DEBIT SUIVANT X",), "M2/S", "Unit discharge X", STYLE_VELOCITY, "QX"),
    _V("flowrate_y", ("FLOWRATE ALONG Y",), ("DEBIT SUIVANT Y",), "M2/S", "Unit discharge Y", STYLE_VELOCITY, "QY"),
    _V("scalar_velocity", ("SCALAR VELOCITY",), ("VITESSE SCALAIRE",), "M/S", "Velocity magnitude", STYLE_VELOCITY, "M"),
    _V("wind_x", ("WIND ALONG X",), ("VENT X",), "M/S", "Wind X", STYLE_GENERIC, "WX"),
    _V("wind_y", ("WIND ALONG Y",), ("VENT Y",), "M/S", "Wind Y", STYLE_GENERIC, "WY"),
    _V("air_pressure", ("AIR PRESSURE",), ("PRESSION ATMOS.",), "PASCAL", "Air pressure", STYLE_GENERIC, "P"),
    _V("bottom_friction", ("BOTTOM FRICTION",), ("FROTTEMENT",), "", "Bottom friction", STYLE_GENERIC, "N"),
    _V("drift_x", ("DRIFT ALONG X",), ("DERIVE EN X",), "M", "Drift X", STYLE_GENERIC, ""),
    _V("drift_y", ("DRIFT ALONG Y",), ("DERIVE EN Y",), "M", "Drift Y", STYLE_GENERIC, ""),
    _V("courant", ("COURANT NUMBER",), ("NBRE DE COURANT",), "", "Courant number", STYLE_GENERIC, ""),
    _V("viscosity", ("VISCOSITY",), ("VISCOSITE TURB.",), "M2/S", "Turbulent viscosity", STYLE_GENERIC, ""),
    _V("tke", ("TURBULENT ENERG.",), ("ENERGIE TURBUL.",), "J/KG", "Turbulent energy", STYLE_GENERIC, ""),
    _V("dissipation", ("DISSIPATION",), ("DISSIPATION",), "WATT/KG", "Dissipation", STYLE_GENERIC, ""),
    _V("friction_velocity", ("FRICTION VEL.",), ("VITESSE DE FROT.",), "M/S", "Friction velocity", STYLE_GENERIC, ""),
    _V("high_water_mark", ("HIGH WATER MARK",), ("COTE MAXIMUM",), "M", "Maximum free surface", STYLE_ELEVATION, ""),
    _V("high_water_time", ("HIGH WATER TIME",), ("TEMPS COTE MAXI",), "S", "Time of maximum free surface", STYLE_TIME, ""),
    _V("highest_velocity", ("HIGHEST VELOCITY",), ("VITESSE MAXIMUM",), "M/S", "Maximum velocity", STYLE_VELOCITY, ""),
    _V("time_high_velocity", ("TIME OF HIGH VEL",), ("T VITESSE MAXI",), "S", "Time of maximum velocity", STYLE_TIME, ""),
    _V("max_depth", ("MAXIMUM DEPTH",), ("HAUTEUR MAXIMUM",), "M", "Maximum depth", STYLE_DEPTH, ""),
    _V("time_max_depth", ("TIME OF MAX H",), ("TEMPS H MAXI",), "S", "Time of maximum depth", STYLE_TIME, ""),
    # --- TELEMAC-3D ---------------------------------------------------------
    _V("elevation_z", ("ELEVATION Z",), ("COTE Z",), "M", "Elevation Z", STYLE_ELEVATION, "Z"),
    _V("nuz", ("NUZ FOR VELOCITY",), ("NUZ POUR VITESSE",), "M2/S", "Vertical viscosity", STYLE_GENERIC, ""),
    # --- GAIA / sediment -----------------------------------------------------
    _V("evolution", ("EVOLUTION",), ("EVOLUTION",), "M", "Bed evolution", STYLE_GENERIC, ""),
    _V("rigid_bed", ("RIGID BED",), ("FOND RIGIDE",), "M", "Rigid bed", STYLE_BOTTOM, ""),
)

_BY_NAME: Dict[str, VarInfo] = {}
for _info in CATALOGUE:
    for _n in _info.en + _info.fr:
        _BY_NAME[_n.upper()] = _info


def _norm(name: str) -> str:
    return " ".join(str(name).strip().upper().split())


def lookup(name: str) -> Optional[VarInfo]:
    """Catalogue entry for a TELEMAC variable name, English or French."""
    return _BY_NAME.get(_norm(name))


def canonical(name: str) -> str:
    """Canonical id; unknown names become a slug of the name itself."""
    info = lookup(name)
    if info:
        return info.cid
    return _norm(name).lower().replace(" ", "_")


def display_name(name: str, unit: str = "") -> str:
    """Readable label with unit, e.g. 'Water depth (m)'."""
    info = lookup(name)
    label = info.label if info else str(name).strip().capitalize()
    u = (unit or (info.unit if info else "")).strip()
    return "%s (%s)" % (label, u.lower()) if u else label


# Explicit component pairs, then a generic rule for "<prefix> U/V" and "<prefix> X/Y".
_PAIRS = (
    ("velocity_u", "velocity_v", "Velocity"),
    ("flowrate_x", "flowrate_y", "Unit discharge"),
    ("wind_x", "wind_y", "Wind"),
    ("drift_x", "drift_y", "Drift"),
)


def vector_pairs(names: Iterable[str]) -> List[Tuple[int, int, str]]:
    """Index pairs (i_x, i_y, label) of vector components among variable names."""
    names = list(names)
    cids = [canonical(n) for n in names]
    used = set()
    out = []
    for cx, cy, label in _PAIRS:
        if cx in cids and cy in cids:
            i, j = cids.index(cx), cids.index(cy)
            out.append((i, j, label))
            used.update((i, j))
    norm = [_norm(n) for n in names]
    for i, n in enumerate(norm):
        if i in used:
            continue
        for a, b in ((" U", " V"), (" X", " Y"), (" ALONG X", " ALONG Y")):
            if n.endswith(a):
                partner = n[: -len(a)] + b
                if partner in norm:
                    j = norm.index(partner)
                    if j not in used and j != i:
                        out.append((i, j, n[: -len(a)].strip().capitalize() or "Vector"))
                        used.update((i, j))
                        break
    return out
