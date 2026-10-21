"""Derived SELAFIN files: 2D extraction from 3D results, frame subsets.

QGIS's MDAL reader (QGIS 3.44 and 4.2) does not open 3D SELAFIN files, so the
plugin shows 3D results through a 2D file written here: one horizontal plane
or the depth average, for every frame.
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence, Union

import numpy as np

from .selafin import SelafinFile, SelafinWriter

__all__ = ["extract_2d", "subset_frames"]

Progress = Optional[Callable[[float], bool]]


def _frames(slf: SelafinFile, frames: Optional[Sequence[int]]):
    if frames is None:
        return list(range(slf.nframes))
    return [slf._frame_index(i) for i in frames]


def extract_2d(
    src: Union[str, SelafinFile],
    dst: str,
    plane: Union[int, str] = "average",
    frames: Optional[Sequence[int]] = None,
    variables: Optional[Sequence[str]] = None,
    progress: Progress = None,
) -> str:
    """Write a 2D SELAFIN from a 3D one.

    plane: an index (0 = bottom, -1 = surface), "bottom", "surface" or
    "average" (depth average with ELEVATION Z; Z itself becomes the free
    surface in that case). A 2D source is copied (optionally subset).
    ``progress(fraction)`` may return True to cancel.
    """
    own = isinstance(src, str)
    slf = SelafinFile(src) if own else src
    try:
        h = slf.header
        names = list(h.var_names) if variables is None else [h.var_names[h.var_index(v)] for v in variables]
        units = [h.var_units[h.var_index(n)] for n in names]
        is3d = h.is_3d
        if isinstance(plane, str):
            key = plane.lower()
            if key == "bottom":
                plane = 0
            elif key in ("surface", "top"):
                plane = -1
            elif key != "average":
                raise ValueError("plane must be an index, 'bottom', 'surface' or 'average'")
        zname = next((n for n in h.var_names if n.strip().upper() in ("ELEVATION Z", "COTE Z")), None)

        out_names, out_units = [], []
        for n, u in zip(names, units):
            if is3d and plane == "average" and n == zname:
                out_names.append("FREE SURFACE" if "ELEVATION" in n.upper() else "SURFACE LIBRE")
            else:
                out_names.append(n)
            out_units.append(u)
        # In "average" mode also give the bottom (lowest plane of Z)
        add_bottom = is3d and plane == "average" and zname is not None
        if add_bottom:
            out_names.append("BOTTOM" if "ELEVATION" in zname.upper() else "FOND")
            out_units.append("M")

        title = (h.title + (" - 2D %s" % (plane if isinstance(plane, str) else "plane %d" % plane) if is3d else ""))[:72]
        frs = _frames(slf, frames)
        with SelafinWriter(
            dst, out_names, out_units, h.x2, h.y2, h.ikle2,
            ipobo=h.ipobo[: h.npoin2] if not h.is_partitioned else None,
            title=title, double=h.is_double, endian=h.endian,
            date=h.date if h.start_datetime else None,
            origin=h.origin,
        ) as w:
            times = slf.times
            for k, i in enumerate(frs):
                vals = []
                for n in names:
                    if not is3d:
                        vals.append(slf.read(n, i))
                    elif plane == "average":
                        vals.append(slf.plane(n, i, -1) if n == zname else slf.depth_average(n, i, zname))
                    else:
                        vals.append(slf.plane(n, i, int(plane)))
                if add_bottom:
                    vals.append(slf.plane(zname, i, 0))
                w.write_frame(times[i], vals)
                if progress and progress((k + 1) / len(frs)):
                    break
        return dst
    finally:
        if own:
            slf.close()


def subset_frames(src: str, dst: str, frames: Sequence[int], variables: Optional[Sequence[str]] = None) -> str:
    """Copy selected frames (and variables) of a 2D or 3D file."""
    with SelafinFile(src) as slf:
        h = slf.header
        names = list(h.var_names) if variables is None else [h.var_names[h.var_index(v)] for v in variables]
        units = [h.var_units[h.var_index(n)] for n in names]
        with SelafinWriter(
            dst, names, units, h.x, h.y, h.ikle, ipobo=h.ipobo, title=h.title,
            double=h.is_double, endian=h.endian, date=h.date if h.start_datetime else None,
            origin=h.origin, iparam=h.iparam,
        ) as w:
            for i in _frames(slf, frames):
                w.write_frame(slf.times[i], [slf.read(n, i) for n in names])
    return dst


def writer_like(h, dst: str, names, units, title: Optional[str] = None) -> SelafinWriter:
    """A writer with the same mesh as header ``h``: nodes, elements, boundary numbering (so an existing .cli
    still fits), origin, start date, precision and byte order; only the variables change."""
    return SelafinWriter(
        dst, list(names), list(units), h.x, h.y, h.ikle,
        ipobo=h.ipobo if (h.ipobo is not None and len(h.ipobo) == h.npoin and not h.is_partitioned) else None,
        title=(h.title if title is None else title)[:72], double=h.is_double, endian=h.endian,
        date=h.date if h.start_datetime else None, origin=h.origin,
    )
