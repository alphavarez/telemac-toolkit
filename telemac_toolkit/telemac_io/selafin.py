"""Read and write TELEMAC SELAFIN (Serafin) files.

SELAFIN is a Fortran sequential unformatted file. Every record is framed by
two 4-byte integers holding the record length. Layout:

    TITLE        80 chars, the last 8 give the format: "SERAFIN " (float32)
                 or "SERAFIND" (float64)
    NBV1, NBV2   number of linear / quadratic variables
    NAMES        NBV1 records of 32 chars (16 name + 16 unit)
    IPARAM       10 integers; (3, 4) mesh origin, (7) number of planes,
                 (8, 9) parallel counts, (10) = 1 when a date record follows
    DATE         6 integers (year, month, day, hour, minute, second), optional
    NELEM, NPOIN, NDP, 1
    IKLE         NELEM * NDP connectivity, 1-based
    IPOBO        NPOIN boundary numbering (KNOLG for partitioned files)
    X, Y         NPOIN coordinates (stored relative to the origin)
    then for every frame: TIME, followed by one record of NPOIN values per
    variable.

Every frame has the same size, so the reader exposes all values as a single
zero-copy strided memory map of shape (nframes, nvar, npoin). Reading a time
series at a few nodes touches only those bytes.

This module has no QGIS dependency.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import math
import os
import struct
import warnings
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence, Union

import numpy as np

__all__ = [
    "SelafinError",
    "SelafinHeader",
    "SelafinFile",
    "SelafinWriter",
    "write_selafin",
    "FLOAT32_SAFE_LIMIT",
]

#: Above this absolute coordinate, float32 resolution is worse than ~6 cm.
FLOAT32_SAFE_LIMIT = 500_000.0

_SINGLE_TAGS = ("SERAFIN ", "SERAFINS", "SELAFIN ", "SERAPHIN", "SERAFIN")
_DOUBLE_TAGS = ("SERAFIND", "SELAFIND", "       D")


class SelafinError(Exception):
    """Raised when a file is not a valid SELAFIN file."""


# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
@dataclass
class SelafinHeader:
    """Everything stored before the first time frame."""

    title: str = ""
    file_format: str = "SERAFIN "
    endian: str = ">"
    float_size: int = 4
    var_names: list = field(default_factory=list)
    var_units: list = field(default_factory=list)
    iparam: np.ndarray = field(default_factory=lambda: np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int64))
    date: Optional[tuple] = None
    nelem: int = 0
    npoin: int = 0
    ndp: int = 3
    ikle: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=np.int64))
    ipobo: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    x_stored: np.ndarray = field(default_factory=lambda: np.zeros(0))
    y_stored: np.ndarray = field(default_factory=lambda: np.zeros(0))
    header_size: int = 0

    # ---- derived properties -------------------------------------------
    @property
    def nvar(self) -> int:
        return len(self.var_names)

    @property
    def is_double(self) -> bool:
        return self.float_size == 8

    @property
    def float_dtype(self) -> np.dtype:
        return np.dtype("f%d" % self.float_size).newbyteorder(self.endian)

    @property
    def int_dtype(self) -> np.dtype:
        return np.dtype("i4").newbyteorder(self.endian)

    @property
    def origin(self) -> tuple:
        """Mesh origin from IPARAM(3), IPARAM(4); added to stored coordinates."""
        return float(self.iparam[2]), float(self.iparam[3])

    @property
    def x(self) -> np.ndarray:
        return self.x_stored + self.origin[0]

    @property
    def y(self) -> np.ndarray:
        return self.y_stored + self.origin[1]

    @property
    def nplan(self) -> int:
        """Number of horizontal planes (0 for a 2D file)."""
        return int(self.iparam[6])

    @property
    def is_3d(self) -> bool:
        return self.nplan > 1

    @property
    def npoin2(self) -> int:
        return self.npoin // self.nplan if self.is_3d else self.npoin

    @property
    def nelem2(self) -> int:
        return self.nelem // (self.nplan - 1) if self.is_3d else self.nelem

    @property
    def ikle2(self) -> np.ndarray:
        """2D triangles (0-based). For a 3D file, the bottom face of layer 1."""
        if self.is_3d:
            return self.ikle[: self.nelem2, :3]
        return self.ikle

    @property
    def x2(self) -> np.ndarray:
        return self.x[: self.npoin2]

    @property
    def y2(self) -> np.ndarray:
        return self.y[: self.npoin2]

    @property
    def is_partitioned(self) -> bool:
        """True for a sub-domain file of a parallel run (IPOBO holds KNOLG).

        IPARAM(8) holds NPTFR (boundary nodes) even in serial runs; only
        IPARAM(9) = NPTIR (interface nodes) is specific to sub-domains.
        """
        return int(self.iparam[8]) != 0

    @property
    def has_calendar_date(self) -> bool:
        """False when the date record is absent or TELEMAC's 1900-01-01 default."""
        d = self.start_datetime
        return d is not None and d.year > 1900

    @property
    def start_datetime(self) -> Optional[_dt.datetime]:
        if not self.date:
            return None
        y, mo, d, h, mi, s = (int(v) for v in self.date)
        if y <= 0 or not (1 <= mo <= 12) or not (1 <= d <= 31):
            return None
        try:
            return _dt.datetime(y, mo, d, h, mi, s)
        except ValueError:
            return None

    @property
    def frame_size(self) -> int:
        """Bytes per frame: time record + one record per variable."""
        return (8 + self.float_size) + self.nvar * (8 + self.npoin * self.float_size)

    def var_index(self, name_or_index: Union[str, int]) -> int:
        if isinstance(name_or_index, (int, np.integer)):
            idx = int(name_or_index)
            if not 0 <= idx < self.nvar:
                raise IndexError("variable index %d out of range (0..%d)" % (idx, self.nvar - 1))
            return idx
        key = str(name_or_index).strip().upper()
        for i, n in enumerate(self.var_names):
            if n.strip().upper() == key:
                return i
        raise KeyError("variable %r not in file; available: %s" % (name_or_index, ", ".join(self.var_names)))

    def has_var(self, name: str) -> bool:
        try:
            self.var_index(name)
            return True
        except KeyError:
            return False


# --------------------------------------------------------------------------
# Low-level record reading
# --------------------------------------------------------------------------
def _read_record(f, endian: str) -> bytes:
    head = f.read(4)
    if len(head) < 4:
        raise SelafinError("unexpected end of file while reading a record marker")
    (n,) = struct.unpack(endian + "i", head)
    if n < 0:
        raise SelafinError("negative record length %d" % n)
    data = f.read(n)
    tail = f.read(4)
    if len(data) < n or len(tail) < 4:
        raise SelafinError("unexpected end of file inside a record")
    (n2,) = struct.unpack(endian + "i", tail)
    if n2 != n:
        raise SelafinError("record markers disagree (%d vs %d): file is corrupted" % (n, n2))
    return data


def _detect_endian(first4: bytes) -> str:
    if struct.unpack(">i", first4)[0] == 80:
        return ">"
    if struct.unpack("<i", first4)[0] == 80:
        return "<"
    raise SelafinError("not a SELAFIN file: the first record must be 80 bytes long")


def _decode(b: bytes) -> str:
    return b.decode("latin-1")


def read_header(path: str) -> SelafinHeader:
    """Parse the header of a SELAFIN file."""
    h = SelafinHeader()
    with open(path, "rb") as f:
        first = f.read(4)
        if len(first) < 4:
            raise SelafinError("file is empty: %s" % path)
        h.endian = _detect_endian(first)
        f.seek(0)
        e = h.endian

        title = _decode(_read_record(f, e))
        h.title = title[:72].rstrip()
        fmt = title[72:80]
        if fmt in _DOUBLE_TAGS or fmt.strip() in ("SERAFIND", "SELAFIND"):
            h.file_format, h.float_size = "SERAFIND", 8
        else:
            if fmt not in _SINGLE_TAGS:
                warnings.warn("unknown SELAFIN format tag %r, reading as single precision" % fmt)
            h.file_format, h.float_size = "SERAFIN ", 4

        nbv = struct.unpack(e + "2i", _read_record(f, e))
        if nbv[1] != 0:
            raise SelafinError("quadratic variables (NBV2=%d) are not supported" % nbv[1])
        for _ in range(nbv[0]):
            rec = _decode(_read_record(f, e))
            h.var_names.append(rec[:16].strip())
            h.var_units.append(rec[16:32].strip())

        h.iparam = np.array(struct.unpack(e + "10i", _read_record(f, e)), dtype=np.int64)
        if h.iparam[9] == 1:
            h.date = tuple(int(v) for v in struct.unpack(e + "6i", _read_record(f, e)))

        nelem, npoin, ndp, one = struct.unpack(e + "4i", _read_record(f, e))
        if one != 1:
            # Some writers leave this at another value; it carries no information.
            warnings.warn("4th integer of the dimension record is %d instead of 1" % one)
        if nelem < 0 or npoin < 0 or ndp < 0:
            raise SelafinError("invalid dimensions NELEM=%d NPOIN=%d NDP=%d" % (nelem, npoin, ndp))
        h.nelem, h.npoin, h.ndp = nelem, npoin, ndp

        rec = _read_record(f, e)
        ikle = np.frombuffer(rec, dtype=h.int_dtype).astype(np.int64)
        if ikle.size != nelem * ndp:
            raise SelafinError("IKLE has %d values, expected %d" % (ikle.size, nelem * ndp))
        h.ikle = (ikle.reshape(nelem, ndp) - 1) if ndp else ikle.reshape(nelem, 0)

        rec = _read_record(f, e)
        h.ipobo = np.frombuffer(rec, dtype=h.int_dtype).astype(np.int64)

        rec_x = _read_record(f, e)
        rec_y = _read_record(f, e)
        # Coordinates are stored with the file precision. Some old writers store
        # them as float32 in double files; detect by record length.
        xsize = len(rec_x) // max(npoin, 1) if npoin else h.float_size
        cdtype = np.dtype("f%d" % xsize).newbyteorder(e) if xsize in (4, 8) else h.float_dtype
        h.x_stored = np.frombuffer(rec_x, dtype=cdtype).astype(np.float64)
        h.y_stored = np.frombuffer(rec_y, dtype=cdtype).astype(np.float64)
        h.header_size = f.tell()

    if h.ikle.size and (h.ikle.min() < 0 or h.ikle.max() >= h.npoin):
        raise SelafinError("IKLE references nodes outside 1..NPOIN")
    if h.is_3d and (h.npoin % h.nplan or h.nelem % (h.nplan - 1)):
        raise SelafinError("3D file with %d planes but NPOIN=%d, NELEM=%d" % (h.nplan, h.npoin, h.nelem))
    return h


# --------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------
class SelafinFile:
    """Lazy reader for a SELAFIN file.

    >>> slf = SelafinFile("results.slf")
    >>> slf.var_names
    ['VELOCITY U', 'VELOCITY V', 'WATER DEPTH', ...]
    >>> h = slf.read("WATER DEPTH", frame=-1)        # one frame, float64
    >>> ts = slf.series("WATER DEPTH", nodes=[10, 11]) # (nframes, 2)

    Values are never loaded as a whole: ``values`` is a read-only memory map
    of shape (nframes, nvar, npoin). Call :meth:`refresh` to pick up frames
    written by a simulation that is still running.
    """

    def __init__(self, path: str):
        self.path = os.fspath(path)
        self.header = read_header(self.path)
        self._mm = None
        self._values = None
        self._times = None
        self.nframes = 0
        self.partial_bytes = 0
        self.refresh()

    # ---- life cycle -------------------------------------------------------
    def refresh(self) -> int:
        """Re-scan the file size; returns the number of complete frames."""
        h = self.header
        size = os.path.getsize(self.path)
        body = size - h.header_size
        fs = h.frame_size
        n = body // fs if fs > 0 else 0
        rest = body - n * fs
        # One trailing byte (old Unix line feed) is tolerated silently.
        self.partial_bytes = rest if rest > 1 else 0
        if self.partial_bytes:
            warnings.warn(
                "%s: last frame is incomplete (%d bytes); it is ignored. "
                "The simulation may still be running." % (os.path.basename(self.path), rest)
            )
        self.nframes = int(n)
        self._close_map()
        if self.nframes > 0:
            self._mm = np.memmap(self.path, dtype=np.uint8, mode="r", offset=0, shape=(size,))
            fsz = h.float_size
            fdt = h.float_dtype
            # time values: one per frame, after the 4-byte record marker
            self._times = np.ndarray(
                shape=(self.nframes,),
                dtype=fdt,
                buffer=self._mm,
                offset=h.header_size + 4,
                strides=(fs,),
            )
            if h.nvar and h.npoin:
                self._values = np.ndarray(
                    shape=(self.nframes, h.nvar, h.npoin),
                    dtype=fdt,
                    buffer=self._mm,
                    offset=h.header_size + (8 + fsz) + 4,
                    strides=(fs, 8 + h.npoin * fsz, fsz),
                )
            else:
                self._values = np.zeros((self.nframes, h.nvar, h.npoin), dtype=np.float64)
        else:
            self._times = np.zeros(0)
            self._values = np.zeros((0, h.nvar, h.npoin))
        return self.nframes

    def _close_map(self):
        self._values = None
        self._times = None
        if self._mm is not None:
            mm = self._mm
            self._mm = None
            # release the file handle (Windows locks it); numpy may already have closed it
            with contextlib.suppress(AttributeError, ValueError, BufferError, OSError):
                mm._mmap.close()

    def close(self):
        self._close_map()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self):
        h = self.header
        return "<SelafinFile %s: %s, %d nodes, %d elements, %d variables, %d frames%s>" % (
            os.path.basename(self.path),
            "3D (%d planes)" % h.nplan if h.is_3d else "2D",
            h.npoin,
            h.nelem,
            h.nvar,
            self.nframes,
            ", double" if h.is_double else "",
        )

    # ---- shortcuts ----------------------------------------------------------
    @property
    def var_names(self) -> list:
        return list(self.header.var_names)

    @property
    def var_units(self) -> list:
        return list(self.header.var_units)

    @property
    def x(self) -> np.ndarray:
        return self.header.x

    @property
    def y(self) -> np.ndarray:
        return self.header.y

    @property
    def ikle(self) -> np.ndarray:
        return self.header.ikle

    @property
    def values(self) -> np.ndarray:
        """Read-only (nframes, nvar, npoin) memory map in file byte order."""
        return self._values

    @property
    def times(self) -> np.ndarray:
        """Time of every frame in seconds since the start (float64 copy)."""
        return np.asarray(self._times, dtype=np.float64)

    @property
    def datetimes(self) -> Optional[list]:
        start = self.header.start_datetime
        if start is None:
            return None
        return [start + _dt.timedelta(seconds=float(t)) for t in self.times]

    def _frame_index(self, frame: int) -> int:
        n = self.nframes
        if n == 0:
            raise IndexError("file has no complete frame")
        i = int(frame)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("frame %d out of range (0..%d)" % (frame, n - 1))
        return i

    # ---- data access --------------------------------------------------------
    def read(self, var: Union[str, int], frame: int = -1) -> np.ndarray:
        """All node values of one variable at one frame (float64 copy)."""
        v = self.header.var_index(var)
        return np.asarray(self._values[self._frame_index(frame), v], dtype=np.float64)

    def read_frame(self, frame: int = -1) -> np.ndarray:
        """All variables at one frame, shape (nvar, npoin), float64 copy."""
        return np.asarray(self._values[self._frame_index(frame)], dtype=np.float64)

    def series(self, var: Union[str, int], nodes: Union[int, Sequence[int]]) -> np.ndarray:
        """Time series at one or more nodes, shape (nframes,) or (nframes, k)."""
        v = self.header.var_index(var)
        if np.isscalar(nodes):
            return np.asarray(self._values[:, v, int(nodes)], dtype=np.float64)
        idx = np.asarray(nodes, dtype=np.int64)
        return np.asarray(self._values[:, v, :][:, idx], dtype=np.float64)

    def series_weighted(self, var: Union[str, int], nodes: np.ndarray, weights: np.ndarray) -> np.ndarray:
        """Interpolated time series: sum_k w_k * value(node_k) for every frame."""
        s = self.series(var, np.asarray(nodes).ravel())
        return s @ np.asarray(weights, dtype=np.float64).ravel()

    def find_frame(self, time_s: float) -> int:
        """Index of the frame nearest to a time in seconds."""
        t = self.times
        if t.size == 0:
            raise IndexError("file has no complete frame")
        return int(np.argmin(np.abs(t - float(time_s))))

    # ---- 3D helpers ---------------------------------------------------------
    def plane(self, var: Union[str, int], frame: int, k: int) -> np.ndarray:
        """Values of horizontal plane k (0 = bottom, -1 = surface)."""
        h = self.header
        if not h.is_3d:
            if k not in (0, -1):
                raise IndexError("2D file has a single plane")
            return self.read(var, frame)
        if k < 0:
            k += h.nplan
        if not 0 <= k < h.nplan:
            raise IndexError("plane %d out of range (0..%d)" % (k, h.nplan - 1))
        n2 = h.npoin2
        v = h.var_index(var)
        return np.asarray(self._values[self._frame_index(frame), v, k * n2:(k + 1) * n2], dtype=np.float64)

    def depth_average(self, var: Union[str, int], frame: int, z_var: Optional[str] = None) -> np.ndarray:
        """Depth-averaged values of a 3D variable (trapezoidal rule along z)."""
        h = self.header
        if not h.is_3d:
            return self.read(var, frame)
        if z_var is None:
            z_var = next((n for n in h.var_names if n.strip().upper() in ("ELEVATION Z", "COTE Z")), None)
            if z_var is None:
                raise KeyError("3D file has no ELEVATION Z / COTE Z variable")
        n2, npl = h.npoin2, h.nplan
        z = self.read(z_var, frame).reshape(npl, n2)
        f = self.read(var, frame).reshape(npl, n2)
        dz = np.diff(z, axis=0)
        num = np.sum(0.5 * (f[1:] + f[:-1]) * dz, axis=0)
        depth = z[-1] - z[0]
        out = np.full(n2, np.nan)
        ok = depth > 1e-9
        out[ok] = num[ok] / depth[ok]
        out[~ok] = f[0, ~ok]
        return out

    # ---- summaries ------------------------------------------------------------
    def info(self) -> dict:
        h = self.header
        t = self.times
        d = {
            "path": self.path,
            "title": h.title,
            "format": h.file_format.strip(),
            "precision": "double" if h.is_double else "single",
            "endian": "big" if h.endian == ">" else "little",
            "dimension": "3D" if h.is_3d else "2D",
            "planes": h.nplan if h.is_3d else 1,
            "nodes": h.npoin,
            "nodes_2d": h.npoin2,
            "elements": h.nelem,
            "nodes_per_element": h.ndp,
            "variables": [{"name": n, "unit": u} for n, u in zip(h.var_names, h.var_units)],
            "frames": self.nframes,
            "first_time_s": float(t[0]) if t.size else None,
            "last_time_s": float(t[-1]) if t.size else None,
            "time_step_s": float(np.median(np.diff(t))) if t.size > 1 else None,
            "start_date": h.start_datetime.isoformat() if h.start_datetime else None,
            "origin": h.origin,
            "partitioned": h.is_partitioned,
            "incomplete_last_frame_bytes": self.partial_bytes,
        }
        if h.npoin:
            d["extent"] = (float(h.x.min()), float(h.y.min()), float(h.x.max()), float(h.y.max()))
            d["float32_precision_warning"] = precision_warning(h) is not None
        return d


def precision_warning(h: SelafinHeader) -> Optional[str]:
    """Message when float32 storage cannot resolve the coordinates well."""
    if h.is_double or not h.npoin:
        return None
    m = max(float(np.abs(h.x_stored).max()), float(np.abs(h.y_stored).max()))
    if m <= FLOAT32_SAFE_LIMIT:
        return None
    resolution = math.ldexp(1.0, math.frexp(m)[1] - 24)
    return (
        "Coordinates are stored as float32 with magnitudes up to %.0f m, so node "
        "positions are only resolved to about %.2f m. Write the file with a mesh "
        "origin offset or in double precision." % (m, resolution)
    )


# --------------------------------------------------------------------------
# Writer
# --------------------------------------------------------------------------
class SelafinWriter:
    """Write a SELAFIN file frame by frame.

    >>> with SelafinWriter("geo.slf", ["BOTTOM"], ["M"], x, y, ikle) as w:
    ...     w.write_frame(0.0, [bottom])

    ``ikle`` is 0-based with shape (nelem, ndp). Coordinates are absolute;
    when ``origin`` is "auto" and the file is single precision, an integer
    origin is subtracted so float32 keeps sub-centimetre resolution, and stored
    in IPARAM(3), IPARAM(4) where TELEMAC expects it.
    """

    def __init__(
        self,
        path: str,
        var_names: Sequence[str],
        var_units: Sequence[str],
        x: np.ndarray,
        y: np.ndarray,
        ikle: np.ndarray,
        ipobo: Optional[np.ndarray] = None,
        title: str = "",
        double: bool = False,
        endian: str = ">",
        date: Optional[Union[_dt.datetime, Sequence[int]]] = None,
        origin: Union[str, tuple, None] = "auto",
        nplan: int = 0,
        iparam: Optional[Sequence[int]] = None,
    ):
        if endian not in (">", "<"):
            raise ValueError("endian must be '>' or '<'")
        if len(var_names) != len(var_units):
            raise ValueError("var_names and var_units differ in length")
        for n in var_names:
            if len(n) > 16:
                raise ValueError("variable name %r is longer than 16 characters" % n)
        x = np.asarray(x, dtype=np.float64).ravel()
        y = np.asarray(y, dtype=np.float64).ravel()
        if x.shape != y.shape:
            raise ValueError("x and y differ in length")
        ikle = np.asarray(ikle, dtype=np.int64)
        if ikle.ndim != 2:
            raise ValueError("ikle must have shape (nelem, ndp)")
        npoin = x.size
        if ikle.size and (ikle.min() < 0 or ikle.max() >= npoin):
            raise ValueError("ikle references nodes outside 0..npoin-1")

        self.path = os.fspath(path)
        self.endian = endian
        self.float_size = 8 if double else 4
        self.fdt = np.dtype("f%d" % self.float_size).newbyteorder(endian)
        self.idt = np.dtype("i4").newbyteorder(endian)
        self.nvar = len(var_names)
        self.npoin = npoin

        ip = np.zeros(10, dtype=np.int64) if iparam is None else np.array(iparam, dtype=np.int64).copy()
        ip[0] = 1
        if nplan:
            ip[6] = int(nplan)
        if isinstance(origin, str) and origin == "auto":
            if not double and npoin and max(np.abs(x).max(), np.abs(y).max()) > FLOAT32_SAFE_LIMIT:
                ox, oy = math.floor(x.min()), math.floor(y.min())
            else:
                ox, oy = int(ip[2]), int(ip[3])
        elif origin is None:
            ox, oy = 0, 0
        else:
            ox, oy = int(origin[0]), int(origin[1])
        ip[2], ip[3] = ox, oy

        if date is not None:
            if isinstance(date, _dt.datetime):
                date = (date.year, date.month, date.day, date.hour, date.minute, date.second)
            date = tuple(int(v) for v in date)
            if len(date) != 6:
                raise ValueError("date must have 6 integers")
            ip[9] = 1
        else:
            ip[9] = 0
        self.iparam = ip
        self.origin = (ox, oy)

        if ipobo is None:
            ipobo = boundary_numbering(ikle, npoin, x, y) if ikle.shape[1] == 3 and ikle.shape[0] else np.zeros(npoin, dtype=np.int64)
        ipobo = np.asarray(ipobo, dtype=np.int64)
        if ipobo.size != npoin:
            raise ValueError("ipobo must have npoin values")

        self._f = open(self.path, "wb")
        try:
            tag = "SERAFIND" if double else "SERAFIN "
            self._rec(title[:72].ljust(72).encode("latin-1", "replace") + tag.encode("ascii"))
            self._rec(struct.pack(endian + "2i", self.nvar, 0))
            for n, u in zip(var_names, var_units):
                self._rec(n.ljust(16)[:16].encode("latin-1", "replace") + u.ljust(16)[:16].encode("latin-1", "replace"))
            self._rec(struct.pack(endian + "10i", *(int(v) for v in ip)))
            if date is not None:
                self._rec(struct.pack(endian + "6i", *date))
            self._rec(struct.pack(endian + "4i", ikle.shape[0], npoin, ikle.shape[1], 1))
            self._rec((ikle + 1).astype(self.idt).tobytes())
            self._rec(ipobo.astype(self.idt).tobytes())
            self._rec((x - ox).astype(self.fdt).tobytes())
            self._rec((y - oy).astype(self.fdt).tobytes())
        except Exception:
            self._f.close()
            raise
        self.nframes = 0

    def _rec(self, payload: bytes):
        n = struct.pack(self.endian + "i", len(payload))
        self._f.write(n)
        self._f.write(payload)
        self._f.write(n)

    def write_frame(self, time_s: float, values: Union[np.ndarray, Iterable[np.ndarray]]):
        """Write one frame. ``values`` has one array of npoin values per variable."""
        arr = np.asarray(values, dtype=np.float64)
        if arr.ndim == 1 and self.nvar == 1:
            arr = arr[None, :]
        if arr.shape != (self.nvar, self.npoin):
            raise ValueError("values must have shape (%d, %d), got %s" % (self.nvar, self.npoin, arr.shape))
        self._rec(np.array([time_s], dtype=self.fdt).tobytes())
        for row in arr:
            self._rec(row.astype(self.fdt).tobytes())
        self.nframes += 1

    def close(self):
        if self._f and not self._f.closed:
            self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def write_selafin(
    path: str,
    x: np.ndarray,
    y: np.ndarray,
    ikle: np.ndarray,
    variables: dict,
    times: Optional[Sequence[float]] = None,
    units: Optional[dict] = None,
    **kwargs,
) -> str:
    """Write a whole file at once.

    ``variables`` maps a name to an array of shape (npoin,) for a single frame
    or (nframes, npoin) for several frames.
    """
    names = list(variables)
    units = units or {}
    data = [np.atleast_2d(np.asarray(variables[n], dtype=np.float64)) for n in names]
    nfr = data[0].shape[0] if data else 0
    if any(d.shape[0] != nfr for d in data):
        raise ValueError("all variables need the same number of frames")
    if times is None:
        times = np.arange(nfr, dtype=np.float64)
    if len(times) != nfr:
        raise ValueError("times has %d values for %d frames" % (len(times), nfr))
    with SelafinWriter(path, names, [units.get(n, "") for n in names], x, y, ikle, **kwargs) as w:
        for i in range(nfr):
            w.write_frame(times[i], [d[i] for d in data])
    return path


def boundary_numbering(ikle: np.ndarray, npoin: int, x: Optional[np.ndarray] = None, y: Optional[np.ndarray] = None) -> np.ndarray:
    """IPOBO for a triangle mesh: boundary nodes numbered 1..NPTFR, others 0.

    The outer boundary is walked counter-clockwise first, from the node with
    the smallest index on it, then each island in turn, as TELEMAC numbers them.
    """
    from .mesh import boundary_loops  # local import keeps module import light

    ipobo = np.zeros(npoin, dtype=np.int64)
    k = 1
    for loop in boundary_loops(ikle, x, y):
        for n in loop:
            if ipobo[n] == 0:
                ipobo[n] = k
                k += 1
    return ipobo
