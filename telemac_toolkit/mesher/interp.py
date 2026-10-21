"""Interpolation of node values (bottom, friction, ...) from rasters and points.

Several sources are applied in priority order: a node takes its value from the
first source that covers it (for example channel survey, then LiDAR, then a
regional DEM). Nodes no source covers are reported, then filled from the
nearest node that has a value, so the geometry file never holds gaps.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from ..telemac_io.mesh import TriangleLocator

__all__ = ["RasterSource", "PointSource", "PlaneSource", "PrismaticChannelSource", "FormulaSource", "channel_formula", "bilinear", "interpolate_nodes", "fill_nearest"]


@dataclass
class RasterSource:
    """A raster window: values[row, col], top-left corner, pixel sizes (dy negative)."""

    values: np.ndarray
    x0: float
    y0: float
    dx: float
    dy: float
    nodata: Optional[float] = None
    scale: float = 1.0  # vertical unit factor, e.g. 1200/3937 for US survey feet
    name: str = "raster"

    def sample(self, px, py) -> np.ndarray:
        return bilinear(self.values, self.x0, self.y0, self.dx, self.dy, px, py, self.nodata) * self.scale


@dataclass
class PointSource:
    """Survey points interpolated linearly on their Delaunay triangulation (a TIN)."""

    xy: np.ndarray
    z: np.ndarray
    ikle: np.ndarray  # triangulation of xy
    scale: float = 1.0
    name: str = "points"

    def sample(self, px, py) -> np.ndarray:
        loc = TriangleLocator(self.xy[:, 0], self.xy[:, 1], self.ikle)
        return loc.interpolate(self.z, px, py) * self.scale


def bilinear(values, x0, y0, dx, dy, px, py, nodata=None) -> np.ndarray:
    """Bilinear interpolation between pixel centres; NaN outside or on nodata."""
    v = np.asarray(values, dtype=np.float64)
    if nodata is not None:
        v = np.where(v == nodata, np.nan, v)
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    c = (px - x0) / dx - 0.5
    r = (py - y0) / dy - 0.5
    nr, nc = v.shape
    out = np.full(px.shape, np.nan)
    inside = (c >= -0.5) & (c <= nc - 0.5) & (r >= -0.5) & (r <= nr - 0.5)
    if not inside.any():
        return out
    cc = np.clip(c[inside], 0, nc - 1)
    rr = np.clip(r[inside], 0, nr - 1)
    c0 = np.clip(np.floor(cc).astype(np.int64), 0, max(nc - 2, 0))
    r0 = np.clip(np.floor(rr).astype(np.int64), 0, max(nr - 2, 0))
    c1 = np.minimum(c0 + 1, nc - 1)
    r1 = np.minimum(r0 + 1, nr - 1)
    fc = cc - c0
    fr = rr - r0
    z00, z01, z10, z11 = v[r0, c0], v[r0, c1], v[r1, c0], v[r1, c1]
    w = np.stack([(1 - fr) * (1 - fc), (1 - fr) * fc, fr * (1 - fc), fr * fc], axis=1)
    z = np.stack([z00, z01, z10, z11], axis=1)
    ok = np.isfinite(z)
    wsum = np.sum(np.where(ok, w, 0.0), axis=1)
    val = np.sum(np.where(ok, w * np.nan_to_num(z), 0.0), axis=1)
    res = np.where(wsum > 1e-9, val / np.maximum(wsum, 1e-12), np.nan)
    out[inside] = res
    return out


def interpolate_nodes(px, py, sources: Sequence, progress: Optional[Callable[[float], bool]] = None):
    """Values at nodes from prioritised sources. Returns (values, source_index) (-1 = none)."""
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    vals = np.full(px.shape, np.nan)
    src = np.full(px.shape, -1, dtype=np.int64)
    for k, s in enumerate(sources):
        todo = np.flatnonzero(np.isnan(vals))
        if todo.size == 0:
            break
        v = s.sample(px[todo], py[todo])
        ok = np.isfinite(v)
        vals[todo[ok]] = v[ok]
        src[todo[ok]] = k
        if progress and progress((k + 1) / len(sources)):
            break
    return vals, src


def fill_nearest(px, py, vals: np.ndarray) -> np.ndarray:
    """Fill NaN values with the value of the nearest node that has one."""
    vals = np.asarray(vals, dtype=np.float64).copy()
    bad = np.flatnonzero(~np.isfinite(vals))
    good = np.flatnonzero(np.isfinite(vals))
    if bad.size == 0 or good.size == 0:
        return vals
    gx, gy = px[good], py[good]
    for chunk in np.array_split(bad, max(1, bad.size // 2000)):
        d2 = (px[chunk, None] - gx[None, :]) ** 2 + (py[chunk, None] - gy[None, :]) ** 2
        vals[chunk] = vals[good[np.argmin(d2, axis=1)]]
    return vals


@dataclass
class PlaneSource:
    """A sloping plane: z = z0 + gx (x - x0) + gy (y - y0). For flumes and schematic channels.

    ``from_slope`` builds it from a bed slope S (positive downhill) along a
    direction (degrees, counter-clockwise from the x axis): z falls by S per
    metre in that direction.
    """

    z0: float
    gx: float = 0.0
    gy: float = 0.0
    x0: float = 0.0
    y0: float = 0.0
    name: str = "plane"

    @classmethod
    def from_slope(cls, z0: float, slope: float, direction_deg: float = 0.0, x0: float = 0.0, y0: float = 0.0):
        import math

        a = math.radians(direction_deg)
        return cls(z0, -slope * math.cos(a), -slope * math.sin(a), x0, y0, "plane z0=%g slope=%g" % (z0, slope))

    def sample(self, px, py) -> np.ndarray:
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        return self.z0 + self.gx * (px - self.x0) + self.gy * (py - self.y0)


@dataclass
class PrismaticChannelSource:
    """Bed of a straight prismatic channel (flume or schematic river), along a local axis.

    The axis starts at (x0, y0) and points ``direction_deg`` counter-clockwise
    from the x axis. Along it the bed falls with ``slope``; across it the
    section is trapezoidal: a flat bed of ``bottom_width``, side slopes of
    ``side_slope`` horizontal to 1 vertical up to ``bank_height``, then flat
    banks. ``side_slope = 0`` gives a rectangular flume (flat bed everywhere).
    """

    z0: float
    slope: float = 0.0
    bottom_width: float = 1.0
    side_slope: float = 0.0
    bank_height: float = 0.0
    x0: float = 0.0
    y0: float = 0.0
    direction_deg: float = 0.0
    name: str = "channel section"

    def sample(self, px, py) -> np.ndarray:
        import math

        a = math.radians(self.direction_deg)
        dx = np.asarray(px, dtype=np.float64) - self.x0
        dy = np.asarray(py, dtype=np.float64) - self.y0
        s = dx * math.cos(a) + dy * math.sin(a)  # along the axis
        n = -dx * math.sin(a) + dy * math.cos(a)  # across, from the axis
        z = self.z0 - self.slope * s
        if self.side_slope > 0 and self.bank_height > 0:
            rise = np.clip((np.abs(n) - 0.5 * self.bottom_width) / self.side_slope, 0.0, self.bank_height)
            z = z + rise
        return z


# ---------------------------------------------------------------------------
# bottom from a formula of x and y (flumes, schematic channels, tests)
# ---------------------------------------------------------------------------
_FUNCS = {
    "abs": np.abs, "sqrt": np.sqrt, "exp": np.exp, "log": np.log, "log10": np.log10,
    "sin": np.sin, "cos": np.cos, "tan": np.tan, "arctan2": np.arctan2, "hypot": np.hypot,
    "minimum": np.minimum, "maximum": np.maximum, "min": np.minimum, "max": np.maximum,
    "clip": np.clip, "where": np.where, "floor": np.floor, "ceil": np.ceil,
}
_CONSTS = {"pi": np.pi, "e": np.e}


#: variables of a bottom formula: map coordinates, and channel coordinates along an axis line
FORMULA_VARIABLES = ("x", "y", "s", "n")


def formula_names(text: str) -> set:
    """Variables a formula uses (x, y, s, n)."""
    import ast

    try:
        tree = ast.parse(text.strip(), mode="eval")
    except SyntaxError:
        return set()
    return {node.id for node in ast.walk(tree) if isinstance(node, ast.Name) and node.id in FORMULA_VARIABLES}


def _check_formula(text: str):
    import ast

    tree = ast.parse(text.strip(), mode="eval")
    ok = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load, ast.Call, ast.Compare,
          ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.FloorDiv, ast.USub, ast.UAdd,
          ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.Eq, ast.NotEq, ast.BoolOp, ast.And, ast.Or, ast.IfExp)
    for node in ast.walk(tree):
        if not isinstance(node, ok):
            raise ValueError("not allowed in a bottom formula: %s" % type(node).__name__)
        if isinstance(node, ast.Name) and node.id not in FORMULA_VARIABLES and node.id not in _FUNCS and node.id not in _CONSTS:
            raise ValueError("unknown name %r in the formula (use x, y, s and n (along and across a channel axis), "
                             "numbers and %s)" % (node.id, ", ".join(sorted(_FUNCS))))
        if isinstance(node, ast.Call) and not (isinstance(node.func, ast.Name) and node.func.id in _FUNCS):
            raise ValueError("only these functions are allowed: %s" % ", ".join(sorted(_FUNCS)))
        if isinstance(node, ast.Constant) and not isinstance(node.value, (int, float)):
            raise ValueError("only numbers are allowed as constants")
    return tree.body


_BINOPS = {"Add": np.add, "Sub": np.subtract, "Mult": np.multiply, "Div": np.true_divide, "Pow": np.power,
           "Mod": np.mod, "FloorDiv": np.floor_divide}
_CMPOPS = {"Lt": np.less, "LtE": np.less_equal, "Gt": np.greater, "GtE": np.greater_equal, "Eq": np.equal,
           "NotEq": np.not_equal}


def _evaluate(node, env):
    """Walk a checked formula tree (no eval): numbers, names from ``env``, arithmetic, comparisons, calls."""
    import ast

    if isinstance(node, ast.Constant):
        return float(node.value)
    if isinstance(node, ast.Name):
        return env[node.id]
    if isinstance(node, ast.UnaryOp):
        v = _evaluate(node.operand, env)
        return -v if isinstance(node.op, ast.USub) else +v
    if isinstance(node, ast.BinOp):
        return _BINOPS[type(node.op).__name__](_evaluate(node.left, env), _evaluate(node.right, env))
    if isinstance(node, ast.Compare):
        left, out = _evaluate(node.left, env), None
        for op, comp in zip(node.ops, node.comparators):
            right = _evaluate(comp, env)
            r = _CMPOPS[type(op).__name__](left, right)
            out = r if out is None else np.logical_and(out, r)
            left = right
        return out
    if isinstance(node, ast.BoolOp):
        vals = [_evaluate(v, env) for v in node.values]
        fn = np.logical_and if isinstance(node.op, ast.And) else np.logical_or
        out = vals[0]
        for v in vals[1:]:
            out = fn(out, v)
        return out
    if isinstance(node, ast.IfExp):
        return np.where(_evaluate(node.test, env), _evaluate(node.body, env), _evaluate(node.orelse, env))
    if isinstance(node, ast.Call):
        return _FUNCS[node.func.id](*[_evaluate(a, env) for a in node.args])
    raise ValueError("not allowed in a bottom formula: %s" % type(node).__name__)


@dataclass
class FormulaSource:
    """Values from a formula, e.g. ``10 - 0.001 * x`` (NumPy syntax, a safe subset).

    Besides x and y, a formula may use channel coordinates when an ``axis``
    line is given: s, the distance along the axis from its start, and n, the
    signed distance from it (left positive looking downstream). The bed of a
    curved channel is then ``z0 - S*s + section(n)``, e.g.
    ``0.30 - 0.001*s + clip((abs(n) - 0.25)/1.5, 0, 0.2)``.
    """

    formula: str
    name: str = "formula"
    axis: Optional[np.ndarray] = None

    def __post_init__(self):
        self._tree = _check_formula(self.formula)
        self.uses = formula_names(self.formula)
        if self.uses & {"s", "n"} and (self.axis is None or len(self.axis) < 2):
            raise ValueError("the formula uses s or n (channel coordinates) but no channel axis line is given")

    def sample(self, px, py) -> np.ndarray:
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        env = {"x": px, "y": py}
        if self.uses & {"s", "n"}:
            from .planform import axis_coordinates

            sv, nv = axis_coordinates(self.axis, px, py)
            env["s"] = sv.reshape(px.shape)
            env["n"] = nv.reshape(px.shape)
        env.update(_FUNCS)
        env.update(_CONSTS)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            v = _evaluate(self._tree, env)
        return np.broadcast_to(np.asarray(v, dtype=np.float64), px.shape).copy()


def channel_formula_sn(z0: float, slope: float, bottom_width: float = 0.0, side_slope: float = 0.0,
                       bank_height: float = 0.0) -> str:
    """Bed of a prismatic channel of any planform, in channel coordinates s (along) and n (across)."""
    f = "%g - %g*s" % (z0, slope)
    if side_slope > 0 and bank_height > 0:
        f += " + clip((abs(n) - %g)/%g, 0, %g)" % (0.5 * bottom_width, side_slope, bank_height)
    return f


def channel_formula(z0: float, slope: float, bottom_width: float = 0.0, side_slope: float = 0.0,
                    bank_height: float = 0.0, x0: float = 0.0, y0: float = 0.0, direction_deg: float = 0.0) -> str:
    """The formula of a straight prismatic channel (see PrismaticChannelSource), readable and editable."""
    import math

    a = math.radians(direction_deg)
    c, s = math.cos(a), math.sin(a)
    xs = "(x - %g)" % x0 if x0 else "x"
    ys = "(y - %g)" % y0 if y0 else "y"
    if abs(s) < 1e-12:
        along = xs if c > 0 else "-%s" % xs
        across = ys if c > 0 else "-%s" % ys
    else:
        along = "(%s*%.9g + %s*%.9g)" % (xs, c, ys, s)
        across = "(-%s*%.9g + %s*%.9g)" % (xs, s, ys, c)
    f = "%g - %g*%s" % (z0, slope, along)
    if side_slope > 0 and bank_height > 0:
        f += " + clip((abs(%s) - %g)/%g, 0, %g)" % (across, 0.5 * bottom_width, side_slope, bank_height)
    return f
