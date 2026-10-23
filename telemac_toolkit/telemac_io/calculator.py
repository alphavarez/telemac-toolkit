"""Field calculator for SELAFIN files: algebra, spatial calculus on the mesh, time reductions.

Expressions use NumPy syntax on the variables of the file::

    VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)
    FROUDE = sqrt(U**2 + V**2) / sqrt(g * maximum(H, 1e-6))
    HMAX [M] = tmax(H)
    SLOPE X = ddx(S)
    SLOPE Y = ddy(S)

Variables: TELEMAC symbols (U, V, H, S, B, F, Q, M, ...), names with
underscores (WATER_DEPTH, VELOCITY_U) and canonical ids (water_depth), plus x,
y, t (time in s), area (area around each node) and g = 9.81. Every earlier
output line can be used by its name in the later ones.

Spatial operators are those of linear (P1) finite elements on triangles, the
interpolation TELEMAC itself uses: the gradient is constant in each triangle
and exact for a linear field; values at nodes are the area-weighted average of
the triangles around the node (the lumped L2 projection). ``laplacian`` uses
the weak form (stiffness matrix over lumped mass); like any finite-element
Laplacian it lacks the boundary flux term on the boundary nodes.

Time operators reduce over all frames (or differentiate): tmax, tmin, tmean,
tsum, tint (trapezoidal time integral), tpeak (time of the maximum), tfirst,
tlast, ddt (time derivative, central differences).

No ``eval``: the expression tree is walked by this module, and only the
names and node types listed here are accepted.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .variables import CATALOGUE, canonical, lookup

__all__ = ["MeshOperators", "Calculator", "Output", "parse_outputs", "FUNCTIONS_HELP", "CalculatorError"]

G = 9.81


class CalculatorError(ValueError):
    pass


# ----------------------------------------------------------------------------
# P1 operators
# ----------------------------------------------------------------------------
class MeshOperators:
    """Gradient, divergence, curl, Laplacian and integrals of nodal fields on a triangle mesh."""

    def __init__(self, x, y, ikle):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        ik = np.asarray(ikle, dtype=np.int64)
        if ik.ndim != 2 or ik.shape[1] != 3:
            raise CalculatorError("spatial operators need a triangular mesh")
        self.x, self.y, self.ikle = x, y, ik
        x1, x2, x3 = x[ik[:, 0]], x[ik[:, 1]], x[ik[:, 2]]
        y1, y2, y3 = y[ik[:, 0]], y[ik[:, 1]], y[ik[:, 2]]
        a2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)  # twice the signed area
        bad = np.abs(a2) <= 1e-300
        a2 = np.where(bad, 1.0, a2)
        # d(phi_i)/dx and d(phi_i)/dy of the three linear shape functions
        self.bx = np.stack([y2 - y3, y3 - y1, y1 - y2], axis=1) / a2[:, None]
        self.by = np.stack([x3 - x2, x1 - x3, x2 - x1], axis=1) / a2[:, None]
        self.bx[bad] = 0.0
        self.by[bad] = 0.0
        self.area = np.where(bad, 0.0, 0.5 * np.abs(a2))
        n = x.size
        self.node_area = np.zeros(n)
        np.add.at(self.node_area, ik.ravel(), np.repeat(self.area / 3.0, 3))
        self._wsum = np.zeros(n)
        np.add.at(self._wsum, ik.ravel(), np.repeat(self.area, 3))
        self._wsum[self._wsum == 0] = 1.0

    # element -> node ----------------------------------------------------------
    def to_nodes(self, ve: np.ndarray) -> np.ndarray:
        """Area-weighted average of element values around each node."""
        out = np.zeros(self.x.size)
        np.add.at(out, self.ikle.ravel(), np.repeat(self.area * ve, 3))
        return out / self._wsum

    def element_gradient(self, f) -> Tuple[np.ndarray, np.ndarray]:
        fe = np.asarray(f, dtype=np.float64)[self.ikle]
        return np.einsum("ij,ij->i", fe, self.bx), np.einsum("ij,ij->i", fe, self.by)

    def ddx(self, f):
        return self.to_nodes(self.element_gradient(f)[0])

    def ddy(self, f):
        return self.to_nodes(self.element_gradient(f)[1])

    def grad(self, f):
        gx, gy = self.element_gradient(f)
        return self.to_nodes(gx), self.to_nodes(gy)

    def div(self, u, v):
        ux, _ = self.element_gradient(u)
        _, vy = self.element_gradient(v)
        return self.to_nodes(ux + vy)

    def curl(self, u, v):
        """Vertical vorticity dv/dx - du/dy."""
        _, uy = self.element_gradient(u)
        vx, _ = self.element_gradient(v)
        return self.to_nodes(vx - uy)

    def strain(self, u, v):
        """Magnitude of the strain rate: sqrt(2 ux² + 2 vy² + (uy + vx)²)."""
        ux, uy = self.element_gradient(u)
        vx, vy = self.element_gradient(v)
        return self.to_nodes(np.sqrt(2 * ux ** 2 + 2 * vy ** 2 + (uy + vx) ** 2))

    def okubo_weiss(self, u, v):
        """sn² + ss² - w²: negative where rotation dominates (eddies), positive where strain does."""
        ux, uy = self.element_gradient(u)
        vx, vy = self.element_gradient(v)
        return self.to_nodes((ux - vy) ** 2 + (vx + uy) ** 2 - (vx - uy) ** 2)

    def laplacian(self, f):
        """div(grad f): the divergence of the nodal gradient (exact for quadratic fields on symmetric patches,
        first order on irregular meshes, one-sided on the boundary)."""
        gx, gy = self.grad(f)
        return self.div(gx, gy)

    def smooth(self, f):
        """Average of the element means around each node (one smoothing pass)."""
        fe = np.asarray(f, dtype=np.float64)[self.ikle].mean(axis=1)
        return self.to_nodes(fe)

    def integral(self, f, element_mask: Optional[np.ndarray] = None) -> float:
        """Integral over the mesh (exact for the linear interpolation of f)."""
        fe = np.asarray(f, dtype=np.float64)[self.ikle].mean(axis=1) * self.area
        if element_mask is not None:
            fe = fe[element_mask]
        return float(np.nansum(fe))


# ----------------------------------------------------------------------------
# expressions
# ----------------------------------------------------------------------------
_NP = {
    "abs": np.abs, "sqrt": np.sqrt, "exp": np.exp, "log": np.log, "log10": np.log10, "sin": np.sin, "cos": np.cos,
    "tan": np.tan, "arcsin": np.arcsin, "arccos": np.arccos, "arctan": np.arctan, "arctan2": np.arctan2,
    "hypot": np.hypot, "minimum": np.minimum, "maximum": np.maximum, "min": np.minimum, "max": np.maximum,
    "clip": np.clip, "where": np.where, "floor": np.floor, "ceil": np.ceil, "sign": np.sign, "degrees": np.degrees,
    "radians": np.radians, "isnan": np.isnan,
}
_SPATIAL1 = ("ddx", "ddy", "grad", "laplacian", "smooth")
_SPATIAL2 = ("div", "curl", "vorticity", "strain", "okubo_weiss", "mag", "bearing", "direction")
_TIME = ("tmax", "tmin", "tmean", "tsum", "tint", "tpeak", "tfirst", "tlast", "tstd")
_SPECIAL = ("ddt", "wet")
FUNCTIONS = set(_NP) | set(_SPATIAL1) | set(_SPATIAL2) | set(_TIME) | set(_SPECIAL)

FUNCTIONS_HELP = [
    ("ddx(f), ddy(f)", "spatial derivatives (P1 gradient, averaged to the nodes)"),
    ("grad(f)", "magnitude of the gradient, e.g. grad(S) = water surface slope"),
    ("div(u, v)", "divergence du/dx + dv/dy"),
    ("curl(u, v) or vorticity(u, v)", "vertical vorticity dv/dx - du/dy (1/s)"),
    ("strain(u, v)", "strain rate magnitude sqrt(2ux² + 2vy² + (uy+vx)²)"),
    ("okubo_weiss(u, v)", "sn² + ss² - w², negative inside eddies"),
    ("laplacian(f)", "Laplacian, div(grad f)"),
    ("smooth(f)", "one pass of nodal smoothing"),
    ("mag(u, v)", "vector magnitude"),
    ("bearing(u, v)", "direction towards which the vector points, degrees from north, clockwise"),
    ("tmax(f), tmin(f), tmean(f), tstd(f)", "maximum, minimum, mean, standard deviation over all frames"),
    ("tsum(f), tint(f)", "sum over frames; time integral (trapezoidal, in f x s)"),
    ("tpeak(f)", "time (s) of the maximum of f"),
    ("tfirst(f), tlast(f)", "f at the first or the last frame"),
    ("ddt(f)", "time derivative (central differences), e.g. ddt(S) = rate of rise"),
    ("wet(eps)", "1 where the water depth is above eps, else 0"),
    ("where(c, a, b), clip, minimum, maximum, abs, sqrt, exp, log, sin, cos, arctan2, ...", "NumPy functions"),
]

_ALLOWED_NODES = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load, ast.Call, ast.Compare,
                  ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.FloorDiv, ast.USub, ast.UAdd,
                  ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.Eq, ast.NotEq, ast.BoolOp, ast.And, ast.Or, ast.IfExp,
                  ast.BitAnd, ast.BitOr, ast.Invert, ast.Not)

_LINE = re.compile(r"^\s*(?P<name>[^=\[\]#]+?)\s*(\[(?P<unit>[^\]]*)\])?\s*=(?!=)\s*(?P<expr>.+?)\s*$")


@dataclass
class Output:
    name: str          # SELAFIN name, at most 16 characters, upper case
    unit: str
    expr: str
    tree: ast.AST = None
    time_dependent: bool = True


def _slug(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip().upper())


def parse_outputs(text: str) -> List[Output]:
    """One ``NAME [UNIT] = expression`` per line (or ';'); a bare expression is called RESULT."""
    outs: List[Output] = []
    lines = []
    for raw in (text or "").replace("\r", "\n").split("\n"):
        line = raw.split("#", 1)[0].strip()
        if line:
            lines.extend(p.strip() for p in line.split(";") if p.strip())
    for k, line in enumerate(lines):
        m = _LINE.match(line)
        if m and not re.match(r"^\s*[\w.]+\s*[<>!]$", m.group("name")):
            name, unit, expr = _slug(m.group("name")), (m.group("unit") or "").strip().upper(), m.group("expr")
        else:
            name, unit, expr = "RESULT" if len(lines) == 1 else "RESULT %d" % (k + 1), "", line
        if len(name) > 16:
            raise CalculatorError("output name %r is longer than 16 characters (a SELAFIN limit)" % name)
        if len(unit) > 16:
            raise CalculatorError("unit %r is longer than 16 characters" % unit)
        outs.append(Output(name, unit, expr))
    if not outs:
        raise CalculatorError("no expression given")
    names = [o.name for o in outs]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise CalculatorError("output names used twice: %s" % ", ".join(sorted(dup)))
    return outs


def variable_names(names: Sequence[str]) -> Dict[str, int]:
    """Every way an expression can refer to the file's variables -> variable index."""
    out: Dict[str, int] = {}
    for i, n in enumerate(names):
        info = lookup(n)
        keys = [n.strip().upper().replace(" ", "_").replace("'", "_").replace(".", ""), canonical(n)]
        if info is not None:
            if info.symbol:
                keys.insert(0, info.symbol)
            keys.append(info.cid)
        for k in keys:
            k = k.strip()
            if k and k.isidentifier() and k not in out:
                out[k] = i
    return out


def _check(tree: ast.AST, known: set):
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise CalculatorError("not allowed in an expression: %s" % type(node).__name__)
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS:
                raise CalculatorError("unknown function %r" % getattr(node.func, "id", "?"))
            if node.keywords:
                raise CalculatorError("keyword arguments are not supported")
        elif isinstance(node, ast.Name) and node.id not in known and node.id not in FUNCTIONS:
            raise CalculatorError("unknown name %r" % node.id)
        if isinstance(node, ast.Constant) and not isinstance(node.value, (int, float, bool)):
            raise CalculatorError("only numbers are allowed as constants")


class _Inline(ast.NodeTransformer):
    def __init__(self, defs: Dict[str, ast.AST]):
        self.defs = defs

    def visit_Name(self, node):
        if node.id in self.defs:
            import copy

            return copy.deepcopy(self.defs[node.id])
        return node


def _ident(name: str) -> str:
    return name.strip().upper().replace(" ", "_")


class Calculator:
    """Evaluate output expressions frame by frame on one SELAFIN file (a ``SelafinFile``)."""

    def __init__(self, slf, outputs: Sequence[Output], dry: float = 0.01):
        h = slf.header
        if h.is_3d:
            raise CalculatorError("3D files: extract a 2D file first (depth average or a plane)")
        self.slf = slf
        self.h = h
        self.ops = MeshOperators(h.x, h.y, h.ikle) if h.ndp == 3 else None
        self.var_index = variable_names(h.var_names)
        self.times = np.asarray(slf.times, dtype=np.float64)
        self.nframes = len(self.times)
        self.dry = dry
        self._frame_cache: Dict[int, np.ndarray] = {}
        self._red_cache: Dict[int, np.ndarray] = {}
        self.depth_index = next((i for i, n in enumerate(h.var_names) if canonical(n) == "water_depth"), None)
        consts = {"x", "y", "t", "area", "g", "pi", "e"}
        defs: Dict[str, ast.AST] = {}
        self.outputs = []
        for o in outputs:
            try:
                tree = ast.parse(o.expr.strip(), mode="eval").body
            except SyntaxError as exc:
                raise CalculatorError("%s: %s" % (o.name, exc.msg)) from None
            tree = _Inline(defs).visit(tree)
            _check(ast.Expression(tree), set(self.var_index) | consts | set(defs))
            o.tree = tree
            o.time_dependent = self._time_dependent(tree)
            self.outputs.append(o)
            for alias in {_ident(o.name), o.name.replace(" ", "_"), o.name.replace(" ", "")}:
                if alias.isidentifier():
                    defs[alias] = tree

    # ---- analysis ------------------------------------------------------------
    def _time_dependent(self, node) -> bool:
        if isinstance(node, ast.Call) and node.func.id in _TIME:
            return False
        if isinstance(node, ast.Name):
            return node.id in self.var_index or node.id == "t"
        return any(self._time_dependent(c) for c in ast.iter_child_nodes(node))

    @property
    def time_dependent(self) -> bool:
        return any(o.time_dependent for o in self.outputs)

    # ---- data ----------------------------------------------------------------
    def _frame(self, f: int) -> np.ndarray:
        arr = self._frame_cache.get(f)
        if arr is None:
            if len(self._frame_cache) > 4:
                self._frame_cache.pop(next(iter(self._frame_cache)))
            arr = np.asarray(self.slf.read_frame(f), dtype=np.float64)
            self._frame_cache[f] = arr
        return arr

    def _need_ops(self, fname):
        if self.ops is None:
            raise CalculatorError("%s needs a triangular mesh" % fname)
        return self.ops

    # ---- evaluation ----------------------------------------------------------
    def evaluate(self, o: Output, frame: int) -> np.ndarray:
        v = self._eval(o.tree, frame)
        return np.broadcast_to(np.asarray(v, dtype=np.float64), (self.h.npoin,)).copy()

    def _eval(self, node, f: int):
        ev = self._eval
        if isinstance(node, ast.Constant):
            return float(node.value)
        if isinstance(node, ast.Name):
            n = node.id
            if n in self.var_index:
                return self._frame(f)[self.var_index[n]]
            if n == "x":
                return self.h.x
            if n == "y":
                return self.h.y
            if n == "t":
                return float(self.times[f]) if self.nframes else 0.0
            if n == "area":
                return self._need_ops("area").node_area
            if n == "g":
                return G
            if n == "pi":
                return np.pi
            if n == "e":
                return np.e
            raise CalculatorError("unknown name %r" % n)
        if isinstance(node, ast.UnaryOp):
            a = ev(node.operand, f)
            if isinstance(node.op, ast.USub):
                return -a
            if isinstance(node.op, ast.UAdd):
                return +a
            return np.logical_not(a)
        if isinstance(node, ast.BinOp):
            a, b = ev(node.left, f), ev(node.right, f)
            op = node.op
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                if isinstance(op, ast.Add):
                    return a + b
                if isinstance(op, ast.Sub):
                    return a - b
                if isinstance(op, ast.Mult):
                    return a * b
                if isinstance(op, ast.Div):
                    return np.true_divide(a, b)
                if isinstance(op, ast.Pow):
                    return np.power(a, b)
                if isinstance(op, ast.Mod):
                    return np.mod(a, b)
                if isinstance(op, ast.FloorDiv):
                    return np.floor_divide(a, b)
                if isinstance(op, ast.BitAnd):
                    return np.logical_and(a, b)
                if isinstance(op, ast.BitOr):
                    return np.logical_or(a, b)
        if isinstance(node, ast.Compare):
            left = ev(node.left, f)
            res = None
            for op, comp in zip(node.ops, node.comparators):
                right = ev(comp, f)
                r = {ast.Lt: np.less, ast.LtE: np.less_equal, ast.Gt: np.greater, ast.GtE: np.greater_equal,
                     ast.Eq: np.equal, ast.NotEq: np.not_equal}[type(op)](left, right)
                res = r if res is None else np.logical_and(res, r)
                left = right
            return res
        if isinstance(node, ast.BoolOp):
            vals = [ev(v, f) for v in node.values]
            fn = np.logical_and if isinstance(node.op, ast.And) else np.logical_or
            out = vals[0]
            for v in vals[1:]:
                out = fn(out, v)
            return out
        if isinstance(node, ast.IfExp):
            return np.where(ev(node.test, f), ev(node.body, f), ev(node.orelse, f))
        if isinstance(node, ast.Call):
            return self._call(node, f)
        raise CalculatorError("cannot evaluate %s" % type(node).__name__)

    def _args(self, node, count):
        if len(node.args) not in (count if isinstance(count, tuple) else (count,)):
            raise CalculatorError("%s takes %s argument(s)" % (node.func.id, count))
        return node.args

    def _nodal(self, v):
        return np.broadcast_to(np.asarray(v, dtype=np.float64), (self.h.npoin,))

    def _call(self, node, f: int):
        fn = node.func.id
        if fn in _NP:
            args = [self._eval(a, f) for a in node.args]
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                return _NP[fn](*args)
        if fn in _SPATIAL1:
            (a,) = self._args(node, 1)
            ops = self._need_ops(fn)
            v = self._nodal(self._eval(a, f))
            if fn == "ddx":
                return ops.ddx(v)
            if fn == "ddy":
                return ops.ddy(v)
            if fn == "grad":
                gx, gy = ops.grad(v)
                return np.hypot(gx, gy)
            if fn == "laplacian":
                return ops.laplacian(v)
            return ops.smooth(v)
        if fn in _SPATIAL2:
            a, b = self._args(node, 2)
            u, v = self._nodal(self._eval(a, f)), self._nodal(self._eval(b, f))
            if fn == "mag":
                return np.hypot(u, v)
            if fn in ("bearing", "direction"):
                return (np.degrees(np.arctan2(u, v)) + 360.0) % 360.0
            ops = self._need_ops(fn)
            if fn == "div":
                return ops.div(u, v)
            if fn in ("curl", "vorticity"):
                return ops.curl(u, v)
            if fn == "strain":
                return ops.strain(u, v)
            return ops.okubo_weiss(u, v)
        if fn in _TIME:
            (a,) = self._args(node, 1)
            key = id(node)
            if key not in self._red_cache:
                self._red_cache[key] = self._reduce(fn, a)
            return self._red_cache[key]
        if fn == "ddt":
            (a,) = self._args(node, 1)
            n = self.nframes
            if n < 2:
                return np.zeros(self.h.npoin)
            i0, i1 = max(f - 1, 0), min(f + 1, n - 1)
            dt = self.times[i1] - self.times[i0]
            if dt == 0:
                return np.zeros(self.h.npoin)
            return (self._nodal(self._eval(a, i1)) - self._nodal(self._eval(a, i0))) / dt
        if fn == "wet":
            args = self._args(node, (0, 1))
            eps = float(self._eval(args[0], f)) if args else self.dry
            if self.depth_index is None:
                raise CalculatorError("wet() needs a WATER DEPTH variable in the file")
            return (self._frame(f)[self.depth_index] > eps).astype(np.float64)
        raise CalculatorError("unknown function %r" % fn)

    def _reduce(self, fn, arg) -> np.ndarray:
        n = self.nframes
        if n == 0:
            raise CalculatorError("the file has no frame")
        if fn == "tfirst":
            return self._nodal(self._eval(arg, 0)).copy()
        if fn == "tlast":
            return self._nodal(self._eval(arg, n - 1)).copy()
        acc = None
        acc2 = None
        tpk = None
        prev = None
        for i in range(n):
            v = self._nodal(self._eval(arg, i)).astype(np.float64)
            if fn == "tmax" or fn == "tpeak":
                if acc is None:
                    acc, tpk = v.copy(), np.full(v.shape, self.times[i])
                else:
                    up = v > acc
                    acc = np.where(up, v, acc)
                    tpk = np.where(up, self.times[i], tpk)
            elif fn == "tmin":
                acc = v.copy() if acc is None else np.minimum(acc, v)
            elif fn in ("tmean", "tsum", "tstd"):
                acc = v.copy() if acc is None else acc + v
                if fn == "tstd":
                    acc2 = v * v if acc2 is None else acc2 + v * v
            elif fn == "tint":
                if acc is None:
                    acc = np.zeros_like(v)
                else:
                    acc = acc + 0.5 * (v + prev) * (self.times[i] - self.times[i - 1])
                prev = v
        if fn == "tpeak":
            return tpk
        if fn == "tmean":
            return acc / n
        if fn == "tstd":
            m = acc / n
            return np.sqrt(np.maximum(acc2 / n - m * m, 0.0))
        return acc


def describe_variables(names: Sequence[str]) -> List[Tuple[str, str]]:
    """(symbol to type, variable name) for the help text of a file."""
    vi = variable_names(names)
    out = []
    for i, n in enumerate(names):
        keys = [k for k, j in vi.items() if j == i]
        short = min(keys, key=len) if keys else n
        out.append((short, n.strip()))
    return out


def catalogue_symbols() -> List[Tuple[str, str]]:
    return [(v.symbol, v.label) for v in CATALOGUE if v.symbol]


@dataclass
class RunReport:
    outputs: List[Output] = field(default_factory=list)
    frames: List[int] = field(default_factory=list)
    ranges: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    renamed: Dict[str, str] = field(default_factory=dict)


def vector_file_names(names: Sequence[str], taken: Sequence[str] = ()) -> Dict[str, str]:
    """New names for output pairs "<P> X" / "<P> Y", so that QGIS reads them as one vector.

    QGIS (MDAL) pairs only "... U / ... V", "... ALONG X / Y" and "... SUIVANT X / Y";
    "SLOPE X" and "SLOPE Y" would stay two separate maps without arrows. They are
    written as "SLOPE ALONG X/Y", or "SLOPE U/V" when that is longer than the 16
    characters of a SELAFIN name.
    """
    upper = [n.strip().upper() for n in names]
    used = set(upper) | {t.strip().upper() for t in taken}
    out: Dict[str, str] = {}
    for n in upper:
        if not n.endswith(" X") or n.endswith((" ALONG X", " SUIVANT X")):
            continue
        p = n[:-2].rstrip()
        if not p or p + " Y" not in upper:
            continue
        for a, b in ((p + " ALONG X", p + " ALONG Y"), (p + " U", p + " V")):
            if len(a) <= 16 and a not in used and b not in used:
                out[n], out[p + " Y"] = a, b
                used.update((a, b))
                break
    return out


def run_calculator(src: str, dst: str, text: str, keep: bool = True, frames: Optional[Sequence[int]] = None,
                   progress: Optional[Callable[[float], bool]] = None, dry: float = 0.01) -> RunReport:
    """Write ``dst``: the variables of ``src`` (when ``keep``) plus one variable per output line."""
    from .convert import writer_like
    from .selafin import SelafinFile

    outs = parse_outputs(text)
    with SelafinFile(src) as slf:
        calc = Calculator(slf, outs, dry=dry)
        h = slf.header
        existing = [n.strip().upper() for n in h.var_names] if keep else []
        clash = [o.name for o in outs if o.name in existing]
        if clash:
            raise CalculatorError("these names already exist in the file: %s (choose other names or do not keep "
                                  "the original variables)" % ", ".join(clash))
        n = slf.nframes
        if frames is None or not len(frames):
            frames = list(range(n)) if (calc.time_dependent or keep) else [n - 1]
        if not calc.time_dependent and not keep:
            frames = [frames[-1]]
        renamed = vector_file_names([o.name for o in outs], existing)
        names = (list(h.var_names) if keep else []) + [renamed.get(o.name, o.name) for o in outs]
        units = (list(h.var_units) if keep else []) + [o.unit for o in outs]
        rep = RunReport(outputs=outs, frames=list(frames), renamed=renamed)
        lo = {o.name: np.inf for o in outs}
        hi = {o.name: -np.inf for o in outs}
        with writer_like(h, dst, names, units, title=(h.title.strip() + " + calculator")[:72]) as w:
            for k, fi in enumerate(frames):
                vals = [np.asarray(row, dtype=np.float64) for row in slf.read_frame(fi)] if keep else []
                for o in outs:
                    v = calc.evaluate(o, fi)
                    fin = v[np.isfinite(v)]
                    if fin.size:
                        lo[o.name] = min(lo[o.name], float(fin.min()))
                        hi[o.name] = max(hi[o.name], float(fin.max()))
                    vals.append(np.where(np.isfinite(v), v, 0.0))
                w.write_frame(float(slf.times[fi]) if n else 0.0, vals)
                if progress and progress((k + 1) / len(frames)):
                    break
        rep.ranges = {o.name: (lo[o.name], hi[o.name]) for o in outs}
    return rep
