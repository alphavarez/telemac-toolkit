"""Mesh generation, built for TELEMAC and usable for other 2D solvers.

Unstructured triangles with refinement and constraints, structured channel
strips and regular grids, interpolation, boundary conditions and writers.

NumPy only; the Delaunay step comes from GEOS (through QGIS) or SciPy.
"""

from .channel import channel_mesh  # noqa: F401
from .engine import MeshInput, MeshResult, NODE_TYPES, Submesh, generate  # noqa: F401
from .sizefield import LineDensity, PointDensity, PolygonDensity, SizeField  # noqa: F401
