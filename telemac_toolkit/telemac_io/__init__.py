"""telemac_io: read and write TELEMAC files with NumPy only (no QGIS).

Part of the TELEMAC Toolkit for QGIS. Usable on its own:

    from telemac_io import SelafinFile
    with SelafinFile("r2d.slf") as slf:
        depth = slf.read("WATER DEPTH", frame=-1)
"""

from .selafin import (  # noqa: F401
    SelafinError,
    SelafinFile,
    SelafinHeader,
    SelafinWriter,
    precision_warning,
    read_header,
    write_selafin,
)
from .convert import extract_2d, subset_frames  # noqa: F401
from .cli import BoundaryTable, read_cli, write_cli  # noqa: F401
from .mesh import TriangleLocator, boundary_loops, mesh_statistics, sample_polyline  # noqa: F401
from . import variables  # noqa: F401

__version__ = "0.4.1"
