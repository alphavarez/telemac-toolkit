"""Cross-check against TELEMAC's own Python reader (data_manip), when TELEMAC is installed.

Run inside a TELEMAC environment, e.g. ``micromamba run -n telemac python -m pytest tests/core``.
"""
import glob
import os

import numpy as np
import pytest

from telemac_io import SelafinFile

Selafin = pytest.importorskip("data_manip.formats.selafin").Selafin

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
FILES = sorted(glob.glob(os.path.join(DATA, "*.slf")) + glob.glob(os.path.join(DATA, "demo_valley", "*.slf")))


@pytest.mark.parametrize("path", FILES, ids=[os.path.basename(f) for f in FILES])
def test_same_as_telemac_reader(path):
    ref = Selafin(path)
    with SelafinFile(path) as ours:
        assert [n.strip() for n in ref.varnames] == ours.var_names
        np.testing.assert_allclose(ref.tags["times"], ours.times)
        for fr in sorted({0, ours.nframes // 2, ours.nframes - 1}):
            np.testing.assert_array_equal(np.asarray(ref.get_values(fr), float), ours.read_frame(fr))
        # TELEMAC's tools keep coordinates relative to IPARAM(3,4); QGIS and this library add the origin
        mx = np.asarray(ref.meshx, float)  # one plane for 3D files
        np.testing.assert_allclose(mx, ours.header.x_stored[: mx.size], atol=1e-6)
