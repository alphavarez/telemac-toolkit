"""Cross-check against an independent reader (xarray-selafin), when installed."""
import os

import numpy as np
import pytest

from telemac_io import SelafinFile

xr = pytest.importorskip("xarray")
pytest.importorskip("xarray_selafin")

DATA = os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.mark.parametrize("name", ["r2d_tidal_flats.slf", "r3d_tidal_flats.slf", "r3d_bump.slf", "r1d_tomsail.slf"])
def test_same_values_as_xarray_selafin(name):
    path = os.path.join(DATA, name)
    ds = xr.open_dataset(path, engine="selafin")
    with SelafinFile(path) as s:
        assert ds.sizes["time"] == s.nframes
        ours = np.asarray(s.values, dtype=np.float64)          # (t, v, n)
        name_to_key = {v[0].strip().upper(): k for k, v in ds.attrs["variables"].items()}
        for iv, vname in enumerate(s.var_names):
            theirs = np.asarray(ds[name_to_key[vname.upper()]].values, dtype=np.float64).reshape(s.nframes, -1)
            np.testing.assert_array_equal(ours[:, iv, :], theirs, err_msg=vname)
        np.testing.assert_allclose(s.header.x2, np.asarray(ds["x"].values, float).ravel(), atol=1e-6)
        np.testing.assert_array_equal(s.header.ikle2 + 1, np.asarray(ds.attrs["ikle2"]))
