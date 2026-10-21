"""Core SELAFIN tests: real TELEMAC files, round trips, precision, edge cases."""

import datetime as dt
import os
import struct
import warnings

import numpy as np
import pytest

from telemac_io import SelafinError, SelafinFile, SelafinWriter, extract_2d, subset_frames, write_selafin
from telemac_io.selafin import precision_warning, read_header

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
R2D = os.path.join(DATA, "r2d_tidal_flats.slf")
R3D = os.path.join(DATA, "r3d_tidal_flats.slf")
BUMP = os.path.join(DATA, "r3d_bump.slf")
HIST = os.path.join(DATA, "r1d_tomsail.slf")


def tiny_mesh():
    # two triangles forming a unit square
    x = np.array([0.0, 1.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    ikle = np.array([[0, 1, 2], [0, 2, 3]])
    return x, y, ikle


# ---------------------------------------------------------------- real files
def test_read_2d_header():
    with SelafinFile(R2D) as s:
        h = s.header
        assert h.npoin == 648 and h.nelem == 1030 and h.ndp == 3
        assert s.var_names == ["VELOCITY U", "VELOCITY V", "WATER DEPTH", "FREE SURFACE", "BOTTOM"]
        assert s.nframes == 17
        assert not h.is_3d and not h.is_double and h.endian == ">"
        assert not h.is_partitioned
        assert not h.has_calendar_date  # TELEMAC default 1900-01-01
        np.testing.assert_allclose(s.times[:3], [0, 10000, 20000])


def test_free_surface_equals_depth_plus_bottom():
    with SelafinFile(R2D) as s:
        for fr in (0, 8, -1):
            h = s.read("WATER DEPTH", fr)
            b = s.read("BOTTOM", fr)
            z = s.read("FREE SURFACE", fr)
            np.testing.assert_allclose(z, h + b, atol=1e-4)


def test_read_3d_and_planes():
    with SelafinFile(R3D) as s:
        h = s.header
        assert h.is_3d and h.nplan == 21 and h.ndp == 6
        assert h.npoin2 * 21 == h.npoin
        assert h.ikle2.shape == (h.nelem // 20, 3)
        z0 = s.plane("ELEVATION Z", -1, 0)
        z1 = s.plane("ELEVATION Z", -1, -1)
        assert np.all(z1 >= z0 - 1e-6)
        avg = s.depth_average("VELOCITY U", -1)
        assert avg.shape == (h.npoin2,)
        u = s.read("VELOCITY U", -1).reshape(21, -1)
        lo, hi = u.min(axis=0), u.max(axis=0)
        ok = ~np.isnan(avg)
        assert np.all((avg[ok] >= lo[ok] - 1e-6) & (avg[ok] <= hi[ok] + 1e-6))


def test_history_file_with_point_elements():
    with SelafinFile(HIST) as s:
        assert s.header.ndp == 1
        assert s.header.has_calendar_date
        assert s.header.start_datetime == dt.datetime(2023, 9, 1)
        assert s.datetimes[-1] == dt.datetime(2023, 9, 1, 6)


def test_series_matches_frames():
    with SelafinFile(R2D) as s:
        ts = s.series("WATER DEPTH", [3, 400])
        for fr in range(s.nframes):
            np.testing.assert_allclose(ts[fr], s.read("WATER DEPTH", fr)[[3, 400]])
        w = np.array([0.2, 0.3, 0.5])
        nodes = np.array([1, 2, 3])
        ref = s.series("BOTTOM", nodes) @ w
        np.testing.assert_allclose(s.series_weighted("BOTTOM", nodes, w), ref)


def test_var_lookup_errors():
    with SelafinFile(R2D) as s:
        assert s.header.var_index("water depth") == 2
        with pytest.raises(KeyError):
            s.read("NOT A VARIABLE")
        with pytest.raises(IndexError):
            s.read(0, frame=99)


# ---------------------------------------------------------------- round trips
@pytest.mark.parametrize("double", [False, True])
@pytest.mark.parametrize("endian", [">", "<"])
def test_round_trip_is_byte_exact(tmp_path, double, endian):
    """Re-writing a TELEMAC file reproduces it byte for byte (same precision/endianness)."""
    with SelafinFile(R2D) as s:
        h = s.header
        out = tmp_path / "copy.slf"
        with SelafinWriter(out, h.var_names, h.var_units, h.x, h.y, h.ikle, ipobo=h.ipobo, title=h.title,
                           double=double, endian=endian, date=h.date, origin=h.origin, iparam=h.iparam) as w:
            for i in range(s.nframes):
                w.write_frame(s.times[i], s.read_frame(i))
        with SelafinFile(out) as c:
            assert c.header.is_double == double and c.header.endian == endian
            np.testing.assert_array_equal(c.header.ikle, h.ikle)
            np.testing.assert_array_equal(c.header.ipobo, h.ipobo)
            np.testing.assert_array_equal(c.header.iparam, h.iparam)
            np.testing.assert_allclose(c.x, h.x)
            np.testing.assert_allclose(np.asarray(c.values, float), np.asarray(s.values, float))
        if not double and endian == ">":
            assert open(out, "rb").read() == open(R2D, "rb").read()


def test_write_selafin_helper_and_dates(tmp_path):
    x, y, ikle = tiny_mesh()
    p = tmp_path / "t.slf"
    write_selafin(p, x, y, ikle, {"BOTTOM": [1, 2, 3, 4], "WATER DEPTH": np.ones((3, 4))[0]},
                  units={"BOTTOM": "M"}, date=dt.datetime(2026, 10, 20, 12, 0, 0))
    with SelafinFile(p) as s:
        assert s.header.start_datetime == dt.datetime(2026, 10, 20, 12)
        assert s.var_units == ["M", ""]
        np.testing.assert_allclose(s.read("BOTTOM", 0), [1, 2, 3, 4])
        assert (s.header.ipobo > 0).sum() == 4


def test_origin_offset_keeps_utm_precision(tmp_path):
    x, y, ikle = tiny_mesh()
    X = 512345.123 + x * 10.0
    Y = 2654321.987 + y * 10.0
    p = tmp_path / "utm.slf"
    write_selafin(p, X, Y, ikle, {"BOTTOM": np.zeros(4)})
    h = read_header(p)
    assert h.origin == (512345.0, 2654321.0)
    assert precision_warning(h) is None
    np.testing.assert_allclose(h.x, X, atol=1e-3)
    np.testing.assert_allclose(h.y, Y, atol=1e-3)
    # without an offset float32 loses ~0.25 m at this northing
    q = tmp_path / "noorigin.slf"
    write_selafin(q, X, Y, ikle, {"BOTTOM": np.zeros(4)}, origin=None)
    hq = read_header(q)
    assert precision_warning(hq) is not None
    assert np.abs(hq.y - Y).max() > 0.01


def test_incomplete_last_frame_is_ignored(tmp_path):
    p = tmp_path / "running.slf"
    data = open(R2D, "rb").read()
    open(p, "wb").write(data[:-500])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        s = SelafinFile(p)
    assert s.nframes == 16
    assert any("incomplete" in str(m.message) for m in w)
    # simulation writes the rest
    open(p, "ab").write(data[-500:])
    assert s.refresh() == 17
    s.close()


def test_not_a_selafin(tmp_path):
    p = tmp_path / "bad.slf"
    p.write_bytes(b"hello world, not a selafin file at all")
    with pytest.raises(SelafinError):
        SelafinFile(p)


def test_corrupted_marker(tmp_path):
    data = bytearray(open(R2D, "rb").read())
    data[84:88] = struct.pack(">i", 79)  # closing marker of title record
    p = tmp_path / "corrupt.slf"
    p.write_bytes(bytes(data))
    with pytest.raises(SelafinError):
        SelafinFile(p)


def test_extract_2d_average_and_plane(tmp_path):
    a = extract_2d(BUMP, str(tmp_path / "avg.slf"), "average")
    b = extract_2d(BUMP, str(tmp_path / "bot.slf"), "bottom")
    with SelafinFile(BUMP) as s3, SelafinFile(a) as sa, SelafinFile(b) as sb:
        assert not sa.header.is_3d and sa.header.npoin == s3.header.npoin2
        assert "FREE SURFACE" in sa.var_names and "BOTTOM" in sa.var_names
        np.testing.assert_allclose(sa.read("BOTTOM", -1), s3.plane("ELEVATION Z", -1, 0), atol=1e-5)
        np.testing.assert_allclose(sa.read("FREE SURFACE", -1), s3.plane("ELEVATION Z", -1, -1), atol=1e-5)
        np.testing.assert_allclose(sb.read("VELOCITY U", 1), s3.plane("VELOCITY U", 1, 0), atol=1e-6)
        assert sa.nframes == s3.nframes


def test_subset_frames(tmp_path):
    p = subset_frames(R2D, str(tmp_path / "sub.slf"), [0, -1], ["WATER DEPTH"])
    with SelafinFile(p) as s, SelafinFile(R2D) as r:
        assert s.nframes == 2 and s.var_names == ["WATER DEPTH"]
        np.testing.assert_allclose(s.read(0, 1), r.read("WATER DEPTH", -1))
