"""Channels composed from segments, drawn centrelines, imported files and cross-sections."""

import math

import numpy as np
import pytest

from telemac_toolkit.mesher.composer import compose, from_centreline, parse_section, read_xy_file
from telemac_toolkit.mesher.interp import FormulaSource


def _bed(comp, x, y):
    return FormulaSource(comp.formula, "bed", comp.axis, comp.half_width).sample(np.asarray(x, float), np.asarray(y, float))


def test_s_bend_ends_where_the_geometry_says():
    c = compose("straight 3\nbend left 2 90\nstraight 1\nbend right 2 90\nstraight 3", 1.0, 0.05)
    assert len(c.pieces) == 1 and c.length == pytest.approx(3 + math.pi + 1 + math.pi + 3, abs=1e-6)
    end = c.axis[-1]
    # left 90 then right 90: heading back to +x, shifted 2 + 1 + 2 = 5 m up and 2 + 2 = 4 m along
    assert end[0] == pytest.approx(3 + 2 + 2 + 3, abs=1e-3) and end[1] == pytest.approx(5.0, abs=1e-3)
    assert c.pieces[0].heading[-1] == pytest.approx(0.0, abs=1e-6)


def test_expansion_and_elbow_make_pieces():
    c = compose("width 0.5\nstraight 2\nexpand 1.5\nstraight 4", 1.0, 0.05)
    assert [j["kind"] for j in c.joints] == ["step"] and len(c.pieces) == 2
    assert c.pieces[0].widths.max() == pytest.approx(0.5) and c.pieces[1].widths.min() == pytest.approx(1.5)
    c = compose("straight 3; elbow right 90; straight 3", 1.0, 0.05)
    assert c.joints[0]["kind"] == "elbow" and c.joints[0]["angle"] == pytest.approx(-90)
    assert c.axis[-1] == pytest.approx([3.0, -3.0], abs=1e-9)


def test_taper_changes_the_width_linearly():
    c = compose("width 1\nstraight 1\ntaper 2 2\nstraight 1", 1.0, 0.1)
    p = c.pieces[0]
    k = np.argmin(np.abs(p.s - 2.0))
    assert p.widths[k] == pytest.approx(1.5)
    assert np.interp(2.0, *c.half_width) == pytest.approx(0.75)


def test_bed_slope_drop_and_sections():
    c = compose("width 1\nslope 0.01\nstraight 2\ndrop 0.05\nsection v 0.2\nstraight 2", 1.0, 0.05, z0=1.0)
    z = _bed(c, [1.0, 3.0, 3.0], [0.0, 0.0, 0.5])
    assert z[0] == pytest.approx(0.99)  # 1 - 0.01 x 1
    assert z[1] == pytest.approx(1.0 - 0.02 - 0.05 - 0.01)  # the drop at s = 2, then the slope
    assert z[2] - z[1] == pytest.approx(0.2)  # the V section: 0.2 m higher at the wall (n = b = 0.5)


@pytest.mark.parametrize("words,n,expected", [
    ("rect", 0.4, 0.0),
    ("trapezoid 0.4 2 0.1", 0.4, 0.1),
    ("compound 0.4 1 0.1 0.1", 0.45, 0.1 + 0.1 * (0.45 - 0.3)),
    ("parabolic 0.1", 0.25, 0.1 * 0.25),
    ("asym 0.2 1 2 0.3", -0.3, 0.1),
    ("custom 0.1*(n/b)**4", 0.5, 0.1),
])
def test_section_formulas(words, n, expected):
    sec = parse_section(words.split() if not words.startswith("custom") else ["custom", words[7:]])
    c = compose("width 1; straight 4", 1.0, 0.1, section=sec)
    assert _bed(c, [2.0], [n])[0] == pytest.approx(expected, abs=1e-9)


def test_errors_say_what_is_wrong():
    with pytest.raises(ValueError, match="radius"):
        compose("bend left 0.4 90", 1.0, 0.05)
    with pytest.raises(ValueError, match="unknown segment"):
        compose("straigth 3", 1.0, 0.05)
    with pytest.raises(ValueError, match="missing numbers"):
        compose("bend left 2", 1.0, 0.05)
    with pytest.raises(ValueError, match="section"):
        parse_section(["trapezoid", "1"])


def test_drawn_centreline_corners_become_elbows(tmp_path):
    c = from_centreline([[0, 0], [5, 0], [5, 5]], 0.1, width=1.0)
    assert len(c.pieces) == 2 and c.joints[0]["angle"] == pytest.approx(90)
    # a smooth curve stays one piece
    t = np.linspace(0, math.pi / 2, 60)
    c = from_centreline(np.stack([3 * np.sin(t), 3 - 3 * np.cos(t)], 1), 0.1, width=1.0)
    assert len(c.pieces) == 1
    with pytest.raises(ValueError, match="too tightly"):
        from_centreline(np.stack([0.4 * np.sin(t), 0.4 - 0.4 * np.cos(t)], 1), 0.05, width=1.0)


def test_read_xy_file_parts_and_widths(tmp_path):
    p = tmp_path / "shape.csv"
    p.write_text("x,y,width\n0,0,1\n5,0,1\n10,0,2\n\nEND\n0 2\n1 2\n1 3\n0 2\n")
    parts = read_xy_file(str(p))
    assert len(parts) == 2 and parts[0].shape == (3, 3) and parts[0][2, 2] == 2.0
    c = from_centreline(parts[0][:, :2], 0.5, widths=parts[0][:, 2])
    assert c.pieces[0].widths[-1] == pytest.approx(2.0) and c.pieces[0].widths[0] == pytest.approx(1.0)
