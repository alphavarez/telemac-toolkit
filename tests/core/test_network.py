"""River networks (junctions), strips whose nodes across follow the width, bends and strip quality."""

import math

import numpy as np
import pytest

from telemac_toolkit.mesher.channel import channel_mesh, channel_strip, strip_quality
from telemac_toolkit.mesher.engine import min_angles
from telemac_toolkit.mesher.geometry import ring_area
from telemac_toolkit.mesher.meshdata import pairs_to_quads
from telemac_toolkit.mesher.network import _seg_hits_ring, river_network, split_network
from telemac_toolkit.mesher.river import _normals, reach_from_centreline
from telemac_toolkit.telemac_io.mesh import boundary_loops, signed_areas

ANG = math.radians(60)
TRIB = [[1000 - 700 * math.cos(ANG), 700 * math.sin(ANG)], [1003.0, 4.0]]  # ends 5 m off the main line


def test_split_network_cuts_the_main_line_and_joins_lines_end_to_end():
    pieces, points, _notes = split_network([[[0, 0], [600, 0]], [[600, 0], [2000, 0]], TRIB], tol=10.0)
    assert len(points) == 2 and len(pieces) == 3  # the end-to-end join is undone, the tributary splits the main
    j = [p for p in points if abs(p[0] - 1000) < 5]
    assert len(j) == 1 and abs(j[0][1]) < 1e-9  # the junction lies on the main line
    ends = sorted((tuple(np.round(c[0])), tuple(np.round(c[-1]))) for c, *_r in pieces)
    assert ((0.0, 0.0), (1003.0, 0.0)) in ends or ((0.0, 0.0), (1000.0, 0.0)) in ends
    with pytest.raises(ValueError, match="separate groups"):
        split_network([[[0, 0], [500, 0]], [[0, 300], [500, 300]]], tol=10.0)


def test_confluence_strips_stop_short_and_the_junction_keeps_its_corners():
    net = river_network([[[0, 0], [2000, 0]], TRIB], 10.0, widths=[40.0, 25.0])
    assert len(net.reaches) == 3 and len(net.junctions) == 1
    assert sorted(k for _r, _e, k in net.ends) == ["inflow", "inflow", "outflow"]
    J = net.junctions[0]
    # no strip section touches another branch's channel
    for r, e in J.members:
        st = net.reaches[r]
        k = 0 if e == 0 else -1
        for o, _e in J.members:
            if o != r:
                assert not _seg_hits_ring(st.left[k], st.right[k], net.full[o].outline())
    # the channel area: a simple ring bounded by the strip ends and the banks
    area = abs(ring_area(J.ring[:-1]))
    assert 40 * 25 < area < 40 * 200
    corners = {tuple(np.round(p, 6)) for r, e in J.members for p in
               (net.reaches[r].left[0 if e == 0 else -1], net.reaches[r].right[0 if e == 0 else -1])}
    ends = {tuple(np.round(b[i], 6)) for b in J.banks for i in (0, -1)}
    assert ends == corners  # every bank runs from a strip corner to a strip corner
    # the banks of the tributary meet the main river's bank (y = 20) at its corners
    pts = np.vstack(J.banks)
    assert np.sum(np.abs(pts[:, 1] - 20.0) < 1e-6) >= 2


def test_bifurcation_has_two_outflows():
    lines = [[[0, 0], [1000, 0]], [[1000, 0], [1800, 400]], [[1000, 0], [1800, -300]]]
    net = river_network(lines, 10.0, widths=[40.0, 30.0, 30.0])
    assert sorted(k for _r, _e, k in net.ends) == ["inflow", "outflow", "outflow"]
    assert len(net.junctions[0].banks) == 3


def test_short_reach_between_junctions_is_reported():
    lines = [[[0, 0], [2000, 0]], [[1000, 500], [1000, 0]], [[1015, -500], [1015, 0]]]
    with pytest.raises(ValueError, match="too short"):
        river_network(lines, 10.0, widths=[40.0, 30.0, 30.0], tol=5.0)


def test_strip_follows_the_width_with_rows_of_triangles():
    xs = np.linspace(0, 1000, 101)
    w = 20 + 60 * xs / 1000
    L, R, T = np.stack([xs, w / 2], 1), np.stack([xs, -w / 2], 1), np.stack([xs, 0 * xs], 1)
    s = channel_strip(L, R, T, interval=10.0, across_spacing=8.0)
    assert s.counts[0] == 3 and s.counts[-1] == 11 and np.all(s.counts % 2 == 1)
    assert np.max(np.abs(np.diff(s.counts))) == 2  # one interval on each side of the axis at a time
    assert min_angles(s.x, s.y, s.ikle).min() > 28.0
    assert len(boundary_loops(s.ikle, s.x, s.y)) == 1
    assert abs(np.abs(signed_areas(s.x, s.y, s.ikle)).sum() - 50000.0) < 1e-6
    assert len(s.ikle) > 2 * len(s.quads)  # the transition rows are triangles
    across = np.hypot(*(R - L).T) / (s.counts - 1)
    assert np.all(np.abs(across / 8.0 - 1.0) < 0.45)  # 1 to 2 intervals a side halves the spacing
    assert np.all(np.abs(across[s.counts >= 7] / 8.0 - 1.0) < 0.3)
    # constant counts: the same strip as before, cells as consecutive pairs
    x, y, ik = channel_mesh(L, R, T, cross_nodes=7, interval=10.0)
    c = channel_strip(L, R, T, cross_nodes=7, interval=10.0)
    assert np.array_equal(ik, c.ikle) and np.array_equal(pairs_to_quads(x, y, ik), c.quads)


def test_strip_quality_of_a_straight_channel():
    xs = np.linspace(0, 100, 11)
    q = strip_quality(np.stack([xs, xs * 0 + 10], 1), np.stack([xs, xs * 0 - 10], 1), np.full(11, 5))
    assert q["ratio_min"] == pytest.approx(2.0) and q["ratio_max"] == pytest.approx(2.0)
    assert q["skew_max"] == pytest.approx(0.0, abs=1e-9)


def test_tight_bend_turns_only_the_sections_there():
    r = reach_from_centreline([[0, 0], [1000, 0], [1000, 1000]], 5.0, width=60.0, smoothing=1.0)
    assert r.crossings().size == 0
    raw = _normals(r.axis)
    dev = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", r.normal, raw), -1, 1)))
    corner = int(np.argmin(np.hypot(*(r.axis - [1000, 0]).T)))
    far = np.abs(r.s - r.s[corner]) > 300
    assert dev[far].max() < 0.01 and dev.max() > 10.0
    gaps = np.minimum(np.hypot(*np.diff(r.left, axis=0).T), np.hypot(*np.diff(r.right, axis=0).T))
    assert gaps.min() >= 5.0 / 4.0 - 1e-9  # inner-bank cells at least the spacing / 4
