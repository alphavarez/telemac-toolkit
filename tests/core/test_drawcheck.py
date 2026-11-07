"""Checks of a drawing against the element size: duplicates, spikes, close corners, narrow gaps, unsnapped ends."""

import numpy as np

from telemac_toolkit.mesher.drawcheck import clean_vertices, narrow_gaps, near_ends, vertex_problems


def const(h):
    return lambda x, y: np.full(np.shape(x), float(h))


SQUARE = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 50.0], [0.0, 50.0]])


def codes(problems):
    return sorted(p.code for p in problems)


def test_a_clean_outline_has_no_problem():
    t = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    circle = np.stack([50 + 10 * np.cos(t), 25 + 10 * np.sin(t)], axis=1)
    assert vertex_problems(SQUARE, True, const(5.0)) == []
    assert vertex_problems(circle, True, const(1.0)) == []
    assert narrow_gaps([SQUARE, circle], ["domain", "island"], const(5.0)) == []


def test_near_duplicates_and_spikes_found_and_cleaned():
    ring = np.array([[0.0, 0.0], [50.0, 0.0], [50.02, 0.0], [100.0, 0.0], [100.0, 25.0], [100.4, 25.02], [100.0, 25.04],
                     [100.0, 50.0], [0.0, 50.0]])
    p = vertex_problems(ring, True, const(5.0))
    assert codes(p) == ["duplicate", "spike"]
    assert all(q.fixable for q in p)
    dup = [q for q in p if q.code == "duplicate"][0]
    assert abs(dup.distance - 0.02) < 1e-9 and dup.size == 5.0
    clean = clean_vertices(ring, True, const(5.0))
    assert vertex_problems(clean, True, const(5.0)) == []
    assert len(clean) == 6  # the duplicate, the spike, and the vertex its removal left 4 cm from (100, 25)
    # an open line keeps its end vertices
    line = np.array([[0.0, 0.0], [10.0, 0.0], [10.001, 0.0]])
    out = clean_vertices(line, False, const(1.0))
    assert len(out) == 2 and np.allclose(out[-1], [10.001, 0.0])


def test_close_corners_and_dense_breaklines():
    # a notch 1 m wide and 2 m deep: corners the mesher keeps, 1 m apart, with 5 m elements
    notch = np.array([[0.0, 0.0], [50.0, 0.0], [50.0, 2.0], [51.0, 2.0], [51.0, 0.0], [100.0, 0.0], [100.0, 50.0],
                      [0.0, 50.0]])
    p = vertex_problems(notch, True, const(5.0))
    assert codes(p) == ["corners"] and not p[0].fixable and abs(p[0].distance - 1.0) < 1e-12
    # a 0.3 m notch is smaller than half an element: the mesher smooths it over, no corners there
    small = np.array([[0.0, 0.0], [50.0, 0.0], [50.0, 0.3], [50.3, 0.3], [50.3, 0.0], [100.0, 0.0], [100.0, 50.0],
                      [0.0, 50.0]])
    assert vertex_problems(small, True, const(5.0)) == []
    brk = np.array([[0.0, 0.0], [10.0, 0.0], [10.2, 0.3], [20.0, 0.0]])
    assert codes(vertex_problems(brk, False, const(5.0), kind="breakline")) == ["breakline"]
    assert vertex_problems(brk, False, const(5.0), kind="line") == []


def test_islands_close_to_the_outline_and_necks():
    near = np.array([[10.0, 1.0], [20.0, 1.0], [20.0, 10.0], [10.0, 10.0]])  # 1 m from the bottom edge
    p = narrow_gaps([SQUARE, near], ["domain", "island A"], const(5.0))
    assert codes(p) == ["gap"] and abs(p[0].distance - 1.0) < 0.3 and "island A" in p[0].message
    far = near + [0.0, 10.0]
    assert narrow_gaps([SQUARE, far], ["domain", "island A"], const(5.0)) == []
    # a slot 1 m wide cut 30 m into the domain: the outline comes within 1 m of itself
    slot = np.array([[0.0, 0.0], [49.5, 0.0], [49.5, 30.0], [50.5, 30.0], [50.5, 0.0], [100.0, 0.0], [100.0, 50.0],
                     [0.0, 50.0]])
    p = narrow_gaps([slot], ["domain"], const(5.0))
    assert p and all(q.code == "neck" for q in p) and min(q.distance for q in p) <= 1.0 + 1e-9
    # the same slot with elements of 0.5 m is wide enough
    assert narrow_gaps([slot], ["domain"], const(0.5)) == []


def test_line_ends_short_of_the_outline_get_a_snap_target():
    ends = np.array([[30.0, 0.4], [60.0, 0.0], [70.0, 25.0]])
    p = near_ends(ends, ["breakline 1 start", "breakline 1 end", "guide line"], [SQUARE], const(5.0))
    assert len(p) == 1 and p[0].index == 0 and p[0].fixable
    assert np.allclose(p[0].target, [30.0, 0.0]) and abs(p[0].distance - 0.4) < 1e-12
