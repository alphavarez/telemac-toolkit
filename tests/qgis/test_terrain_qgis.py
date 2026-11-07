"""Terrain-driven meshing in a real (headless) QGIS: size raster, Generate mesh with a terrain tolerance,
breaklines from the DEM into the project, MMG optimisation and adaptation (skipped without mmg2d).

    QT_QPA_PLATFORM=offscreen MMG2D=/path/to/mmg2d_O3 python -m pytest tests/qgis/test_terrain_qgis.py
"""

import os

import numpy as np
import pytest
from qgis.core import QgsGeometry, QgsPointXY, QgsRasterLayer, QgsVectorLayer

from test_bc_qgis import _numbers, case  # noqa: F401  (fixture)
from test_mesher_qgis import CRS, X0, Y0, _add
from test_river_qgis import dem, qgis, terrain  # noqa: F401  (fixtures)


def _mmg():
    from telemac_toolkit.mesher.mmg import find_mmg2d

    return find_mmg2d()


def _min_angle(path):
    from telemac_toolkit.mesher.engine import min_angles
    from telemac_toolkit.telemac_io import SelafinFile

    with SelafinFile(path) as s:
        h = s.header
        return float(min_angles(h.x, h.y, h.ikle).min()), h.npoin


def test_size_raster_and_terrain_tolerance(case, dem, tmp_path):  # noqa: F811
    import processing

    d, path, mesh, _cli = case
    r = processing.run("telemactoolkit:terrainsizes", {"DEM": dem, "PROJECT": path, "TOLERANCE": 0.05, "MIN_SIZE": 2.0,
                                                       "MAX_SIZE": 30.0, "OUTPUT": str(tmp_path / "sizes.tif")})
    lyr = QgsRasterLayer(r["OUTPUT"], "sizes")
    assert lyr.isValid()
    st = lyr.dataProvider().bandStatistics(1)
    assert 1.9 < st.minimumValue < 4.0 and st.maximumValue == pytest.approx(30.0, rel=1e-3)
    base_angle, base_n = _min_angle(mesh)
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "TERRAIN_TOLERANCE": 0.05,
                                                       "TERRAIN_MIN_SIZE": 2.0, "OUTPUT": str(d / "terrain.slf"), "LOAD": False})
    a, n = _min_angle(g["OUTPUT"])
    assert n > 1.3 * base_n  # refined along the banks and the valley sides
    assert a > 20.0  # the floor next to the structured channel keeps the elements there well shaped
    assert set(_numbers(g["OUTPUT_CLI"], g["OUTPUT"])) == {"inflow", "outflow"}


def test_breaklines_skip_structured_channels(case, dem):  # noqa: F811
    import processing

    _d, path, _mesh, _cli = case
    # the river project's structured strip already covers the banks: no line is drawn on it or next to it
    r = processing.run("telemactoolkit:dembreaklines", {"DEM": dem, "PROJECT": path, "TOLERANCE": 0.05, "SIZE": 30.0,
                                                        "OUTPUT": "memory:"})
    for f in r["OUTPUT"].getFeatures():
        pts = np.array([[p.x() - X0, p.y() - Y0] for p in f.geometry().vertices()])
        assert np.min(np.abs(pts[:, 1] - 70 * np.sin(2 * np.pi * pts[:, 0] / 900))) > 25.0


def test_breaklines_into_the_project(dem, tmp_path):  # noqa: F811
    import processing
    from qgis.core import QgsCoordinateReferenceSystem
    from telemac_toolkit.qgis_bridge import meshing as M

    # the same valley without a structured channel: the lines follow the bank tops and feet
    path = str(tmp_path / "valley.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), settings={"default_size": "30"})
    _add(lay["domain"], QgsGeometry.fromPolygonXY([[QgsPointXY(X0 + a, Y0 + b) for a, b in
                                                    [(0, -250), (2000, -250), (2000, 250), (0, 250), (0, -250)]]]))
    d = tmp_path
    before = 0
    r = processing.run("telemactoolkit:dembreaklines", {"DEM": dem, "PROJECT": path, "TOLERANCE": 0.05, "SIZE": 30.0,
                                                        "ADD": True, "MODE": 0, "OUTPUT": "memory:"})
    kinds = {f["kind"] for f in r["OUTPUT"].getFeatures()}
    assert kinds == {"crest", "toe"}
    # the bank tops of the synthetic valley are 26 m from the axis, the bank feet 18 m
    for f in r["OUTPUT"].getFeatures():
        pts = np.array([[p.x() - X0, p.y() - Y0] for p in f.geometry().vertices()])
        off = np.median(np.abs(pts[:, 1] - 70 * np.sin(2 * np.pi * pts[:, 0] / 900)))
        if f["length"] > 100:
            assert (abs(off - 26) < 3) if f["kind"] == "crest" else (abs(off - 18) < 3), (f["kind"], off)
    lay = M.open_mesh_project(path)
    assert lay["constraint_lines"].featureCount() > before
    sizes = [f["size"] for f in lay["constraint_lines"].getFeatures()]
    assert any(s is not None and 6.0 < s < 10.0 for s in sizes)  # bank tops and feet are 8 m apart
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "OUTPUT": str(d / "brk.slf"), "LOAD": False})
    a, n = _min_angle(g["OUTPUT"])
    assert a > 20.0, a


@pytest.mark.skipif(_mmg() is None, reason="mmg2d not found (set MMG2D)")
def test_mmg_modes(case, dem, tmp_path):  # noqa: F811
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile
    from telemac_toolkit.telemac_io.selafin import SelafinWriter

    d, path, mesh, cli = case
    lay = M.open_mesh_project(path)
    left = next(f for f in lay["channels"].getFeatures() if f["role"] == "left").geometry().asPolyline()
    right = next(f for f in lay["channels"].getFeatures() if f["role"] == "right").geometry().asPolyline()
    keep = QgsVectorLayer("Polygon?crs=%s" % CRS, "keep", "memory")
    _add(keep, QgsGeometry.fromPolygonXY([list(left) + list(reversed(right)) + [left[0]]]))
    base = _numbers(cli, mesh)
    with SelafinFile(mesh) as s:
        n0 = s.header.npoin
    # optimise
    o = processing.run("telemactoolkit:mmgadapt", {"INPUT": mesh, "MODE": 0, "KEEP": keep, "MESH_CRS": CRS,
                                                   "OUTPUT": str(tmp_path / "opt.slf")})
    assert set(_numbers(o["OUTPUT_CLI"], o["OUTPUT"])) == set(base)
    # terrain, with the channel kept and the bottom from the DEM
    t = processing.run("telemactoolkit:mmgadapt", {"INPUT": mesh, "MODE": 1, "DEM": dem, "TOLERANCE": 0.05, "HMIN": 2.0, "HMAX": 30.0,
                                                   "KEEP": keep, "MESH_CRS": CRS, "OUTPUT": str(tmp_path / "terrain.slf")})
    with SelafinFile(t["OUTPUT"]) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        assert h.npoin > 1.3 * n0
        assert np.max(np.abs(z - terrain(h.x - X0, h.y - Y0))) < 0.05  # sampled from the DEM at every node
    a, _n = _min_angle(t["OUTPUT"])
    assert a > 15.0
    assert set(_numbers(t["OUTPUT_CLI"], t["OUTPUT"])) == set(base)
    # a result with a sharp front across the valley at x = 1000 m: refined there
    with SelafinFile(mesh) as s:
        h = s.header
        x, y, ik = h.x.copy(), h.y.copy(), h.ikle.copy()
    res = str(tmp_path / "res.slf")
    fs = 5.0 + np.tanh((x - X0 - 1000.0) / 20.0)
    with SelafinWriter(res, ["FREE SURFACE"], ["M"], x, y, ik) as w:
        w.write_frame(0.0, fs[None, :])
        w.write_frame(60.0, (fs + 0.1)[None, :])
    r = processing.run("telemactoolkit:mmgadapt", {"INPUT": mesh, "MODE": 2, "RESULT": res, "VARIABLE": "FREE SURFACE", "WHEN": 1,
                                                   "TOLERANCE": 0.005, "HMIN": 2.0, "HMAX": 30.0, "KEEP": keep, "MESH_CRS": CRS,
                                                   "OUTPUT": str(tmp_path / "adapt.slf")})
    with SelafinFile(r["OUTPUT"]) as s:
        h = s.header
        near = np.abs(h.x - X0 - 1000.0) < 60
        far = np.abs(h.x - X0 - 400.0) < 60
        assert near.sum() > 1.5 * far.sum()
    assert os.path.exists(r["OUTPUT_CLI"])
