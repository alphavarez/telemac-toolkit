"""River reach helper in a real (headless) QGIS: DEM banks, domain, boundary lines, Generate mesh.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_river_qgis.py
"""

import numpy as np
import pytest
from qgis.core import QgsCoordinateReferenceSystem, QgsGeometry, QgsPointXY, QgsRasterLayer, QgsVectorLayer

from test_mesher_qgis import CRS, X0, Y0, _add, _app


def terrain(x, y):
    """A meandering river 3 m deep in a valley sloping 2 per mil downstream (+x); relative coordinates."""
    yc = 70.0 * np.sin(2 * np.pi * x / 900.0)
    d = np.abs(y - yc)
    floor = 10.0 - 0.002 * x + 0.006 * np.maximum(d - 25.0, 0.0) + 2.5e-5 * np.maximum(d - 25.0, 0.0) ** 2
    bank = np.clip((d - 18.0) / 8.0, 0.0, 1.0)
    return floor - 3.0 * (1.0 - bank)


@pytest.fixture(scope="module")
def qgis():
    return _app()


@pytest.fixture(scope="module")
def dem(qgis, tmp_path_factory):
    from osgeo import gdal, osr

    res = 2.0
    nx, ny = 1150, 450  # -150..2150 x -450..450
    xs = -150 + (np.arange(nx) + 0.5) * res
    ys = 450 - (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    path = str(tmp_path_factory.mktemp("river") / "valley.tif")
    ds = gdal.GetDriverByName("GTiff").Create(path, nx, ny, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((X0 - 150, res, 0, Y0 + 450, 0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32614)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(terrain(X, Y).astype(np.float32))
    ds = None
    lyr = QgsRasterLayer(path, "valley dem")
    assert lyr.isValid()
    return lyr


def _line_layer(lines):
    lyr = QgsVectorLayer("LineString?crs=%s" % CRS, "lines", "memory")
    for c in lines:
        _add(lyr, QgsGeometry.fromPolylineXY([QgsPointXY(X0 + a, Y0 + b) for a, b in c]))
    return lyr


def _liquid(cli_path):
    from telemac_toolkit.telemac_io import read_cli
    from telemac_toolkit.telemac_io.cli import liquid_boundaries

    t = read_cli(cli_path)
    return {t.comments[b["start"]]: (len(b["nodes"]), b["codes"]) for b in liquid_boundaries(t)}


def test_centreline_with_dem_to_mesh(dem, tmp_path):
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile

    rng = np.random.default_rng(3)
    xs = np.linspace(2000, 0, 60)  # drawn upstream: the DEM turns it round
    cl = np.stack([xs, 70 * np.sin(2 * np.pi * xs / 900) + rng.normal(0, 3, xs.size)], axis=1)
    path = str(tmp_path / "river.gpkg")
    r = processing.run("telemactoolkit:riverreach", {
        "NEW_PROJECT": path, "SOURCE": 0, "CENTRELINE": _line_layer([cl]), "WIDTH_MODE": 1, "DEM": dem,
        "BANKFULL_DEPTH": 2.5, "MAX_WIDTH": 120.0, "SPACING": 8.0, "FLOODPLAIN": 150.0, "Q": 20.0, "H": 7.0})
    lay = M.open_mesh_project(r["OUTPUT"])
    assert lay["domain"].crs() == QgsCoordinateReferenceSystem(CRS)
    roles = sorted(f["role"] for f in lay["channels"].getFeatures())
    assert roles == ["axis", "left", "right"]
    across = {f["nodes_across"] for f in lay["channels"].getFeatures()}.pop()
    assert across % 2 == 1 and across >= 5
    bnd = {f["name"]: (f["type"], f["value"]) for f in lay["boundaries"].getFeatures()}
    assert bnd == {"river inflow": ("Q", 20.0), "river outflow": ("H", 7.0)}
    # the axis now runs downstream (+x) and follows the river within a few metres
    ax = next(f for f in lay["channels"].getFeatures() if f["role"] == "axis").geometry().asPolyline()
    ax = np.array([[p.x() - X0, p.y() - Y0] for p in ax])
    assert ax[0, 0] < ax[-1, 0]
    assert np.max(np.abs(ax[:, 1] - 70 * np.sin(2 * np.pi * ax[:, 0] / 900))) < 3.0
    assert M.read_settings(path)["default_size"] == "32"

    out = str(tmp_path / "river.slf")
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "OUTPUT": out,
                                                       "OUTPUT_LIQ": str(tmp_path / "river.liq"), "LOAD": False})
    with SelafinFile(g["OUTPUT"]) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
        assert 3000 < h.npoin < 20000
        zt = terrain(h.x - X0, h.y - Y0)
        assert np.max(np.abs(z - zt)) < 0.35  # bilinear on a 2 m DEM, steep banks
    liq = _liquid(g["OUTPUT_CLI"])
    assert set(liq) == {"river inflow", "river outflow"}
    n_in, codes_in = liq["river inflow"]
    assert n_in == across - 2 and codes_in == (4, 5, 5)  # the channel without its bank nodes
    n_out, codes_out = liq["river outflow"]
    assert n_out > across and codes_out == (5, 4, 4)  # the whole downstream end of the valley
    text = open(str(tmp_path / "river.liq")).read()
    assert "20" in text and "7" in text

    # the quality-guaranteed method on the same river
    from telemac_toolkit.mesher.engine import min_angles

    out_q = str(tmp_path / "river_q.slf")
    gq = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "OUTPUT": out_q,
                                                        "METHOD": 1, "MIN_ANGLE": 30.0, "LOAD": False})
    with SelafinFile(gq["OUTPUT"]) as s:
        h = s.header
        a = min_angles(h.x, h.y, h.ikle)
    with SelafinFile(g["OUTPUT"]) as s:
        a0 = min_angles(s.header.x, s.header.y, s.header.ikle)
    # every triangle reaches 30 degrees except, possibly, one with an edge on the structured channel (fixed)
    assert a.min() >= 29.0 and a.min() >= a0.min() and np.sum(a < 30.0) <= 3
    assert set(_liquid(gq["OUTPUT_CLI"])) == {"river inflow", "river outflow"}


def test_bank_lines_in_an_existing_domain(dem, tmp_path):
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M

    path = str(tmp_path / "banks.gpkg")
    lay = M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS), settings={"default_size": "40"})
    ring = [(100, -300), (1900, -300), (1900, 300), (100, 300), (100, -300)]
    _add(lay["domain"], QgsGeometry.fromPolygonXY([[QgsPointXY(X0 + a, Y0 + b) for a, b in ring]]), name="valley")
    x = np.linspace(-50, 2050, 300)
    yc = 70 * np.sin(2 * np.pi * x / 900)
    sl = np.gradient(yc, x)
    nx_, ny_ = -sl / np.hypot(1, sl), 1 / np.hypot(1, sl)
    left = np.stack([x + 24 * nx_, yc + 24 * ny_], 1)
    right = np.stack([x - 24 * nx_, yc - 24 * ny_], 1)
    # the right bank drawn first and backwards: the helper sorts them out
    banks = _line_layer([right[::-1], left])
    processing.run("telemactoolkit:riverreach", {
        "PROJECT": path, "SOURCE": 1, "BANKS": banks, "DEM": dem, "SPACING": 10.0, "DOMAIN_MODE": 1, "NAME": "main"})
    lay = M.open_mesh_project(path)
    assert lay["domain"].featureCount() == 1  # kept, not replaced
    L = next(f for f in lay["channels"].getFeatures() if f["role"] == "left").geometry().asPolyline()
    xs = [p.x() - X0 for p in L]
    assert min(xs) > 100 and max(xs) < 1900  # cut to the cross-sections inside the domain
    assert {f["name"] for f in lay["boundaries"].getFeatures()} == {"main inflow", "main outflow"}
    out = str(tmp_path / "banks.slf")
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "RASTERS": [dem], "OUTPUT": out, "LOAD": False})
    liq = _liquid(g["OUTPUT_CLI"])
    assert set(liq) == {"main inflow", "main outflow"}
    assert liq["main inflow"][0] >= 3


def test_river_reach_errors(dem, tmp_path):
    import processing
    from qgis.core import QgsProcessingException

    two = _line_layer([[(0, 0), (500, 0)], [(0, 50), (500, 60)]])
    with pytest.raises(QgsProcessingException, match="separate groups"):
        processing.run("telemactoolkit:riverreach", {"NEW_PROJECT": str(tmp_path / "e1.gpkg"), "CENTRELINE": two, "WIDTH": 20.0})
    one = _line_layer([[(0, 0), (500, 0)]])
    with pytest.raises(QgsProcessingException, match="need a DEM"):
        processing.run("telemactoolkit:riverreach", {"NEW_PROJECT": str(tmp_path / "e2.gpkg"), "CENTRELINE": one, "WIDTH_MODE": 1})
