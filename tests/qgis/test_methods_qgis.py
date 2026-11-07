"""Meshing methods in a real (headless) QGIS: Gmsh through Generate mesh, Compare meshing methods (with a Blue
Kenue mesh), Blue Kenue line import; with TELEMAC_CMD a Gmsh mesh is also run in TELEMAC-2D.

    QT_QPA_PLATFORM=offscreen python -m pytest tests/qgis/test_methods_qgis.py
"""

import os

import numpy as np
import pytest
from qgis.core import QgsCoordinateReferenceSystem

from test_bc_qgis import case, _numbers  # noqa: F401  (fixture)
from test_mesher_qgis import CRS
from test_river_qgis import dem, qgis  # noqa: F401  (fixtures)


def _gmsh():
    from telemac_toolkit.mesher import gmsh_engine

    return gmsh_engine.available()


def _write_t3s(path, x, y, z, ik):
    with open(path, "w") as fh:
        fh.write(":FileType t3s  ASCII  EnSim 1.0\n:NodeCount %d\n:ElementCount %d\n:EndHeader\n" % (len(x), len(ik)))
        for a, b, c in zip(x, y, z):
            fh.write("%.4f %.4f %.4f\n" % (a, b, c))
        for t in ik + 1:
            fh.write("%d %d %d\n" % tuple(t))


@pytest.mark.skipif(not _gmsh(), reason="gmsh is not installed")
def test_gmsh_methods_through_generate_mesh(case):  # noqa: F811
    import processing
    from telemac_toolkit.mesher.engine import min_angles
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, mesh, cli = case
    base = _numbers(cli, mesh)
    for k, key in enumerate(M.METHOD_KEYS):
        if not key.startswith("gmsh_"):
            continue
        out = str(d / ("%s.slf" % key))
        g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "METHOD": k, "OUTPUT": out, "LOAD": False})
        with SelafinFile(g["OUTPUT"]) as s:
            h = s.header
            a = min_angles(h.x, h.y, h.ikle)
        assert set(_numbers(g["OUTPUT_CLI"], g["OUTPUT"])) == set(base), key
        assert a.min() > (10.0 if key == "gmsh_quads" else 25.0), key


def test_compare_methods_and_a_blue_kenue_mesh(case, dem, tmp_path):  # noqa: F811
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io import SelafinFile

    d, path, mesh, _cli = case
    with SelafinFile(mesh) as s:
        h = s.header
        z = s.read("BOTTOM", 0)
    bk = str(tmp_path / "from_blue_kenue.t3s")
    _write_t3s(bk, h.x, h.y, z, h.ikle)
    keys = ["balanced", "quality"] + (["gmsh_frontal"] if _gmsh() else [])
    from telemac_toolkit.tt_processing.compare_algorithms import method_choices

    idx = [i for i, (k, _l) in enumerate(method_choices()) if k in keys]
    report = str(tmp_path / "methods.html")
    r = processing.run("telemactoolkit:comparemeshing", {
        "PROJECT": path, "METHODS": idx, "MESHES": [bk], "DEM": dem, "LEVEL": 6.0, "VELOCITY": 1.5,
        "TABLE": "memory:", "REPORT": report, "MESH_FOLDER": str(tmp_path / "meshes")})
    rows = {f["mesh"]: f for f in r["TABLE"].getFeatures()}
    assert len(rows) == len(keys) + 1 and "from_blue_kenue.t3s" in rows
    for name, f in rows.items():
        assert f["nodes"] > 3000 and 0 < f["bed_rms"] < 0.5 and f["dt_courant1"] > 0.1, name
        assert f["outline_dev_max"] < 1.0, name  # every mesh follows the drawn domain
    q = next(f for n, f in rows.items() if n.startswith("Quality"))
    assert q["angle_min"] >= 29.0
    html = open(report, encoding="utf-8").read()
    assert "class='best'" in html and "Bed error vs DEM, RMS" in html and "from_blue_kenue.t3s" in html
    assert os.path.exists(str(tmp_path / "meshes" / "balanced.slf"))
    # the Blue Kenue mesh is the plugin's own balanced mesh: same measures
    b = next(f for n, f in rows.items() if n.startswith("Balanced"))
    assert abs(b["nodes"] - rows["from_blue_kenue.t3s"]["nodes"]) == 0
    assert M.read_mesh_file(bk)[0].z is not None


def test_blue_kenue_lines_into_a_project(qgis, tmp_path):  # noqa: F811
    import processing
    from telemac_toolkit.qgis_bridge import meshing as M
    from telemac_toolkit.telemac_io.bluekenue import write_i2s

    path = str(tmp_path / "bk.gpkg")
    M.create_mesh_project(path, QgsCoordinateReferenceSystem(CRS))
    i2s = str(tmp_path / "lines.i2s")
    write_i2s(i2s, [np.array([[512000, 2155000], [512100, 2155000]], float),
                    np.array([[512000, 2155050], [512050, 2155050], [512050, 2155100], [512000, 2155050]], float)], [0, 4.0])
    r = processing.run("telemactoolkit:bluekenuelines", {"INPUT": i2s, "CRS": CRS, "PROJECT": path, "TARGET": 1, "OUTPUT": "memory:"})
    assert r["OUTPUT"].featureCount() == 2
    assert sorted(f["closed"] for f in r["OUTPUT"].getFeatures()) == [0, 1]
    lay = M.open_mesh_project(path)
    assert lay["constraint_lines"].featureCount() == 2
    assert {f["mode"] for f in lay["constraint_lines"].getFeatures()} == {"breakline"}
    processing.run("telemactoolkit:bluekenuelines", {"INPUT": i2s, "CRS": CRS, "PROJECT": path, "TARGET": 5, "OUTPUT": "memory:"})
    lay = M.open_mesh_project(path)
    assert [f["size"] for f in lay["refinement_areas"].getFeatures()] == [4.0]  # only the closed line is an area


@pytest.mark.skipif(not (_gmsh() and os.environ.get("TELEMAC_CMD")), reason="needs gmsh and TELEMAC_CMD")
def test_gmsh_mesh_runs_in_telemac(case, dem):  # noqa: F811
    import processing

    from test_edit_qgis import _run_chain

    d, path, _mesh, _cli = case
    g = processing.run("telemactoolkit:generatemesh", {"PROJECT": path, "METHOD": 3, "RASTERS": [dem],
                                                       "OUTPUT": str(d / "gmshfrontal.slf"), "LOAD": False})
    flux, _log = _run_chain(d, path, g["OUTPUT"], g["OUTPUT_CLI"], "gmshfrontal")
    assert abs(flux[0] - 60) < 1.5 and abs(flux[1] - 60) < 1.5, flux
