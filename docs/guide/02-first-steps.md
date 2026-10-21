# 2. First steps

## A result in five minutes

1. *Open TELEMAC result…* on the toolbar, or drag a `.slf` file onto the QGIS window.
2. If the plugin cannot tell what the coordinates are, it asks: **local coordinates** (a laboratory flume
   or a schematic channel) or a **map CRS**. Leave *Save the choice next to the file* ticked so the
   question is not asked again ([details](03-results.md#coordinates-of-a-result)).
3. The **TELEMAC Results** panel opens and the map shows the water depth with velocity arrows.
4. Press **play**, or drag the slider. The QGIS Temporal Controller moves with it.
5. Change *Map shows* to *Velocity magnitude*; open *Arrow settings* to make the arrows denser, sparser,
   longer or shorter.
6. Turn on the **probe** (Time series tab) and click two points on the mesh: their time series are drawn.
   Turn on the **profile** tool, click along a line and right-click to finish: the bed and the water
   surface along it follow the clock.

## A laboratory flume in ten minutes

1. Toolbar → **Mesh workflow panel** → *New…*, save `flume.gpkg`, answer **Yes** to local coordinates.
2. Step 1 → *Channel or flume template…*:
   - Planform *Straight*, width 1 m, length 20 m, element size 0.05 m, *Structured*.
   - Bed level at the inlet 0.30 m, slope 0.001, rectangular section.
   - Inflow 0.05 m³/s, outlet level 0.42 m. *Run*.
3. The domain, the inlet (Q) and outlet (H) lines and a structured grid appear on the map; the bed formula
   is saved in the project.
4. Step 4 → format *TELEMAC*, mesh file `geo_flume.slf`, **Generate mesh**. The log lists the mesh, the
   boundary numbers TELEMAC will use and the lines to paste in the steering file
   (`PRESCRIBED FLOWRATES`, `PRESCRIBED ELEVATIONS`).
5. The plugin writes `geo_flume.slf`, `geo_flume.cli` and (if you ask for it) a `.liq` template. Write
   a steering file around them and run `telemac2d.py t2d_flume.cas`.
6. Open the result `r2d_flume.slf` as above.

The whole case, with its steering file, is in `scripts/make_flume_demo.py` and
`tests/data/demo_flume/`. A curved version (180° bend) is in `scripts/make_bend_demo.py` and
`tests/data/demo_bend/`.

[← Install](01-install.md) · [Contents](README.md) · [Next: Viewing results →](03-results.md)
