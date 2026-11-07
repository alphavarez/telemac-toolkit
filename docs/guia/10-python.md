# 10. Scripts con Python

## Herramientas de Procesos desde la consola de Python de QGIS

Todas las herramientas se ejecutan con `processing.run`; las claves están en la
[referencia](09-referencia-procesos.md).

```python
import processing

# un canal con curva de 180° en un proyecto nuevo, mallado para TELEMAC
processing.run("telemactoolkit:channeltemplate", {
    "NEW_PROJECT": "C:/work/bend.gpkg", "PLANFORM": 1, "WIDTH": 1.0, "INLET": 3, "OUTLET": 3,
    "RADIUS": 2.0, "BEND_ANGLE": 180, "TURN": 0, "CELL": 0.05, "STRUCTURED": True,
    "BED_Z0": 0.30, "SLOPE": 0.001, "Q": 0.05, "H": 0.43})
processing.run("telemactoolkit:generatemesh", {
    "PROJECT": "C:/work/bend.gpkg", "FRICTION": 0.025, "OUTPUT": "C:/work/geo_bend.slf",
    "OUTPUT_LIQ": "C:/work/bend.liq"})

# máximos de un resultado
processing.run("telemactoolkit:calculator", {
    "INPUT": "C:/work/r2d_bend.slf", "EXPRESSIONS": "HMAX [M] = tmax(H)\nVMAX [M/S] = tmax(mag(U, V))",
    "KEEP": False, "OUTPUT": "C:/work/r2d_bend_max.slf"})

# caudal a través de las líneas de una capa, para cada paso de tiempo
res = processing.run("telemactoolkit:integrals", {
    "INPUT": "C:/work/r2d_bend.slf", "MODE": 0, "LINES": "C:/work/sections.gpkg", "OUTPUT": "memory:"})
for f in res["OUTPUT"].getFeatures():
    print(f["name"], f["time_s"], f["discharge"])
```

Los casos de demostración completos son scripts: `scripts/make_flume_demo.py`, `scripts/make_bend_demo.py`
(+ `check_bend_demo.py`) y `scripts/make_mesh_demo.py` (un valle a partir de capas y un MDE).

## Desde la línea de comandos, sin abrir QGIS

`qgis_process` ejecuta las mismas herramientas desde una terminal o un archivo por lotes. Viene con QGIS (en
Windows, `qgis_process-qgis.bat` o `qgis_process-qgis-ltr.bat` en la carpeta `bin` de la instalación de
QGIS) y encuentra el plugin en tu perfil de QGIS una vez que está habilitado:

```bat
qgis_process plugins enable telemac_toolkit
qgis_process run telemactoolkit:calculator -- INPUT=C:/work/r2d_bend.slf "EXPRESSIONS=HMAX [M] = tmax(H)" OUTPUT=C:/work/r2d_bend_max.slf
qgis_process help telemactoolkit:generatemesh
```

Con muchos parámetros, pásalos como JSON por la entrada estándar:

```bat
echo {"inputs": {"PROJECT": "C:/work/bend.gpkg", "FRICTION": 0.025, "OUTPUT": "C:/work/geo_bend.slf"}} | qgis_process run telemactoolkit:generatemesh -
```

Esta es la forma de encadenar el plugin con una corrida de TELEMAC en un script (construir la malla,
ejecutar `telemac2d.py`, calcular máximos y caudales), por ejemplo con TELEMAC en un contenedor Docker y los
archivos en una carpeta compartida.

## Las bibliotecas sin QGIS

`telemac_toolkit/telemac_io` (SELAFIN, `.cli`, calculadora, integrales) y `telemac_toolkit/mesher`
(generador, plantas, escritores) solo necesitan Python y NumPy (el generador de mallas también SciPy, para
la triangulación), así que también se ejecutan dentro de un contenedor de TELEMAC o en un clúster, junto al
solver:

```python
import sys; sys.path.insert(0, "path/to/telemac_toolkit")
from telemac_io import SelafinFile, write_selafin
from telemac_io.calculator import run_calculator, MeshOperators
from telemac_io.integrals import Section

with SelafinFile("r2d.slf") as slf:
    print(slf)                                     # dimensiones, variables, cuadros
    depth = slf.read("WATER DEPTH", frame=-1)      # último cuadro
    ts = slf.series("FREE SURFACE", nodes=[10])    # serie de tiempo, leída bajo demanda
    h = slf.header
    ops = MeshOperators(h.x, h.y, h.ikle)
    w = ops.curl(slf.read("VELOCITY U", -1), slf.read("VELOCITY V", -1))
    sec = Section.build([(1.5, 0.5), (1.5, -0.5)], h.x, h.y, h.ikle)
    q = sec.discharge(depth, slf.read("VELOCITY U", -1), slf.read("VELOCITY V", -1))

run_calculator("r2d.slf", "r2d_calc.slf", "FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))")
write_selafin("geo.slf", x, y, ikle, {"BOTTOM": z}, units={"BOTTOM": "M"})
```

El lector acepta precisión simple y doble, big endian y little endian, archivos 2D y 3D, archivos con
registro de fecha y archivos incompletos de simulaciones en ejecución.

[← Referencia de Procesos](09-referencia-procesos.md) · [Contenido](README.md) · [Siguiente: Solución de problemas →](11-solucion-de-problemas.md)
