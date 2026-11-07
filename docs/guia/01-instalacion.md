# 1. Instalación

## Requisitos

- QGIS 3.44 LTR (Qt5) o QGIS 4.x (Qt6), en Windows, Linux o macOS.
- Nada más. El plugin usa solamente Python y NumPy, que vienen con QGIS. No necesitas TELEMAC para construir
  mallas ni para ver resultados; solo lo necesitas para ejecutar las simulaciones.

## Desde el repositorio de plugins de QGIS

*Complementos → Administrar e instalar complementos → Todos*, busca **TELEMAC Toolkit** y pulsa *Instalar
complemento*. Mientras el plugin esté marcado como experimental, primero activa *Mostrar también los
complementos experimentales* en la pestaña *Configuración* del administrador. QGIS te avisa cuando hay una
versión nueva (pestaña *Actualizables*).

## Desde un archivo zip

1. Descarga `telemac_toolkit-<versión>.zip` (el zip instalable, no el zip del código fuente).
2. *Complementos → Administrar e instalar complementos → Instalar a partir de ZIP*, elige el archivo y pulsa
   *Instalar complemento*.
3. Si había una versión anterior instalada, se reemplaza. Reinicia QGIS para que los paneles y menús se
   carguen bien.

## Opcional: Gmsh

Los métodos de mallado con Gmsh necesitan el paquete de Python **gmsh** (GPL-2.0 o posterior, compatible con la
licencia de este plugin) en el propio Python de QGIS. Todo lo demás funciona sin él, y el panel muestra los
métodos de Gmsh desactivados hasta que lo instales.

- **Windows**: abre el *OSGeo4W Shell* que viene con QGIS y ejecuta `python -m pip install gmsh`.
- **macOS**: `/Applications/QGIS.app/Contents/MacOS/bin/python3 -m pip install gmsh` (la ruta de tu aplicación
  de QGIS).
- **Linux**: `python3 -m pip install gmsh` con el Python que usa QGIS; además la rueda necesita la biblioteca
  del sistema libGLU (`libglu1-mesa` en Debian y Ubuntu).

Después reinicia QGIS. *Compare meshing methods* indica si encontró Gmsh.

## Opcional: MMG

*Optimise or adapt a mesh (MMG)* ejecuta el programa independiente de MMG, `mmg2d_O3` (LGPL); no se instala nada
en el Python de QGIS. Descárgalo de [mmgtools.org](https://www.mmgtools.org) (versiones para Windows, macOS y
Linux) o de conda-forge (`mmgsuite`), y la primera vez que uses la herramienta indica la ruta de `mmg2d_O3`; se
recuerda. También sirve ponerlo en el PATH o en una variable de entorno `MMG2D`. No se usa el paquete de Python
de MMG (mmgpy): necesita un NumPy más reciente que el de algunas instalaciones de QGIS.

## Qué aparece

- Una barra de herramientas **TELEMAC Toolkit**: abrir resultado, panel de resultados, cursor en vivo,
  sonda, perfil, panel *Mesh workflow* y nuevo proyecto de malla.
- Un menú **Complementos → TELEMAC Toolkit** con las mismas herramientas, la calculadora, integrales y
  flujos, información SELAFIN, extracción 3D y un submenú **Mesh** con todas las herramientas de malla.
- Un proveedor **TELEMAC Toolkit** en la caja de herramientas de Procesos con todas las herramientas
  ([referencia](09-referencia-procesos.md)).

## Actualizar y desinstalar

Las actualizaciones llegan por el administrador de complementos (o instalando un zip más nuevo). Los
proyectos de malla hechos con la 0.2.0 y posteriores siguen funcionando: los nombres de capa antiguos
(`outline`, `density`, `hard_lines`, `soft_lines`, `hard_points`, `grids`) se leen como sus equivalentes
actuales. Para desinstalar: administrador de complementos → *Instalados* → *Desinstalar complemento*. Tus
archivos no se tocan.

[← Contenido](README.md) · [Siguiente: Primeros pasos →](02-primeros-pasos.md)
