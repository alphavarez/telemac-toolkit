# 8. Formatos, coordenadas y condiciones de frontera

## Formatos

El generador de mallas está hecho para TELEMAC; la misma malla se puede escribir para otros programas desde
el paso 4 del flujo de trabajo de la malla o con *Mesh → Export or convert mesh* (que además lee SELAFIN,
2DM, UGRID, MIKE y cualquier malla que QGIS abra).

| Formato | Para | Fronteras | Rugosidad | Cuadriláteros |
| --- | --- | --- | --- | --- |
| SELAFIN `.slf` + `.cli` + `.liq` | TELEMAC-2D/3D y la suite TELEMAC | `.cli` en el orden de TELEMAC, plantilla `.liq` | variable de nodo `BOTTOM FRICTION` | no |
| SMS 2DM `.2dm` | BASEMENT, TUFLOW FV, SMS, QGIS | *nodestrings* (+ `_nodestrings.csv` con sus nombres) | identificadores `MAT` | sí |
| Gmsh `.msh` (ASCII 2.2) | Gmsh y los solvers que lo leen | líneas físicas por frontera | superficies físicas por material | sí |
| Malla GiD `_gid.msh` | GiD y programas basados en GiD, como Iber | — | material por elemento | sí |
| UGRID `.nc` | Delft3D FM / D-HYDRO | — | — | sí |
| MIKE 21 FM `.mesh` | MIKE 21 FM | código de frontera por nodo | — | sí |
| ADCIRC `fort.14` (`.14`) | ADCIRC, SCHISM (`hgrid.gr3`) | listas de nodos de frontera abierta y de tierra | — | no |
| HEC-RAS `.gpkg` | HEC-RAS (RAS Mapper) | líneas de condiciones de frontera | — | — |

- **HEC-RAS** construye sus propias celdas a partir de un perímetro, líneas de quiebre y tamaños de celda,
  así que recibe el dibujo (perímetro, líneas de quiebre, regiones de refinamiento y líneas de frontera en
  un GeoPackage para RAS Mapper), no los triángulos. Usa *Generate mesh* para él, no *Export*.
- **Cuadriláteros**: los canales y las cuadrículas estructurados pueden quedar como cuadriláteros en los
  formatos que los aceptan (*Keep structured cells as quadrilaterals*).
- **Elevación**: se escribe como elevación del nodo (ADCIRC: profundidad = −elevación, que es su
  convención).
- **Comprobado**: los archivos SELAFIN, 2DM, UGRID y MIKE se vuelven a leer con QGIS en las pruebas
  automáticas, los de Gmsh con meshio, y los de TELEMAC se ejecutan en TELEMAC-2D. Los archivos GiD y ADCIRC
  siguen sus formatos publicados, pero todavía no se han cargado en Iber ni en ADCIRC: revisa el primero que
  hagas.

## Coordenadas en el archivo

| Opción | El archivo contiene |
| --- | --- |
| **Map coordinates as drawn** | las coordenadas del SRC del proyecto |
| **Lower-left corner at (0, 0)** | coordenadas desplazadas para que la esquina inferior izquierda de la malla sea el origen |
| **A chosen point at (0, 0)** | coordenadas desplazadas para que ese punto sea el origen |
| **Rotation** | ejes girados en sentido antihorario según el ángulo (por ejemplo, x a lo largo de un canal dibujado con un ángulo) |

Para TELEMAC, un desplazamiento sin rotación se puede **conservar en el encabezado del archivo** (IPARAM 3 y
4): TELEMAC calcula cerca de (0, 0), con toda la precisión de la precisión simple, y QGIS sigue dibujando la
malla en su lugar. En otro caso, una nota `<nombre>_frame.txt` junto al archivo da el desplazamiento y la
rotación para volver a las coordenadas del mapa, y la malla se abre en coordenadas locales.

Las coordenadas grandes (UTM, cuadrículas nacionales) en archivos SELAFIN de precisión simple se guardan
automáticamente con un desplazamiento de origen entero, de modo que float32 conserva una resolución
submétrica de menos de un centímetro.

### Coordenadas locales

La cuadrícula local del plugin ("Local metric grid (laboratory)") es un plano con (0, 0) en el origen y
longitudes en metros. Se comporta como un SRC proyectado normal en QGIS (medir, ensamblar, etiquetas,
todos los formatos), pero no es un lugar de la Tierra: no se puede mostrar junto con capas de mapa en un
SRC real. La barra de estado de QGIS la llama "SRC desconocido" (*Unknown CRS*). Junto a cada malla se
escribe un archivo `.prj` con su definición.

## Condiciones de frontera

*Generate mesh* (a partir de la capa de fronteras) y *Boundary conditions from lines* (para una malla
existente) escriben el `.cli` de TELEMAC:

- Todo nodo de frontera es una **pared** (2 2 2), salvo los nodos cercanos a una línea de frontera, que
  toman su tipo: `Q` (4 5 5), `H` (5 4 4), `QH` (5 5 5), `UV` (4 6 6), `UVH` (5 6 6), `free` (4 4 4),
  `waves` (1 1 1), o un código de tres dígitos.
- "Cercano" = a menos de la mitad de la arista de frontera más corta de ese nodo (o la distancia de
  ensamblado que indiques), así que una línea termina limpiamente donde las aristas gruesas de la planicie
  se encuentran con las aristas finas del río.
- Los nodos se enumeran en el orden de TELEMAC: la frontera exterior en sentido antihorario y después las
  islas.

### Números de las fronteras líquidas

TELEMAC numera las fronteras líquidas a su manera: empieza en el **punto de la frontera con el menor x + y**
(la esquina más baja a la izquierda) y camina en sentido antihorario, numerando las fronteras líquidas a
medida que las encuentra; dos tramos líquidos contiguos con tipos distintos son dos fronteras. Los números
importan: el archivo de configuración (`PRESCRIBED FLOWRATES = 0.05;0.`) y las columnas del `.liq` (`Q(1)`,
`SL(2)`) se refieren a ellos.

El plugin numera las fronteras exactamente como lo hace TELEMAC y las informa en el registro ("Liquid
boundary 1: outlet, codes 5 4 4…"), junto con las líneas del archivo de configuración que debes pegar. La
plantilla `.liq` que escribe tiene las columnas correctas. En un canal que corre a lo largo de +x desde
(0, 0), la **salida** es la frontera 1, porque el recorrido empieza en (0, 0) y avanza por la margen
derecha hasta la salida primero.

[← Editar una malla existente](07-editar-malla.md) · [Contenido](README.md) · [Siguiente: Referencia de las herramientas de Procesos →](09-referencia-procesos.md)
