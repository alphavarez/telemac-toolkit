# 7. Editar una malla existente

*Complementos → TELEMAC Toolkit → Mesh → Edit node values of a mesh* cambia una variable de nodo de una
geometría SELAFIN existente (el fondo, la fricción o cualquier otra) sin volver a mallar. El archivo nuevo
tiene los **mismos nodos, elementos, numeración de fronteras, origen y fecha**, así que el archivo de
condiciones de frontera (`.cli`) sigue sirviendo; se copia junto al archivo nuevo. Las demás variables se
conservan.

## De dónde salen los valores nuevos

En orden de prioridad (la primera fuente que cubre un nodo da su valor):

1. **Rásteres** (bilineal), en el orden en que los enumeres;
2. **Puntos de levantamiento** (lineal sobre su triangulación), con un campo de valor o la Z del punto;
3. una **fórmula** de `x` e `y`, o de las coordenadas de canal `s` y `n` a lo largo de una **línea de eje**
   que indiques;
4. una **constante**, cuando no hay otra fuente (o, con *Fill nodes no source covers* desactivada en las
   opciones avanzadas, para los nodos que las otras fuentes no cubren; con ella activada, toman el valor más
   cercano).

*Unit factor* convierte rásteres y puntos en pies (0.3048) o en centímetros (0.01).

## Qué pasa con los valores anteriores

| Operación | Valor nuevo | Para |
| --- | --- | --- |
| **Replace** (reemplazar) | el valor nuevo | un levantamiento nuevo, un MDE corregido, un mapa de fricción |
| **Add** (sumar) | anterior + valor nuevo | subir o bajar un espesor (una constante negativa baja) |
| **Raise to** (subir hasta) | el mayor entre el anterior y el nuevo | rellenos, diques, terraplenes, coronas de camino |
| **Lower to** (bajar hasta) | el menor entre el anterior y el nuevo | dragado, canales excavados, excavaciones |

## Dónde

- **La malla completa**, o **solo dentro de polígonos**.
- **Transition width** (ancho de transición): dentro del polígono, el cambio pasa de nulo en el borde a
  completo a esa distancia de él, de modo que una zona dragada no termina en un escalón vertical.

## Recetas

**Levantamiento nuevo sobre un tramo.** Rasters = el MDE nuevo; operación *Replace*; área = el polígono del
levantamiento; ancho de transición de 2 a 3 elementos.

**Canal dragado con taludes.** Dibuja el eje del canal como una línea. Fórmula (5 m de ancho de plantilla,
taludes 2:1, fondo 3 m por debajo de la referencia de la fórmula, pendiente 0.0005):
`20.0 - 0.0005*s - 3 + clip((abs(n) - 2.5)/2, 0, 3)`; operación *Lower to*; área = un polígono un poco más
ancho que el corte.

**Corona de un dique o bordo.** Constante = cota de la corona; operación *Raise to*; área = un polígono
delgado a lo largo de la corona; ancho de transición = la mitad del ancho de la base del dique. Para una
corona definida, la malla también necesita una línea de quiebre (*breakline*) sobre ella (flujo de trabajo
de la malla, líneas de restricción).

**Zona de fricción.** Variable `BOTTOM FRICTION`; constante = el valor; operación *Replace*; área = la zona.
Si el archivo no tiene esa variable, se agrega (con 0 en el resto: primero asigna un valor a toda la malla
y luego a las zonas).

## Cambiar la propia malla

Las herramientas anteriores cambian valores en los nodos. Estas cambian los nodos y los elementos, y
conservan todo lo demás: cada variable de cada instante (los nodos que no cambian conservan sus valores
exactamente, los nodos nuevos reciben valores interpolados en los elementos originales, o del nodo original
más cercano fuera de ellos), el `.cli` reconstruido para la nueva numeración, y el `.prj`.

### Editar a mano en QGIS

1. Flujo de trabajo de la malla, paso 5 → **Edit in QGIS (Mesh Digitizing)**. Una copia de la malla
   (`nombre_editing.slf`) se abre en modo de edición.
2. Usa la barra **Mesh Digitizing** de QGIS (en QGIS en español, *Digitalización de malla*; *Ver → Barras de
   herramientas*): agregar, quitar o mover vértices, voltear aristas, dividir caras, forzar la malla a lo
   largo de líneas.
3. **Rebuild after editing** guarda las ediciones y escribe `nombre_edited.slf` y `nombre_edited.cli`, y luego
   agrega el resultado al mapa.

QGIS guarda solo la geometría de una malla editada (descarta el fondo y las demás variables) y no sabe nada
del `.cli`; la reconstrucción devuelve ambos. En la caja de herramientas el mismo paso es *Mesh → Rebuild a
TELEMAC mesh after editing it in QGIS*, que también acepta una malla que QGIS guardó en otro formato.

**Condiciones de frontera después de editar**: los nodos de frontera que ya existían conservan sus
condiciones. Un nodo de frontera nuevo toma el tipo de la arista de frontera original sobre la que queda
cuando ambos extremos de esa arista lo comparten (un nodo agregado en medio de una entrada es un nodo de
entrada); en otro caso es una pared. El registro informa cuántos nodos conservaron sus condiciones y
cuántos son nuevos.

### Re-mallar dentro de polígonos

*Mesh → Re-mesh inside polygons* (paso 5 → **Re-mesh inside polygons…**) reemplaza los elementos cuyo
centro cae dentro de los polígonos por una malla nueva de la misma área: para refinar alrededor de la pila de
un puente o de un vertedor, reparar un parche de elementos malos o hacer más gruesa una zona.

- **Tamaño**: un número, un campo de los polígonos, o 0 para la mediana de las aristas del contorno del área.
  El tamaño crece desde las aristas del contorno hasta el tamaño interior con la tasa de crecimiento.
- **Método**: *Balanced* o *Quality-guaranteed* (ángulo mínimo de 30°), como en el paso 1.
- Donde el área colinda con elementos viejos, su contorno conserva sus nodos, así que los elementos nuevos se
  unen a los viejos nodo a nodo. Donde llega a la frontera de la malla pueden aparecer nodos de frontera
  nuevos; toman las condiciones de frontera como se explicó arriba.

### Líneas de frontera a partir de un .cli

*Mesh → Boundary lines from a .cli* (paso 5 → **Boundary lines from the .cli**) escribe una línea por cada
frontera líquida de una malla existente, con los campos `name`, `type` y `value` de la capa de fronteras de
un proyecto de malla (los valores salen de la primera fila de un archivo `.liq` si lo indicas). Edita las
líneas en el mapa (mueve sus extremos a lo largo de la frontera, cambia un tipo, agrega una frontera) y
escríbelas de vuelta con *Boundary conditions from lines*. Ensambla a la frontera de la malla mientras
editas (*Mesh boundary as lines* la da como capa).

**Comprobado en TELEMAC-2D**: una malla de río re-mallada en su extremo aguas abajo (nodos nuevos en la
frontera de salida) y una malla de río editada en QGIS corren ambas en TELEMAC-2D v9.1 con su `.cli`
reconstruido, y el caudal prescrito pasa (60 m³/s de entrada y de salida).

## Otros cambios a una malla existente

| Cambio | Herramienta |
| --- | --- |
| Condiciones de frontera (`.cli`, `.liq`) | *Mesh → Boundary conditions from lines*; las actuales como líneas: *Boundary lines from a .cli* |
| Nodos y elementos a mano | Paso 5 → *Edit in QGIS*, luego *Rebuild after editing* |
| Elementos más finos, más gruesos o mejores en una zona | *Mesh → Re-mesh inside polygons* |
| Coordenadas: desplazar a (0, 0), rotar | *Mesh → Export or convert mesh*, opciones de coordenadas ([detalles](08-formatos-coordenadas-fronteras.md#coordenadas-en-el-archivo)) |
| Recortar una parte | *Mesh → Extract part of a mesh* |
| Otro formato | *Mesh → Export or convert mesh* |
| Variables nuevas a partir de fórmulas (por ejemplo, fricción según clases de tirante) | [Calculadora](04-calculadora-integrales.md) |
| Otra malla (una línea de quiebre nueva, un refinamiento nuevo) | Flujo de trabajo de la malla: edita las capas del proyecto y genera de nuevo |

## Para quienes vienen de Blue Kenue

En Blue Kenue las mismas tareas pasan por proyectar un conjunto de datos sobre la malla (*Map Object*), por
la calculadora y por construir un objeto SELAFIN nuevo a partir de la malla y sus variables antes de
guardarlo. Aquí los pasos son una sola herramienta que escribe un archivo de geometría listo, y cada paso se
puede repetir o ejecutar por lotes desde la caja de herramientas de Procesos o desde un modelo. La edición de
nodos de Blue Kenue (mover, agregar y borrar nodos, intercambiar aristas) es la barra Mesh Digitizing de QGIS
seguida de *Rebuild after editing*, que además conserva las variables y el `.cli`.

**Traer el trabajo de Blue Kenue**:

- *Mesh → Import Blue Kenue lines (.i2s, .i3s)* lee un conjunto de líneas en una capa de líneas, con el valor de
  cada línea, y puede agregarlas directamente a un proyecto de malla: como líneas de quiebre o líneas guía, como
  el dominio o islas (las líneas cerradas se vuelven polígonos), o como áreas de refinamiento con el valor de la
  línea como su tamaño. Los archivos de Blue Kenue no traen SRC, así que indícalo.
- Las mallas T3 de Blue Kenue (`.t3s`) se leen en cualquier lugar donde el plugin lee una malla: *Compare
  meshing methods*, *Export or convert mesh*, *Rebuild a TELEMAC mesh after editing it in QGIS*. Con las mismas
  líneas de entrada, una malla de Blue Kenue y los métodos del plugin se pueden comparar medida por medida.

[← Plantillas de canal](06-plantillas-de-canal.md) · [Contenido](README.md) · [Siguiente: Formatos, coordenadas y condiciones de frontera →](08-formatos-coordenadas-fronteras.md)
