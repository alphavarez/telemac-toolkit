# 12. Condiciones de frontera y archivos de corrida

Después de *Generate mesh*, el paso **6 · Run files** del panel *Mesh workflow* prepara el resto de una
corrida de TELEMAC-2D. Cada botón abre una herramienta de Procesos con la malla generada ya indicada. Las
herramientas funcionan con cualquier malla SELAFIN con su `.cli` y también están en la caja de herramientas,
en *Boundary conditions and run*.

| Botón | Escribe |
| --- | --- |
| *Rating curve (outflow)…* | el archivo de curvas de descarga |
| *Hydrographs (.liq) from a table…* | el archivo de fronteras líquidas |
| *Initial state…* | un archivo de condiciones iniciales (SELAFIN) |
| *Steering file (.cas)…* | el archivo de configuración de TELEMAC-2D |
| *Check boundary conditions…* | una capa de las fronteras líquidas con sus números de TELEMAC, y un informe |

Las fronteras líquidas siempre se numeran como las numera TELEMAC, así que las columnas del archivo `.liq`,
las curvas y las listas del archivo de configuración coinciden con lo que lee TELEMAC.

## Curva de descarga para la salida

Rara vez se conoce el nivel del agua en el extremo aguas abajo. Una **curva de descarga** (curva
elevación-caudal) permite que TELEMAC fije el nivel a partir del caudal que llega a la frontera. La
herramienta calcula la curva de tirante normal sobre la propia frontera, con sus nodos y su fondo,
exactamente como TELEMAC ve la sección:

- Fórmula de Manning Q = √S · Σ Aᵢ Rᵢ^(2/3) / nᵢ, con la sección dividida donde cambia la rugosidad (el
  método habitual de sección dividida).
- *Pendiente de energía*: en general la pendiente del fondo o de la superficie del agua del tramo aguas
  abajo.
- *n de Manning*: un valor, o la variable `BOTTOM FRICTION` de la malla.
- Por encima del extremo más bajo de la sección, sus extremos actúan como paredes verticales, igual que las
  fronteras cerradas de la malla.
- Se advierte si el flujo en la sección se vuelve supercrítico, donde una curva de descarga no controla el
  nivel.

El archivo sigue el formato de TELEMAC, y un archivo existente conserva sus otras curvas. La frontera debe ser
de elevación prescrita (`H`, 5 4 4). El registro da las líneas del archivo de configuración:
`STAGE-DISCHARGE CURVES = 1;0` (un valor por frontera líquida) y `STAGE-DISCHARGE CURVES FILE`.

## Hidrogramas desde una tabla

Carga en QGIS una hoja CSV o Excel (por ejemplo, la salida de HEC-HMS) y elige:

- el **campo de tiempo** y su unidad: segundos, minutos, horas, días, o fechas y horas;
- los **campos de valores**, en el orden de las columnas que espera TELEMAC (el registro las enumera: `Q(n)`
  para una frontera de caudal, `SL(n)` para un nivel, `U(n)` y `V(n)` para velocidades), o las columnas por
  nombre, por ejemplo `Q(2)=caudal_m3s; SL(1)=nivel`.

Una frontera con curva de descarga no necesita columna `SL`; indica el archivo de curvas y se omite. Los
tiempos se escriben en segundos desde la primera fila. Con fechas, la primera se indica para
`ORIGINAL DATE OF TIME` y `ORIGINAL HOUR OF TIME`.

## Estado inicial

TELEMAC arranca con más suavidad desde un estado cercano al flujo esperado. La herramienta escribe
velocidades nulas, el tirante, la superficie libre y el fondo, a partir de:

- un **nivel de agua constante**;
- un **nivel a lo largo del eje de un canal**, desde un nivel aguas arriba hasta uno aguas abajo (cada nodo
  toma el nivel de su posición a lo largo del eje), lo que conviene para un río;
- un **tirante constante sobre el fondo**;
- un **tirante sobre el fondo de los ejes de los canales**, para ríos y redes de ríos: cada nodo toma el fondo
  del punto más cercano de cualquier eje más el tirante, así que las márgenes y la planicie por encima de ese
  nivel quedan secas y todos los brazos empiezan con agua. La capa de ejes puede ser la capa de canales del
  proyecto de malla (solo se leen sus ejes).

Los nodos donde el nivel queda por debajo del fondo quedan secos. Para un río, usa aguas abajo el nivel de
la curva de descarga para el primer caudal, y aguas arriba ese nivel más pendiente × longitud. Un estado
lejano de los primeros valores de frontera hace que TELEMAC vacíe o llene el modelo al arrancar, lo que puede
dejar seca una frontera de caudal.

## Archivo de configuración

El archivo de configuración nombra los archivos anteriores con rutas relativas a sí mismo, y contiene:

- las fronteras líquidas en el orden de TELEMAC, con `PRESCRIBED FLOWRATES` y `PRESCRIBED ELEVATIONS` tomados
  de la primera fila del `.liq`;
- `STAGE-DISCHARGE CURVES` para las fronteras que tienen curva;
- `VELOCITY PROFILES = 4` (velocidad proporcional a √h) en las fronteras de caudal;
- el paso de tiempo: indicado, o estimado para un número de Courant de 1 en el elemento más pequeño, a partir
  de un tirante y una velocidad típicos;
- la ley y el coeficiente de fricción; cuando la malla tiene la variable `BOTTOM FRICTION`, TELEMAC la lee de
  ahí;
- ajustes numéricos para mojado y secado con tirantes positivos, los de los casos de demostración de este
  plugin que corrieron en TELEMAC-2D v9.1, incluido `MAXIMUM NUMBER OF ITERATIONS FOR ADVECTION SCHEMES = 200` y,
  desde la 0.5.0, el esquema NERD para las velocidades (`SCHEME FOR ADVECTION OF VELOCITIES = 14`).

Revisa cada valor antes de correr; los comentarios del inicio enumeran las fronteras líquidas.

### Probado en TELEMAC-2D

La cadena completa se corrió en TELEMAC-2D v9.1.1 sobre un río con meandros de 2 km: asistente de río con el
MDE, *Generate mesh*, curva de descarga en la salida, una entrada de 60 m³/s desde un archivo CSV, nivel
inicial a lo largo del eje y archivo de configuración. Después de 3 horas la entrada era de 59.9 m³/s, la
salida de 60.0 m³/s, el error de volumen de 10⁻¹⁵, y el nivel en la salida coincidía con la curva de descarga
(4.212 m) con menos de un milímetro de diferencia. Esta corrida forma parte de las pruebas.

Cuatro hallazgos de estas corridas ya están incorporados:

- **Perfil de velocidad en las fronteras de caudal.** Con una velocidad constante a lo ancho de la entrada (la
  opción por defecto de TELEMAC), los nodos junto a las márgenes, donde el agua es menos profunda, se vaciaban
  hasta que la corrida se detenía con `DEBIMP: PROBLEM ON BOUNDARY`. El perfil en √h los mantiene mojados.
- **Diagonales de las celdas estructuradas.** Las franjas de canal y las zonas de cuadrícula cortaban sus
  celdas con diagonales alternadas, de modo que los nodos interiores tenían cuatro y ocho triángulos
  alternadamente. Los tirantes de TELEMAC oscilaban entonces de nodo a nodo. Ahora las celdas se cortan por su
  diagonal más corta (en una sola dirección en las cuadrículas), y cada nodo interior tiene seis triángulos.
- **Advección de las velocidades (0.5.0).** Con el método de las características, el mismo río en una malla
  más fina (paso de tiempo de 0.24 s en lugar de 0.55 s) se detuvo a los 616 s: la entrada se volvió
  supercrítica en un nodo poco profundo, se secó, y TELEMAC se detuvo con `DEBIMP: PROBLEM ON BOUNDARY`. La
  malla original falló de la misma manera, así que era la numérica y no la malla. Se descartaron otros dos
  remedios: el método de Thompson en la entrada corrió pero dejó entrar 67 m³/s en lugar de los 60 prescritos,
  y el perfil de velocidad con tirante virtual (5) falló igual. Con el esquema NERD (14) ambos pasos de tiempo
  corrieron las 3 horas, con 60.000 m³/s de entrada y de salida. Es unas dos a tres veces más lento que las
  características; vuelve a poner `SCHEME FOR ADVECTION OF VELOCITIES = 1` para ganar velocidad cuando ya sepas
  que un modelo es estable.
- **Un paso de tiempo que divide la duración (0.5.0).** TELEMAC corre round(duración / paso de tiempo)
  pasos, así que un paso de 0.61 s durante 3 horas terminaba su último paso en 10 800.05 s, después del final de
  un hidrograma que acababa en 10 800 s, y la corrida se detenía con `OUT OF RANGE OF THE FILE OF LIQUID
  BOUNDARIES` (TELEMAC tolera 1 ms). Ahora el archivo de configuración toma el mayor paso por debajo del límite
  de Courant que divide exactamente la duración, con cuatro decimales como máximo (0.61 s pasa a 0.60 s, 18 000
  pasos). Un paso de tiempo que indiques se ajusta de la misma forma, con un aviso, y la herramienta avisa si
  el archivo de fronteras líquidas termina antes de la duración.

## Revisar las condiciones de frontera

La revisión enumera cada frontera líquida con su número, tipo, nodos, longitud y rango del fondo, la dibuja en
el mapa etiquetada con su número, y advierte sobre:

- un valor prescrito sin columna en el archivo `.liq`;
- un nivel que nunca sube por encima del fondo de su frontera;
- columnas para fronteras que la malla no tiene;
- tiempos que no aumentan;
- una curva de descarga en una frontera que no es de elevación prescrita;
- listas del archivo de configuración con un número equivocado de valores;
- dos fronteras líquidas de tipos distintos que se tocan sin un nodo de pared entre ellas;
- fronteras con menos de cinco nodos.

Si se indica el archivo de configuración, encuentra los demás archivos a través de él.

[← Solución de problemas](11-solucion-de-problemas.md) · [Contenido](README.md)
