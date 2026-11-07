# 4. Calculadora, integrales y flujos

## Calculadora

Ábrela desde el **botón de la calculadora** del panel de resultados (el cuadro de expresiones entonces
enumera las variables del resultado actual) o desde *Complementos → TELEMAC Toolkit → Calculator*. Escribe
un **archivo SELAFIN nuevo** con las variables que definas, más las originales si está marcada la opción
*Keep the variables of the file*, y lo abre en el panel de resultados, donde los colores, las flechas, las
sondas y los perfiles funcionan con las variables nuevas.

### Escribir expresiones

Una salida por línea:

```
NOMBRE [UNIDAD] = expresión
```

- `NOMBRE`: hasta 16 caracteres (un límite de SELAFIN); se permiten espacios. `[UNIDAD]` es opcional.
- Las líneas que empiezan con `#` son comentarios. Varias salidas pueden compartir una línea separadas por
  `;`.
- Cada salida se puede usar por su nombre en las líneas siguientes (los espacios pasan a `_`:
  `SLOPE X` → `SLOPE_X`).
- Dos salidas llamadas `… X` y `… Y` (o `… U` y `… V`) se muestran como un **vector con flechas**. QGIS solo
  lee como vectores los pares `… U / V` y `… ALONG X / Y`, así que `SLOPE X` y `SLOPE Y` se escriben en el
  archivo como `SLOPE ALONG X` y `SLOPE ALONG Y` (o `… U / V` cuando pasarían de 16 caracteres); el informe
  lo indica. En las líneas siguientes del mismo cálculo conservan sus nombres (`SLOPE_X`).
- Una línea sin `NOMBRE =` se llama `RESULT`.

### Variables

| Escribe | Significa |
| --- | --- |
| `U`, `V`, `H`, `S`, `B`, `F`, `Q`, `M`, `W` … | Símbolos de TELEMAC: velocidad U y V, tirante, superficie libre, fondo, Froude, caudal unitario, magnitud de la velocidad, … (archivos en inglés o en francés) |
| `WATER_DEPTH`, `VELOCITY_U`, `MY_TRACER` | Cualquier variable por su nombre, con los espacios como `_` |
| `water_depth`, `velocity_u` | Nombres canónicos, iguales en archivos en inglés y en francés |
| `x`, `y` | Coordenadas de los nodos |
| `t` | Tiempo del cuadro, en segundos |
| `area` | Área alrededor de cada nodo (un tercio de los triángulos que lo tocan) |
| `g`, `pi`, `e` | 9.81, π, e |

### Operadores y funciones

**Aritmética**: `+ - * / **` (potencia) `%` (residuo) y paréntesis. **Comparaciones** `> < >= <= == !=`,
`and`, `or`, `where(condición, a, b)` (el SI de otras calculadoras) y `a if condición else b`.

**Matemáticas**: `abs sqrt exp log log10 sin cos tan arcsin arccos arctan arctan2 hypot min max minimum
maximum clip floor ceil sign degrees radians isnan`.

**Cálculo espacial sobre la malla**:

| Función | Resultado |
| --- | --- |
| `ddx(f)`, `ddy(f)` | ∂f/∂x, ∂f/∂y |
| `grad(f)` | magnitud del gradiente; por ejemplo `grad(S)` = pendiente de la superficie del agua |
| `div(u, v)` | divergencia ∂u/∂x + ∂v/∂y |
| `curl(u, v)` o `vorticity(u, v)` | vorticidad vertical ∂v/∂x − ∂u/∂y (1/s) |
| `strain(u, v)` | tasa de deformación √(2u<sub>x</sub>² + 2v<sub>y</sub>² + (u<sub>y</sub> + v<sub>x</sub>)²) |
| `okubo_weiss(u, v)` | s<sub>n</sub>² + s<sub>s</sub>² − ω²; negativo dentro de los remolinos |
| `laplacian(f)` | ∇²f, calculado como div(grad f) |
| `smooth(f)` | una pasada de suavizado nodal |
| `mag(u, v)`, `bearing(u, v)` | magnitud; dirección hacia la que apunta el vector, en grados desde el norte, en sentido horario |

**Tiempo**:

| Función | Resultado |
| --- | --- |
| `tmax(f)`, `tmin(f)`, `tmean(f)`, `tstd(f)` | máximo, mínimo, media y desviación estándar sobre todos los cuadros |
| `tsum(f)`, `tint(f)` | suma sobre los cuadros; integral en el tiempo (trapecial, unidad de f × s) |
| `tpeak(f)` | tiempo (s) en que f es máxima |
| `tfirst(f)`, `tlast(f)` | f en el primer o el último cuadro |
| `ddt(f)` | tasa de cambio en el tiempo (diferencias centrales) |
| `wet(eps)` | 1 donde el tirante supera eps (por defecto 0.01 m), 0 en otro caso |

Una salida que usa solo reducciones en el tiempo (`HMAX = tmax(H)`) es igual en todos los cuadros; cuando no
se conservan las variables originales, se escribe una sola vez.

### Ejemplos

```
# hidráulica
FROUDE = mag(U, V) / sqrt(g * maximum(H, 1e-6))
SHEAR [PA] = 1000 * g * 0.03**2 * (U**2 + V**2) / maximum(H, 0.01)**(1/3)
HV [M2/S] = H * mag(U, V)
SPECIFIC ENERGY [M] = H + (U**2 + V**2) / (2 * g)

# mapas de extremos
HMAX [M] = tmax(H)
TIME HMAX [S] = tpeak(H)
VMAX [M/S] = tmax(mag(U, V))
WET TIME [S] = tint(wet(0.05))
DEPTH CHANGE [M] = tlast(H) - tfirst(H)

# estructura del flujo
VORTICITY [1/S] = where(H > 0.01, curl(U, V), 0)
SLOPE X = ddx(S)
SLOPE Y = ddy(S)
RISE [M/S] = ddt(S)
```

`SHEAR` usa n de Manning = 0.03 (τ = ρ g n² |U|² / h<sup>1/3</sup>); cambia el 0.03 o usa en su lugar la
variable `BOTTOM FRICTION` del archivo (`N`).

### Precisión

Las derivadas espaciales usan la interpolación lineal sobre triángulos que usa el propio TELEMAC: el
gradiente es constante en cada triángulo y se promedia hacia los nodos, ponderado por el área. Es exacto
para campos lineales, de segundo orden en mallas regulares y de primer orden en mallas irregulares. En los
bordes húmedo/seco la velocidad cae a cero, así que enmascara las derivadas con `where(H > 0.01, …, 0)`. En
la frontera las derivadas son unilaterales. Los operadores se comprueban en las pruebas automáticas contra
campos exactos (campos lineales, rotación de cuerpo rígido, deformación pura, campos cuadráticos) y en su
convergencia.

![Vorticidad del flujo en una curva de 180°](../img/bend_vorticity.png)

### Opciones

- **Frames**: vacío = todos; `-1` = el último; `0:10`; `0:100:5` (uno de cada 5).
- **Keep the variables of the file**: desactivada, genera un archivo pequeño solo con las variables nuevas.
- *Advanced*: el tirante que usa `wet()`, y el SRC que se guarda como `.prj` junto a la salida (por defecto
  se copia el `.prj` de la entrada).

La calculadora trabaja con un solo archivo. Las diferencias entre dos escenarios sobre la misma malla están
planeadas.

## Integrales y flujos

*Complementos → TELEMAC Toolkit → Integrals and fluxes*. La salida es una tabla con una fila por línea o
polígono y por paso de tiempo, que puedes graficar, unir o exportar a CSV/Excel.

| Integral | Qué das | Columnas |
| --- | --- | --- |
| **Discharge through section lines** (caudal a través de líneas de sección) | una capa de líneas | caudal (m³/s), longitud, longitud mojada |
| **Flux of an expression** (flujo de una expresión) | líneas y una expresión, por ejemplo una concentración de trazador `T1` | flujo (expresión × m³/s) |
| **Line integral and mean** (integral de línea y media) | líneas y una expresión, por ejemplo `mag(U, V)` | integral, longitud, media |
| **Area integral and mean** (integral de área y media) | polígonos (o nada = toda la malla) y una expresión: `H` da el volumen, `1` el área, `wet(0.01)` el área mojada | integral, área, media, área mojada |

**Signo de un caudal**: positivo de la derecha de la línea hacia su izquierda, mirando en el sentido en que
se dibujó la línea. Dibuja las secciones de la **margen izquierda a la margen derecha, mirando aguas abajo**
(como las secciones transversales de HEC-RAS) y el flujo hacia aguas abajo sale positivo.

**Precisión**: la línea se corta en cada arista de la malla que cruza; entre dos cortes el tirante y la
velocidad son lineales y su producto se integra de forma exacta. En una curva de 180° corrida en TELEMAC-2D
con 0.05 m³/s de entrada, cinco secciones dan de 0.04999 a 0.05001 m³/s
([resumen de la corrida](../../tests/data/demo_bend/run_summary.txt)). Las integrales de área cuentan los
elementos cuyo centro está dentro del polígono, así que el borde del polígono corta la malla con la
resolución de los elementos; la tabla indica el área que realmente se contó.

Las líneas y los polígonos en otro SRC se transforman al SRC de la malla (el del `.prj`, o el que indiques).

[← Ver resultados](03-resultados.md) · [Contenido](README.md) · [Siguiente: El flujo de trabajo de la malla →](05-flujo-de-malla.md)
