# 11. Solución de problemas y preguntas frecuentes

**Generate mesh se queda en 6 % (o en cualquier porcentaje bajo) y nunca termina.**
Antes de la 0.5.0 casi siempre era un tamaño de elemento demasiado pequeño para el tamaño del dominio: la
barra de progreso no puede avanzar mientras el programa coloca millones de nodos. Los tamaños (tamaño por
defecto y los de las áreas, puntos y líneas de refinamiento) están en unidades del mapa, así que 0.05 o 0.2
en un dominio de 100 km de ancho es imposible. Desde la 0.5.0 el plugin estima primero el número de nodos
y se detiene con un mensaje que da la estimación, el tamaño del dominio y el tamaño mínimo pedido. Sube los
tamaños (un tramo de río de 5 km suele estar bien con 10–25 m en el cauce y 50–100 m en la planicie) o
reduce las áreas refinadas. Si de verdad necesitas más de 2 millones de nodos (un tramo largo de un río grande con resolución fina),
sube el **Node limit** del paso 4.
Revisa también la unidad: las coordenadas en un SRC desconocido se leen como metros.

**Abrí un resultado y el mapa está en blanco.**
El resultado está fuera de la vista actual o en un SRC que no coincide con el del proyecto. Pulsa el **botón
de zoom** del panel de resultados (o clic derecho en la capa → *Zoom a la(s) capa(s)*). Si todavía no se ve,
revisa sus coordenadas con el **botón de coordenadas**: un modelo de laboratorio debe estar en coordenadas
locales y un modelo de campo en su SRC de mapa. Desde la 0.4.0 el plugin ajusta el mapa por sí mismo al
abrir un resultado.

**La malla aparece en un lugar equivocado, o la barra de escala muestra miles de kilómetros.**
Al resultado se le asignó un SRC equivocado (por ejemplo, metros leídos como grados). Usa el botón de
coordenadas y marca *Save the choice next to the file*.

**La barra de estado dice "SRC desconocido" (*Unknown CRS*).**
Es la etiqueta de QGIS para un SRC del que no tiene nombre ni código. Desde la 0.5.0 el plugin registra su
cuadrícula local de laboratorio como SRC de usuario de QGIS la primera vez que la usa, así que QGIS muestra
`USER:100000 - Local metric grid (laboratory)` (el número puede variar) y lista las capas con ese código. A un
SRC de mapa sin código EPSG se le ofrece uno al abrir el proyecto (véase el paso 1 del
[flujo de malla](05-flujo-de-malla.md)). En cualquier caso las coordenadas son correctas.

**Generate mesh se detiene con "OSError: [Errno 22] Invalid argument" (o "cannot be written").**
En Windows no se puede sobrescribir un archivo que QGIS tiene cargado (una capa de malla de una corrida
anterior). Desde la 0.5.0 el panel quita esa capa antes de mallar y carga la malla nueva en su lugar; la
herramienta de Processing hace lo mismo, o escribe la malla al lado (`nombre_2.slf`) cuando corre en segundo
plano, y lo indica. Antes de la 0.5.0: quita la capa de la malla anterior del panel de capas y vuelve a correr.

**El borde de la malla no coincide exactamente con el contorno del dominio: un poco adentro aquí, un poco afuera allá.**
El contorno se remuestrea al tamaño de elemento, y cada arista del borde es una cuerda recta entre dos puntos
del dibujo; donde el dibujo se curva suavemente entre ellos (un vértice que gira menos de 30°), la arista corta
la curva: hacia adentro en una curva convexa, hacia afuera en una cóncava, sobre todo a lo largo de los lados,
lejos de las esquinas. Antes de la 0.5.0 la diferencia podía llegar a cerca del 10 % del tamaño de elemento.
Desde la 0.5.0 esos vértices quedan como nodos donde hace falta, así que la diferencia queda por debajo de
**Follow the outline within** en el paso 4 (*auto*: 5 % del tamaño de elemento, p. ej. 1 m para elementos de
20 m); el registro da la mayor que queda. Pon una distancia en metros (0.01 m conserva prácticamente todos los
vértices) para un ajuste más fiel, o dibuja ese tramo como línea de quiebre (breakline), que siempre conserva
todos sus vértices. La 0.5.0 también corrige los contornos digitalizados con muchos vértices muy juntos (uno
cada metro o dos, en zigzag): el generador podía quedarse colgado en *Triangulating*.

**Re-mesh inside polygons dice que "no pudo unir limpiamente" los elementos anteriores.**
Los elementos nuevos de un área deben coincidir nodo por nodo con los anteriores; cuando no pueden (un
polígono que corta elementos muy delgados, un tamaño muy distinto de las aristas del contorno del área), la
herramienta dice dónde y no cambia nada. Usa un tamaño más cercano al de las aristas del contorno (el
registro da su mediana) o mueve un poco el polígono.

**A escala 1:2 las esquinas de la malla quedan a unos milímetros de las esquinas dibujadas.**
SELAFIN guarda las coordenadas en precisión simple (7 cifras significativas): a 125 000 m una coordenada se
redondea a 7.8 mm, a 189 000 m a 15.6 mm, así que antes de la 0.5.0 un nodo podía quedar hasta a 8 mm del
dibujo. Desde la 0.5.0 una malla a más de 4 km del origen se escribe respecto a un origen entero (IPARAM 3 y 4,
que TELEMAC y QGIS leen), y los nodos conservan su posición a 0.01 mm. Vuelve a generar la malla.

**Las flechas son enormes, están amontonadas o se traslapan.**
*Arrow settings*: colocación *Regular grid*, una separación que convenga al canal (15–25 px para canales
angostos), longitud *Proportional* con *Longest arrow* en auto. Sube *Hide below* para quitar las flechas en
agua en reposo.

**No hay flechas dentro de un canal angosto.**
Los puntos de la cuadrícula caen fuera del canal con este zoom: baja la separación, acércate o usa *Every
mesh node*.

**Las flechas muestran siempre la velocidad, aunque elija otra variable.**
En las versiones anteriores a la 0.4.3 solo se podía dibujar como flechas la velocidad, y el cambio de
variable en *Map shows* no movía las flechas. Desde la 0.4.3, la lista *Arrows* enumera todos los vectores
del archivo (velocidad, caudal unitario, viento) y las flechas siguen al vector que eliges en *Map shows*. Si
eliges un escalar (tirante, superficie libre), no tiene dirección y las flechas se quedan como estaban. Si
tu archivo trae un vector con nombres que QGIS no reconoce como par (por ejemplo `VENT X` / `VENT Y`), se
muestra como dos mapas separados, sin flechas.

**TELEMAC se detiene con "SL: MORE PRESCRIBED ELEVATIONS ARE REQUIRED" (o lo mismo para caudales).**
Las fronteras líquidas del archivo de configuración o del archivo `.liq` están en el orden equivocado. Usa
los números y las líneas del archivo de configuración que el plugin imprime al escribir el `.cli`, y la
plantilla `.liq` que escribe (véanse los
[números de las fronteras líquidas](08-formatos-coordenadas-fronteras.md#números-de-las-fronteras-líquidas)).
Las versiones anteriores a la 0.4.0 empezaban la numeración en el primer nodo de frontera en lugar de en el
punto de inicio de TELEMAC, así que algunas mallas (una curva en forma de U, por ejemplo) quedaban con el
orden equivocado: escribe de nuevo el `.cli` y el `.liq` con la 0.4.0 o posterior.

**TELEMAC ignora mi FRICTION COEFFICIENT.**
El archivo de geometría tiene una variable `BOTTOM FRICTION` (zonas de rugosidad o una fricción por
defecto) y TELEMAC-2D la lee en su lugar (el listado dice `FONSTR : FRICTION COEFFICIENTS READ IN THE
GEOMETRY FILE`). Quita las zonas y la fricción por defecto para usar el valor del archivo de configuración.

**"The formula uses s or n but no channel axis line is given."**
El proyecto necesita una línea de eje en su capa de canales (rol `axis`); con varios ejes, ajusta
`bottom_axis` en los ajustes del proyecto con el nombre del canal. En *Edit node values*, indica la capa del
eje.

**La calculadora dice que un nombre es desconocido.**
Usa el símbolo de TELEMAC (`H`, `U`…) o el nombre de la variable con `_` en lugar de espacios
(`WATER_DEPTH`). El cuadro de expresiones que se abre desde el panel de resultados enumera los nombres del
archivo actual.

**La calculadora dice que el nombre de la salida es demasiado largo.**
Los nombres de variable de SELAFIN tienen como máximo 16 caracteres.

**Las derivadas se ven ruidosas en el borde del agua.**
La velocidad cae a cero en los nodos secos. Enmascara: `where(H > 0.01, curl(U, V), 0)`.

**Un resultado 3D se abre como 2D.**
QGIS no puede leer archivos SELAFIN 3D; el plugin muestra un promedio en la profundidad o un plano escrito
en `<nombre>_2d_<opción>.slf`. Abre el resultado de nuevo para elegir otro plano.

**"This is a sub-domain file of a parallel run".**
Primero une las particiones con `gretel` de TELEMAC.

**El generador de mallas informa de una línea de restricción que no toca ningún nodo de frontera / una línea
de frontera que no toca ningún nodo.**
Ensambla la línea al borde del dominio (barra de autoensamblado) y comprueba que esté a menos de medio
elemento de la frontera.

**¿Dónde reporto un problema?**
En los *Issues* del repositorio de GitHub, con la versión de QGIS, el mensaje de *Ver → Paneles → Mensajes
de registro → TELEMAC Toolkit* y, si es posible, un archivo pequeño que muestre el problema.

[← Python](10-python.md) · [Contenido](README.md) · [Siguiente: Condiciones de frontera y archivos de corrida →](12-archivos-de-corrida.md)
