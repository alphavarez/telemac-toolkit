# 2. Primeros pasos

## Un resultado en cinco minutos

1. *Open TELEMAC result…* en la barra de herramientas, o arrastra un archivo `.slf` a la ventana de QGIS.
2. Si el plugin no puede saber cuáles son las coordenadas, te pregunta: **coordenadas locales** (un canal de
   laboratorio o un canal esquemático) o un **SRC de mapa**. Deja marcada la opción *Save the choice next to
   the file* (guardar la elección junto al archivo) para que no vuelva a preguntar
   ([detalles](03-resultados.md#coordenadas-de-un-resultado)).
3. Se abre el panel **TELEMAC Results** y el mapa muestra el tirante (profundidad del agua) con flechas de
   velocidad.
4. Pulsa **play**, o arrastra el deslizador. El Controlador temporal de QGIS se mueve con él.
5. Cambia *Map shows* (qué muestra el mapa) a *Velocity magnitude*; abre *Arrow settings* para hacer las
   flechas más densas, más separadas, más largas o más cortas.
6. Activa la **sonda** (pestaña *Time series*) y haz clic en dos puntos de la malla: se dibujan sus series de
   tiempo. Activa la herramienta de **perfil**, haz clic a lo largo de una línea y termina con clic derecho:
   el fondo y la superficie del agua a lo largo de la línea siguen al reloj.

## Un canal de laboratorio en diez minutos

1. Barra de herramientas → **Mesh workflow panel** → *New…*, guarda `flume.gpkg` y responde **Sí** a
   coordenadas locales.
2. Paso 1 → *Channel or flume template…*:
   - Planta *Straight* (recta), ancho 1 m, largo 20 m, tamaño de elemento 0.05 m, *Structured*.
   - Cota del fondo a la entrada 0.30 m, pendiente 0.001, sección rectangular.
   - Caudal de entrada 0.05 m³/s, nivel de salida 0.42 m. *Run*.
3. En el mapa aparecen el dominio, las líneas de entrada (Q) y de salida (H) y una malla estructurada; la
   fórmula del fondo se guarda en el proyecto.
4. Paso 4 → formato *TELEMAC*, archivo de malla `geo_flume.slf`, **Generate mesh**. El registro enumera la
   malla, los números de frontera que usará TELEMAC y las líneas que debes pegar en el archivo de
   configuración (`PRESCRIBED FLOWRATES`, `PRESCRIBED ELEVATIONS`).
5. El plugin escribe `geo_flume.slf`, `geo_flume.cli` y (si lo pides) una plantilla `.liq`. Escribe un
   archivo de configuración (*steering file*) con ellos y ejecuta `telemac2d.py t2d_flume.cas`.
6. Abre el resultado `r2d_flume.slf` como en la sección anterior.

El caso completo, con su archivo de configuración, está en `scripts/make_flume_demo.py` y
`tests/data/demo_flume/`. Una versión curva (curva de 180°) está en `scripts/make_bend_demo.py` y
`tests/data/demo_bend/`.

[← Instalación](01-instalacion.md) · [Contenido](README.md) · [Siguiente: Ver resultados →](03-resultados.md)
