"""Write docs/guide/09-processing-reference.md (English) and docs/guia/09-referencia-procesos.md (Spanish)
from the algorithms themselves (names, help, parameters). The help texts and parameter names are the ones QGIS
shows, in English, in both pages; the headings, introduction and table headers are translated.

    QT_QPA_PLATFORM=offscreen <qgis-python> scripts/make_algorithm_reference.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
for _p in (os.path.join(sys.prefix, "share", "qgis", "python", "plugins"),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.append(_p)

from qgis.core import Qgis, QgsApplication  # noqa: E402

KIND = {
    "file": "file", "crs": "CRS", "enum": "choice", "boolean": "yes/no", "number": "number", "string": "text",
    "source": "vector layer", "sink": "output layer", "fileDestination": "output file", "multilayer": "layers",
    "field": "field", "point": "point",
}


def help_markdown(text):
    """Help strings use one line per paragraph and '• ' bullets; Markdown needs blank lines and '- '."""
    md, in_list = [], False
    for line in text.split("\n"):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith("• "):
            if not in_list:
                md.append("")
            md.append("- " + line[2:])
            in_list = True
        else:
            md += ["", line]
            in_list = False
    return md[1:] + [""] if md and md[0] == "" else md + [""]


KIND_ES = {
    "file": "archivo", "crs": "SRC", "enum": "opción", "boolean": "sí/no", "number": "número", "string": "texto",
    "source": "capa vectorial", "sink": "capa de salida", "fileDestination": "archivo de salida",
    "multilayer": "capas", "field": "campo", "point": "punto",
}

TEXT = {
    "en": {
        "kind": KIND, "yes": "yes", "no": "no", "example": "(example text)", "choice": "choice: ",
        "optional": "optional", "advanced": "advanced",
        "header": ["Parameter", "Key", "Type", "Default"],
        "results": "Results (group SELAFIN)", "mesh": "Mesh (group Mesh)",
        "title": "# 9. Processing tools reference",
        "intro": ("Every tool of the plugin is in the Processing toolbox under **TELEMAC Toolkit**, so it can run in "
                  "batch mode, in the Model Designer or from Python (`processing.run(\"telemactoolkit:<id>\", {...})`). "
                  "This page is generated from the tools themselves (`scripts/make_algorithm_reference.py`)."),
        "id": "Id",
        "footer": ("[← Formats, coordinates and boundary conditions](08-formats-coordinates-boundaries.md) · "
                   "[Contents](README.md) · [Next: Scripting with Python →](10-python.md)"),
        "out": ("docs", "guide", "09-processing-reference.md"),
    },
    "es": {
        "kind": KIND_ES, "yes": "sí", "no": "no", "example": "(texto de ejemplo)", "choice": "opción: ",
        "optional": "opcional", "advanced": "avanzado",
        "header": ["Parámetro", "Clave", "Tipo", "Valor por defecto"],
        "results": "Resultados (grupo SELAFIN)", "mesh": "Malla (grupo Mesh)",
        "title": "# 9. Referencia de las herramientas de Procesos",
        "intro": ("Todas las herramientas del plugin están en la caja de herramientas de Procesos, bajo **TELEMAC "
                  "Toolkit**, así que se pueden ejecutar por lotes, en el Diseñador de modelos o desde Python "
                  "(`processing.run(\"telemactoolkit:<id>\", {...})`). Esta página se genera a partir de las propias "
                  "herramientas (`scripts/make_algorithm_reference.py`). **Los nombres de las herramientas, sus "
                  "textos de ayuda y los nombres de los parámetros están en inglés, tal como los muestra QGIS**, "
                  "para que los encuentres igual en pantalla; los títulos, las columnas y los tipos están en "
                  "español. La explicación de cada función está en las páginas de esta guía."),
        "id": "Id",
        "footer": ("[← Formatos, coordenadas y condiciones de frontera](08-formatos-coordenadas-fronteras.md) · "
                   "[Contenido](README.md) · [Siguiente: Scripts con Python →](10-python.md)"),
        "out": ("docs", "guia", "09-referencia-procesos.md"),
    },
}


def default_text(p, t):
    d = p.defaultValue()
    if d is None or d == "":
        return ""
    if p.type() == "enum":
        try:
            opts = p.options()
            return opts[int(d)] if str(d).isdigit() and int(d) < len(opts) else str(d)
        except Exception:
            return str(d)
    if p.type() == "boolean":
        return t["yes"] if d else t["no"]
    if isinstance(d, str) and "\n" in d:
        return t["example"]
    return str(d)


def main():
    app = QgsApplication([], False)
    app.initQgis()
    from processing.core.Processing import Processing

    Processing.initialize()
    from telemac_toolkit.tt_processing.provider import TelemacProvider

    prov = TelemacProvider()
    QgsApplication.processingRegistry().addProvider(prov)
    groups = {}
    for alg in prov.algorithms():
        groups.setdefault(alg.group(), []).append(alg)
    order = ["SELAFIN", "Mesh"]
    for lang, t in TEXT.items():
        out = [t["title"], "", t["intro"], ""]
        for g in sorted(groups, key=lambda k: order.index(k) if k in order else 9):
            out += ["## %s" % (t["results"] if g == "SELAFIN" else t["mesh"]), ""]
            for alg in sorted(groups[g], key=lambda a: a.displayName()):
                out += ["### %s" % alg.displayName(), "", "%s: `telemactoolkit:%s`" % (t["id"], alg.name()), ""]
                help_ = alg.shortHelpString().strip()
                if help_:
                    out += help_markdown(help_)
                out += ["| %s | %s | %s | %s |" % tuple(t["header"]), "| --- | --- | --- | --- |"]
                for p in alg.parameterDefinitions():
                    adv = " *(%s)*" % t["advanced"] if p.flags() & Qgis.ProcessingParameterFlag.Advanced else ""
                    opt = " *(%s)*" % t["optional"] if p.flags() & Qgis.ProcessingParameterFlag.Optional else ""
                    kind = t["kind"].get(p.type(), p.type())
                    if p.type() == "enum":
                        kind = t["choice"] + "; ".join(p.options())
                    desc = p.description().replace("|", "/")
                    out.append("| %s%s%s | `%s` | %s | %s |" % (desc, opt, adv, p.name(), kind.replace("|", "/"),
                                                                default_text(p, t).replace("|", "/")))
                out.append("")
        path = os.path.join(ROOT, *t["out"])
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(out) + "\n" + t["footer"] + "\n")
        print("wrote", path, sum(len(v) for v in groups.values()), "tools")


if __name__ == "__main__":
    main()
