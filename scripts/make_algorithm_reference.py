"""Write docs/guide/09-processing-reference.md from the algorithms themselves (names, help, parameters).

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


def default_text(p):
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
        return "yes" if d else "no"
    if isinstance(d, str) and "\n" in d:
        return "(example text)"
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
    out = ["# 9. Processing tools reference", "",
           "Every tool of the plugin is in the Processing toolbox under **TELEMAC Toolkit**, so it can run in batch "
           "mode, in the Model Designer or from Python (`processing.run(\"telemactoolkit:<id>\", {...})`). This page "
           "is generated from the tools themselves (`scripts/make_algorithm_reference.py`).", ""]
    order = ["SELAFIN", "Mesh"]
    for g in sorted(groups, key=lambda k: order.index(k) if k in order else 9):
        out += ["## %s" % ("Results (group SELAFIN)" if g == "SELAFIN" else "Mesh (group Mesh)"), ""]
        for alg in sorted(groups[g], key=lambda a: a.displayName()):
            out += ["### %s" % alg.displayName(), "", "Id: `telemactoolkit:%s`" % alg.name(), ""]
            help_ = alg.shortHelpString().strip()
            if help_:
                out += help_markdown(help_)
            out += ["| Parameter | Key | Type | Default |", "| --- | --- | --- | --- |"]
            for p in alg.parameterDefinitions():
                adv = " *(advanced)*" if p.flags() & Qgis.ProcessingParameterFlag.Advanced else ""
                opt = " *(optional)*" if p.flags() & Qgis.ProcessingParameterFlag.Optional else ""
                kind = KIND.get(p.type(), p.type())
                if p.type() == "enum":
                    kind = "choice: " + "; ".join(p.options())
                desc = p.description().replace("|", "/")
                out.append("| %s%s%s | `%s` | %s | %s |" % (desc, opt, adv, p.name(), kind.replace("|", "/"),
                                                            default_text(p).replace("|", "/")))
            out.append("")
    path = os.path.join(ROOT, "docs", "guide", "09-processing-reference.md")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n[← Formats, coordinates and boundary conditions](08-formats-coordinates-boundaries.md) · "
                "[Contents](README.md) · [Next: Scripting with Python →](10-python.md)\n")
    print("wrote", path, sum(len(v) for v in groups.values()), "tools")


if __name__ == "__main__":
    main()
