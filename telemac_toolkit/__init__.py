"""TELEMAC Toolkit for QGIS."""


def classFactory(iface):  # noqa: N802 - QGIS API name
    from .plugin import TelemacToolkitPlugin

    return TelemacToolkitPlugin(iface)
