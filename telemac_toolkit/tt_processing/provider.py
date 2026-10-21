import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithms import ALGORITHMS
from .mesh_algorithms import MESH_ALGORITHMS
from .calc_algorithms import CALC_ALGORITHMS


class TelemacProvider(QgsProcessingProvider):
    def id(self):
        return "telemactoolkit"

    def name(self):
        return "TELEMAC Toolkit"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(os.path.dirname(__file__)), "icons", "telemac.svg"))

    def loadAlgorithms(self):
        for cls in ALGORITHMS + CALC_ALGORITHMS + MESH_ALGORITHMS:
            self.addAlgorithm(cls())
