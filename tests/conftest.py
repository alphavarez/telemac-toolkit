import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# core library lives inside the plugin folder
sys.path.insert(0, os.path.join(ROOT, "telemac_toolkit"))
sys.path.insert(0, ROOT)
