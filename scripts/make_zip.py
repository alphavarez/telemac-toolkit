"""Build the installable plugin zip: dist/telemac_toolkit-<version>.zip

The zip holds one top-level folder, telemac_toolkit/, as QGIS expects.
"""

import configparser
import os
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLUGIN = os.path.join(ROOT, "telemac_toolkit")
SKIP_DIRS = {"__pycache__", ".pytest_cache"}
SKIP_EXT = {".pyc", ".pyo"}


def main():
    meta = configparser.ConfigParser()
    meta.read(os.path.join(PLUGIN, "metadata.txt"), encoding="utf-8")
    version = meta["general"]["version"]
    os.makedirs(os.path.join(ROOT, "dist"), exist_ok=True)
    out = os.path.join(ROOT, "dist", "telemac_toolkit-%s.zip" % version)
    n = 0
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for base, dirs, files in os.walk(PLUGIN):
            dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
            for f in sorted(files):
                if os.path.splitext(f)[1] in SKIP_EXT:
                    continue
                full = os.path.join(base, f)
                z.write(full, os.path.relpath(full, ROOT))
                n += 1
    print("%s (%d files, %.0f kB)" % (out, n, os.path.getsize(out) / 1024))
    return out


if __name__ == "__main__":
    main()
