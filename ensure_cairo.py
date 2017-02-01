# Called by tox to add Homebrew's cairo to the testing virtualenv
# without otherwise relying on the global site-packages.

import os
import subprocess
import sys

if sys.platform == "darwin":
    cairo_path = subprocess.check_output(
        ["brew", "--prefix", "py2cairo"]
    )
    cairo_path = cairo_path.strip()
    cairo_path += "/lib/python2.7/site-packages/cairo"
else:
    # debuntu
    cairo_path = "/usr/lib/python2.7/dist-packages/cairo"

cairo_dest = os.path.join(sys.argv[1], "cairo")
if not os.path.exists(cairo_dest):
    os.symlink(cairo_path, cairo_dest)
