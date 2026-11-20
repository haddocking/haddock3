"""Static assets

Static assets for analysis report.

The Mol* viewer (MIT license, https://molstar.org) files were downloaded from
https://cdn.jsdelivr.net/npm/molstar@5.13.1/build/viewer/molstar.js
and
https://cdn.jsdelivr.net/npm/molstar@5.13.1/build/viewer/molstar.css
"""

import importlib.resources as importlib_resources

haddock_ui_path = importlib_resources.files("haddock.libs.assets")

MOLSTAR_VERSION = "5.13.1"
"""Version of the Mol* files in this folder, also used for the CDN links."""
