"""Static assets

Static assets for analysis report.

Files where downloaded from
https://cdn.jsdelivr.net/npm/@i-vresse/haddock3-ui@~0.3.0/dist/index.css
and
https://cdn.jsdelivr.net/npm/@i-vresse/haddock3-ui@~0.3.0/dist/report.bundle.js
"""
import importlib.resources as importlib_resources

haddock_ui_path = importlib_resources.files("haddock.libs.assets")
