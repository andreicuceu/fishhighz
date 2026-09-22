"""FishHighz documentation configuration; imports never run forecasts."""

from importlib.metadata import version as package_version

project = "FishHighz"
author = "FishHighz contributors"
release = package_version("fishhighz")
version = release
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
nitpicky = True
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "FishHighz documentation"
autosummary_generate = True
napoleon_numpy_docstring = True
napoleon_google_docstring = False
# Keep NumPy type descriptions (shapes, choices, optional) as text, not class links.
napoleon_use_param = False
napoleon_use_rtype = False
autodoc_typehints = "none"
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_heading_anchors = 6
