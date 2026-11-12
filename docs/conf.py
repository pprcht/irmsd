# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import subprocess
import sys
from pathlib import Path

# -- Project information -----------------------------------------------------

project = "irmsd"
copyright = "2025, Philipp Pracht, Tobias Kaczun"
author = "Philipp Pracht, Tobias Kaczun"

# -- General configuration ---------------------------------------------------

# autodoc imports the installed package, which loads the compiled library, so
# the docs are built against `pip install .` and not against the source tree.

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "sphinxarg.ext",
]

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]
# anchors for "#section" links between the Markdown pages
myst_heading_anchors = 3

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

napoleon_numpy_docstring = True
napoleon_preprocess_types = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "flowchart"]

# autosummary: automatically create stub pages
autosummary_generate = True
autosummary_imported_members = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
    "show-inheritance": False,
}

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Options for HTML output -------------------------------------------------

html_theme = "shibuya"
html_title = "iRMSD"
html_static_path = ["_static"]

# the interactive flowchart (algorithm.md) is a standalone page, generated
# here and copied to the output root
_flowchart = Path(__file__).parent / "flowchart"
subprocess.run([sys.executable, str(_flowchart / "build.py")], check=True)
html_extra_path = [str(_flowchart / "irmsd-flowchart.html")]
