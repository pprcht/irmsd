# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "irmsd"
copyright = "2025, Philipp Pracht, Tobias Kaczun"
author = "Philipp Pracht, Tobias Kaczun"


import os

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import sys

sys.path.insert(0, os.path.abspath("../../src"))  # src layout
# sys.path.insert(0, os.path.abspath("../.."))  # project root

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
]

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]


hintersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# Sphinx defaults to automatically resolve *unresolved* labels using all your Intersphinx mappings.
# This behavior has unintended side-effects, namely that documentations local references can
# suddenly resolve to an external location.
# See also:
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#confval-intersphinx_disabled_reftypes
# If this should become a problem, uncomment the following line:
# intersphinx_disabled_reftypes = ["*"]

napoleon_numpy_docstring = True
napoleon_preprocess_types = True

# # Hide module names in description unit titles (such as .. function::).
# add_module_names = Falsetml_static_path = ["_static"]
templates_path = ["_templates"]
exclude_patterns = []

# autosummary: automatically create stub pages
autosummary_generate = True
autosummary_imported_members = True
typehints_fully_qualified = True

autodoc_default_options = {
    "members": True,  # document class/method/function members
    "undoc-members": True,  # even if no docstring
    "private-members": True,  # set to True if needed
    "show-inheritance": False,
}


def setup(app):
    """Add custom CSS overrides."""
    app.add_css_file("theme_override.css")
    app.add_css_file("svg_scaling.css")


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "sphinx_book_theme"

html_theme = "shibuya"

html_static_path = [
    "_static",
]
#

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
