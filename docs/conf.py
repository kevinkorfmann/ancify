"""Sphinx configuration for ancify documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "ancify"
copyright = "2025, ancify contributors"
author = "ancify contributors"
release = "1.0.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "ancify — Ancestral Allele Polarization"
html_logo = "_static/gpu.png"
html_css_files = ["custom.css"]

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "logo_only": True,
    "prev_next_buttons_location": "both",
    "style_external_links": True,
    "style_nav_header_background": "#2c3e50",
    "titles_only": False,
}

html_context = {
    "display_github": True,
    "github_user": "kevinkorfmann",
    "github_repo": "ancify",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "substitution",
    "tasklist",
]

myst_heading_anchors = 3
