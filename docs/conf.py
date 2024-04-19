"""Configuration file for the Sphinx documentation builder."""

# Project information
project = "mzSpecLib"
author = "HUPO-PSI"
copyright = "Creative Commons Attribution 4.0 International"
github_project_url = "https://github.com/HUPO-PSI/mzspeclib-py"
github_doc_root = "https://github.com/HUPO-PSI/mzspeclib-py/tree/main/docs"

# General configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_click.ext",
    "myst_parser",
]
source_suffix = [".rst"]
master_doc = "index"
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

# Autodoc options
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autoclass_content = "init"

# Intersphinx options
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "psims": ("https://mobiusklein.github.io/psims/docs/build/html/", None),
    "pyteomics": ("https://pyteomics.readthedocs.io/en/stable/", None),
}


def setup(app):
    """Set up the Sphinx application."""
    config = {"enable_eval_rst": True}  # noqa: F841
