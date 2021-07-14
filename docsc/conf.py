import os
import sys

sys.path.insert(0, os.path.abspath(".."))

extensions = [
    "breathe",
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx_click.ext",
]
# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Breathe Configuration
breathe_default_project = "IBS"

# General information about the project.
project = u"ibs"
copyright = u"2021, Tom Mertens"
author = u"Tom Mertens"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
try:
    import sphinx_rtd_theme

    html_theme = "sphinx_rtd_theme"
except ImportError:
    print("Sphinx html theme 'sphinx_rtd_theme' not found. Using 'classic' instead.")
    html_theme = "classic"
