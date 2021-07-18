import os
import sys

sys.path.insert(0, os.path.abspath(".."))

extensions = [
    "matplotlib.sphinxext.mathmpl",
    "matplotlib.sphinxext.plot_directive",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx.ext.mathjax",
    "breathe",
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx_click.ext",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.napoleon",
    "nbsphinx",
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
exclude_patterns = ["build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

plot_include_source = True
nofigs = False

# -- Options for HTML output -------------------------------------------
plot_formats = ["png", "hires.png", "pdf", "svg"]

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
try:
    import sphinx_rtd_theme

    html_theme = "sphinx_rtd_theme"
except ImportError:
    print("Sphinx html theme 'sphinx_rtd_theme' not found. Using 'classic' instead.")
    html_theme = "classic"

import os

# Runnig Doxygen on readthedocs
import subprocess


# can be used to use more general doxygen input files with auto
# addition of input and output dirs.
def configureDoxyfile(input_dir, output_dir):
    with open("Doxyfile", "r") as file:
        filedata = file.read()

    # filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    # filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    # with open('Doxyfile', 'w') as file:
    #     file.write(filedata)


# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

breathe_projects = {}

if read_the_docs_build:
    subprocess.call("pip install poetry", shell=True)
    subprocess.call("pwd", shell=True)
    subprocess.call("cd .. && bash rtd_build_all.sh", shell=True)
    import ibs

    output_dir = ""
    # configureDoxyfile(input_dir, output_dir)
    subprocess.call("doxygen", shell=True)
    breathe_projects["ibs"] = output_dir + "xml"
