# Configuration file for the Sphinx documentation builder.
#
# For a full list of configuration options, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html.

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../../build/'))


# -- Project information -----------------------------------------------------

project = 'PROFESS'
copyright = '2019-2021, William C. Witt'
author = 'William C. Witt'
release = '4-dev'


# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosectionlabel',
              'sphinx.ext.autosummary',
              'numpydoc',
              'myst_nb']

autosummary_generate = True

templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 1,
    "show_prev_next": False,
}
html_static_path = ['_static']
html_css_files = ['custom.css']

html_title = ''
html_logo = 'profess.svg'

nb_execution_mode = 'cache'
nb_execution_timeout = -1
nb_execution_excludepatterns = []

# -- Options for math --------------------------------------------------------

math_number_all = True
