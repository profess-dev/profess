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
#              'sphinx.ext.napoleon',
              'numpydoc',
              'myst_nb']

autosummary_generate = True

templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_title = ''
html_logo = 'profess.svg'

jupyter_execute_notebooks = 'auto'
execution_timeout = -1
execution_excludepatterns = ['HydrogenAtom.ipynb',
                             'HydrogenDimer.ipynb', 
                             'Aluminium.ipynb',
                             'Magnesium.ipynb']

# -- Options for math --------------------------------------------------------

math_number_all = True
