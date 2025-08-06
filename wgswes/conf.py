# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------
project = 'WGS/WES Pipeline'
copyright = '2025, Akhilesh Kaushal'
author = 'Akhilesh Kaushal'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.githubpages',
    "sphinx.ext.graphviz",
    'sphinx_copybutton'  # Optional but recommended
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = '.rst'
source_encoding = 'utf-8-sig'
todo_include_todos = True
html_show_sphinx = False

# -- HTML output -------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['css/custom.css']

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "style_external_links": True,
    "display_version": True,
    "sticky_navigation": True,
    "prev_next_buttons_location": "both",
    "navigation_with_keys": True
}

html_baseurl = "https://akhileshkaushal.github.io/wgswes/"

# Optional Logo/Favicon
# html_logo = "_static/images/logo.png"
# html_favicon = "_static/images/favicon.ico"