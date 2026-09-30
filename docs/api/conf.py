# Sphinx Configuration for PhyloDater API Documentation

project = 'PhyloDater'
copyright = '2025, Zichao Zeng'
author = 'Zichao Zeng'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'numpydoc',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

napoleon_google_docstring = True
napoleon_numpy_docstring = True
numpydoc_show_class_members = False
