# -*- coding: utf-8 -*-
"""Sphinx config file."""
from __future__ import unicode_literals

import os
import mock
import sys


mock_modules = [
    'biopython',
    'jsonpickle',
    'numpy',
    'pdbtools',
    'pdbtools.pdb_segxchain',
    'pdbtools.pdb_splitchain',
    'pdbtools.pdb_splitmodel',
    'pdbtools.pdb_tidy',
    'pyyaml',
    ]

for modulename in mock_modules:
    sys.modules[modulename] = mock.Mock()

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
    'sphinx.ext.autosectionlabel',
    'myst_parser',
    ]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
    }

master_doc = 'index'
project = 'haddock3'
year = '2022'
author = 'BonvinLab'
copyright = '{0}, {1}'.format(year, author)
version = release = '3.0.0'

todo_include_todos = True
pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/haddocking/haddock3/issues/%s', '#'),
    'pr': ('https://github.com/haddocking/haddock3/pull/%s', 'PR #'),
    }
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

linkcheck_ignore = [r'https://codecov.io/*']

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

# html_logo = 'img/taurenmd_logo_black.png'
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
    '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
    }
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
