#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
paper_figures.py
----------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import search
import numpy as np
np.random.seed(206265)

# Search
search.Search()