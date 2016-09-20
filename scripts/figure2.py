#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
figure2.py
----------

Script to reproduce Figure 2 in Luger et al. (2016)

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import search
import numpy as np
np.random.seed(123)

# Build up FAP
#search.Compute(plot = False, fap_iter = 5000)

# Injection detection plot (8 sigma)
#search.Plot(frame = 'planet', inject_contrast = 8e-3, crop_outliers = True, fap_iter = 0,
#             figname = 'strong.pdf', suptitle = 'Injected')

# Search
search.Search(orb_iter = 100, inclination = [55.])
quit()

planet = search.ProxCenB()
planet.inclination = 55
planet.mass = 1.55
planet.period = 11.18349
planet.stellar_mass = 0.11
planet.mean_longitude = 117.26
search.Plot(frame = 'planet', inject_contrast = 0, crop_outliers = True, fap_iter = 0,
            figname = 'fourth_quarter.pdf', planet = planet, load_fap = False)