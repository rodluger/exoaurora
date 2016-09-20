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

# Build up FAP, then plot the injection detection plot (8 sigma)
search.Plot(inject_contrast = 8e-3, figname = 'strong.pdf', suptitle = 'Injected')

# The actual "detection" (~5 sigma)
planet = search.ProxCenB()
planet.inclination = 55
planet.mass = 1.55
planet.period = 11.18349
planet.stellar_mass = 0.11
planet.mean_longitude = 117.26
search.Plot(planet = planet, figname = 'real.pdf', suptitle = 'Real')