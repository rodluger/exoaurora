# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 18:17:41 2016

@author: dflemin3 [David P. Fleming, University of Washington]

@email: dflemin3 (at) uw (dot) edu

Compute the relative blackbody spectral intensity for Proxima Centauri and
TRAPPIST-1 and the minimum IWA required to observed the OI green line at
the TRAPPIST-1-Earth distance.
"""

from __future__ import print_function
import contrast_utils as cu

# OI line in cm
lam_0 = cu.OGREEN_LINE*1.0e-8

# TRAPPIST, Prox Cen effective temperatures
T_trap = 2550. # from Gillon et al 2016
T_pc = 3050. # From Anglada-Escude et al 2016

lam_0 = 2000.*1.0e-8

# Spectral intensity ratio
frac = cu.planck(T_pc,lam_0)/cu.planck(T_trap,lam_0)

# Output!
print(frac)

# IWA calculation

# Separation in radians
theta = ((cu.AUCM*0.01)/(cu.DIST*10.))

# Green line in microns
lam_0 = cu.OGREEN_LINE*1.0e-4

# Output!
print(theta * (10. * 100.) / (lam_0 * 1.0e-4))