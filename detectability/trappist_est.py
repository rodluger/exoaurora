# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 18:17:41 2016

@author: dflemin3

Compute the relative blackbody spectral intensity for Proxima Centauri and
TRAPPIST.
"""

from __future__ import print_function
import contrast_utils as cu

# At the OI line
lam_0 = cu.OGREEN_LINE*1.0e-8

# TRAPPIST, Prox Cen temperatures
T_sun = 5800.
T_trap = 2550. # from Gillon et al 2016
T_pc = 3050. # From Anglada-Escude et al 2016

lam_0 = 2000.*1.0e-8

frac = cu.planck(T_pc,lam_0)/cu.planck(T_trap,lam_0)

print(frac)

# IWA calculation

# Separation in radians
theta = ((cu.AUCM*0.01)/(cu.DIST*10.))

lam_0 = cu.OGREEN_LINE*1.0e-4
print(theta * (10. * 100.) / (lam_0 * 1.0e-4))
