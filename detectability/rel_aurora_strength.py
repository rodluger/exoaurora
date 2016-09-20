# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 13:43:00 2016

@author: dflemin3

This script computes the relative aurora strength between Earth and
Prox Cen b assuming Prox Cen b has an Earth-like magnetic field and
atmosphere.

Note: EV Lac has a magnetic field that is ~ 3x that of Proxima Centauri
so the final estimate should be scaled down by a factor of 3 or so for
a conservative estimate.
"""

# Imports
from __future__ import print_function, division
import numpy as np

########################################
#
# Cohen et al 2014 values
#
########################################

# For Sub-Alfvenic planet B
n_e = 433. # cm^-3
T_p = 3.42e5 # K
T_e = T_p/3. # K where T_p/T_e ~ 3 in slow solar wind
vel = [-630.,-1.,30.] # km/s
B_field = [-804.,-173.,63.] # nT, B = [B_x,B_y,B_z]
B_mag = np.sqrt(np.dot(B_field,B_field))
R_mag = 0.73*5 # Magnetopause radius in Earth Units

# Solar Quantities, same units but from Wang et al 2014
n_sol = 5.
Tp_sol = 1.0e4
Te_sol = Tp_sol/3.
vel_sol = [-400.,0.,0.]
B_sol = [0.,0.,-5.]

########################################
#
# Wang et al 2014 Functions
#
########################################

def E_aurora(n,vel,B):
    """
    Auroral energy from Wang et al 2014.

    Parameters
    ----------
    n : float
        Solar wind number density [cm^-3]
    vel : array
        Cartesian components of solar wind at planet
    B : array
        Cartesian components of magnetic field at planet

    Returns
    -------
    Esw : float
        Energy of solar wind injected into aurora -> aurora power
    """

    # Correction here: 0.12 is the fraction of Esw
    # driving particle precipitation. Wang 2014
    # From this point, we should include the specific line production
    # efficiencies for the brightest ones. I am working on this, and will
    # correct the code here as well as write it up in the next day or so

    A_aur = 0.36 # FIXME Fraction of stellar wind energy fed into aurora from Matt Tilley

    # Compute intermediate quantities (see Wang et al. 2014)
    theta = np.arctan2(B[1],B[2]) # IMF Clock angle
    Bt = np.sqrt(B[1]**2 + B[2]**2) # Transverse magnetic field
    v = np.sqrt(np.dot(vel,vel)) # Solar wind velocity

    # Convert theta range to [0,2pi)
    while theta < 0:
        theta += 2.0*np.pi
    while theta > 2.0*np.pi:
        theta -= 2.0*np.pi

    # For simplicity, assume same clock angle
    theta = 0.0

    # Returns in W
    return A_aur*3.78e7*np.power(n,0.24)*np.power(v,1.47)*np.power(Bt,0.86)* \
    (np.power(np.sin(theta/2.),2.7)+0.25)
# End function

########################################
#
# Compute!
#
########################################

# Estimate auroral power of Prox b to Earth using EV Lac conditions
Proxb = E_aurora(n_e,vel,B_field)/E_aurora(n_sol,vel_sol,B_sol)
print("Assuming Proxima Centauri is EV Lac-like:")
print("Aurora on Proxb are about %.2lf times stronger than on Earth" % Proxb)

print("Assuming B_Proxima_Centauri ~ (1/3)B_EV_Lac:")
print("Aurora on Proxb are about %.2lf times stronger than on Earth" % (Proxb/3.))
