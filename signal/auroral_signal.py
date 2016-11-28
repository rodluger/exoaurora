# -*- coding: utf-8 -*-
"""
this script computes the expected power coupling between the
stellar wind and an earth-like magnetosphere, as predicted by
wang et al 2014.

@author: mtilley [matt a. tilley, university of washington]
@email: mtilley (at) uw (dot) edu

"""

# imports
from __future__ import print_function, division
import numpy as np
import scipy as sp
import scipy.constants as spcon
from scipy import integrate
'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
relevant parameters
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
pi = spcon.pi          # pi
mu_0 = spcon.mu_0      # vacuum permeability
k_b = spcon.k          # boltzmann's constant
m_earth = 8.01e15   # earth's dipole moment [T m^3]
m_neptune = 2.2e17  # neptune's dipole moment [T m^3]
r_earth = 6.371e6   # earth's radius [m]

'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
power_calc
    calculates the power delivered to the magnetosphere
    given stellar wind conditions

inputs:
      n_sw  (float) - number density of stellar wind [cm^-3]
      v_sw  (float) - velocity of the stellar wind [km s^-1]
      b_t   (float) - transverse imf in stellar wind [nt]
      theta (float) - imf clock angle [rad]
      m_p   (float) - magnetic dipole moment [T m^3]

output:
      (float) power delivered to auroral regions [w]
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
def power_calc( n_sw=4., v_sw=400., b_t=4., theta=pi/2., m_p=m_earth ):

    aurora_frac = 0.12  # fraction of power delivered to auroral regions
    k_c = 3.78e7        # coupling constant for earth from wang et al 2014
    scale_m = (m_p/m_earth)**0.666667       # scale planetary dipole

    const = aurora_frac*scale_m*k_c
    wind = n_sw**0.24*v_sw**1.47*b_t**0.86
    imf_comp = np.sin(theta/2.)**2.7+0.25

    return const*wind*imf_comp

# end power_calc


'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
mpause_dist
    calculates magnetopause distance given
    stellar wind parameters

inputs:
      n_sw  (float) - number density of stellar wind [cm^-3]
      v_sw  (float) - velocity of the stellar wind [km s^-1]
      t_p   (float) - proton temperature [10^5 K]
      b_imf (float) - magnitude of imf [nT]
      m_p   (float) - magnetic dipole moment [T m^3]

output:
      (float) magnetopause sub-stellar distance [m]
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
def mpause_dist( n_sw=4., v_sw=400., m_p=m_earth ):

    # stellar wind ram pressure
    ram = 1.67e-27*(n_sw*1.e6)*(v_sw*1.e3)**2.

    # magnetopause distance from pressure balance
    dist = ( 4*m_p**2./(2*mu_0*(ram) ) )**0.166667

    return dist

# end mpause_dist

'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
auroral_oval
    calculates auroral oval coverage area given
    magnetopause distance and planetary radius

inputs:
      mpause (float) - magnetopause distance [m]
      r_p    (float) - planetary radius [m]

output:
      (float) area of auroral oval [m^2]
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
def auroral_oval( mpause=10.5*r_earth, r_p=r_earth ):

    # get colatitude of magnetopause foot point
    colat = np.arcsin( 1./np.sqrt( mpause/r_p ) )
    width = pi*5./180.    # assume 5 deg auroral oval

    # integrate over spherical coords to obtain
    # steradians of coverage
    sinx = lambda x: np.sin(x)
    frac_sr,err = integrate.quad(sinx,colat,colat+width)

    area = 2*pi*r_p**2.*frac_sr

    return area

# end auroral_oval
