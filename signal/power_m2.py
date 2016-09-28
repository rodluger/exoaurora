# -*- coding: utf-8 -*-
"""
this script computes the expected auroral power output for the
oi 5577 angstrom line for proxima b, given stellar wind conditions
for planet 'b' from cohen et al 2014

@author: mtilley [matt a. tilley, university of washington]
@email: mtilley (at) uw (dot) edu

"""

# imports
from __future__ import print_function, division
import numpy as np
from numpy import linalg as la
import scipy as sp
import scipy.constants as spcon
import auroral_signal as asig
'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
relevant parameters
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
# sub-alfvenic stellar wind
nsub = 433.                 # [cm^-3]
vsub = [-630.,-1,30.]       # [km s^-1]
bsub = [-804.,-173.,63.]    # [nt]
tsub = 3.42                 # [10^5 k]

# super-alfvenic stellar wind
nsup = 12895.               # [cm^-3]
vsup = [-202.,102.,22.]     # [km s^-1]
bsup = [-57.,-223.,92.]     # [nt]
tsup = 4.77                 # [10^5 k]

# cme scaling using wang et al formula
# and nominal increases in density, velocity, and imf [arb]
cme_scale = 10**0.24 * 3.**1.47 * 15**0.86

# energy per photon [j photon^-1]
e_5577 = spcon.h*spcon.c/5.577e-7

# electron fraction of auroral precip [arb]
e_frac = 0.8

# steele & mcewen conversion efficiency [photons (erg cm^-2 s^-1)^-1]
oi_eff = 1.55e9

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

'''
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
power estimations
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
'''
# power estimated at earth
print( '\n-=-=-=-=-=-=-=-=-=-=-' )
print( '\n\tearth\n' )
print( '-=-=-=-=-=-=-=-=-=-=-\n' )

# magnetopause distance
dist = asig.mpause_dist()
print( ' estimated sub-stellar magnetopause distance: %.2e' % dist + ' m\n' )

# auroral oval
oval = asig.auroral_oval()*1.e4
print( ' estimated auroral oval area: %.2e' % oval + ' cm^2\n' )

print( ' estimated auroral energetic particle power delivered to auroral regions:\n' )

# quiet magnetosphere
power_out_q = asig.power_calc()
print( ' quiet\t\t-\t%.2e' % power_out_q + ' w' )

# stormy magnetosphere
power_out_s = asig.power_calc( theta=spcon.pi )
print( ' substorm\t-\t%.2e' % power_out_s + ' w' )

# cme wind conditions
power_out_c = asig.power_calc()*cme_scale
print( ' cme\t\t-\t%.2e' % power_out_c + ' w' )

# cme wind conditions + stormy magnetosphere
power_out_cs = asig.power_calc( theta=spcon.pi )*cme_scale
print( ' cme+substorm\t-\t%.2e' % power_out_cs + ' w' )

print( '\n estimated auroral power for the oi 5577 a line:\n' )

# 5577 power out for quiet msphere, both hemispheres
out_5577_q = oi_eff*power_out_q*1.e7*e_frac*2*e_5577
print( ' quiet 5577\t-\t%.2e' % out_5577_q + ' w' )

# 5577 power out for stormy msphere, both hemispheres
out_5577_s = oi_eff*power_out_s*1.e7*e_frac*2*e_5577
print( ' substorm 5577\t-\t%.2e' % out_5577_s + ' w' )

# 5577 power out for cme winds, both hemispheres
out_5577_c = oi_eff*power_out_c*1.e7*e_frac*2*e_5577
print( ' cme 5577\t-\t%.2e' % out_5577_c + ' w' )

# 5577 power out for cme winds + stormy msphere, both hemispheres
out_5577_cs = oi_eff*power_out_cs*1.e7*e_frac*2*e_5577
print( ' cme+ss 5577\t-\t%.2e' % out_5577_cs + ' w' )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# power estimated at
print( '\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-' )
print( '\n proxima b - sub-alfvenic,earth-like dipole\n' )
print( '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n' )

# magnetopause distance
dist = asig.mpause_dist( nsub, la.norm(vsub), tsub, la.norm(bsub), asig.m_earth )
print( ' estimated sub-stellar magnetopause distance: %.2e' % dist + ' m\n' )

# auroral oval
oval = asig.auroral_oval( dist )*1.e4
print( ' estimated auroral oval area: %.2e' % oval + ' cm^2\n' )

# imf clock angle and transverse imf
imf_clock = np.arctan2( abs(bsub[1]), bsub[2] )
b_t = np.sqrt( bsub[1]**2. + bsub[2]**2. )

print( ' estimated auroral energetic particle power delivered to auroral regions:\n' )

# quiet magnetosphere
power_out_q = asig.power_calc( nsub, la.norm(vsub), b_t, imf_clock, asig.m_earth )
print( ' quiet\t\t-\t%.2e' % power_out_q + ' w' )

# stormy magnetosphere
power_out_s = asig.power_calc( nsub, la.norm(vsub), b_t, spcon.pi, asig.m_earth )
print( ' substorm\t-\t%.2e' % power_out_s + ' w' )

# cme wind conditions
power_out_c = asig.power_calc( nsub, la.norm(vsub), b_t, imf_clock, asig.m_earth )*cme_scale
print( ' cme\t\t-\t%.2e' % power_out_c + ' w' )

# cme wind conditions + stormy magnetosphere
power_out_cs = asig.power_calc( nsub, la.norm(vsub), b_t, spcon.pi, asig.m_earth )*cme_scale
print( ' cme+substorm\t-\t%.2e' % power_out_cs + ' w' )

print( '\n estimated auroral power for the oi 5577 a line:\n' )

# 5577 power out for quiet msphere, both hemispheres
out_5577_q = oi_eff*power_out_q*1.e7*e_frac*2*e_5577
print( ' quiet 5577\t-\t%.2e' % out_5577_q + ' w' )

# 5577 power out for stormy msphere, both hemispheres
out_5577_s = oi_eff*power_out_s*1.e7*e_frac*2*e_5577
print( ' substorm 5577\t-\t%.2e' % out_5577_s + ' w' )

# 5577 power out for cme winds, both hemispheres
out_5577_c = oi_eff*power_out_c*1.e7*e_frac*2*e_5577
print( ' cme 5577\t-\t%.2e' % out_5577_c + ' w' )

# 5577 power out for cme winds + stormy msphere, both hemispheres
out_5577_cs = oi_eff*power_out_cs*1.e7*e_frac*2*e_5577
print( ' cme+ss 5577\t-\t%.2e' % out_5577_cs + ' w' )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# power estimated at
print( '\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-' )
print( '\n proxima b - super-alfvenic,earth-like dipole\n' )
print( '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n' )

# magnetopause distance
dist = asig.mpause_dist( nsup, la.norm(vsup), tsup, la.norm(bsup), asig.m_earth )
print( ' estimated sub-stellar magnetopause distance: %.2e' % dist + ' m\n' )

# auroral oval
oval = asig.auroral_oval( dist )*1.e4
print( ' estimated auroral oval area: %.2e' % oval + ' cm^2\n' )

# imf clock angle and transverse imf
imf_clock = np.arctan2( abs(bsup[1]), bsup[2] )
b_t = np.sqrt( bsup[1]**2. + bsup[2]**2. )

print( ' estimated auroral energetic particle power delivered to auroral regions:\n' )

# quiet magnetosphere
power_out_q = asig.power_calc( nsup, la.norm(vsup), b_t, imf_clock, asig.m_earth )
print( ' quiet\t\t-\t%.2e' % power_out_q + ' w' )

# stormy magnetosphere
power_out_s = asig.power_calc( nsup, la.norm(vsup), b_t, spcon.pi, asig.m_earth )
print( ' substorm\t-\t%.2e' % power_out_s + ' w' )

# cme wind conditions
power_out_c = asig.power_calc( nsup, la.norm(vsup), b_t, imf_clock, asig.m_earth )*cme_scale
print( ' cme\t\t-\t%.2e' % power_out_c + ' w' )

# cme wind conditions + stormy magnetosphere
power_out_cs = asig.power_calc( nsup, la.norm(vsup), b_t, spcon.pi, asig.m_earth )*cme_scale
print( ' cme+substorm\t-\t%.2e' % power_out_cs + ' w' )

print( '\n estimated auroral power for the oi 5577 a line:\n' )

# 5577 power out for quiet msphere, both hemispheres
out_5577_q = oi_eff*power_out_q*1.e7*e_frac*2*e_5577
print( ' quiet 5577\t-\t%.2e' % out_5577_q + ' w' )

# 5577 power out for stormy msphere, both hemispheres
out_5577_s = oi_eff*power_out_s*1.e7*e_frac*2*e_5577
print( ' substorm 5577\t-\t%.2e' % out_5577_s + ' w' )

# 5577 power out for cme winds, both hemispheres
out_5577_c = oi_eff*power_out_c*1.e7*e_frac*2*e_5577
print( ' cme 5577\t-\t%.2e' % out_5577_c + ' w' )

# 5577 power out for cme winds + stormy msphere, both hemispheres
out_5577_cs = oi_eff*power_out_cs*1.e7*e_frac*2*e_5577
print( ' cme+ss 5577\t-\t%.2e' % out_5577_cs + ' w' )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# power estimated at
print( '\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-' )
print( '\n proxima b - sub-alfvenic,neptune-like dipole\n' )
print( '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n' )

# magnetopause distance
dist = asig.mpause_dist( nsub, la.norm(vsub), tsub, la.norm(bsub), asig.m_neptune )
print( ' estimated sub-stellar magnetopause distance: %.2e' % dist + ' m\n' )

# auroral oval
oval = asig.auroral_oval( dist )*1.e4
print( ' estimated auroral oval area: %.2e' % oval + ' cm^2\n' )

# imf clock angle and transverse imf
imf_clock = np.arctan2( abs(bsub[1]), bsub[2] )
b_t = np.sqrt( bsub[1]**2. + bsub[2]**2. )

print( ' estimated auroral energetic particle power delivered to auroral regions:\n' )

# quiet magnetosphere
power_out_q = asig.power_calc( nsub, la.norm(vsub), b_t, imf_clock, asig.m_neptune )
print( ' quiet\t\t-\t%.2e' % power_out_q + ' w' )

# stormy magnetosphere
power_out_s = asig.power_calc( nsub, la.norm(vsub), b_t, spcon.pi, asig.m_neptune )
print( ' substorm\t-\t%.2e' % power_out_s + ' w' )

# cme wind conditions
power_out_c = asig.power_calc( nsub, la.norm(vsub), b_t, imf_clock, asig.m_neptune )*cme_scale
print( ' cme\t\t-\t%.2e' % power_out_c + ' w' )

# cme wind conditions + stormy magnetosphere
power_out_cs = asig.power_calc( nsub, la.norm(vsub), b_t, spcon.pi, asig.m_neptune )*cme_scale
print( ' cme+substorm\t-\t%.2e' % power_out_cs + ' w' )

print( '\n estimated auroral power for the oi 5577 a line:\n' )

# 5577 power out for quiet msphere, both hemispheres
out_5577_q = oi_eff*power_out_q*1.e7*e_frac*2*e_5577
print( ' quiet 5577\t-\t%.2e' % out_5577_q + ' w' )

# 5577 power out for stormy msphere, both hemispheres
out_5577_s = oi_eff*power_out_s*1.e7*e_frac*2*e_5577
print( ' substorm 5577\t-\t%.2e' % out_5577_s + ' w' )

# 5577 power out for cme winds, both hemispheres
out_5577_c = oi_eff*power_out_c*1.e7*e_frac*2*e_5577
print( ' cme 5577\t-\t%.2e' % out_5577_c + ' w' )

# 5577 power out for cme winds + stormy msphere, both hemispheres
out_5577_cs = oi_eff*power_out_cs*1.e7*e_frac*2*e_5577
print( ' cme+ss 5577\t-\t%.2e' % out_5577_cs + ' w' )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# power estimated at
print( '\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-' )
print( '\n proxima b - super-alfvenic,neptune-like dipole\n' )
print( '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n' )

# magnetopause distance
dist = asig.mpause_dist( nsup, la.norm(vsup), tsup, la.norm(bsup), asig.m_neptune )
print( ' estimated sub-stellar magnetopause distance: %.2e' % dist + ' m\n' )

# auroral oval
oval = asig.auroral_oval( dist )*1.e4
print( ' estimated auroral oval area: %.2e' % oval + ' cm^2\n' )

# imf clock angle and transverse imf
imf_clock = np.arctan2( abs(bsup[1]), bsup[2] )
b_t = np.sqrt( bsup[1]**2. + bsup[2]**2. )

print( ' estimated auroral energetic particle power delivered to auroral regions:\n' )

# quiet magnetosphere
power_out_q = asig.power_calc( nsup, la.norm(vsup), b_t, imf_clock, asig.m_neptune )
print( ' quiet\t\t-\t%.2e' % power_out_q + ' w' )

# stormy magnetosphere
power_out_s = asig.power_calc( nsup, la.norm(vsup), b_t, spcon.pi, asig.m_neptune )
print( ' substorm\t-\t%.2e' % power_out_s + ' w' )

# cme wind conditions
power_out_c = asig.power_calc( nsup, la.norm(vsup), b_t, imf_clock, asig.m_neptune )*cme_scale
print( ' cme\t\t-\t%.2e' % power_out_c + ' w' )

# cme wind conditions + stormy magnetosphere
power_out_cs = asig.power_calc( nsup, la.norm(vsup), b_t, spcon.pi, asig.m_neptune )*cme_scale
print( ' cme+substorm\t-\t%.2e' % power_out_cs + ' w' )

print( '\n estimated auroral power for the oi 5577 a line:\n' )

# 5577 power out for quiet msphere, both hemispheres
out_5577_q = oi_eff*power_out_q*1.e7*e_frac*2*e_5577
print( ' quiet 5577\t-\t%.2e' % out_5577_q + ' w' )

# 5577 power out for stormy msphere, both hemispheres
out_5577_s = oi_eff*power_out_s*1.e7*e_frac*2*e_5577
print( ' substorm 5577\t-\t%.2e' % out_5577_s + ' w' )

# 5577 power out for cme winds, both hemispheres
out_5577_c = oi_eff*power_out_c*1.e7*e_frac*2*e_5577
print( ' cme 5577\t-\t%.2e' % out_5577_c + ' w' )

# 5577 power out for cme winds + stormy msphere, both hemispheres
out_5577_cs = oi_eff*power_out_cs*1.e7*e_frac*2*e_5577
print( ' cme+ss 5577\t-\t%.2e' % out_5577_cs + ' w' )

