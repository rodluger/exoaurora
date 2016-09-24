#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler.py
---------

Keplerian orbit solver, for RV calcs.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np

G = 6.67428e-11
DAYSEC = 86400.0
MSUN = 1.988416e30
MEARTH = 5.972186e24

def solve_kepler(manom, e):
  '''
  solve kepler's equation for eccentric anomaly, using mean anomaly in radians
  and eccentricity 
  
  '''
  
  # make initial guess of ecc anom = mean anom
  eanom = manom + np.sign(np.sin(manom))*0.85*e
  di_3 = 1.0
  
  while di_3 > 1e-15: #tolerance 
    fi = eanom - e*np.sin(eanom) - manom
    fi_1 = 1.0 - e*np.cos(eanom)
    fi_2 = e*np.sin(eanom)
    fi_3 = e*np.cos(eanom)
    di_1 = -fi / fi_1
    di_2 = -fi / (fi_1 + 0.5*di_1*fi_2)
    di_3 = -fi / (fi_1 + 0.5*di_2*fi_2 + 1./6.*di_2**2.0*fi_3)
    eanom = eanom + di_3
  return eanom

def ecc2true_anom(eanom, e):
  '''
  calculate true anomaly from eccentric anomaly (radians) and eccentricity
  will not work if e >= 1 !!!!!
  
  '''
  
  tanf2 = np.sqrt((1+e)/(1-e))*np.tan(eanom/2.0)
  
  f = 2.0*np.arctan(tanf2)
  
  if f < 0:
    f += 2*np.pi
    
  return f

def meanl2truea(e, argp, meanl):
  '''
  input argp and meanl in degrees
  
  '''
  
  meana = (meanl - argp) % 360
  ecca = solve_kepler(meana*np.pi/180.,e)  #ecca will be in radians
  f = ecc2true_anom(ecca, e)*180./np.pi  #get true anom and convert to degrees
  
  #return true anomaly in degrees
  return f

def RadialVelocity(jd, planet): 
  '''
  Velocity in m/s of planet RELATIVE to star, along line of sight
  
  '''
  
  e = planet.eccentricity
  P = planet.period
  t0 = planet.t0
  l0 = planet.mean_longitude
  m1 = planet.stellar_mass * MSUN
  m2 = planet.mass * MEARTH
  inc = planet.inclination
  argp = planet.arg_periastron
  meanl = (l0 + 360. * (jd - t0) / P) % 360.
  truea = meanl2truea(e,argp,meanl)
  a = ((P*DAYSEC)**2 * G *(m1+m2)/(4*np.pi**2))**(1./3)
  return -np.sqrt(G*(m1+m2)/(a*(1-e**2)))*np.sin(inc*np.pi/180.)\
         *(np.cos((argp+truea)*np.pi/180.)+e*np.cos(argp*np.pi/180.))