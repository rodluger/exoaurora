# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:06:49 2016

@author: dflemin3 [David P. Fleming, University of Washington]

@email: dflemin3 (at) uw (dot) edu

Produces figure that compares/contrasts blackbody
flux from Prox Cen to actual (Meadows et al 2016)

Note: Flux spectrum is what an observer on Earth would see
from the Proxima Centauri system in units of W/m^2/micron
"""

# Imports
from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import contrast_utils as cu

#Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

############################################
#
# Load in Prox Cen Spectra (from Meadows et al 2016)
#
############################################

# Default star object loads in Meadows spectrum
proxcen = cu.Star()

############################################
#
# Create blackbody flux spectrum
#
############################################

# Conv goes from cgs -> W/m^2/micron and scale from 1 AU -> distance
# Also do intensity -> flux where F = PI * B * (R/d)^2
conv = np.pi*1.0e-7*100*100*((cu.RSTAR/cu.DIST)**2)*1.0e-4

bb = cu.planck(proxcen.Teff,proxcen.wave*1.0e-4)*conv

############################################
#
# Make figure, overplot lines of interest
# Note: goes from microns -> Angstroms
#
############################################

fig, ax = plt.subplots()

# Ignore that huge initial jump (likely unphysical)
cut = 3

# Actual Spectrum
prox_spec = proxcen.spectrum(proxcen.wave)

ax.plot(proxcen.wave[cut:]*1.0e4,prox_spec[cut:],lw=3,label="Meadows et al. 2016")

# Blackbody
ax.plot(proxcen.wave[cut:]*1.0e4, bb[cut:], lw=3, color="black", label=r"T = %.0lf K Blackbody"
        % proxcen.Teff)

# Oxygen lines of interest
ax.axvline(x=cu.OGREEN_LINE,lw=3,ls="--",color="green",label=r"OI %.0lf $\AA$"
            % cu.OGREEN_LINE)
ax.axvline(x=cu.OREDI_LINE,lw=3,ls="--",color="red",label=r"OI %.0lf, %.0lf $\AA$"
           % (cu.OREDI_LINE, cu.OREDII_LINE))
ax.axvline(x=cu.OREDII_LINE,lw=3,ls="--",color="red")

# Format (and go from microns -> Angstroms)
ax.set_xlim(0.3*1.0e4,1*1.0e4)
ax.set_ylim(0,3.0e-11)
leg = ax.legend(loc="upper left")

ax.set_xlabel(r"Wavelength [$\AA$]")
ax.set_ylabel(r"Flux [W/m$^2$/$\mu$m]")

plt.show()

# Uncomment if you want to save
#fig.savefig("prox_cen_spec.pdf")