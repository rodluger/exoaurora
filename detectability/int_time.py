# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:19:19 2016

@author: dflemin3 [David P. Fleming, University of Washington]

@email: dflemin3 (at) uw (dot) edu

This script computes integration time required to achieve a
given signal-to-noise in the Oxygen 5577 Angstroms auroral
line band for various telescope configurations.  All calculations
assume photon-limited noise or a coronagraph with noise from Robinson et al 2016
and perfect subtraction of telluric features (excellent optical AO) and rigorous
data reduction (i.e. how theorists hope we'll be able to observe).
"""

# Imports
from __future__ import print_function, division
import numpy as np
import contrast_utils as cu
import coronagraph_inputs as ci
import pandas as pd

#################################################
#
# Initialize system, telescopes
#
#################################################

# For Oxygen green line (lambda = 5577 Angstroms)
lam_0 = cu.OGREEN_LINE * 1.0e-4 # -> microns

# Nitrogen N_2+ UV line
#lam_0 = cu.NUV_LINE * 1.0e-4 # -> microns

# Default Proxima Centauri system
proxcen = cu.System()

# HARPS telescope
harps = cu.Telescope(D=3.6)

# 8.2m SPHERE + ESPRESSO (Lovis et al 2016) with R = 120,000 spectrograph
vlt = cu.Telescope(D=8.2,R=120000)

# 30m TMT-like
tmt = cu.Telescope(D=30.)

# LUVOIR-concept space-based telescope
luvoir = cu.Telescope(D=16.)

# HabEX-concept space-based telescope
habex = cu.Telescope(D=6.5)

# Define a desired signal-to-noise
SN = 6.0

#################################################
#
# Compute coronagraph inner working angles
#
#################################################

# Print auroral line to observe
print("Auroral emission line wavelength [microns]: %.6lf" % lam_0)

print("VLT IWA: %.1lf" % vlt.calc_IWA(proxcen, lam_0))
print("TMT IWA: %.1lf" % tmt.calc_IWA(proxcen, lam_0))
print("LUVOIR IWA: %.1lf" % luvoir.calc_IWA(proxcen, lam_0))
print("HABEX IWA: %.1lf" % habex.calc_IWA(proxcen, lam_0))

# Print design contrast
print("Coronagraph design contrast : %.1e" % ci.C)

#################################################
#
# Compute auroral strengths, integration times
#
#################################################

# Define auroral values

# Create array of total auroral power in Watts in the given line
aurora = np.array([1.0e9,1.0e10,1.0e11,1.0e12,1.0e13,1.0e14,1.0e15])

# Init arrays to hold results
FpFs = np.zeros(len(aurora))
HARPS_Time = np.zeros(len(aurora))
VLT_Time = np.zeros(len(aurora))
VLT_Time_Cor = np.zeros(len(aurora))
TMT_Time = np.zeros(len(aurora))
TMT_Time_Cor = np.zeros(len(aurora))
LUVOIR_TIME = np.zeros(len(aurora))
HABEX_TIME = np.zeros(len(aurora))

# Loop because the objects hate vectorized calculations (?? but why)
for ii in range(len(aurora)):

    # Compute reflected photon count rate to add to star-planet contrast
    cs = luvoir.cstar(proxcen, lam_0) # Star photon counts
    cref = luvoir.cref(cs, proxcen, alpha=90.) # Take at quadrature

    # HARPS
    HARPS_Time[ii] = harps.int_time(SN, proxcen, aurora[ii], lam_0)

    # VLT
    VLT_Time[ii] = vlt.int_time(SN, proxcen, aurora[ii], lam_0)

    # VLT +  coronagraph
    VLT_Time_Cor[ii] = vlt.int_time(SN, proxcen, aurora[ii], lam_0,
                                          coronagraph=True)

    # TMT
    TMT_Time[ii] = tmt.int_time(SN, proxcen, aurora[ii], lam_0)

    # TMT +  coronagraph
    TMT_Time_Cor[ii] = tmt.int_time(SN, proxcen, aurora[ii], lam_0,
                                          coronagraph=True)

    # LUVOIR +  coronagraph
    LUVOIR_TIME[ii] = luvoir.int_time(SN, proxcen, aurora[ii], lam_0,
                                          coronagraph=True)

    # HABEX +  coronagraph
    HABEX_TIME[ii] = habex.int_time(SN, proxcen, aurora[ii], lam_0,
                                          coronagraph=True)

    # Auroral Planet - Star contrast in line's spectral element
    # Use LUVOIR because it doesn't matter as long as all count rates
    # are computed using same telescope (count rates were scaled by mirror D)
    FpFs[ii] = (cref+luvoir.cphot(proxcen, aurora[ii]))/luvoir.cstar(proxcen,lam_0)


# Plop into dataframe
columns = ["Auroral Power [W]",
           "Planet-Star Contrast",
           "HARPS",
           "VLT",
           "VLT + Coronagraph",
           "TMT",
           "TMT + Coronagraph",
           "HaBEX",
           "LUVOIR"]
df = pd.DataFrame(np.array([aurora,FpFs,HARPS_Time,VLT_Time,VLT_Time_Cor,
                            TMT_Time,TMT_Time_Cor,HABEX_TIME,LUVOIR_TIME
                            ]).T,columns=columns)

# Function to format floats for dataframe
def fmt(flo):
    return "%.0e" % flo

# Print dataframe into latex table!
print(df.to_latex(float_format=fmt,index_names=False))