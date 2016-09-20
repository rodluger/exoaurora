# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:19:19 2016

@author: dflemin3

This script computes integration time required to achieve a
given signal-to-noise in the Oxygen 5577 Angstroms auroral
line band for various telescope configurations.  All calculations
assume photon-limited noise, an ideal coronagraph (when used), and
perfect subtraction of Earth's auroral features.
"""

# Imports
from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import contrast_utils as cu
import pandas as pd

#Parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 25.0
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)


# Set to True if you want to save the plots
save_plots = False

# If you want to see the plots
show_plots = False

#################################################
#
# Initialize system, telescopes
#
#################################################

# For Oxygen green line (lambda = 5577 Angstroms)
lam_0 = cu.OGREEN_LINE*1.0e-4 # -> microns

# Default Proxima Centauri system
proxcen = cu.System()

# Ideal HARPS telescope
harps = cu.Telescope(D=3.6)

# 8.2m SPHERE + ESPRESSO R = 120,000 telescope
vlt = cu.Telescope(D=8.2,R=120000)

# 30m TMT-like with/without coronagraph
tmt = cu.Telescope(D=30.)

# Define a desired signal-to-noise
SN = 3.0

#################################################
#
# Compute coronagraph inner working angle
#
#################################################

print("VLT IWA: %.1lf" % vlt.calc_IWA(proxcen, lam_0))
print("TMT IWA: %.1lf" % tmt.calc_IWA(proxcen, lam_0))


#################################################
#
# Compute auroral strengths, integration times
#
#################################################

# Define auroral values

# Create array of total auroral power in Watts
tot_aurora = np.array([1.0e9,1.0e10,1.0e11,1.0e12,1.0e13,1.0e14,1.0e15])

# Scale to oxygen green aurora from branching ratio
aurora = tot_aurora * cu.EPS_EARTH_OGREEN

# Init arrays
FpFs = np.zeros(len(aurora))
HARPS_Time = np.zeros(len(aurora))
VLT_Time = np.zeros(len(aurora))
VLT_Time_Cor = np.zeros(len(aurora))
TMT_Time = np.zeros(len(aurora))
TMT_Time_Cor = np.zeros(len(aurora))

# Loop because the objects hate vectorized calculations (?? but why)
for ii in range(len(aurora)):

    # Compute reflected photon count rate to add to star-planet contrast
    cs = harps.cstar(proxcen, lam_0) # Star photon counts
    cref = harps.cref(cs, proxcen, alpha=90.)

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

    # Auroral Planet - Star contrast in line's spectral element
    FpFs[ii] = (cref+harps.cphot(proxcen, aurora[ii]))/harps.cstar(proxcen,lam_0)


# Plop into dataframe
columns = ["Auroral Power [W]",
           "Planet-Star Contrast",
           "HARPS",
           "VLT",
           "VLT + Coronagraph",
           "TMT",
           "TMT + Coronagraph"]
df = pd.DataFrame(np.array([tot_aurora,FpFs,HARPS_Time,VLT_Time,VLT_Time_Cor,
                            TMT_Time,TMT_Time_Cor]).T,columns=columns)

# Function to format floats for dataframe
def fmt(flo):
    return "%.1e" % flo

# Print dataframe into latex table!
print(df.to_latex(float_format=fmt,index_names=False))

#################################################
#
# Generate plots
#
#################################################

# Photon count rates as a function of auroral power
if show_plots:
    fig, ax = plt.subplots()

    cs = harps.cstar(proxcen, lam_0) # Star photon counts
    cp = harps.cphot(proxcen, aurora) # Planet auroral phot counts
    cref = harps.cref(cs, proxcen, alpha=90.)

    # Plot auroral count rates
    ax.plot(aurora,cp,lw=3,color="k",label="Aurora")

    # Plot reflected light count rates
    ax.axhline(y=cref,lw=3,color="blue",label="Reflected")

    # Format plot
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_ylabel("Photon Count Rates")
    ax.set_xlabel("Auroral Power [W]")

    ax.legend(loc="upper left")

    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("../plots/count_rates.pdf")

# Integration time as a function of auroral power
# !!!
# Note: This plot probably doesn't make sense and should
# likely be ignored
# !!!
if False:

    fig, ax = plt.subplots()

    # Make multiple axes
    ax1 = ax.twinx()
    ax2 = ax.twinx()

    # Move the last y-axis spine over to the right by 10% of the width of the axes
    ax2.spines['right'].set_position(('axes', 1.125))

    # To make the border of the right-most axis visible, we need to turn the frame
    # on. This hides the other plots, however, so we need to turn its fill off.
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)

    fig.subplots_adjust(right=0.75)

    # Plot contrast ratio on left y axis
    ax.plot(aurora,FpFs,lw=3)

    # Format
    ax.set_ylim(FpFs.min(),FpFs.max())
    ax.grid(True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Auroral Power [W]")
    ax.set_ylabel(r"F$_p$/F$_* \vert _{\lambda=5577 \AA}$")

    # Plot Integration times on right y axes

    # Plot integration time as is
    ax1.plot(aurora,Ten_Time,lw=3, color = "black")

    # Format
    ax1.set_ylim(Ten_Time.min(),Ten_Time.max())
    ax1.set_yscale("log")
    ax1.invert_yaxis()

    # Plot integration time with coronagraph
    ax2.plot(aurora,Ten_Time_Cor, lw=3, color="blue")

    # Format
    # Color axes
    color = "blue"
    ax2.tick_params(axis="y", colors=color)
    ax2.spines["right"].set_color(color)
    ax2.set_ylim(Ten_Time_Cor.min(),Ten_Time_Cor.max())
    ax2.set_yscale("log")
    ax2.invert_yaxis()
    ax2.set_ylabel(r"Integration Time [hr]", rotation=270, labelpad=25,
                   color="black")

    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("../plots/planet_star.pdf")