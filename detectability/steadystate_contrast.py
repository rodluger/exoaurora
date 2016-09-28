"""
Created on Thu Sep 27 09:01:33 2016

@author: Jave

"""

from __future__ import print_function, division
import numpy as np
import contrast_utils as cu
import coronagraph_inputs as ci

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes

#Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 30.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

# Import from npz

nfile = np.load("TMTC_int_time_C1e-05_NOnoise2.npz")
resolver, watts, CN5 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]
nfile = np.load("TMTC_int_time_C1e-06_NOnoise2.npz")
resolver, watts, CN6 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]
nfile = np.load("TMTC_int_time_C1e-07_NOnoise2.npz")
resolver, watts, CN7 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]

nfile = np.load("TMTC_int_time_C1e-05_noise2.npz")
resolver, watts, C5 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]
nfile = np.load("TMTC_int_time_C1e-06_noise2.npz")
resolver, watts, C6 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]
nfile = np.load("TMTC_int_time_C1e-07_noise2.npz")
resolver, watts, C7 = nfile["arr_0"], nfile["arr_1"], nfile["arr_2"]

lw = 3.0
bins = 200
aurora = 0.1

watts = np.logspace(9,15,bins) / 1e12   # Auroral power [TW]
resolver = np.logspace(1,6,bins) # Resolving power

wi = np.fabs(watts - aurora).argmin()

# Create figure
fig, ax = plt.subplots()

# With noise
ax.plot(resolver, C5[wi,:], lw=lw, color="red", label=r"$C=10^{-5}$")
ax.plot(resolver, C6[wi,:], lw=lw, color="blue", label=r"$C=10^{-6}$")
ax.plot(resolver, C7[wi,:], lw=lw, color="green", label=r"$C=10^{-7}$")

# Without noise
ax.plot(resolver, CN5[wi,:], lw=lw, ls="--", color="red")
ax.plot(resolver, CN6[wi,:], lw=lw, ls="--", color="blue")
ax.plot(resolver, CN7[wi,:], lw=lw, ls="--", color="green")

ax.text(0.05, 0.05, r"Auroral Power = 0.1 TW",\
         verticalalignment='bottom', horizontalalignment='left',\
         transform=ax.transAxes,\
         color='black', fontsize=30)
# Axis Format
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(3e2,resolver.max())
ax.set_ylim(10,1e4)

# Axis lables
ax.set_ylabel(r"Integration Time [hrs]")
ax.set_xlabel(r"Resolving Power ($\lambda / \Delta \lambda$)")

# Legend
leg=ax.legend(loc=1, fontsize=20)
leg.get_frame().set_alpha(0.0)
for text in leg.get_texts():
    text.set_color("black")

plt.show()

if True:
    fig.savefig("steadystate_contrast2.pdf", bbox_inches='tight')
