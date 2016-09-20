# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:23:34 2016

@author: dflemin3

This script produces star-planet contrast ratios as a function of phase for
reflected light (when phase matters) and auroral features.  Note that for
contrast ratios, the real stellar spectrum (Meadows et al 2016) is used
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


# Set to True if you want to see the plots
show_plots = True

# Set to True if you want to save the plots
save_plots = False

############################################
#
# Init system
#
############################################

# Use defaults (see contrast_utils.py for default params)
# Defaults are best fit observed parameters
proxcen = cu.System()
tel = cu.Telescope() # Fiducial telescope
cor = cu.Telescope(D=30, R=100) # Future class ground-based telescope

# Echo defaults
print(proxcen)
print(tel)
print(cor)

############################################
#
# Plot Planet-star reflected light contrast ratio
# Note that this assumes coplanar orbit
#
############################################

if False:

    plt.clf()

    fig, ax = plt.subplots()

    al = np.linspace(0,180,380)

    ax.plot(al,proxcen.contrast_ratio(al),lw=3,label="Reflected Light")

    # Format
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$ [$^{\circ}$]")
    ax.set_ylabel(r"F$_p$/F$_*$")
    ax.legend(loc="lower left")
    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("planet_star.pdf")

############################################
#
# Plot Planet-star reflected light contrast ratio
# as a function of phase and albedo
# Note that this assumes coplanar orbit
#
############################################

num = 100

# true anomaly grid
alpha = np.linspace(0.,160.0,num)

# Albedo grid
A = np.linspace(0.05,0.5,num)

ALPHA, A = np.meshgrid(alpha, A) # grid of point
Z = proxcen.contrast_ratio(ALPHA) # evaluation of the function on the grid

if False:

    plt.clf()

    fig, ax = plt.subplots()

    extent = [alpha.min(),alpha.max(),A.min(),A.max()]

    im = ax.imshow(Z,cmap="magma",extent=extent,aspect="auto",
                  origin="lower",interpolation="nearest",
                  norm=mpl.colors.LogNorm()) # drawing the function

    # Plot contours on top
    cset = ax.contour(Z, np.arange(Z.min(), Z.max(), 5.0e-8), linewidths=2,
                       cmap="jet",
                       extent=extent)
    ax.clabel(cset, inline=True, fmt='%.1e', fontsize=15)

    # Format plot
    ax.set_xlabel(r"$\alpha$ [$^{\circ}$]")
    ax.set_ylabel("Albedo")

    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel('Planet-Star Flux Ratio', rotation=270, labelpad=30)

    plt.show()

    if save_plots:
        fig.savefig("contast_alpha_albedo.pdf")

############################################
#
# Plot Reflected photon flux at O Green Line
# in 0.1 Angstrom bin
# as a function of phase
# Note that this assumes coplanar orbit
#
############################################

if False:

    plt.clf()

    fig, ax = plt.subplots()

    al = np.linspace(0,180,380)
    dl = 0.1 * 1.0e-4 # Wavelength bin in microns

    # Plot the 3 lines
    lam_0 = [cu.OGREEN_LINE, cu.OREDI_LINE, cu.OREDII_LINE] # Angstroms
    colors = ["green", "red", "red"]
    labels = ["O Green", "O Red I", "O Red II"]

    for ii in range(3):
        ax.plot(al,proxcen.ref_phot_flux(al,lam_0[ii]*1.0e-4,dl),lw=3, color=colors[ii],
                label=labels[ii])

    # Format
    ax.set_yscale("log")
    ax.set_xlabel(r"$\alpha$ [$^{\circ}$]")
    ax.set_ylabel(r"Reflected Photon Flux [photons/cm$^2$/s")
    ax.legend(loc="lower left")
    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("planet_star.pdf")


############################################
#
# Compute auroral, reflected light flux
# ratios at given lines at quadrature
#
############################################

# Central wavelength of aurora
lam_0 = [cu.OGREEN_LINE, cu.OREDI_LINE, cu.OREDII_LINE] # Angstroms
dl = 0.1 * 1.0e-4 # Wavelength bin in microns

# Scale aurora strength to Prox Cen b
aurora = np.array([cu.OGREEN,cu.OREDI,cu.OREDII])*cu.REL_AURORA_STRENGTH

# Aurora name
labels = ["O Green", "O Red I", "O Red II"]

print("Aurora Photon Flux / Reflected photon flux at quadrature:")
for ii in range(3):
    tmp = proxcen.aurora_phot_flux(aurora[ii])/proxcen.ref_phot_flux(90.,lam_0[ii]*1.0e-4,dl)
    print("%s: %.1e." % (labels[ii], tmp))

print("")

print("Aurora Photon Flux / Stellar photon flux:")
for ii in range(3):
    tmp = proxcen.aurora_phot_flux(aurora[ii])/proxcen.star.flux(lam_0[ii]*1.0e-4,dl)
    print("%s: %.1e." % (labels[ii], tmp))

print("")

print("Reflected Photon Flux / Stellar photon flux at quadrature:")
for ii in range(3):
    tmp = proxcen.ref_phot_flux(90.,lam_0[ii]*1.0e-4,dl)/proxcen.star.flux(lam_0[ii]*1.0e-4,dl)
    print("%s: %.1e." % (labels[ii], tmp))

print("")
print("Stellar Photon Flux:")
for ii in range(3):
    tmp = proxcen.star.flux(lam_0[ii]*1.0e-4,dl)
    print("%s: %.1e" % (labels[ii], tmp))

############################################
#
# Plot aurora - star contrast as a function
# of aurora strength for our fiducial telescope
# and the future ground-based 30m coronagraph 'scope
#
############################################

if show_plots:

    plt.clf()

    fig, ax = plt.subplots()

    # Make multiple axes
    ax1 = ax.twinx()
    ax2 = ax.twinx()

    # Move the last y-axis spine over to the right by 10% of the width of the axes
    ax2.spines['right'].set_position(('axes', 1.1))

    # To make the border of the right-most axis visible, we need to turn the frame
    # on. This hides the other plots, however, so we need to turn its fill off.
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)

    fig.subplots_adjust(right=0.75)

    # Define auroral values
    lam_0 = cu.OGREEN_LINE*1.0e-4 # microns
    ray = np.linspace(1.0e3,1.0e8,10000)
    aurora = ray * cu.REL_AURORA_STRENGTH * cu.RUNIT
    dl = 0.1 * 1.0e-4 # Wavelength bin in microns

    # Estimate integration times
    SN = 3

    times = np.zeros(len(aurora))
    cont = np.zeros(len(aurora))
    times_cor = np.zeros(len(aurora))

    for ii in range(len(times)):
        times[ii] = tel.int_time(SN, proxcen, aurora[ii], lam_0)
        cont[ii] = tel.cphot(proxcen, aurora[ii])/tel.cstar(proxcen,lam_0)
        times_cor[ii] = cor.int_time(SN, proxcen, aurora[ii], lam_0,
                                    coronagraph=True)

    # Plot contrast ratio on left y axis
    ax.plot(ray/1000.,cont,lw=3)

    # Format
    ax.set_ylim(cont.min(),cont.max())
    ax.grid(True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"log(kR)")
    ax.set_ylabel(r"F$_p$/F$_* \vert _{\lambda=5577 \AA}$")

    # Plot Integration times on right y axes

    # Plot integration time as is
    ax1.plot(ray/1000.,times,lw=3, color = "black")

    # Format
    ax1.set_ylim(times.min(),times.max())
    ax1.set_yscale("log")
    ax1.invert_yaxis()
    #ax1.set_ylabel(r"Integration Time [hr]", rotation=270, labelpad=20)

    # Plot integration time with coronagraph
    ax2.plot(ray/1000.,times_cor, lw=3, color="black")

    # Format
    # Color axes
    color = "blue"
    ax2.tick_params(axis="y", colors=color)
    ax2.spines["right"].set_color(color)
    ax2.set_ylim(times_cor.min(),times_cor.max())
    ax2.set_yscale("log")
    ax2.invert_yaxis()
    ax2.set_ylabel(r"Integration Time [hr]", rotation=270, labelpad=20,
                   color="black")


    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("../plots/planet_star.pdf")