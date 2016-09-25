# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:01:33 2016

@author: dflemin3

Visualize the line's equivalent width as a function of contrast ratio.
"""

from __future__ import print_function, division
import numpy as np
import contrast_utils as cu

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes

#Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

show_spec = True

show_eqw = False

save_plots = False

show_contour = False

# Init system object with Proxima Centauri defaults
proxcen = cu.System()

# Specify an auroral emission line
lam_0 = cu.OGREEN_LINE*1.0e-4 # microns

# Grid over auroral powers
watts = np.linspace(1.0e10,1.0e15,20)
equivalent_width = np.zeros_like(watts)

for ii in range(len(watts)):

    # FWHM from Eric
    FWHM = 0.0186*1.0e-4
    dl = 50*FWHM

    wave_hires = cu.make_wave_array(lam_0,dl,FWHM)

    # Compute auroral flux line profile and hi-res wavelength array
    aurora = cu.create_auroral_line(lam_0, watts[ii], 1.0, proxcen, wave_hires)

    # Using star's spectrum, get stellar flux over the hi-res array
    # Note: stellar spectrum is R~10,000, so flux will be interpolated
    stellar_flux = proxcen.star.spectrum(wave_hires)

    # Now compute the reflected flux at quadrature assuming a grey planet with
    # geometric albedo = 0.3
    ref_flux = stellar_flux * proxcen.contrast_ratio(90.0)

    # Compute total planetary flux (reflected + auroral)
    planet_flux = aurora + ref_flux

    # Compute continuum
    mask = (wave_hires >= lam_0 - FWHM) & (wave_hires <= lam_0 + FWHM)
    continuum = np.median(ref_flux[mask])

    # Compute equivalent width in angstroms
    equivalent_width[ii] = 1.0e4*cu.equivalent_width(planet_flux[mask],
                                                     wave_hires[mask],
                                                     continuum)

# Plot Planet flux
if show_eqw:

    # Plot equivalent widths as a function of auroral power
    fig, ax = plt.subplots()

    ax.plot(watts,equivalent_width, lw=3)

    # Format plot
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_ylabel(r"Equivalent Width [\AA]")
    ax.set_xlabel(r"Auroral Power [W]")

    plt.show()

    if save_plots:
        fig.savefig("eq_vs_power.pdf")

if show_spec:

    # Plot spectrum with injected aurora
    fig, ax = plt.subplots()

    # FWHM from Eric
    vFWHM = 1.0
    FWHM = 0.0186*1.0e-4
    dl = 70000*FWHM

    # Auroral power
    apow = 6.5e10

    wave_hires = cu.make_wave_array(lam_0,dl,FWHM)

     #Compute auroral flux line profile and hi-res wavelength array
    aurora = cu.create_auroral_line(lam_0, apow, vFWHM, proxcen,
                                    wave_hires)

    # Using star's spectrum, get stellar flux over the hi-res array
    # Note: stellar spectrum is R~10,000, so flux will be interpolated
    stellar_flux = proxcen.star.spectrum(wave_hires)

    # Now compute the reflected flux at quadrature assuming a grey planet with
    # geometric albedo = 0.3
    ref_flux = stellar_flux * proxcen.contrast_ratio(90.0)

    # Compute total planetary flux (reflected + auroral)
    planet_flux = aurora + ref_flux

    ax.plot((wave_hires-lam_0)*1.0e4,planet_flux,color="black",lw=3,
            label="Net Planet Flux")

    mask = (wave_hires >= lam_0 - FWHM) & (wave_hires <= lam_0 + FWHM)

    # Compute continuum
    continuum = np.median(ref_flux[mask])

    # Plot continuum
    #ax.axhline(y=continuum, color = "red", ls = "--", lw=3, label="Continuum")

    # Compute equivalent width
    ew = cu.equivalent_width(planet_flux[mask],wave_hires[mask],continuum)

    # Plot equivalent width box
    ax.fill_between((wave_hires-lam_0)*1.0e4, continuum, 2.0*continuum,
                    where=((wave_hires < lam_0 + ew/2) & (wave_hires > lam_0 - ew/2.)),
                    facecolor='green', alpha = 0.2)

    # Format
    ax.set_ylabel(r"Flux [W/m$^2$/$\mu$m]")
    ax.set_xlabel(r"Wavelength [\AA]$- 5577$\AA")
    ax.set_xlim(((wave_hires-lam_0)*1.0e4).min(),
                ((wave_hires-lam_0)*1.0e4).max())
    ax.set_yscale("log")

    # Inset axis
    xmin = lam_0 - ew/2
    xmax = lam_0 + ew/2
    mask2 = (wave_hires >= xmin) & (wave_hires <= xmax)
    ax2 = plt.axes([.61, .6, .32, .3])
    ax2.plot((wave_hires[mask2])*1.0e4,planet_flux[mask2],color="black",lw=3)
    ax2.fill_between((wave_hires[mask2])*1.0e4, continuum, 2.0*continuum,
                    where=((wave_hires[mask2] < lam_0 + ew/2) & (wave_hires[mask2] > lam_0 - ew/2.)),
                    facecolor='green', alpha = 0.2)
    ax2.set_xlim((xmin*1.0e4, xmax*1.0e4))
    xticks = [xmin*1.0e4, lam_0*1.0e4, xmax*1.0e4]
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(["%.2f" % x for x in xticks])
    plt.setp(ax2.get_xticklabels(), fontsize=12, rotation=0)
    plt.setp(ax2.get_yticklabels(), fontsize=12, rotation=0)
    ax2.set_yscale("log")

    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("ew_test_case.pdf")

if show_contour:

    # Contrast, Resolving power, auroral power grid
    bins = 100

    watts = np.logspace(9,15,bins)
    resolver = np.logspace(1,6,bins)

    contrast = np.zeros((bins,bins))
    ew = np.zeros_like(watts)

    R_limit = lam_0/FWHM

    for ii in range(len(watts)):
        for jj in range(len(resolver)):

            # Compute auroral flux line profile and hi-res wavelength array
            # Compute delta_lam
            dl = lam_0/resolver[jj]
            FWHM = 0.0186*1.0e-4
            wave_hires = cu.make_wave_array(lam_0,dl,FWHM)

            aurora = cu.create_auroral_line(lam_0, watts[ii], 1.0, proxcen,
                                            wave_hires)

            # Using star's spectrum, get stellar flux over the hi-res array
            # Note: stellar spectrum is R~10,000, so flux will be interpolated
            stellar_flux = proxcen.star.spectrum(wave_hires)

            # Now compute the reflected flux at quadrature assuming a grey planet with
            # geometric albedo = 0.3
            ref_flux = stellar_flux * proxcen.contrast_ratio(90.0)

            # Compute total planetary flux (reflected + auroral)
            planet_flux = aurora + ref_flux

            # Integrate planet and star flux across spectral element;
            # divide for contrast ratio
            num = np.trapz(planet_flux,wave_hires)
            denom = np.trapz(stellar_flux,wave_hires)
            contrast[ii,jj] = num/denom
            # End inner loop

        # Compute equivalent width for each power step

        # Compute continuum
        continuum = np.median(ref_flux)
        ew[ii] = cu.equivalent_width(planet_flux,wave_hires,continuum)

    # Compute equivlent resolving power (?)
    ew_r = lam_0/ew

    fig, ax = plt.subplots()

    vmin = contrast.min()
    vmax = contrast.max()

    # Color contour
    cax = ax.pcolor(watts, resolver, contrast.T, cmap="viridis",
                    norm=colors.LogNorm(vmin=vmin,vmax=vmax))

    # Contour Lines
    contour_levels = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
    # Automate "manual locations" for contour labels
    ia = (np.abs(watts-6.5e14)).argmin() # A. Power index
    ic = [np.abs(contrast[ia, :] - clvl).argmin() for clvl in contour_levels]
    manual_locations = [(watts[ia], ic[i]) for i in range(len(contour_levels))]
    # Log levels
    fmt = ticker.LogFormatterMathtext()
    cln = ax.contour(watts, resolver, contrast.T, contour_levels,
                     colors=["black","black", "black", "black", "black"])
    plt.clabel(cln, inline=1, fontsize=15, fmt=fmt, inline_spacing=20,
               use_clabeltext=True)#, manual=manual_locations)
    # Thicken the contour lines
    plt.setp(cln.collections, linewidth=2)

    # Colorbar
    cbar = fig.colorbar(cax)
    cbar.set_label(r"Planet-Star Contrast", rotation=270, labelpad=25)

    # Plot equivalent resolving power (equivalent width converted to resolving
    # power at 5577 angstroms)
    ax.plot(watts,ew_r, color="white", lw=3, ls="--",
            label=r"Equiv. Resolving Power ($\lambda / W_{\lambda}$)")

    # Plot limiting resolving power
    ax.axhline(R_limit,color="orange", ls="--", lw=3,
               label=r"Limiting Resolving Power")

    # Axis Format
    ax.set_ylim(resolver.min(),resolver.max())
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Axis lables
    ax.set_xlabel(r"OI Auroral Power [W]")
    ax.set_ylabel(r"Resolving Power ($\lambda / \Delta \lambda$)")

    # Legend
    leg=ax.legend(loc=0, fontsize=15)
    leg.get_frame().set_alpha(0.0)
    for text in leg.get_texts():
        text.set_color("white")
        text.set_va('center') # va is alias for "verticalalignment"

    plt.show()

    if save_plots:
        fig.savefig("contrast_vs_R_vs_watts.pdf")
