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

#Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

show_plots = False

save_plots = False

show_contour = True

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
if show_plots:

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

    #plt.clf()

    # Plot 1 test case

    fig, ax = plt.subplots()

    # FWHM from Eric
    FWHM = 0.0186*1.0e-4
    dl = 50000*FWHM

    wave_hires = cu.make_wave_array(lam_0,dl,FWHM)

     #Compute auroral flux line profile and hi-res wavelength array
    aurora = cu.create_auroral_line(lam_0, 1e11, 1.0, proxcen,
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
    ax.axhline(y=continuum, color = "red", ls = "--", lw=3, label="Continuum")

    # Compute equivalent width
    ew = cu.equivalent_width(planet_flux[mask],wave_hires[mask],continuum)

    # Plot equivalent width box
    ax.fill_between((wave_hires-lam_0)*1.0e4, continuum, 2.0*continuum,
                    where=((wave_hires < lam_0 + ew/2) & (wave_hires > lam_0 - ew/2.)),
                    facecolor='green', alpha = 0.2)

    # Format
    ax.set_ylabel(r"Flux [W/m$^2$/$\mu$m]")
    ax.set_xlabel(r"Wavelength [\AA]")
    ax.set_xlim(((wave_hires-lam_0)*1.0e4).min(),
                ((wave_hires-lam_0)*1.0e4).max())
    ax.set_yscale("log")

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
    #plt.clabel(cln, inline=1, fontsize=15, fmt=fmt, manual=manual_locations)
    # Thicken the contour lines
    zc = cln.collections; plt.setp(zc, linewidth=2)

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

    # Format
    ax.set_ylim(resolver.min(),resolver.max())
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Axis lables
    ax.set_xlabel(r"OI Auroral Power [W]")
    ax.set_ylabel(r"Resolving Power ($\lambda / \Delta \lambda$)")

    # legend
    leg=ax.legend(loc=0, fontsize=14)
    leg.get_frame().set_alpha(0.7)

    plt.show()

    if save_plots:
        fig.savefig("contrast_vs_R_vs_watts.pdf")
