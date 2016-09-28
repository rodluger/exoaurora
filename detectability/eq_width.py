# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:01:33 2016

@author: dflemin3

Visualize the line's equivalent width as a function of contrast ratio.
"""

from __future__ import print_function, division
import numpy as np
import contrast_utils as cu
import coronagraph_inputs as ci

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes

#Typical plot parameters that make for pretty plots
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 30.0

## for Palatino and other serif fonts use:
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rc('text', usetex=True)

show_spec = False

show_eqw = False

show_contrast_contour = False

show_time_contour = True

save_plots = False

# Init system object with Proxima Centauri defaults
proxcen = cu.System()

# Specify an auroral emission line
lam_0 = cu.OGREEN_LINE*1.0e-4 # microns

# Grid over auroral powers
watts = np.linspace(1.0e10,1.0e15,20)
equivalent_width = np.zeros_like(watts)

# FWHM from Eric
Temp = 200.0 # K
FWHM = 0.014*1.0e-4 * np.power(Temp / 200, 0.5)
dl = 50*FWHM
print("FWHM = "+str(FWHM*1.0e4)+" A")

for ii in range(len(watts)):

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
    Temp = 200.0 # K
    FWHM = 0.014*1.0e-4 * np.power(Temp / 200, 0.5)
    dl = 6.33e5*FWHM # Arbitrary scaling to make plot extend to 1 micron

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

    xscale = 0 #lam_0
    yscale = 1e18 #1.0

    ax.plot((wave_hires - xscale)*1.0e4,planet_flux*yscale,color="black",lw=3,
            label="Net Planet Flux")

    mask = (wave_hires >= lam_0 - FWHM) & (wave_hires <= lam_0 + FWHM)

    # Compute continuum
    continuum = np.median(ref_flux[mask])

    # Plot continuum
    #ax.axhline(y=continuum, color = "red", ls = "--", lw=3, label="Continuum")

    # Compute equivalent width
    ew = cu.equivalent_width(planet_flux[mask],wave_hires[mask],continuum)

    print("Equiv. Width = "+str(ew*1.0e4)+" A")

    # Plot equivalent width box
    ax.fill_between((wave_hires - xscale)*1.0e4, continuum*yscale, 2.0*continuum*yscale,
                    where=((wave_hires < lam_0 + ew/2) & (wave_hires > lam_0 - ew/2.)),
                    facecolor='green', alpha = 0.2)

    # Format
    ax.set_ylabel(r"Flux Density $\times 10^{18}$ [W/m$^2$/$\mu$m]")
    ax.set_xlabel(r"Wavelength [\AA]")
    xmin = (0.48 - xscale) * 1.0e4
    xmax = ((wave_hires - xscale)*1.0e4).max()
    xmask = (wave_hires * 1.0e4 >= xmin) & (wave_hires * 1.0e4 <= xmax)
    ymin = 0.0
    ymax = np.nanmax(planet_flux[xmask]*yscale) * 1.05
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin, ymax)
    #ax.set_yscale("log")

    # Inset axis
    """
    pad = ew/20
    xmin = lam_0 - ew/2 - pad
    xmax = lam_0 + ew/2 + pad
    mask2 = (wave_hires >= xmin) & (wave_hires <= xmax)
    ax2 = plt.axes([.27, .64, .32, .3])
    ax2.plot((wave_hires[mask2])*1.0e4,planet_flux[mask2]*yscale,color="black",lw=3)
    ax2.fill_between((wave_hires[mask2])*1.0e4, 0, 1.0*continuum*yscale,
                    where=((wave_hires[mask2] < lam_0 + ew/2) & (wave_hires[mask2] > lam_0 - ew/2.)),
                    facecolor='green', alpha = 0.5)
    ax2.set_xlim((xmin*1.0e4, xmax*1.0e4))
    #xticks = [xmin*1.0e4, lam_0*1.0e4, xmax*1.0e4]
    xticks = np.linspace(xmin*1.0e4, xmax*1.0e4, 5)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(["%.2f" % x for x in xticks])
    plt.setp(ax2.get_xticklabels(), fontsize=14, rotation=0)
    plt.setp(ax2.get_yticklabels(), fontsize=14, rotation=0)
    mark_inset(ax, ax2, loc1=2, loc2=4, fc="none", ec="black", zorder=0, alpha=0.3)
    #ax2.set_yscale("log")
    """

    fig.tight_layout()

    plt.show()

    if save_plots:
        fig.savefig("OI_ref_spec.pdf", bbox_inches='tight')

if show_contrast_contour:

    mpl.rcParams['font.size'] = 28.0

    infig_font_size = 18.0

    # Contrast, Resolving power, auroral power grid
    bins = 100

    # FWHM from Eric
    Temp = 200.0 # K
    FWHM = 0.014*1.0e-4 * np.power(Temp / 200, 0.5)

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

    mpl.rcParams['figure.figsize'] = (10,10)

    # Create figure
    #fig, ax = plt.subplots()
    fig = plt.figure()
    gs = gridspec.GridSpec(2,1, height_ratios=[0.05, 1.0])
    plt.subplots_adjust(wspace=0, hspace=0.16)
    ax = plt.subplot(gs[1])
    cbar_ax = plt.subplot(gs[0])


    vmin = contrast.min()
    vmax = contrast.max()

    # Convert watts to TW
    watts = watts / 1e12

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
    plt.clabel(cln, inline=1, fontsize=15.5,
               fmt=fmt, inline_spacing=30,
               use_clabeltext=True)#, manual=manual_locations)
    # Thicken the contour lines
    plt.setp(cln.collections, linewidth=2)

    # Colorbar
    cbar = fig.colorbar(cax, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(r"Planet-Star Contrast", rotation=0, labelpad=-85)
    cbar_ax.tick_params(labelsize=26.0)

    # Plot equivalent resolving power (equivalent width converted to resolving
    # power at 5577 angstroms)
    ax.plot(watts,ew_r, color="white", lw=3, ls="--",
            label=r"Equiv. Resolving Power ($\lambda / W_{\lambda}$)")

    # Plot limiting resolving power
    ax.axhline(R_limit,color="orange", ls="--", lw=3,
               label=r"FWHM Resolving Power ($\lambda / FWHM$)")

    # Axis Format
    ax.set_ylim(resolver.min(),resolver.max())
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Axis lables
    ax.set_xlabel(r"OI Auroral Power [TW]")
    ax.set_ylabel(r"Resolving Power ($\lambda / \Delta \lambda$)")

    # Legend
    leg=ax.legend(loc=0, fontsize=infig_font_size)
    leg.get_frame().set_alpha(0.0)
    for text in leg.get_texts():
        text.set_color("white")
        text.set_va('center') # va is alias for "verticalalignment"

    plt.show()

    if save_plots:
        fig.savefig("contrast_vs_R_vs_watts_new.pdf", bbox_inches='tight')

if show_time_contour:

    mpl.rcParams['font.size'] = 28.0

    infig_font_size = 18.0

    # Contrast, Resolving power, auroral power grid
    bins = 200

    # FWHM from Eric
    Temp = 200.0 # K
    FWHM = 0.014*1.0e-4 * np.power(Temp / 200, 0.5)

    # Define a desired signal-to-noise
    SN = 6.0

    # Default Proxima Centauri system
    proxcen = cu.System()

    # Create arrays
    watts = np.logspace(9,15,bins)   # Auroral power
    resolver = np.logspace(1,6,bins) # Resolving power
    contrast = np.zeros((bins,bins)) # Planet-star contrast ratio
    ew = np.zeros_like(watts)        # Equivalent width of line
    exptime = np.zeros((bins,bins))  # Exposure time to get SN

    R_limit = lam_0/FWHM

    # Loop over auroral power grid
    for ii in range(len(watts)):
        # Loop over resolving power grid
        for jj in range(len(resolver)):

            # Compute auroral flux line profile and hi-res wavelength array
            # Compute delta_lam
            dl = lam_0/resolver[jj]

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

            # BEGIN INT TIME CALC

            # How many resolution elements span the FWHM?
            # Set number of bins to integrate over as the ceiling so there is
            # always one resolution element in the case that the resolution
            # element is much larger than the FWHM
            nbins = int(np.ceil(FWHM / dl))

            # LUVOIR-concept space-based telescope
            tmt = cu.Telescope(D=30., R=resolver[jj], bins=nbins)

            # Compute integration time
            exptime[ii,jj] = tmt.int_time(SN, proxcen, watts[ii], lam_0,
                                             coronagraph=True)
            # End inner loop

        # Compute equivalent width for each power step

        # Compute continuum
        continuum = np.median(ref_flux)
        ew[ii] = cu.equivalent_width(planet_flux,wave_hires,continuum)

    # Compute equivlent resolving power (?)
    ew_r = lam_0/ew

    # Save data?
    if True:
        np.savez("TMTC_int_time_C"+("%.0e" % ci.C)+"_noise2.npz", resolver, watts, exptime)

    # Make plot?
    if False:

        mpl.rcParams['figure.figsize'] = (10,10)

        # Create figure
        #fig, ax = plt.subplots()
        fig = plt.figure()
        gs = gridspec.GridSpec(2,1, height_ratios=[0.05, 1.0])
        plt.subplots_adjust(wspace=0, hspace=0.16)
        ax = plt.subplot(gs[1])
        cbar_ax = plt.subplot(gs[0])

        vmin = exptime.min()
        vmax = 365*24 #exptime.max()

        # Convert watts to TW
        watts = watts / 1e12

        # Color contour
        cax = ax.pcolor(watts, resolver, exptime.T, cmap="viridis_r",
                        norm=colors.LogNorm(vmin=vmin,vmax=vmax))

        # Contour Lines
        contour_levels = np.array([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])
        # Automate "manual locations" for contour labels
        ypos = 2e5
        iy = (np.abs(resolver-ypos)).argmin() # R. Power index
        ix = np.array([np.abs(exptime[:, iy] - clvl).argmin() for clvl in contour_levels])
        manual_locations = ((watts[ix[i]], resolver[iy]) for i in range(len(contour_levels)))
        # Log levels
        fmt = ticker.LogFormatterMathtext()
        cln = ax.contour(watts, resolver, exptime.T, contour_levels,
                         colors="black")
        plt.clabel(cln, inline=1, fontsize=15.5, fmt=fmt, inline_spacing=30,
                   use_clabeltext=True, manual=manual_locations)
        # Thicken the contour lines
        plt.setp(cln.collections, linewidth=2)

        # Colorbar
        cbar = fig.colorbar(cax, cax=cbar_ax, orientation="horizontal")
        cbar.set_label(r"Integration Time [hrs]", rotation=0, labelpad=-85)
        cbar_ax.tick_params(labelsize=26.0)

        # Plot equivalent resolving power (equivalent width converted to resolving
        # power at 5577 angstroms)
        ax.plot(watts,ew_r, color="white", lw=3, ls="--",
                label=r"Equiv. Resolving Power ($\lambda / W_{\lambda}$)")

        # Plot limiting resolving power
        ax.axhline(R_limit,color="orange", ls="--", lw=3,
                   label=r"FWHM Resolving Power ($\lambda / FWHM$)")

        #ax.axhline(115000.,color="red", ls="--", lw=3,
        #           label=r"$R=115,000$")

        # Axis Format
        ax.set_ylim(resolver.min(),resolver.max())
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Axis lables
        ax.set_xlabel(r"OI Auroral Power [TW]")
        ax.set_ylabel(r"Resolving Power ($\lambda / \Delta \lambda$)")

        # Legend
        leg=ax.legend(loc=0, fontsize=infig_font_size)
        leg.get_frame().set_alpha(0.0)
        for text in leg.get_texts():
            text.set_color("white")
            text.set_va('center') # va is alias for "verticalalignment"

        plt.show()

        if save_plots:
            fig.savefig("exptime_B_inline_TMTC-7_noInstrumental_new.pdf", bbox_inches='tight')
