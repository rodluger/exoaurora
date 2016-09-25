#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
search.py
---------

Searching the HARPS data for the OI emission signal.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from pool import Pool
import matplotlib as mpl; mpl.use('Agg')
mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Times New Roman']
from kepler import RadialVelocity
import numpy as np
from sklearn.decomposition import PCA
from wpca import WPCA
from scipy.stats import norm, tukeylambda
from scipy.stats import t as student_t
from scipy.signal import savgol_filter, medfilt
import glob
import os, sys
import itertools
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator, ScalarFormatter, FormatStrFormatter
from matplotlib.patches import Rectangle
import matplotlib.mlab as mlab
import subprocess
import argparse
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
    
# Empirical system offset, see comment below
SYSTEM_OFFSET = 0.999927

class Spectrum(object):
  '''
  Important aurora/airglow lines in angstroms.
  
  '''
  
  NitrogenUVI = 3914.
  NitrogenUVII = 4278.
  NitrogenUVIII = 4709.
  NitrogenBlueI = 5199.
  NitrogenBlueII = 5201.
  OxygenGreen = 5577.345
  OxygenRedI = 6300.308 
  OxygenRedII = 6363.790
  HAlpha = 6562.8
  MysteryLineI = 5515.46 # Prominent in the star, not sure what it is
  SodiumDI = 5889.950
  SodiumDII = 5895.924
  HeliumYellow = 5875.6 # Activity tracer?
  MysteryLineII = 5540.07 # Prominent in both star and earth... Weird!
  
  EarthLines = [MysteryLineII, NitrogenUVI, OxygenGreen, OxygenRedI, OxygenRedII]
  StarLines = [MysteryLineII, OxygenGreen, MysteryLineI, SodiumDI, SodiumDII, HeliumYellow]

class Planet(object):
  '''
  Miscellaneous planet properties.
  
  '''
  
  def __init__(self):
    '''
    
    '''
    
    pass
  
  def phase(self, date):
    '''
    
    '''
    
    return ((date - self.t0) % self.period) / self.period
  
  def earth_phase(self, date):
    '''
    
    '''
    
    return (date % 365.25) / 365.25
  
  def rv(self, date):
    '''
    
    '''
    
    return RadialVelocity(date, self)
  
  def doppler_factor(self, date):
    '''
        
    '''
    
    c = 299792458.
    v = self.rv(date)
    one_plus_z = np.sqrt((1 + v / c) / (1 - v / c))
    # This is > 1 for positive RV (redshift)
    # and < 1 for negative RV (blueshift)
    return np.log10(one_plus_z)

  def earth_doppler_factor(self, earth_rv):
    '''
    
    '''
    
    c = 299792458.
    one_plus_z = np.sqrt((1 + earth_rv / c) / (1 - earth_rv / c))
    return np.log10(one_plus_z) - np.log10(SYSTEM_OFFSET)
    
class ProxCenB(Planet):
  '''
  
  '''
  
  def __init__(self):
    '''
    
    '''
    
    self.t0 = 2451634.73146
    self.stellar_mass = 0.12
    self.mass = 1.27
    self.period = 11.186
    self.eccentricity = 0.
    self.mean_longitude = 110.
    self.arg_periastron = 310.
    self.inclination = 90.
 
def GaussianLine(x, line = Spectrum.OxygenGreen, fwhm = 0.05, A = 1.):
  '''
  
  '''
  
  sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
  return A / np.sqrt(2 * sig ** 2 * np.pi) * np.exp(-(x - line) ** 2. / (2 * sig ** 2))
   
def ApplyDopplerFactor(wav, fac):
  '''
  
  '''

  return 10 ** (np.log10(wav) + fac)

def ReadFITS():
  '''
  
  '''
  
  # Time, wavelength, flux, and flux error arrays
  jd = []
  wav = []
  flx = []
  err = []
  dat = []
  erv = []
  exp = []
  air = []
  
  # Load in all the fits files
  # This is be coded up to allow importing of the UVES dataset as well
  for dataset in ['HARPS']:
    
    files = glob.glob(os.path.join(os.path.dirname(__file__), dataset, '*.fits'))
    for file in files:

      # Read the data
      f = pyfits.open(file)
      x = np.array(f[1].data['wave'][0])
      try:
        y = np.array(f[1].data['flux_reduced'][0])
      except:
        y = np.array(f[1].data['flux'][0])
      try:
        e = np.array(f[1].data['err_reduced'][0])
      except:
        e = None
      
      # Correct from MJD to JD
      date = 0.5 * (f[0].header['MJD-OBS'] + f[0].header['MJD-END'])
      date += 2400000.5

      # Build our arrays
      jd.append(date)
      wav.append(x)
      flx.append(y)
      err.append(e)
      dat.append(dataset)
      erv.append(f[0].header['HIERARCH ESO DRS BERV'] * 1.e3) 
      exp.append(f[0].header['texptime'])
      air.append(0.5 * (f[0].header['HIERARCH ESO TEL AIRM START'] + f[0].header['HIERARCH ESO TEL AIRM END']))
      f.close()
  
  return np.array(jd), wav, flx, err, np.array(dat), np.array(erv), np.array(exp), air

def GetData():
  '''
  
  '''
  
  filename = os.path.join(os.path.dirname(__file__), 'data.npz')
  
  if not os.path.exists(filename):
  
    # Grab the data from the FITS files
    print("Loading FITS files...")
    jd, wavs, flx, err, dataset, earth_rv, exptime, air = ReadFITS()
    
    # Bin all spectra to the first wavelength grid (regularized)
    wav = np.linspace(wavs[0][0], wavs[0][-1], len(wavs[0]))
    flx = np.array([(np.interp(wav, wavs[i], flx[i])) for i in range(len(wavs))])
    
    # Apply the system offset. This was determined empirically by matching the
    # stellar Na D I and II lines. This calibration is accurate to ~ 0.01 A.
    # Note that HARPS spectra are actually all in the same reference frame
    # (heliocentric? I'm not sure), so this offset puts us in the frame of 
    # Proxima Centauri for all spectra.
    wav /= SYSTEM_OFFSET
    
    # Save    
    print("Saving dataset...")
    result = {'wav': wav, 'jd': jd, 'flx': flx, 'err': err, 'dataset': dataset,
              'earth_rv': earth_rv, 'exptime': exptime, 'air': air}
    np.savez(filename, **result)
  
  else:
    
    print("Loading dataset...")
    result = dict(np.load(filename))

  return result
  
def RemoveStellarLines(wav, flx, weights = None, npc = 5, inds = None):
  '''
  
  '''
    
  # Get the principal components
  if npc > 0:
    pca = WPCA(n_components = npc)
    if inds is None:
      xpca = pca.fit_transform(flx.T, weights = weights)
    else:
      xpca = pca.fit_transform(flx[inds].T, weights = weights[:,inds])
    xpca = np.hstack([np.ones(xpca.shape[0]).reshape(-1, 1), xpca])
  else:
    xpca = np.ones_like(wav)
    
  # Fit a linear model to each spectrum and subtract it
  for i in range(len(flx)):
    w = np.linalg.solve(np.dot(xpca.T, xpca), np.dot(xpca.T, flx[i]))
    model = np.dot(xpca, w)
    model -= np.nanmedian(model)
    flx[i] -= model 
       
  return flx

def Compute(planet = ProxCenB(), data = None, line = Spectrum.OxygenGreen, plot = True, 
            plot_sz = 10, npc = 10, frame = 'planet', wpca_sz = 250, wpca_mask_sz = 0.2, 
            bin_sz = 0.05, mask_star = True, med_mask_sz = 5, fwhm = 0.05,
            inject_contrast = 0, inject_planet = ProxCenB(), airmass_correction = False,
            crop_outliers = True, clobber = False, quiet = False, spectrum_filter = 'True', 
            max_frac_noise = 3, wpca_weights = 'exptime', filter_sz = 1.):
  '''
  
  '''
  
  # Check what rest frame we're plotting in
  assert frame in ['star', 'planet', 'earth'], "Argument `frame` must be one of `star`, `planet`, or `earth`."
    
  # Read the raw data
  if data is None:
    data = GetData()
  jd = data['jd']
  wav = data['wav']
  flx = data['flx']
  err = data['err']
  earth_rv = data['earth_rv']
  air = data['air']
  exp = data['exptime']
  dataset = data['dataset']
  total_exp = 0
  
  # Some wavelength-specific params. The noise changes a lot from the UV to the vis,
  # so we need to change our plotting scale.
  if line < 4000:
    spec_amp = 0.001
    stack_lims = [(0.6, 1.4), (0.8, 1.2)]
  elif line < 4500:
    spec_amp = 0.0025
    stack_lims = [(0.9, 1.1), (0.95, 1.05)]
  else:
    spec_amp = 0.01
    stack_lims = [(0.9840, 1.0160), (0.9960, 1.0100)]
  
  # Crop to a smaller window centered on the line
  inds = np.where((wav >= line - wpca_sz / 2.) & (wav <= line + wpca_sz / 2.))[0]
  wav = wav[inds]
  flx = flx[:,inds]

  # Inject a signal for injection/recovery tests?
  if inject_contrast > 0.:
    if not quiet: print("Injecting signal...")
    for i in range(len(flx)):
      # Get planet line wavelength in the stellar frame
      x = ApplyDopplerFactor(line, inject_planet.doppler_factor(jd[i]))
      # Get indices corresponding to the FWHM of the line (TODO: Check this)
      inds = np.where(np.abs(wav - x) <= fwhm / 2.)[0]
      # Compute the flux in the stellar spectrum at those indices
      z = np.trapz(flx[i][inds], wav[inds])
      # Compute the line profile. The amplitude is contrast * z
      flx[i] += GaussianLine(wav, line = x, fwhm = fwhm, A = inject_contrast * z)

  # Get masks in the stellar frame
  if not quiet: print("Applying masks...")
  weights = np.ones((len(wav), len(flx)))
  earth = [[] for i in range(len(flx))]
  star = [[] for i in range(len(flx))]
  for i in range(len(flx)):
  
    # Set weights according to exposure time
    if wpca_weights == 'exptime':
      weights[:,i] = np.sqrt(exp[i])
    elif wpca_weights == 'std':
      weights[:,i] = 1. / np.nanstd(flx[i] / np.nanmedian(flx[i]))
    else:
      raise ValueError('Parameter `wpca_weights` must be one of `exptime` or `std`.')
      
    # Add mask for each of the star's lines
    for star_line in Spectrum.StarLines:
      m = list(np.where((wav > star_line - wpca_mask_sz / 2.) & (wav < star_line + wpca_mask_sz / 2.))[0])
      star[i].append(m)
  
    # Add mask for each of Earth's lines in the stellar frame
    for earth_line in Spectrum.EarthLines:
      x = ApplyDopplerFactor(earth_line, planet.earth_doppler_factor(earth_rv[i]))
      m = list(np.where((wav > x - wpca_mask_sz / 2.) & (wav < x + wpca_mask_sz / 2.))[0])
      earth[i].append(m)
  
    # Combine the masks
    if mask_star:
      m = [item for sublist in earth[i] + star[i] for item in sublist]
    else:
      m = [item for sublist in earth[i] for item in sublist]
    m = np.array(sorted(list(set(m))), dtype = int)
  
    # Set the masked weights to a very small number (~zero)
    weights[m,i] = 1e-10

  # Remove stellar features
  if npc > 0:
    if not quiet: print("Removing stellar features...")
    flx = RemoveStellarLines(wav, flx, weights = weights, npc = npc, inds = None)

  # Remove earth features
  if airmass_correction:
    if not quiet: print("Applying airmass correction...")
  
    # Shift to Earth frame. See note on the NEGATIVE sign below.
    x = [ApplyDopplerFactor(wav, -planet.earth_doppler_factor(e)) for e in earth_rv]
    flxe = np.array([np.interp(wav, xi, flxi) for xi, flxi in zip(x, flx)])
  
    # Subtract a linear function of the airmass at each wavelength
    med = np.nanmedian(flxe, axis = 1)
    w = np.sqrt(exp)
    m = np.zeros_like(wav)
    b = np.zeros_like(wav)
    for j, _ in enumerate(wav):
      m[j], b[j] = np.polyfit(air, flxe[:,j] / med, 1, w = w)
      corr = med * (b[j] + air * m[j])
      corr -= np.median(corr)
      flxe[:,j] -= corr

    # Shift back to stellar frame
    x = [ApplyDopplerFactor(wav, planet.earth_doppler_factor(e)) for e in earth_rv]
    flx = np.array([np.interp(wav, xi, flxei) for xi, flxei in zip(x, flxe)])
  
  # Compute the planet doppler factor now, since we use it a couple times below
  pdf = [planet.doppler_factor(j) for j in jd]
  
  # Crop outliers in planet frame? This is useful only for the FAP calculation,
  # since we don't want an outlier in a single spectrum affecting our FAP. Note
  # that this is somewhat asymmetrical, since we're not removing outliers in the
  # 5577A window: we could then detect a "signal" due to a single outlier in one
  # spectrum. But we must allow for time-variable emission; plus, if we do find
  # such a signal because of an outlier, it would be easy to tell which spectrum 
  # it came from, and we can deal with it then.
  if crop_outliers:

    # Shift to planet frame. If the planet is moving away from us (redshifted),
    # we must BLUESHIFT the wavelength array to move into the planet frame,
    # hence the NEGATIVE sign below.
    x = [ApplyDopplerFactor(wav, -p) for p in pdf]
    flxp = np.array([np.interp(wav, xi, flxi) for xi, flxi in zip(x, flx)])
  
    # Remove 10-sigma outliers in each wavelength bin, except where we're looking
    # for our signal (which could be highly variable).
    med = np.nanmedian(flxp, axis = 1)
    for j, _ in enumerate(wav):
      if (wav[j] > line - wpca_mask_sz / 2.) and (wav[j] < line + wpca_mask_sz / 2.):
        continue
      y = flxp[:,j] / med
      m = np.nanmedian(y)
      MAD = 1.4826 * np.nanmedian(np.abs(y - m))
      bad = np.where((y > m + 10. * MAD) | (y < m - 10. * MAD))[0]
      for b in bad:
        flxp[b,j] = np.nanmedian(flxp[b,max(0,j-50):j+50])
  
    # Shift back to stellar frame. The sign below is the opposite of that above.
    x = [ApplyDopplerFactor(wav, p) for p in pdf]
    flx = np.array([np.interp(wav, xi, flxpi) for xi, flxpi in zip(x, flxp)])
  
  # The stacked flux
  xstack = np.array(wav)
  fstack = np.zeros_like(xstack)
  
  # The figure
  if plot:
    fig = pl.figure(figsize = (10,10))
    ax = [None, None, None, None]
    ax[0] = pl.subplot2grid((8, 8), (0, 0), colspan=8, rowspan=6)
    ax[1] = pl.subplot2grid((8, 8), (6, 0), colspan=8, rowspan=1, sharex = ax[0])
    ax[2] = pl.subplot2grid((8, 8), (7, 0), colspan=8, rowspan=1, sharex = ax[0])
  else:
    fig = None
    ax = None
    
  # Loop over spectra
  if not quiet: print("Stacking spectra...")
  pcb = [[] for i in range(len(flx))]
  for i in range(len(flx)):
    
    # User-defined filter
    if not eval(spectrum_filter):
      continue
    
    # Compute the planet mask in the stellar frame. If the planet is moving away
    # from us (redshifted), we REDSHIFT the line accordingly.
    x = ApplyDopplerFactor(line, pdf[i])
    pcb[i] = np.where((wav > x - wpca_mask_sz / 2.) & (wav < x + wpca_mask_sz / 2.))[0]
    
    # Plot and stack only spectra where there's no overlap between planet and Earth / planet and star
    if len(set(pcb[i]) & set([item for sublist in earth[i] for item in sublist])) + len(set(pcb[i]) & set([item for sublist in star[i] for item in sublist])) == 0:
      
      # The de-trended spectra
      x = np.array(wav)
      y = planet.phase(jd[i]) + spec_amp * flx[i] / np.nanmedian(flx[i])
      y0 = np.array(flx[i])
      
      # Is the spectrum too noisy?
      if np.std(y0 / np.nanmedian(y0)) > max_frac_noise:
        continue
           
      # DOPPLER SHIFTING SANITY CHECK:
      #
      #   if (planet.phase(jd[i]) < 0.22) and (planet.phase(jd[i]) > 0.18):
      #     N = int((2457400 - jd[i]) / planet.period) + 1.
      #     print(jd[i] + N * planet.period)
      # 
      # At a phase of 0.2, the figure in our paper shows that the planet is
      # at a maximum redshift relative to the star (and to Earth). This means that
      # the star should be at a maximum blueshift relative to the barycenter,
      # corresponding to a minimum in the RV curve.
      # The above code snippet yields JD = 2457408.82233, which by inspection of
      # Figure 3 in Anglada-Escude et al. (2016), is correct! 
           
      # Mask problematic stellar and earth lines
      for e in earth[i]:
        if len(e):
          y0[e] = np.nanmedian(y0[np.where(np.abs(x - wav[e][0]) < med_mask_sz)])
      
      for s in star[i]:
        if len(s):
          y0[s] = np.nanmedian(y0[np.where(np.abs(x - wav[s][0]) < med_mask_sz)])
                
      # Switch frames?
      if frame == 'star':
        # This is the default frame the data is stored in
        pass
      elif (frame == 'planet'):
        # Doppler-shift into the frame of the planet. If the planet is moving away
        # from us/from the star (redshifted), we must BLUESHIFT its lines to get into
        # its frame, hence the NEGATIVE sign.
        x = ApplyDopplerFactor(x, -pdf[i])
      else:
        # Doppler-shift into the frame of the Earth
        x = ApplyDopplerFactor(x, -planet.earth_doppler_factor(earth_rv[i]))
    
      # Stacking: interpolate to original grid and add
      fstack += np.interp(xstack, x, y0)
      total_exp += exp[i]
          
      # Plot!
      if plot:
        mask = np.array(sorted(list(set(np.concatenate([pcb[i], [item for sublist in earth[i] for item in sublist], [item for sublist in star[i] for item in sublist]])))), dtype = int)
        ymask = np.array(y)
        ymask[mask] = np.nan
        ax[0].plot(x, ymask, color = 'k', alpha = 0.25)
        ax[0].plot(x[pcb[i]], y[pcb[i]], color = 'b', alpha = 0.75, lw = 1.5)
        for e in earth[i]:
          ax[0].plot(x[e], y[e], color = 'g', alpha = 0.75, lw = 1.5)
        for s in star[i]:
          ax[0].plot(x[s], y[s], color = 'r', alpha = 0.75, lw = 1.5)
  
  # High pass median filter? `filter_sz` angstrom(s) wide
  if filter_sz:
    window = int(filter_sz / np.nanmedian(xstack[1:] - xstack[:-1]))
    if not (window % 2):
      window += 1
    filt = medfilt(fstack, window)
    fstack = fstack - filt + np.nanmedian(fstack)

  # Plot the stacked flux
  if plot:
    a = np.argmax(xstack > line - wpca_mask_sz / 2.)
    b = np.argmax(xstack > line + wpca_mask_sz / 2.)
    z = fstack / np.nanmedian(fstack)
    ax[1].plot(xstack[a:b], z[a:b], 'b-', alpha = 0.75)
    z[a+1:b-1] = np.nan
    ax[1].plot(xstack, z, 'k-', alpha = 0.75)
  
  # Get the stacked, binned flux
  bin_centers = np.append(np.arange(line, line - wpca_sz / 2., -bin_sz)[::-1], np.arange(line, line + wpca_sz / 2., bin_sz)[1:])
  bflx = np.zeros_like(bin_centers)
  bnrm = np.zeros_like(bin_centers)
  for x, z in zip(xstack, fstack / np.nanmedian(fstack)):
    if (x >= bin_centers[0] - bin_sz / 2.) and (x < bin_centers[-1] + bin_sz / 2.):
      i = np.argmin(np.abs(bin_centers - x))
      bflx[i] += z
      bnrm[i] += 1
  bflx /= bnrm
  
  # Pad the binned flux by 5% on each side to get rid of edge effects
  pad = int(0.05 * len(bin_centers))
  bin_centers = bin_centers[pad:-pad]
  bflx = bflx[pad:-pad]
  
  # This is the signal in the bin centered at the line
  signal = bflx[len(bflx) // 2]
  std = np.nanstd(bflx)
  snr = (signal - 1.) / std
        
  # The histogram inset
  if plot:
    ax[3] = pl.axes([0.195 + 0.01, 0.36 + 0.01, 0.2, 0.15], zorder = 2)
    n, bins, _ = ax[3].hist((bflx - 1.) / std, bins = 30, color = 'w')
    d = np.digitize(snr, bins)
    if d == len(bins): 
      d -= 1
        
    ax[3].axvline(bins[d] - 0.5 * (bins[1] - bins[0]), color = 'r', ls = '--')
    ax[3].set_yscale('log')
    ax[3].set_ylim(0.8, 1e3)
    ax[3].margins(0.1, None)
    ax[3].set_ylabel(r'$\log\ N$', fontsize = 14)
    ax[3].set_xlabel(r'SNR', fontsize = 14)
    rect = Rectangle((0.125 + 0.01, 0.36 - 0.055 + 0.01), 0.28, 0.235, facecolor='w', edgecolor='k',
                      transform=fig.transFigure, alpha = 0.85, zorder=1)
    fig.patches.append(rect)

  # The binned flux
  if plot:
    ax[2].plot(bin_centers, bflx, 'k.', alpha = 0.1)
    ax[2].bar(bin_centers - bin_sz / 2., bflx, color = 'None', width = bin_sz)
    ax[2].axhline(1., color = 'r', ls = '--', alpha = 0.5)
    lo, hi = np.min(bflx), np.max(bflx)
    rng = hi - lo
    ax[2].set_ylim(lo - 0.1 * rng, hi + 0.1 * rng)
  
  # Print some info
  if not quiet:
    print('Total integration time: %.1f hours' % (total_exp / 3600.))
    print('SNR: %.2f' % snr)
  
  # Appearance
  if plot:
    ax[0].axvline(line, color = 'k', ls = '--') 
    ax[0].set_xlim(line - plot_sz / 2., line + plot_sz / 2.)
    ax[0].set_ylim(-0.025, 1.025)
    ax[1].axvline(line, color = 'k', ls = '--') 
    ax[2].axvline(line, color = 'k', ls = '--') 
    ax[1].set_xlim(line - plot_sz / 2., line + plot_sz / 2.)
    ax[2].set_xlim(line - plot_sz / 2., line + plot_sz / 2.)
    ax[1].yaxis.set_major_locator(MaxNLocator(nbins = 4))
    ax[2].yaxis.set_major_locator(MaxNLocator(nbins = 4))
    ax[2].set_xlabel('Wavelength (Angstroms)', fontsize = 24)
    ax[0].set_ylabel('Planet Orbital Phase', fontsize = 24)
    ax[1].set_ylabel('Stacked', fontsize = 18)
    [tick.label.set_fontsize(8) for tick in ax[1].yaxis.get_major_ticks()]
    ax[2].set_ylabel('Binned', fontsize = 18)
    [tick.label.set_fontsize(8) for tick in ax[2].yaxis.get_major_ticks()]
    
    # Axes limits for the stacked/binned flux subplots
    if stack_lims is not None:
      ax[1].set_ylim(*stack_lims[0])
      ax[2].set_ylim(*stack_lims[1])
    
    pl.setp(ax[0].get_xticklabels(), visible = False)
    pl.setp(ax[1].get_xticklabels(), visible = False)
    ax[0].ticklabel_format(useOffset=False)
    ax[1].yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    ax[2].yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    ax[0].plot(0, 0, color = 'b', lw = 2, label = 'Planet')
    ax[0].plot(0, 0, color = 'g', lw = 2, label = 'Earth')
    ax[0].plot(0, 0, color = 'r', lw = 2, label = 'Star')
    ax[0].legend(loc = 'upper right', fontsize = 20)

  return {'fig': fig, 'ax': ax, 'signal': bflx[len(bflx) // 2],
          'bins': bin_centers, 'bflx': bflx, 'snr': snr}

class SearchWrap(object):
  '''
  
  '''
  
  def __init__(self, **kwargs):
    '''
    
    '''
    
    self.planet = ProxCenB()
    self.kwargs = kwargs
    self.kwargs.update({'quiet': True, 'plot': False})
  
  def __call__(self, params):
    '''
    
    '''
    
    inclination, period, stellar_mass, mean_longitude = params
    self.planet.inclination = inclination
    self.planet.mass = 1.27 / np.sin(self.planet.inclination * np.pi / 180)
    self.planet.period = period
    self.planet.stellar_mass = stellar_mass
    self.planet.mean_longitude = mean_longitude
    res = Compute(planet = self.planet, **self.kwargs)
    return res['bflx']

def PBSSearch(line = Spectrum.OxygenGreen, nodes = 12, ppn = 12, walltime = 100):
  '''
  Submits a PBS cluster job to do the line search.

  :param int walltime: The number of hours to request. Default `100`
  :param int nodes: The number of nodes to request. Default `5`
  :param int ppn: The number of processors per node to request. Default `12`
  
  '''
  
  # Submit the cluster job      
  pbsfile = os.path.join(os.path.dirname(__file__), 'search.pbs')
  name = '%d' % np.floor(line)
  str_n = 'nodes=%d:ppn=%d,feature=%dcore' % (nodes, ppn, ppn)
  str_w = 'walltime=%d:00:00' % walltime
  str_v = 'LINE=%.3f,NODES=%d,SEARCH_DIR=%s' % (line, nodes, os.path.dirname(os.path.realpath(sys.argv[0])))
  str_out = os.path.join(os.path.dirname(__file__), '%s.log' % name)
  qsub_args = ['qsub', pbsfile, 
               '-v', str_v, 
               '-o', str_out,
               '-j', 'oe', 
               '-N', name,
               '-l', str_n,
               '-l', str_w]          
  print("Submitting the job...")
  subprocess.call(qsub_args)

def Search(inclination = np.arange(30., 90., 2.), 
           period = np.arange(11.186 - 3 * 0.002, 11.186 + 3 * 0.002, 0.002 / 2),
           mean_longitude = np.arange(110. - 3 * 8., 110. + 3 * 8., 8. / 2), 
           stellar_mass = [0.120], clobber = False, 
           period_ticks = [11.182, 11.184, 11.186, 11.188, 11.190],
           mean_longitude_ticks = [90., 100., 110., 120., 130.],
           inclination_ticks = [35, 45, 55, 65, 75, 85], pool = None,
           **kwargs):
  '''
  
  '''
  
  # Begin multiprocessing
  with Pool() as pool:
  
    # Single or multithreaded?
    if pool is None:
      M = map
    else:
      M = pool.map
  
    line = kwargs.get('line', Spectrum.OxygenGreen)
    pref = "%d" % np.floor(line)
    search_file = os.path.join(os.path.dirname(__file__), '%s_search.npz' % pref)
    data = GetData()
    kwargs.update({'data': data})
  
    if clobber or not os.path.exists(search_file):

      # Get the bin array
      wpca_sz = kwargs.get('wpca_sz', 250)
      bin_sz = kwargs.get('bin_sz', 0.05)
      bins = np.append(np.arange(line, line - wpca_sz / 2., -bin_sz)[::-1], np.arange(line, line + wpca_sz / 2., bin_sz)[1:])
      pad = int(0.05 * len(bins))
      bins = bins[pad:-pad]
  
      # Loop over planet params
      print("Running grid search...")
      params = list(itertools.product(inclination, period, mean_longitude, stellar_mass))
      sw = SearchWrap(**kwargs)
      bflx = np.array(list(M(sw, params))).reshape((-1, len(inclination), len(period), len(mean_longitude), len(stellar_mass)))
    
      print("Saving...")
      np.savez(search_file, bins = bins, bflx = bflx, inclination = inclination, period = period,
               mean_longitude = mean_longitude, stellar_mass = stellar_mass)
    else:
      print("Loading saved search...")
      data = np.load(search_file)
      bins = data['bins']
      bflx = data['bflx']
      inclination = data['inclination']
      period = data['period']
      mean_longitude = data['mean_longitude']
      stellar_mass = data['stellar_mass']
    
    # Here we compute the distribution of the values of the
    # maximum signals at each wavelength (to compute significance later)
    bmax = np.max(bflx, axis = (1,2,3,4))
    bmu = np.nanmean(bmax)
    bstd = np.nanstd(bmax)
  
    # The binned flux at the line as a function of all the grid params
    bline = bflx[np.argmin(np.abs(line - bins))]
    blinemax = bmax[np.argmin(np.abs(line - bins))]
  
    # The best-fitting planet params
    planet = ProxCenB()
    i, p, m, s = np.unravel_index(np.nanargmax(bline), bline.shape)
    planet.inclination = inclination[i]
    planet.mass = 1.27 / np.sin(planet.inclination * np.pi / 180)
    planet.period = period[p]
    planet.stellar_mass = stellar_mass[s]
    planet.mean_longitude = mean_longitude[m]

    # --- FIGURE 1: Injection test (~8 sigma). This is our nominal detection threshold.
    # Note that we inject 10 angstroms redward of the line we're interested in, so we
    # don't stack an injected signal on top of an actual signal. Note also that we 
    # purposefully inject assuming the same planet params as those of the peak signal,
    # so that the number of spectra and the noise properties are the same.
    print("Plotting figure 1...")
    inj_kwargs = dict(kwargs)
    inj_kwargs.update({'line': line - 10., 'quiet': True})
    res = Compute(inject_contrast = 2e-2, inject_planet = planet, planet = planet, **inj_kwargs)
    fig1 = res['fig']
    inj_sig = (res['signal'] - bmu) / bstd
    fig1.suptitle('%d$\sigma$ Injected Signal' % inj_sig, fontsize = 30, y = 0.95)
    fig1.savefig('%s_injection_river.pdf' % pref, bbox_inches = 'tight')

    # --- FIGURE 2: Triangle plot
    print("Plotting figure 2...")
    fig2, ax2 = pl.subplots(3,3)
    fig2.subplots_adjust(wspace = 0.08, hspace = 0.1, top = 0.975, bottom = 0.15)
    # The marginalized distributions
    ax2[0,0].plot(inclination, np.max(bline, axis = (1,2,3)) - 1., color = 'k')
    ax2[1,1].plot(period, np.max(bline, axis = (0,2,3)) - 1., color = 'k')
    ax2[2,2].plot(mean_longitude, np.max(bline, axis = (0,1,3)) - 1., color = 'k')
    # Indicate the 1-sigma bounds
    ax2[1,1].axvline(11.186, color = 'k', ls = '--')
    ax2[1,1].axvspan(11.186 - 0.002, 11.186 + 0.002, color = 'k', alpha = 0.075)
    ax2[2,2].axvline(110., color = 'k', ls = '--')
    ax2[2,2].axvspan(110. - 8., 110. + 8., color = 'k', alpha = 0.075)
    # The two-parameter heatmaps
    ax2[1,0].imshow(np.max(bline, axis = (2,3)).T, aspect = 'auto', extent = (np.min(inclination), np.max(inclination), np.min(period), np.max(period)), cmap = pl.get_cmap('Greys'), origin = 'lower')
    ax2[2,0].imshow(np.max(bline, axis = (1,3)).T, aspect = 'auto', extent = (np.min(inclination), np.max(inclination), np.min(mean_longitude), np.max(mean_longitude)), cmap = pl.get_cmap('Greys'), origin = 'lower')
    ax2[2,1].imshow(np.max(bline, axis = (0,3)).T, aspect = 'auto', extent = (np.min(period), np.max(period), np.min(mean_longitude), np.max(mean_longitude)), cmap = pl.get_cmap('Greys'), origin = 'lower')
    # Tweak the appearance
    for axis in ax2.flatten():
      axis.margins(0,0)
      axis.ticklabel_format(useOffset = False)
      for tick in axis.get_xticklabels() + axis.get_yticklabels():
        tick.set_rotation(45)
    for axis in [ax2[0,1], ax2[0,2], ax2[1,2]]:
      axis.set_visible(False)
    for axis in [ax2[0,0], ax2[1,0], ax2[1,1]]:
      axis.xaxis.set_ticklabels([])
    ylims = np.array([axis.get_ylim() for axis in [ax2[0,0], ax2[1,1], ax2[2,2]]])
    ymin = np.min(ylims)
    ymax = np.max(ylims) 
    yrng = ymax - ymin
    ymin -= 0.1 * yrng
    ymax += 0.1 * yrng
    ax2[2,1].yaxis.set_ticklabels([])
    for axis in [ax2[0,0], ax2[1,1], ax2[2,2]]:
      axis.yaxis.tick_right()
      axis.set_ylim(ymin, ymax)
      axis.margins(0, None)
    ax2[0,0].set_xticks(inclination_ticks)
    ax2[1,1].set_xticks(period_ticks)
    ax2[2,2].set_xticks(mean_longitude_ticks)
    ax2[1,0].set_xticks(inclination_ticks)
    ax2[1,0].set_yticks(period_ticks)
    ax2[2,0].set_xticks(inclination_ticks)
    ax2[2,0].set_yticks(mean_longitude_ticks)
    ax2[2,1].set_xticks(period_ticks)
    ax2[2,1].set_yticks(mean_longitude_ticks)
    ax2[1,0].set_ylabel('Period (days)', labelpad = 10, fontsize = 16)
    ax2[2,0].set_ylabel('Mean longitude ($^\circ$)', labelpad = 18, fontsize = 16)
    ax2[2,0].set_xlabel('Inclination ($^\circ$)', labelpad = 21, fontsize = 18)
    ax2[2,1].set_xlabel('Period (days)', labelpad = 8, fontsize = 18)
    ax2[2,2].set_xlabel('Mean longitude ($^\circ$)', labelpad = 17, fontsize = 18)
    fig2.savefig('%s_triangle.pdf' % pref, bbox_inches = 'tight')
  
    #--- FIGURE 3: Plot bmax versus wavelength
    print("Plotting figure 3...")
    fig4, ax4 = pl.subplots(1, figsize = (12,4))
    ax4.plot(bins, bmax - 1., 'k-', lw = 0.5, zorder = -2)
    ax4.plot(line, blinemax - 1., 'ro', markeredgecolor = 'none')
    ax4.axhline(blinemax - 1., color = 'r', ls = '--')
    ax4.margins(0, None)
    ax4.yaxis.set_major_locator(MaxNLocator(nbins = 5))
    ax4.set_xlabel('Wavelength ($\AA$)', fontsize = 28)
    ax4.set_ylabel('Fractional Signal', fontsize = 28)
    [tick.label.set_fontsize(22) for tick in ax4.xaxis.get_major_ticks() + ax4.yaxis.get_major_ticks()]
    fig4.savefig('%s_max_signal_vs_wavelength.pdf' % pref, bbox_inches = 'tight') 

    # --- FIGURE 4: The distribution of signal maxima at each wavelength (this gives us the FAP)
    print("Plotting figure 4...")
    fig3, ax3 = pl.subplots(1)
    nb = len(np.where(bmax >= blinemax)[0])
    fap = nb / len(bmax)
    fap_mult, fap_exp = [int(n) for n in ("%.0e" % fap).split("e")]
    if nb == 1:
      fap_sgn = "<"
    else:
      fap_sgn = r"\approx"
    fap_str = r'$\mathrm{FAP} %s %d \times 10^{%d}$' % (fap_sgn, fap_mult, fap_exp)
    n, b, _ = ax3.hist(bmax - 1., bins = 30, color = 'w')
    d = np.digitize(blinemax - 1., b)
    if d == len(b): 
      d -= 1
    ax3.axvline(b[d] - 0.5 * (b[1] - b[0]), color = 'r', ls = '--')
    ax3.margins(0.1, None)
    ax3.set_ylabel(r'Number of signals', fontsize = 22)
    ax3.set_xlabel(r'Fractional strength', fontsize = 22)
    ax3.annotate(fap_str, xy = (0.975, 0.95), 
                 xycoords = 'axes fraction', ha = 'right', 
                 va= 'top', fontsize = 20)
    [tick.label.set_fontsize(16) for tick in ax3.xaxis.get_major_ticks() + ax3.yaxis.get_major_ticks()]
    fig3.savefig('%s_fap.pdf' % pref, bbox_inches = 'tight')

    #--- FIGURE 5: Plot the "river plot" for the best solution
    print("Plotting figure 5...")
    res = Compute(planet = planet, quiet = True, **kwargs)
    fig5 = res['fig']
    fig5.suptitle('Strongest Signal', fontsize = 30, y = 0.95)
    fig5.savefig('%s_strongest_river.pdf' % pref, bbox_inches = 'tight')

# Reproduce Figures 2-6 in the paper
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--pbs", action = 'store_true', help = 'Run on a PBS cluster')
  parser.add_argument("-l", "--line", default = Spectrum.OxygenGreen, type = float, help = 'Line wavelength (angstroms)')
  args = parser.parse_args()
  
  if args.pbs:
    PBSSearch(line = args.line)
  else:
    Search(line = args.line)