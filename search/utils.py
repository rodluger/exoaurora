#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
utils.py
--------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib as mpl
mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Times New Roman']
from .kepler import RadialVelocity
import numpy as np
from sklearn.decomposition import PCA
from wpca import WPCA
from scipy.stats import norm, tukeylambda
from scipy.stats import t as student_t
from scipy.signal import savgol_filter, medfilt
import glob
import os
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator, ScalarFormatter, FormatStrFormatter
from matplotlib.patches import Rectangle
import matplotlib.mlab as mlab
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
  
  EarthLines = [MysteryLineII, OxygenGreen, OxygenRedI, OxygenRedII]
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
    
    WARNING: In an earlier version of this code, there
    was a minus sign preceding np.log10(...) in the output.
    I changed this to match the detection paper RV figure.
    Check this!
    
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
    
    self._period = lambda: np.random.normal(11.186, 0.002)
    self._stellar_mass = lambda: np.random.normal(0.120, 0.015)
    self._eccentricity = lambda: np.random.rand() * 0.35
    self._mean_longitude = lambda: np.random.normal(110., 8.)
    self._arg_periastron = lambda: np.random.rand() * 360.
    self._inclination = lambda: np.arccos(np.random.rand()) * 180. / np.pi
    self._mass = lambda: np.random.normal(1.27, 0.18) / np.sin(np.pi / 180. * self.inclination)
  
  def randomize(self):
    '''
    
    '''
    
    self.period = self._period()
    self.stellar_mass = self._stellar_mass()
    self.eccentricity = self._eccentricity()
    self.mean_longitude = self._mean_longitude()
    self.arg_periastron = self._arg_periastron()
    self.inclination = self._inclination()
    # Our method doesn't handle low inclinations well,
    # so we won't allow them
    while self.inclination < 20.:
      self.inclination = self._inclination()
    self.mass = self._mass()
 
def GaussianLine(x, line = Spectrum.OxygenGreen, fwhm = 0.1, A = 1.):
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
  for dataset in ['HARPS2016', 'HARPS']:
    
    files = glob.glob(os.path.join(os.path.dirname(__file__), 'fits', dataset, 'ADP*fits'))
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
  
  filename = os.path.join(os.path.dirname(__file__), 'fits', 'data.npz')
  
  if not os.path.exists(filename):
  
    # Grab the data from the FITS files
    print("Loading FITS files...")
    jd, wavs, flx, err, dataset, earth_rv, exptime, air = ReadFITS()
    
    # Bin all spectra to the first wavelength grid (regularized)
    wav = np.linspace(wavs[0][0], wavs[0][-1], len(wavs[0]))
    flx = np.array([(np.interp(wav, wavs[i], flx[i])) for i in range(len(wavs))])
    
    # Apply the system offset. This was determined empirically by matching the
    # stellar Na D I and II lines. This calibration is accurate to ~ 0.01 A.
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
            plot_sz = 10, npc = 25, frame = 'planet', wpca_sz = 250, wpca_mask_sz = 0.2, 
            bin_sz = 0.1, mask_star = True, med_mask_sz = 5,
            inject_contrast = 0, inject_planet = ProxCenB(), airmass_correction = False,
            high_pass_filter = True, crop_outliers = True, stack_lims = None, fap_iter = 0,
            clobber = False, quiet = False, plot_gaussian_fit = False, filter = 'True',
            load_fap = True):
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
  
  # Crop to a smaller window centered on the line
  inds = np.where((wav >= line - wpca_sz / 2.) & (wav <= line + wpca_sz / 2.))[0]
  wav = wav[inds]
  flx = flx[:,inds]

  # Inject a signal for injection/recovery tests?
  if inject_contrast > 0.:
    print("Injecting signal...")
    for i in range(len(flx)):
      # Get planet line wavelength in the stellar frame
      x = ApplyDopplerFactor(line, inject_planet.doppler_factor(jd[i]))
      # Get indices corresponding to the FWHM of the line
      inds = np.where(np.abs(wav - x) < 0.1)[0]
      # Compute the flux in the stellar spectrum at those indices
      z = np.trapz(flx[i][inds], wav[inds])
      # Compute the line profile. The amplitude is contrast * z
      flx[i] += GaussianLine(wav, line = x, A = inject_contrast * z)

  # Get masks in the stellar frame
  if not quiet: print("Applying masks...")
  weights = np.ones((len(wav), len(flx)))
  earth = [[] for i in range(len(flx))]
  star = [[] for i in range(len(flx))]
  for i in range(len(flx)):
  
    # Set weights according to exposure time
    weights[:,i] = np.sqrt(exp[i])
  
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
  
  # Crop outliers in planet frame?
  if crop_outliers:

    # Shift to planet frame. If the planet is moving away from us (redshifted),
    # we must BLUESHIFT the wavelength array to move into the planet frame,
    # hence the NEGATIVE sign below.
    x = [ApplyDopplerFactor(wav, -p) for p in pdf]
    flxp = np.array([np.interp(wav, xi, flxi) for xi, flxi in zip(x, flx)])
  
    # Remove 10-sigma outliers in each wavelength bin
    med = np.nanmedian(flxp, axis = 1)
    for j, _ in enumerate(wav):
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
  fstack_fap = np.empty((fap_iter, len(xstack)))
  planet_fap = ProxCenB()
  
  # The figure
  if plot:
  
    fig = pl.figure(figsize = (10,10))
    ax = [None, None, None, None]
    ax[0] = pl.subplot2grid((8, 8), (0, 0), colspan=8, rowspan=6)
    ax[1] = pl.subplot2grid((8, 8), (6, 0), colspan=8, rowspan=1, sharex = ax[0])
    ax[2] = pl.subplot2grid((8, 8), (7, 0), colspan=8, rowspan=1, sharex = ax[0])
    
  # Loop over spectra
  if not quiet: print("Stacking spectra...")
  pcb = [[] for i in range(len(flx))]
  for i in range(len(flx)):
    
    # User-defined filter
    if not eval(filter):
      continue
    
    # Compute the planet mask in the stellar frame. If the planet is moving away
    # from us (redshifted), we REDSHIFT the line accordingly.
    x = ApplyDopplerFactor(line, pdf[i])
    pcb[i] = np.where((wav > x - wpca_mask_sz / 2.) & (wav < x + wpca_mask_sz / 2.))[0]
    
    # Plot and stack only spectra where there's no overlap between planet and Earth / planet and star
    if len(set(pcb[i]) & set([item for sublist in earth[i] for item in sublist])) + len(set(pcb[i]) & set([item for sublist in star[i] for item in sublist])) == 0:
      
      # The de-trended spectra
      x = np.array(wav)
      y = planet.phase(jd[i]) + 0.01 * flx[i] / np.nanmedian(flx[i])
      y0 = np.array(flx[i])
      
      # DOPPLER SHIFTING SANITY CHECK:
      #
      #   if (planet.phase(jd[i]) < 0.22) and (planet.phase(jd[i]) > 0.18):
      #     N = int((2457400 - jd[i]) / planet.period) + 1.
      #     print(jd[i] + N * planet.period)
      # 
      # At a phase of 0.2, the figure in the paper shows that the planet is
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
      
      # High pass filter? 5 angstroms wide, 2nd order SavGol
      if high_pass_filter:
        window = int(5. / np.nanmedian(x[1:] - x[:-1]))
        if not (window % 2):
          window += 1
        filt = savgol_filter(y0, window, 2)
        y0 = y0 - filt + np.nanmedian(y0)
    
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
      
      # Stack to compute the false alarm probability
      # NOTE: This is a bit hacky, but we purposefully Doppler-shift the spectra in
      # the **opposite** direction here (no negative sign preceding `pdf[i]` below)
      # to prevent stacking up features from the actual planet when computing the FAP.
      for f in range(fap_iter):
        planet_fap.randomize()
        x_fap = ApplyDopplerFactor(x, pdf[i])
        fstack_fap[f] += np.interp(xstack, x_fap, y0)
    
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
    
  # Generate a histogram to compute the false alarm probability
  # We'll save these to disk and add to the array each time we run the script
  fap_file = os.path.join(os.path.dirname(__file__), 'fits', 'fap.npz')
  bflx_fap = np.array(bflx)
  if load_fap:
    bflx_fap = np.append(bflx_fap, np.load(fap_file)['bflx_fap'])
  for f in range(fap_iter):
    bflxf = np.zeros_like(bin_centers)
    bnrm = np.zeros_like(bin_centers)
    for x, z in zip(xstack, fstack_fap[f] / np.nanmedian(fstack_fap[f])):
      if (x >= bin_centers[0] - bin_sz / 2.) and (x < bin_centers[-1] + bin_sz / 2.):
        i = np.argmin(np.abs(bin_centers - x))
        bflxf[i] += z
        bnrm[i] += 1
    bflxf /= bnrm
    bflxf = bflxf[pad:-pad]
    bflx_fap = np.append(bflx_fap, bflxf)
  if inject_contrast == 0:
    np.savez(fap_file, bflx_fap = bflx_fap)
    
  # The detection significance
  std = np.nanstd(bflx_fap)
  sig = (bflx[len(bflx) // 2] - 1) / std
    
  # Print some info
  if not quiet:
    print('Total integration time: %.1f hours' % (total_exp / 3600.))
    print('Number of samples in FAP calculation: %d' % len(bflx_fap))
    print('Significance (sigma): %.2f' % sig)
    
  # Compute the false alarm probability: number of bins with at least the 
  # significance of the detection, divided by total number of bins
  nb = len(np.where(np.abs(bflx_fap - 1) / std >= sig)[0])
  fap = nb / len(bflx_fap)
  
  if (fap < 1e-4) and (sig > 0):
    fap_mult, fap_exp = [int(n) for n in ("%.0e" % fap).split("e")]
    if nb == 1:
      fap_sgn = "<"
    else:
      fap_sgn = "="
    fap_str = r'$\mathrm{FAP} %s %d \times 10^{%d}$' % (fap_sgn, fap_mult, fap_exp)
  else:
    fap_str = r'$\mathrm{No\ detection}$'
  
  # The histogram inset
  if plot:
    ax[3] = pl.axes([0.175 + 0.01, 0.36 + 0.01, 0.2, 0.15], zorder = 2)
    n, bins, _ = ax[3].hist((bflx_fap - 1.) / std, bins = 30, color = 'w')
    d = np.digitize(sig, bins)
    if d == len(bins): 
      d -= 1
      
    if plot_gaussian_fit:
      df = np.linspace(-sig, sig, 1000)
      mu, sigma = norm.fit((bflx - 1) / std)
      N = norm.pdf(df, mu, sigma)
      ax[3].plot(df, N * np.max(n) / np.max(N), 'b-')
    
    ax[3].axvline(bins[d] - 0.5 * (bins[1] - bins[0]), color = 'r', ls = '--')
    ax[3].set_yscale('log')
    ax[3].set_ylim(0.5, 2 * np.max(n))
    ax[3].set_yticks([1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    ax[3].set_yticklabels(['0', '1', '2', '3', '4', '5', '6'])
    ax[3].margins(0.1, None)
    ax[3].set_ylabel(r'$\log\ N$', fontsize = 14)
    ax[3].set_xlabel(r'$\Delta f\ (\sigma$)', fontsize = 14)
    ax[3].set_title(fap_str, fontsize = 16)
    rect = Rectangle((0.125 + 0.01, 0.36 - 0.055 + 0.01), 0.26, 0.245, facecolor='w', edgecolor='k',
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
    return fig, ax
  else:
    return sig

def Search(orb_iter = 50, inclination = np.arange(20., 90., 1.), **kwargs):
  '''
  
  '''
  
  # Get the data
  data = GetData()
  planet = ProxCenB()
  kwargs.update({'data': data, 'quiet': True, 'load_fap': False})
  
  # Loop over planet params
  sig = np.zeros_like(inclination)
  par = [None for i in range(len(sig))]
  for i, inc in enumerate(inclination):
    print("Iteration %d/%d..." % (i + 1, len(inclination)))
    s = []
    p = []
    planet.inclination = inc
    planet.mass = 1.27 / np.sin(planet.inclination * np.pi / 180)

    if orb_iter > 1:
      for j in range(orb_iter): 
        planet.period = np.random.normal(11.186, 0.002)
        planet.stellar_mass = np.random.normal(0.120, 0.015)
        planet.mean_longitude = np.random.normal(110., 8.)
        s.append(Compute(planet = planet, plot = False, **kwargs))
        p.append([planet.inclination, planet.mass, planet.period, planet.stellar_mass, planet.mean_longitude])
    else:
      s.append(Compute(planet = planet, plot = False, **kwargs))
      p.append([planet.inclination, planet.mass, planet.period, planet.stellar_mass, planet.mean_longitude])
    sig[i] = np.nanmax(s)
    par[i] = p[np.nanargmax(s)]

  # Plot the inclination - signal strength plot
  fig = pl.figure()
  pl.plot(inclination, sig, 'bo')
  pl.plot(inclination, sig, 'b-', alpha = 0.2)
  pl.xlabel('Inclination (degrees)')
  pl.ylabel('Max signal strength (sigma)')
  fig.savefig('inclination.pdf')
  
  # Plot best solution
  planet.inclination, planet.mass, planet.period, planet.stellar_mass, planet.mean_longitude = par[np.argmax(sig)]
  print("")
  print("BEST SOLUTION")
  print("-------------")
  print("Signal amp:  %.2f" % sig[np.argmax(sig)])
  print("Inclination: %.1f" % planet.inclination)
  print("Planet mass: %.2f" % planet.mass)
  print("Period:      %.5f" % planet.period)
  print("Star mass:   %.2f" % planet.stellar_mass)
  print("Mean long.:  %.2f" % planet.mean_longitude)
  fig2, ax = Compute(planet = planet, stack_lims = [(0.9840, 1.0160), (0.9960, 1.0090)], **kwargs)
  fig2.suptitle('Real', fontsize = 30, y = 0.95)
  fig2.savefig('real.pdf', bbox_inches = 'tight')

def Plot(figname = 'plot.pdf', suptitle = None, planet = ProxCenB(), **kwargs):
  '''
  
  '''
  
  # Get the data
  data = GetData()
  fig, ax = Compute(planet = planet, data = data, 
                            stack_lims = [(0.9840, 1.0160), (0.9960, 1.0090)], **kwargs)
  if suptitle is not None:
    fig.suptitle(suptitle, fontsize = 30, y = 0.95)
  fig.savefig(figname, bbox_inches = 'tight')