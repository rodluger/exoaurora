# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:38:40 2016

@author: dflemin3

Utility functions, defs for contrast and blackbody-real spectral
comparisons for Proxima Centauri system
"""

# Imports
from __future__ import print_function, division
import numpy as np
from scipy import interpolate
import noise_routines as nr
import coronagraph_inputs as ci

############################################
#
# Physical Constants (in cgs unless stated otherwise)
#
############################################

RSUN = 6.955e10 # radius of sun [cm]
REARTH = 6.3674447e8 # radius of earth [cm]
AUCM = 1.496e13 # cm per AU
ANGCM = 1e-8 # cm per Angstrom
PCCM = 3.086e18 # cm per parsec

# Physical constastns [cgs]
h = 6.6260755e-27 # Planck constant
c = 3.0e10 # Speed of light [cm/s]
k = 1.380658e-16 # Boltzman const

# Prox consts
TSTAR = 3100.0 # Proxima Centauri's effective temperature [K] (fit: 3042 +/- 150)
TP = 288.0 # Planet's effective temperature if Earth-like [K]
DIST = 3.991e18# Distance to proxima cen [cm]

RSTAR = 0.14*RSUN
RP = 1.07*REARTH # Prox Cen b radius if terrestrial, Sotin et al 2007 [cm]
AP = 0.0485*AUCM # Prox Cen b semimajor axis [cm]

# Auroral Lines [Angstroms]
OGREEN_LINE = 5577.345
OREDI_LINE = 6300.308
OREDII_LINE = 6363.790

# Auroral stengths (Chamberlain 1961, Appendix 2)
RUNIT = 1.0e6 # 1 Rayleigh == 1.0e6 photons/(cm^2 s)
REL_AURORA_STRENGTH = 100. #3.0*42.0 # Relative strength of Prox Cen b to Earth aurora
                           # Computed via rel_aurora_strength.py
OGREEN = 1.0e6*RUNIT # Class 4 -> 1000 KiloRayleigh signal
OREDI = 100.0*RUNIT # Note: Oxygen Red line strengths are for night airglow
OREDII = 100.0*RUNIT

# Branching ratios: fraction of total auroral power that goes into line
EPS_EARTH_OGREEN = 1.0
EPS_NEP_OGREEN = 1.0

# Noise scaling derived from comparing integration time required for a detection
# from an injected signal to photon-limited noise estimates
EPS_TEL_NOISE = 1.0

############################################
#
# Function definitions
#
############################################

def lambertPhaseFunction(alpha):
    """
    Calculate the Lambertian Phase Function from the phase angle.
    Parameters
    ----------
        alpha : float
            Planet phase angle (degrees)

    Returns
    -------
        Psi : float
            The Lambertian phase function
    """
    conv = np.pi / 180.
    return np.fabs((np.sin(alpha*conv) + (np.pi - alpha*conv) * np.cos(alpha*conv)) / np.pi)
# end function


def planck(T,lam):
    """
    For a given effective temperature, compute the blackbody intensity
    over the wavelength array.  All units in cgs!

    Parameters
    ----------
    T : float
        effective temperature [K]
    lam : float/array
        wavelength value(s) [cm]
    """
    # Takes cgs, outputs blackbody in cgs
    return (2.0*h*c*c/np.power(lam,5.))/(np.exp(h*c/(lam*k*T)) - 1.0)
# end function


def watts_to_phots(watts,lam):
    """
    Convert from auroral power in watts to photons
    at the given wavelength.

    Parameters
    ----------
    watts : float(s)
        auroral power [W]
    lam : float
        central wavelength of auroral line [microns]

    Returns
    -------
    phots : float
        photons/second
    """

    return watts*1.0e7*(lam*1.0e-4)/(h*c)


############################################
#
# Class definitions
#
############################################

class Planet(object):
    """
    Planet object that contains typical planet parameters like geometric
    albedo, radius, etc. Planet is assumed to lie on a circular orbit.
    """
    def __init__(self,A=0.3,Rp=RP,a=AP,name="Proxima Centauri b"):
        # Assumed circular orbit
        self.name = name # Planet's name
        self.A = A # Geometric albedo (assumed grey)
        self.Rp = Rp # Radius of planet [cm]
        self.a = AP # Orbital semimajor axis [cm]


    def __repr__(self):
        # For pretty prints
        return ("%s: Geometric Albedo: %.2lf. Radius (REARTH): %.2lf. Semimajor Axis (AU): %.2lf"
             % (self.name, self.A, self.Rp/REARTH, self.a/AUCM))

# end class


class Star(object):
    """
    Star object that contains typical stellar parameters (effective temperature,
    radius, etc) and spectrum.
    """
    def __init__(self,Teff=TSTAR,Rs=RSTAR,d=DIST,spectrum=None,wave=None,
                 name="Proxima Centauri"):
        # Init as Prox Cen by default
        self.name = name # Star's name
        self.Teff = Teff # Stellar effective temperature [K]
        self.Rs = Rs # Stellar Radius [cm]
        self.d = d # Earth - star distance [cm]
        self.spectrum = spectrum # Stellar flux spectrum in W/m^2/micron
                                 # self.spectrum is callable.
                                 # Arg is wavelength in microns
                                 # within range(self.wave)
        self.wave = wave # Spectrum wavelength array [microns]

        # Load in Meadows et al 2016 spectrum
        if self.spectrum is None and self.wave is None:
            spec_path = "../airglow/fits/ProxCenHubbleSpectrum.txt"
            spec = np.genfromtxt(spec_path,skip_header=25)

            # Covert from 1 AU normalization to Prox Cen Distance
            # so flux is flux received at Earth
            spec[:,1] *= (AUCM/DIST)**2

            # Linearly interpolate spectrum so it's callable
            # as a function of microns
            # Call like self.spectrum(wavelength_in_microns)
            self.wave = spec[:,0]
            self.spectrum =  interpolate.interp1d(self.wave, spec[:,1],
                                                  kind='linear')

    def __repr__(self):
        # For pretty prints
        return ("%s: Teff (k): %.2lf. Radius (RSUN): %.2lf.  Earth-Star Distance (PC): %.2lf"
             % (self.name, self.Teff, self.Rs/RSUN, self.d/PCCM))


    def flux(self,lam_0,delta_lam):
        """
        Photon flux from star at lam_0 in delta_lam window (assumed const in window)
        where flux = PI * B * (r/d)^2 an is received at Earth

        Uses linear interpolation of Meadows et al 2016 spectrum

        Parameters
        ----------
        lam_0 : float
            central wavelength in microns
        delta_lam : float
            wavelength band width in microns

        Returns
        -------
        flux : float
            Flux at that wavelength in photons/cm^2/s
        """
        conv = 1.0e7 * 1.0e-4 * 1.0e-4 # converts to cgs

        # Integrate over the band
        wav = np.linspace(lam_0-delta_lam/2.,lam_0+delta_lam/2.,100)
        flux = np.trapz(self.spectrum(wav), x=wav)
        return flux * conv * lam_0/(h*c)
#end class


class System(object):
    """
    Class to contain Planet, Star classes and routines that involve
    fluxes/properties from both bodies
    """
    # Defaults to Proxima Centauri
    def __init__(self,star=Star(),planet=Planet()):
        self.planet = planet
        self.star = star

    # Star - Planet interaction caused flux functions


    def contrast_ratio(self,alpha=0.0):
        """
        Planet-star reflected light contrast ratio following Robinson et al 2016

        Parameters
        ----------
        alpha : float
            planet's phase angle [degrees] (0 == full phase, 90 == quadrature)

        Returns
        -------
        contrast ratio : float
        """

        return self.planet.A*lambertPhaseFunction(alpha)*np.power(self.planet.Rp/self.planet.a,2.)

    def ref_phot_flux(self,alpha,lam_0, delta_lam):
        """
        Given phase (alpha), central wavelength (lam_0) and small bandwitdh (delta_lam),
        compute the planet's reflected flux in photons/(cm^2 s) received at Earth

        Parameters
        ----------
        alpha : float
            planet's phase angle [degrees] (0 == full phase, 90 == quadrature)
        lam_0 : float
           central wavelength [microns]
        delta_lam : float
            size of wavelength band [microns]

        Returns
        -------
        reflected photon flux : float
            [photons/cm^2/s]
        """

        return self.contrast_ratio(alpha) * self.star.flux(lam_0, delta_lam)

    def aurora_phot_flux(self,aurora,lam_0=OGREEN_LINE*1.0e-4):
        """
        Given an auroral strength, compute the auroral photon flux received at
        Earth from the planet.
        Assumes aurora is delta function at its central wavelength

        Parameters
        ----------
        aurora : float
           aurora strength [W]

        Returns
        -------
        aurora photon flux at Earth : float
            [photons/cm^2/s]
        """
        #return 2.0*np.pi*aurora*self.planet.aurora_area*np.power(self.planet.Rp/self.star.d,2)/(4.0*np.pi)
        # Convert to photons/s
        return watts_to_phots(aurora,lam_0) / (4.0*np.pi*DIST**2)


    def __repr__(self):
        return "Star: " + self.star.__repr__() + "\nPlanet: " + self.planet.__repr__()

# end class


class Telescope(object):
    """
    Telescope class for observing!
    """
    def __init__(self,D=10.0,eps=0.05,R=115000):
        self.D = D # Telescope diameter in m
        self.eps = eps # Net telescope efficiency (throughput)
        self.R = R # Resolving Power (lambda/delta_lambda)

    def __repr__(self):
        word = "Telescope: Diameter (m): %.2lf. Efficiency: %.2lf." % (self.D, self.eps)
        word += " Resolving Power: %.2lf." % self.R
        return word

    def calc_IWA(self, System, lam_0):
        """
        Compute the maximum inner working angle to reach a given wavelength.

        Parameters
        ----------
        System : System obj
           star-planet system to observe
       lam_0 : float
           central wavelength [microns]

       Returns
       -------
       IWA : float
           Inner working angle
        """

        # Separation in radians
        theta = (System.planet.a/System.star.d)

        return theta * (self.D * 100.) / (lam_0 * 1.0e-4)

    def cphot(self, System, aurora):
        """
        Compute the source photon count rate due to planetary aurora

        Parameters
        ----------
       System : System obj
           star-planet system to observe
       aurora : float
           auroral strength [W]

        Returns
        -------
        cphot : float
            Auroral emission photon count rates
        """
        flux = System.aurora_phot_flux(aurora)
        return self.eps * flux * np.pi * np.power(self.D*100./2.,2.)

    def cref(self,cs, System, alpha=90.):
        """
        Planet reflected light photon count rates.

        Parameters
        ----------
        cs : float
            star photon count rates
        System : System obj
            star-planet system to observe
        aurora : float
            auroral strength [W]

        Returns
        -------
        cref : float
            Reflected photon count rates
        """
        return cs * System.contrast_ratio(alpha=alpha)

    def cstar(self, System, lam_0):
        """
        Compute the stellar photon count rate

        Parameters
        ----------
       System : System obj
           star-planet system to observe
       lam_0 : float
           central wavelength [microns]

        """
        # Compute spectral element
        dl = lam_0/self.R
        flux = System.star.flux(lam_0,dl)
        return self.eps * flux * np.pi * np.power(self.D*100./2.,2.)

    def int_time(self, SN, System, aurora, lam_0, coronagraph=False):
        """
        Compute the integration time required for a photon-limited
        observation to achieve a given signal-to-noise

        Parameters
        ----------
        SN : float
            required signal-to-noise
       System : System obj
           star-planet system to observe
       aurora : float
           aurora strength [Rayleighs]
       lam_0 : float
           central auroral wavelength [microns]
       coronagraph : bool (optional)
           whether or not to null starlight

        Returns
        -------
        dt : float
            Integration time [hrs]
        """
        # Compute spectral bin
        dl = lam_0/self.R

        # Integration time in hours for given signal-to-noise
        cs = self.cstar(System, lam_0) # Star photon counts
        cp = self.cphot(System, aurora) # Planet auroral phot counts
        cref = self.cref(cs, System, alpha=90.) # Planet reflected phot counts

        # Using coronograph + coronagraph noise model (Robinson et al 2016)
        if coronagraph:

            # Scale planet counts by quantum efficiency, airy pattern,
            # and use coronagraph throughput, not telescope!
            cp *= (ci.Q * nr.f_airy(ci.X) / self.eps * ci.T)

            # Compute noise terms

            # Zodaical dust noise
            cz = nr.czodi(ci.Q,ci.X,ci.T,np.array([lam_0]),np.array([dl]),
                          self.D,ci.MZV)

            # Compute Flux of star at 1 AU in W/m^2/micron
            fstar = System.star.spectrum(lam_0) * (System.star.d/AUCM)**2

            # Exozodiacal dust noise
            cez = nr.cezodi(ci.Q,ci.X,ci.T,np.array([lam_0]),np.array([dl]),
                            self.D,System.planet.a/AUCM,
                            np.array([fstar]),ci.NEZ,ci.MEZV)

            # Speckle noise
            csp = nr.cspeck(ci.Q,ci.T,ci.C,lam_0,dl,(fstar*(AUCM/System.star.d)**2),
                            self.D)

            # Lenslet angular diameter assuming sampled at lambda/(2D)
            theta = (0.4/1.0e6)/(self.D/2.)*(180./np.pi)*3600.

            # Dark noise
            cd = nr.cdark(ci.DE,ci.X,lam_0,self.D,theta,ci.DNHPIX)

            # Read noise
            cr = nr.cread(ci.RE,ci.X,lam_0,self.D,theta,ci.DNHPIX,ci.DTMAX)

            # Thermal noise
            cth = nr.ctherm(ci.Q,ci.X,lam_0,dl,self.D,ci.TSYS,ci.EMIS)

            # Sum the noise -- factor of 2 from background subtraction
            # from telescope roll (Brown 2005)
            cn = cp + 2.0*(cz + cez + csp + cd + cr + cth)

        # No corongraph, use photon-limited telescope scheme with empircal
        # correction from HARPS 2016 data fit
        else:
            cn = cp + cs + cref

        # Integration time in hours
        return (SN**2 * cn)/(cp**2)/3600.

# end class