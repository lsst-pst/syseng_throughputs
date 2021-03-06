import os
import numpy as np
import pandas as pd
from lsst.sims.photUtils import Bandpass, Sed, SignalToNoise
from lsst.sims.photUtils import PhotometricParameters, LSSTdefaults
from lsst.utils import getPackageDir

filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
filtercolors = {'u':'b', 'g':'c', 'r':'g',
                'i':'y', 'z':'r', 'y':'m'}

# Fiducial M5 values from the SRD
m5_fid = {'u': 23.9, 'g': 25.0, 'r': 24.7, 'i': 24.0, 'z': 23.3, 'y': 22.1}
m5_min = {'u': 23.4, 'g': 24.6, 'r': 24.3, 'i': 23.6, 'z': 22.9, 'y': 21.7}


def get_effwavelens(system_bandpasses, filterlist):
    """Calculate the effective wavelengths.

    Parameters
    ----------
    system_bandpasses: dict of ~lsst.sims.photUtils.Bandpass

    Returns
    -------
    numpy.ndarray of floats
        The effective wavelength values, in a numpy array, in order of filterlist.
    """
    eff_wavelen = np.zeros(len(filterlist), float)
    for i, f in enumerate(filterlist):
        eff_wavelen[i] = system_bandpasses[f].calcEffWavelen()[1]
    return eff_wavelen


def makeM5(hardware, system, darksky=None, exptime=15, nexp=2,
           readnoise=8.8, othernoise=0, darkcurrent=0.2,
           effarea=np.pi*(6.423/2*100)**2, X=1.0):
    """Calculate values which are related to m5 (basically 'table2' of overview paper).

    Parameters
    ----------
    hardware : dict of ~lsst.sims.photUtils.Bandpass
        The bandpasses for the hardware only.
    system : dict of ~lsst.sims.photUtils.Bandpass
        The bandpasses for the total system (hardware + atmosphere)
    darksky : ~lsst.sims.photUtils.Sed or None
        The Sed of the dark night sky.
        Default None, in which case it will be read from siteProperties/darksky.dat
    exptime : float, opt
        The open-shutter exposure time for one exposure (seconds). Default 15s.
    nexp : int, opt
        The number of exposures in one visit. Default 2.
    readnoise : float, opt
        The readnoise for one exposure (electrons). Default 8.8 e-
    othernoise : float, opt
        Additional noise to be handled like readnoise (electrons). Default 0.
    darkcurrent : float, opt
        Dark current (electrons/second). Default 0.2 e-/s.
    effarea : float, opt
        Effective area of the primary mirror (cm^2). Default pi*(6.432/2 m)**2
    X : float, opt
        The airmass for the system bandpasses. Used to modify FWHMeff and Cm values.
         Default 1.0.

    Returns
    -------
    pd.DataFrame
    """
    # PhotometricParameters object to calculate telescope zeropoint (1s exposure).
    photParams_zp = PhotometricParameters(exptime=1, nexp=1, gain=1, effarea=effarea,
                                          readnoise=readnoise, othernoise=othernoise,
                                          darkcurrent=darkcurrent)
    # PhotometricParameters object for "real" visit.
    photParams_std = PhotometricParameters(exptime=exptime, nexp=nexp,
                                           gain=1.0, effarea=effarea, readnoise=readnoise,
                                           othernoise=othernoise, darkcurrent=darkcurrent)
    photParams_double = PhotometricParameters(exptime=2*exptime, nexp=nexp,
                                              gain=1.0, effarea=effarea, readnoise=readnoise,
                                              othernoise=othernoise, darkcurrent=darkcurrent)
    # PhotometricParameters object for "no noise" visit.
    photParams_infinity = PhotometricParameters(exptime=exptime, nexp=nexp,
                                                gain=1.0, readnoise=0, darkcurrent=0,
                                                othernoise=0, effarea=effarea)
    # lsstDefaults stores default values for the FWHMeff.
    # See https://github.com/lsst/sims_photUtils/blob/master/python/lsst/sims/photUtils/LSSTdefaults.py
    lsstDefaults = LSSTdefaults()
    # Set up dark sky and flat seds.
    if darksky is None:
        darksky = Sed()
        darksky.readSED_flambda(os.path.join(getPackageDir('syseng_throughputs'),
                                             'siteProperties', 'darksky.dat'))
    flatSed = Sed()
    flatSed.setFlatSED()

    # Now set up dataframe. filters x properties.
    properties = ['FWHMeff', 'FWHMgeom', 'skyMag', 'skyCounts', 'Zp_t',
                  'Tb', 'Sb', 'kAtm', 'gamma', 'Cm', 'dCm_infinity', 'dCm_double', 'm5', 'sourceCounts',
                  'm5_fid', 'm5_min']

    # Set filterlist - usually this will be just 'filterlist' above, but could have added extras
    for f in system:
        if f not in filterlist:
            filterlist.append(f)
    # Now set up dataframe to hold outputs
    d = pd.DataFrame(index=filterlist, columns=properties, dtype='float')

    # Calculate effective wavelengths
    eff_wavelengths = get_effwavelens(system, filterlist)
    # Calculate the FWHM at each wavelength, for this airmass.
    fwhm_eff_zenith = np.zeros(len(filterlist), float)
    for i, f in enumerate(filterlist[0:7]):
        # The original LSST filters
        fwhm_eff_zenith[i] = lsstDefaults.FWHMeff(f)
    if len(filterlist) > 6:
        for i, f in enumerate(filterlist[7:]):
            fwhm_eff_zenith[i] = np.interp(system[f].calcEffWavelen()[1],
                                            eff_wavelengths[0:7],
                                            fwhm_eff_zenith[0:7])
    fwhm_eff = {}
    for i, f in enumerate(filterlist):
        fwhm_eff[f] = np.power(X, 0.6) * fwhm_eff_zenith[i]

    for f in system:
        # add any missing m5 fiducials and mininum values - no requirements on non-standard bands
        if f in m5_min.keys():
            d.m5_fid.loc[f] = m5_fid[f]
            d.m5_min.loc[f] = m5_min[f]
        else:
            d.m5_fid.loc[f] = -666
            d.m5_min.loc[f] = -666
        d.Zp_t.loc[f] = system[f].calcZP_t(photParams_zp)
        d.FWHMeff.loc[f] = fwhm_eff[f]
        d.FWHMgeom.loc[f] = 0.822 * d.FWHMeff.loc[f] + 0.052
        d.m5.loc[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f],
                                           photParams_std, FWHMeff=d.FWHMeff.loc[f])
        fNorm = flatSed.calcFluxNorm(d.m5.loc[f], system[f])
        flatSed.multiplyFluxNorm(fNorm)
        d.sourceCounts.loc[f] = flatSed.calcADU(system[f], photParams=photParams_std)
        # Calculate the Skycounts expected in this bandpass.
        d.skyCounts.loc[f] = (darksky.calcADU(hardware[f], photParams=photParams_std)
                              * photParams_std.platescale**2)
        # Calculate the sky surface brightness.
        d.skyMag.loc[f] = darksky.calcMag(hardware[f])
        # Calculate the gamma value.
        d.gamma.loc[f] = SignalToNoise.calcGamma(system[f], d.m5.loc[f], photParams_std)
        # Calculate the "Throughput Integral" (this is the hardware + atmosphere)
        dwavelen = np.mean(np.diff(system[f].wavelen))
        d.Tb.loc[f] = np.sum(system[f].sb / system[f].wavelen) * dwavelen
        # Calculate the "Sigma" 'system integral' (this is the hardware only)
        d.Sb.loc[f] = np.sum(hardware[f].sb / hardware[f].wavelen) * dwavelen
        # Calculate km - atmospheric extinction in a particular bandpass
        d.kAtm.loc[f] = -2.5 * np.log10(d.Tb.loc[f] / d.Sb.loc[f])
        # Calculate the Cm and Cm_Infinity values.
        # m5 = Cm + 0.5*(msky - 21) + 2.5log10(0.7/FWHMeff) + 1.25log10(t/30) - km(X-1.0)
        # Assumes atmosphere used in system throughput is X=1.0
        d.Cm.loc[f] = (d.m5.loc[f] - 0.5 * (d.skyMag.loc[f] - 21) - 2.5 * np.log10(0.7 / d.FWHMeff.loc[f])
                       - 1.25 * np.log10((photParams_std.exptime * photParams_std.nexp) / 30.0)
                       + d.kAtm.loc[f] * (X - 1.0))
        # Calculate Cm_Infinity by setting readout noise to zero.
        m5inf = SignalToNoise.calcM5(darksky, system[f], hardware[f],  photParams_infinity,
                                     FWHMeff=d.FWHMeff.loc[f])
        Cm_infinity = (m5inf - 0.5 * (d.skyMag.loc[f] - 21) - 2.5 * np.log10(0.7 / d.FWHMeff.loc[f])
                       - 1.25 * np.log10((photParams_infinity.exptime * photParams_infinity.nexp) / 30.0)
                       + d.kAtm.loc[f] * (X - 1.0))
        d.dCm_infinity.loc[f] = Cm_infinity - d.Cm.loc[f]
        m5double = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams_double,
                                        FWHMeff=d.FWHMeff.loc[f])
        Cm_double = (m5double - 0.5 * (d.skyMag.loc[f] - 21) - 2.5 * np.log10(0.7 / d.FWHMeff.loc[f])
                     - 1.25 * np.log10(photParams_double.exptime * photParams_double.nexp / 30.0)
                     + d.kAtm.loc[f] * (X - 1.0))
        d.dCm_double.loc[f] = Cm_infinity - Cm_double

    m5_cm = (d.Cm + 0.5*(d.skyMag - 21.0) + 2.5*np.log10(0.7/d.FWHMeff) - d.kAtm*(X-1.0)
             + 1.25 * np.log10((photParams_infinity.exptime * photParams_infinity.nexp) / 30.0))
    if np.any(m5_cm - d.m5 > 0.001):
        raise ValueError('m5 from Cm does not match m5 from photUtils.')
    return d
