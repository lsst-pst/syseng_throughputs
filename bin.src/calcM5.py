# Calculate m5 and table2 values using the SYSENG_THROUGHPUTS files. Saves figures to plots directory automatically.
from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lsst.sims.photUtils import Bandpass, Sed, SignalToNoise
from lsst.sims.photUtils import PhotometricParameters, LSSTdefaults
import bandpassUtils as bu
from lsst.utils import getPackageDir

filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
filtercolors = {'u':'b', 'g':'c', 'r':'g',
                'i':'y', 'z':'r', 'y':'m'}

plotformat = 'png'

def calcM5(hardware, system, atmos, title='m5', X=1.0, return_t2_values=False):
    """
    Calculate m5 values for all filters in hardware and system.
    Prints all values that go into "table 2" of the overview paper.
    Returns dictionary of m5 values.
    """
    # photParams stores default values for the exposure time, nexp, size of the primary,
    #  readnoise, gain, platescale, etc.
    # See https://github.com/lsst/sims_photUtils/blob/master/python/lsst/sims/photUtils/PhotometricParameters.py
    effarea = np.pi * (6.423/2.*100.)**2
    photParams_zp = PhotometricParameters(exptime=1, nexp=1, gain=1, effarea=effarea,
                                          readnoise=8.8, othernoise=0, darkcurrent=0.2)
    photParams = PhotometricParameters(gain=1.0, effarea=effarea, readnoise=8.8, othernoise=0, darkcurrent=0.2)
    photParams_infinity = PhotometricParameters(gain=1.0, readnoise=0, darkcurrent=0,
                                                othernoise=0, effarea=effarea)
    # lsstDefaults stores default values for the FWHMeff.
    # See https://github.com/lsst/sims_photUtils/blob/master/python/lsst/sims/photUtils/LSSTdefaults.py
    lsstDefaults = LSSTdefaults()
    darksky = Sed()
    darksky.readSED_flambda(os.path.join(getPackageDir('syseng_throughputs'),
                                         'siteProperties', 'darksky.dat'))
    flatSed = Sed()
    flatSed.setFlatSED()
    m5 = {}
    Tb = {}
    Sb = {}
    kAtm = {}
    Cm = {}
    dCm_infinity = {}
    sourceCounts = {}
    skyCounts = {}
    skyMag = {}
    gamma = {}
    zpT = {}
    FWHMgeom = {}
    FWHMeff = {}
    for f in system:
        zpT[f] = system[f].calcZP_t(photParams_zp)
        FWHMeff[f] = np.power(X, 0.6) * lsstDefaults.FWHMeff(f)
        m5[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams, FWHMeff=FWHMeff[f])
        fNorm = flatSed.calcFluxNorm(m5[f], system[f])
        flatSed.multiplyFluxNorm(fNorm)
        sourceCounts[f] = flatSed.calcADU(system[f], photParams=photParams)
        # Calculate the Skycounts expected in this bandpass.
        skyCounts[f] = (darksky.calcADU(hardware[f], photParams=photParams)
                        * photParams.platescale**2)
        # Calculate the sky surface brightness.
        skyMag[f] = darksky.calcMag(hardware[f])
        # Calculate the gamma value.
        gamma[f] = SignalToNoise.calcGamma(system[f], m5[f], photParams)
        # Calculate the "Throughput Integral" (this is the hardware + atmosphere)
        dwavelen = np.mean(np.diff(system[f].wavelen))
        Tb[f] = np.sum(system[f].sb / system[f].wavelen) * dwavelen
        # Calculate the "Sigma" 'system integral' (this is the hardware only)
        Sb[f] = np.sum(hardware[f].sb / hardware[f].wavelen) * dwavelen
        # Calculate km - atmospheric extinction in a particular bandpass
        kAtm[f] = -2.5*np.log10(Tb[f] / Sb[f])
        # Calculate the Cm and Cm_Infinity values.
        # m5 = Cm + 0.5*(msky - 21) + 2.5log10(0.7/FWHMeff) + 1.25log10(t/30) - km(X-1.0)
        # Assumes atmosphere used in system throughput is X=1.0
        Cm[f] = (m5[f] - 0.5*(skyMag[f] - 21) - 2.5*np.log10(0.7/FWHMeff[f])
                 - 1.25*np.log10((photParams.exptime*photParams.nexp)/30.0) + kAtm[f]*(X-1.0))
        # Calculate Cm_Infinity by setting readout noise to zero.
        m5inf = SignalToNoise.calcM5(darksky, system[f], hardware[f],  photParams_infinity,
                                     FWHMeff=FWHMeff[f])
        Cm_infinity = (m5inf - 0.5*(skyMag[f] - 21) - 2.5*np.log10(0.7/FWHMeff[f])
                       - 1.25*np.log10((photParams.exptime*photParams.nexp)/30.0) + kAtm[f]*(X-1.0))
        dCm_infinity[f] = Cm_infinity - Cm[f]
    print('Filter FWHMeff FWHMgeom SkyMag SkyCounts Zp_t Tb Sb kAtm Gamma Cm dCm_infinity m5 SourceCounts')
    for f in ('u', 'g' ,'r', 'i', 'z', 'y'):
        FWHMgeom[f] = SignalToNoise.FWHMeff2FWHMgeom(FWHMeff[f])
        print('%s %.2f %.2f %.2f %.1f %.2f %.3f %.3f %.4f %.6f %.2f %.2f %.2f %.2f'\
           % (f, FWHMeff[f], FWHMgeom[f],
              skyMag[f], skyCounts[f], zpT[f], Tb[f], Sb[f], kAtm[f],
              gamma[f], Cm[f], dCm_infinity[f], m5[f], sourceCounts[f]))

    for f in filterlist:
        m5_cm = Cm[f] + 0.5*(skyMag[f] - 21.0) + 2.5*np.log10(0.7/FWHMeff[f]) - kAtm[f]*(X-1.0)
        if m5_cm - m5[f] > 0.001:
            raise ValueError('Cm calculation for %s band is incorrect! m5_cm != m5_snr' %f)

    if return_t2_values:
        return {'FWHMeff': FWHMeff, 'FWHMgeom': FWHMgeom, 'skyMag': skyMag, 'skycounts': skyCounts,
                'zpT': zpT, 'Tb': Tb, 'Sb': Sb, 'kAtm': kAtm,
                'gamma': gamma, 'Cm': Cm, 'dCm_infinity': dCm_infinity,
                'm5': m5, 'sourceCounts': sourceCounts}

    # Show what these look like individually (add sky & m5 limits on throughput curves)
    plt.figure()
    for f in filterlist:
        plt.plot(system[f].wavelen, system[f].sb, color=filtercolors[f], linewidth=2, label=f)
    plt.plot(atmos.wavelen, atmos.sb, 'k:', label='X=%.1f' % X)
    plt.legend(loc='center right', fontsize='smaller')
    plt.xlim(300, 1100)
    plt.ylim(0, 1)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.title('System Throughputs')
    plt.grid(True)
    plt.savefig('../plots/throughputs%s.%s' % (title, plotformat), format=plotformat)

    plt.figure()
    ax = plt.gca()
    # Add dark sky
    ax2 = ax.twinx()
    plt.sca(ax2)
    skyab = np.zeros(len(darksky.fnu))
    condition = np.where(darksky.fnu > 0)
    skyab[condition] = -2.5*np.log10(darksky.fnu[condition]) - darksky.zp
    ax2.plot(darksky.wavelen, skyab,
             'k-', linewidth=0.8, label='Dark sky mags')
    ax2.set_ylabel('AB mags')
    ax2.set_ylim(24, 14)
    plt.sca(ax)
    # end of dark sky
    handles = []
    for f in filterlist:
        plt.plot(system[f].wavelen, system[f].sb, color=filtercolors[f], linewidth=2)
        myline = mlines.Line2D([], [], color=filtercolors[f], linestyle='-', linewidth=2,
                               label = '%s: m5 %.1f (sky %.1f)' %(f, m5[f], skyMag[f]))
        handles.append(myline)
    plt.plot(atmos.wavelen, atmos.sb, 'k:', label='Atmosphere, X=1.0')
    # Add legend for dark sky.
    myline = mlines.Line2D([], [], color='k', linestyle='-', label='Dark sky AB mags/arcsec^2')
    handles.append(myline)
    # end of dark sky legend line
    plt.legend(loc=(0.01, 0.6), handles=handles, fancybox=True, numpoints=1, fontsize='small')
    plt.ylim(0, 1)
    plt.xlim(300, 1100)
    plt.grid(True, alpha=0.5)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')
    plt.title('System total response curves %s' %(title))
    plt.savefig('../plots/system+sky%s.%s' % (title, plotformat), format=plotformat, dpi=600)
    return m5


def get_effwavelens(system_bandpasses, filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Calculate the effective wavelengths.

    Parameters
    ----------
    system_bandpasses: dict of ~lsst.sims.photUtils.Bandpass
        The bandpasses for the system.

    Returns
    -------
    numpy.ndarray
        The effective wavelength values, in a numpy array.
        Order matches filterlist.
    """
    eff_wavelen = np.zeros(len(filterlist), float)
    for i, f in enumerate(filterlist):
        eff_wavelen[i] = system_bandpasses[f].calcEffWavelen()[1]
    return eff_wavelen

def scale_seeing(fwhm_500, eff_wavelen, airmass=1.0):
    """Calculate the FWHMeff and FWHMgeom given FWHM_500.

    Parameters
    ----------
    fwhm_500 : float
        The FWHM_500 (FWHM at 500nm at zenith). Fiducial values is 0.6".
    eff_wavelen : numpy.ndarray
        The effective wavelengths of the system bandpasses.
        Can be calculated using get_effwavelens.
    airmass : float or numpy.ndarray
        The airmass at which to calculate the FWHMeff and FWHMgeom values.
        Default 1.0.
        Can be a single value or a numpy array.

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        FWHMeff, FWHMgeom: both are the same shape numpy.ndarray.
        If airmass is a single value, FWHMeff & FWHMgeom are 1-d arrays,
        with the same order as eff_wavelen (i.e. eff_wavelen[0] = u, then FWHMeff[0] = u).
        If airmass is a numpy array, FWHMeff and FWHMgeom are 2-d arrays,
        in the order of <filter><airmass> (i.e. eff_wavelen[0] = u, 1-d array over airmass range).
    """
    airmass_correction = np.power(airmass, 0.6)
    wavelen_correction = np.power(500.0 / eff_wavelen, 0.3)
    # telSeeing = 0.25”, opticalDesign = 0.08”, and cameraSeeing = 0.30”.
    # Document-20160: total system contribution at zenith = 0.4" (add above, in quadrature)
    if isinstance(airmass, np.ndarray):
        fwhm_system = 0.4 * np.outer(np.ones(len(wavelen_correction)), airmass_correction)
        fwhm_atmo = fwhm_500 * np.outer(wavelen_correction, airmass_correction)
    else:
        fwhm_system = 0.4 * airmass_correction
        fwhm_atmo = fwhm_500 * wavelen_correction * airmass_correction
    # Calculate combined FWHMeff.
    fwhm_eff = 1.16 * np.sqrt(fwhm_system**2 + 1.04 * fwhm_atmo**2)
    # Translate to FWHMgeom.
    fwhm_geom = 0.822 * fwhm_eff + 0.052
    return fwhm_eff, fwhm_geom



if __name__ == '__main__':

    # Set the directories for each component.
    # Note that this sets the detector to be the 'generic detector' (minimum of each vendor).
    defaultDirs = bu.setDefaultDirs(rootDir = '..')
    # To use a particular vendor, uncomment one of the following lines or edit as necessary.
    # defaultDirs['detector'].replace('joint_mininimum', 'itl')
    # defaultDirs['detector'].replace('joint_minimum', 'e2v')

    # Add losses to each component?
    addLosses = True

    # Build the system and hardware throughput curves (with aerosols, with X=1.2).
    atmosphere = bu.readAtmosphere(defaultDirs['atmosphere'], atmosFile='pachonModtranAtm_12_aerosol.dat')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses, atmosphereOverride=atmosphere)
    m5 = calcM5(hardware, system, atmosphere, title='X=1.2', X=1.2)

    atmosphere = bu.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses, atmosphereOverride=atmosphere)
    m5 = calcM5(hardware, system, atmosphere, title='')


    plt.show()
