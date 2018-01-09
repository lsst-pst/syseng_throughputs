from __future__ import print_function
# Calculate m5 and table2 values using the *throughputs* repo files. Does NOT save figures.
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lsst.sims.photUtils import Bandpass, Sed, SignalToNoise
from lsst.sims.photUtils import PhotometricParameters, LSSTdefaults
import bandpassUtils as bu


filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
filtercolors = {'u':'b', 'g':'c', 'r':'g',
                'i':'y', 'z':'r', 'y':'m'}

def calcM5(hardware, system, atmos, title='m5'):
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
    darksky.readSED_flambda(os.path.join('../siteProperties', 'darksky.dat'))
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
    for f in system:
        m5[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams,
                                     FWHMeff=lsstDefaults.FWHMeff(f))
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
        # Exptime should be 30 seconds and X=1.0
        exptime = photParams.exptime * photParams.nexp
        if exptime != 30.0:
            print("Whoa, exposure time was not as expected - got %s not 30 seconds. Please edit Cm calculation." %(exptime))
        # Assumes atmosphere used in system throughput is X=1.0
        X = 1.0
        Cm[f] = (m5[f] - 0.5*(skyMag[f] - 21) - 2.5*np.log10(0.7/lsstDefaults.FWHMeff(f)))
        # Calculate Cm_Infinity by setting readout noise to zero.
        m5inf = SignalToNoise.calcM5(darksky, system[f], hardware[f],  photParams_infinity,
                                     FWHMeff=lsstDefaults.FWHMeff(f))
        Cm_infinity = (m5inf - 0.5*(skyMag[f] - 21)
                       - 2.5*np.log10(0.7/lsstDefaults.FWHMeff(f)))
        dCm_infinity[f] = Cm_infinity - Cm[f]
    print(title)
    print('Filter FWHMeff FWHMgeom SkyMag SkyCounts Tb Sb kAtm Gamma Cm dCm_infinity m5 SourceCounts')
    for f in ('u', 'g' ,'r', 'i', 'z', 'y'):
        print('%s %.2f %.2f %.2f %.1f %.3f %.3f %.4f %.6f %.2f %.2f %.2f %.2f'
              %(f, lsstDefaults.FWHMeff(f),
                SignalToNoise.FWHMeff2FWHMgeom(lsstDefaults.FWHMeff(f)),
                skyMag[f], skyCounts[f], Tb[f], Sb[f], kAtm[f],
                gamma[f], Cm[f], dCm_infinity[f], m5[f], sourceCounts[f]))

    # Show what these look like individually (add sky & m5 limits on throughput curves)
    plt.figure()
    for f in filterlist:
        plt.plot(system[f].wavelen, system[f].sb, color=filtercolors[f], linewidth=2, label=f)
    plt.plot(atmosphere.wavelen, atmosphere.sb, 'k:', label='X=1.0')
    plt.legend(loc='center right', fontsize='smaller')
    plt.xlim(300, 1100)
    plt.ylim(0, 1)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Throughput')
    plt.title('System Throughputs')
    plt.grid(True)

    plt.figure()
    ax = plt.gca()
    # Add dark sky
    ax2 = ax.twinx()
    plt.sca(ax2)
    skyab = -2.5*np.log10(darksky.fnu) - darksky.zp
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
    plt.legend(loc=(0.01, 0.69), handles=handles, fancybox=True, numpoints=1, fontsize='small')
    plt.ylim(0, 1)
    plt.xlim(300, 1100)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')
    plt.title('System total response curves %s' %(title))
    return m5


if __name__ == '__main__':

    # Set the directories for each component.
    # Note that this sets the detector to be the 'generic detector' (minimum of each vendor).
    throughputDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    defaultDirs = bu.setDefaultDirs()

    # Build the system and hardware throughput curves (without aerosols, with X=1.0).
    atmosphere = bu.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
    hardware = {}
    system = {}
    for f in ['u', 'g', 'r', 'i', 'z', 'y']:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        hardware[f].readThroughputList(componentList=['detector.dat', 'filter_'+f+'.dat','lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat'],
                                       rootDir=throughputDir)
        system[f].wavelen, system[f].sb = hardware[f].multiplyThroughputs(atmosphere.wavelen, atmosphere.sb)

    m5 = calcM5(hardware, system, atmosphere, title='')


    plt.show()
