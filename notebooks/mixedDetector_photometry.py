from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, LSSTdefaults, SignalToNoise
import bandpassUtils as bu
import sedUtils as su

filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
filtercolors = {'u':'b', 'g':'c', 'r':'g',
                'i':'y', 'z':'r', 'y':'m'}

def calcM5(hardware, system, atmos, title='m5'):
    effarea = np.pi * (6.423/2.0*100.)**2
    photParams = PhotometricParameters(effarea = effarea)
    lsstDefaults = LSSTdefaults()
    darksky = Sed()
    darksky.readSED_flambda(os.path.join('../siteProperties', 'darksky.dat'))
    flatSed = Sed()
    flatSed.setFlatSED()
    m5 = {}
    sourceCounts = {}
    skyCounts = {}
    skyMag = {}
    gamma = {}
    for f in system:
        m5[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams, FWHMeff=lsstDefaults.FWHMeff(f))
        fNorm = flatSed.calcFluxNorm(m5[f], system[f])
        flatSed.multiplyFluxNorm(fNorm)
        sourceCounts[f] = flatSed.calcADU(system[f], photParams=photParams)
        # Calculate the Skycounts expected in this bandpass.
        skyCounts[f] = darksky.calcADU(hardware[f], photParams=photParams) * photParams.platescale**2
        # Calculate the sky surface brightness.
        skyMag[f] = darksky.calcMag(hardware[f])
        # Calculate the gamma value.
        gamma[f] = SignalToNoise.calcGamma(system[f], m5[f], photParams)
    print(title)
    print('Filter m5 SourceCounts SkyCounts SkyMag Gamma')
    for f in ('u', 'g' ,'r', 'i', 'z', 'y'):
        print('%s %.2f %.1f %.2f %.2f %.6f' %(f, m5[f], sourceCounts[f], skyCounts[f], skyMag[f], gamma[f]))

    # Show what these look like individually (add sky & m5 limits on throughput curves)
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
    plt.plot(atmos.wavelen, atmos.sb, 'k:', label='Atmosphere, X=1.0 with aerosols')
    # Add legend for dark sky.
    myline = mlines.Line2D([], [], color='k', linestyle='-', label='Dark sky AB mags')
    handles.append(myline)
    # end of dark sky legend line
    plt.legend(loc=(0.01, 0.69), handles=handles, fancybox=True, numpoints=1, fontsize='small')
    plt.ylim(0, 1)
    plt.xlim(300, 1100)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')
    if title == 'Vendor combo':
        title = ''
    plt.title('System total response curves %s' %(title))
    plt.savefig('../plots/system+sky' + title + '.png', format='png', dpi=600)
    return m5


if __name__ == '__main__':

    defaultDirs = bu.setDefaultDirs(rootDir = '..')
    addLosses = True

    hardware = {}
    system = {}
    m5 = {}
    atmosphere = bu.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
    hardware['combo'], system['combo'] = bu.buildHardwareAndSystem(defaultDirs, atmosphereOverride=atmosphere)
    m5['combo'] = calcM5(hardware['combo'], system['combo'], atmosphere, title='combo')
    detectors = ['itl', 'e2v']
    for det in detectors:
        defaultDirs['detector'].replace('joint_minimum', det)
        hardware[det], system[det] = bu.buildHardwareAndSystem(defaultDirs, atmosphereOverride=atmosphere)
        m5[det] = calcM5(hardware[det], system[det], atmosphere, title=det)

    # Show what these look like (print m5 limits on throughput curves)
    plt.figure()
    for f in filterlist:
        for det in detectors:
            if det == 'itl':
                linestyle = ':'
                spacer= ' '
            else:
                linestyle = '-'
                spacer = ''
            plt.plot(system[det][f].wavelen, system[det][f].sb, linestyle=linestyle, color=filtercolors[f],
                    label='%s %s:%s m5 %.2f' %(det, f, spacer, m5[det][f]))
    plt.legend(loc=(0.85, 0.5), fancybox=True, numpoints=1, fontsize='small')
    plt.ylim(0, 0.6)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')
    plt.title('System total response curves')

    plt.figure()
    for det in detectors:
        for f in filterlist:
            if det == 'itl':
                linestyle = ':'
            else:
                linestyle = '-'
            plt.plot(system[det][f].wavelen, system[det][f].phi, linestyle=linestyle, color=filtercolors[f], label='%s %s' %(det, f))
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='small')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('$\phi$')
    plt.ylim(ymin=0)
    plt.title('Normalized system response curves ($\phi$)')


    # Calculate color terms (magnitudes)

    sedDir = '../seds'
    seds = su.readPhotSeds(sedDir = sedDir)
    redshifts = {}
    redshifts['galaxies'] = np.array([0.5, 1.0], 'float')
    redshifts['quasar'] = np.array([1.0, 1.5, 2.5], 'float')
    redshifts['sn'] = np.array([0.3, 0.8, 1.2, 1.5], 'float')
    redshifts['photoZ_outliers'] = np.array([0, 0.2, 1.0], 'float')
    seds = su.makeRedshiftedSeds(seds, redshifts)
    seds, system[det] = su.matchSedsBp(seds, system[det])

    mags = {}
    mags[1] = {}
    mags[2] = {}
    for det in detectors:
        mags[det] = su.calcNatMags(system[det], seds)
    dmags = su.calcDeltaMags(mags['itl'], mags['e2v'], mmags=True, matchBlue=False)

    gi = su.calcGiColors(mags['itl'])
    ug = su.calcAnyColor(mags['itl'], 'u', 'g')

    su.plotDmags(gi, dmags, titletext='Color terms between detectors, as function of ITL g-i color')
    su.plotDmagsSingle(gi, dmags, titletext='Color terms between detectors, as function of ITL g-i color')
    plt.ylim(-60, 60)
    plt.xlim(-1, 3)
    plt.grid()

    su.plotDmags(ug, dmags, colorname='u-g', xlim=[-1, 3],
                titletext='Color terms between detectors, as function of ITL u-g color')
    su.plotDmagsSingle(ug, dmags, colorname='u-g',
                    titletext='Color terms between detectors, as function of ITL u-g color')
    plt.ylim(-60, 60)
    plt.xlim(-1, 3)
    plt.grid()

    plt.show()
