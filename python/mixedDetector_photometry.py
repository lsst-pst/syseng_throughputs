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


def calcM5s(hardware, system, atmos, title='m5'):
    photParams = PhotometricParameters()
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
        m5[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams, seeing=lsstDefaults.seeing(f))
        fNorm = flatSed.calcFluxNorm(m5[f], system[f])
        flatSed.multiplyFluxNorm(fNorm)
        sourceCounts[f] = flatSed.calcADU(system[f], photParams=photParams)
        # Calculate the Skycounts expected in this bandpass.
        skyCounts[f] = darksky.calcADU(hardware[f], photParams=photParams) * photParams.platescale**2
        # Calculate the sky surface brightness.
        skyMag[f] = darksky.calcMag(hardware[f])
        # Calculate the gamma value.
        gamma[f] = SignalToNoise.calcGamma(system[f], m5[f], photParams)
    print title
    print 'Filter m5 SourceCounts SkyCounts SkyMag Gamma'
    for f in ('u', 'g' ,'r', 'i', 'z', 'y'):
        print '%s %.2f %.1f %.2f %.2f %.6f' %(f, m5[f], sourceCounts[f], skyCounts[f], skyMag[f], gamma[f])

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
    plt.plot(atmos.wavelen, atmos.sb, 'k:', label='Atmosphere, X=1.2')
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



    photParams = PhotometricParameters()
    lsstDefaults = LSSTdefaults()

    # Build the separate vendor detectors.
    qevendors = {}
    qevendors[1] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor1'), addLosses)
    qevendors[2] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor2'), addLosses)
    qevendors['combo'] = bu.buildGenericDetector(defaultDirs['detector'], addLosses)
    bu.plotBandpasses(qevendors, title='Vendor Detector Responses')

    # Build the other components.
    lens1 = bu.buildLens(defaultDirs['lens1'], addLosses)
    lens2 = bu.buildLens(defaultDirs['lens2'], addLosses)
    lens3 = bu.buildLens(defaultDirs['lens3'], addLosses)
    filters = bu.buildFilters(defaultDirs['filters'], addLosses)
    mirror1 = bu.buildMirror(defaultDirs['mirror1'], addLosses)
    mirror2 = bu.buildMirror(defaultDirs['mirror2'], addLosses)
    mirror3 = bu.buildMirror(defaultDirs['mirror3'], addLosses)
    atmosphere = bu.buildAtmosphere(defaultDirs['atmosphere'])

    # Plot all components.
    plt.figure()
    plt.plot(qevendors['combo'].wavelen, qevendors['combo'].sb, 'k-', linewidth=2, label='Detector')
    plt.plot(lens1.wavelen, lens2.sb, 'g-', linewidth=2, label='L1')
    plt.plot(lens2.wavelen, lens2.sb, 'r-', linewidth=2, label='L2')
    plt.plot(lens3.wavelen, lens3.sb, 'b-', linewidth=2, label='L3')
    for f in ['u', 'g', 'r', 'i', 'z', 'y']:
        plt.plot(filters[f].wavelen, filters[f].sb, linestyle=':', linewidth=5, label=f)
    plt.plot(mirror1.wavelen, mirror1.sb, 'g-.', linewidth=2, label='M1')
    plt.plot(mirror2.wavelen, mirror2.sb, 'r--', linewidth=2, label='M2')
    plt.plot(mirror3.wavelen, mirror3.sb, 'b--', linewidth=2, label='M3')
    plt.plot(atmosphere.wavelen, atmosphere.sb, 'k:', linewidth=2, label='X=1.2')
    plt.legend(loc=(0.96, 0.2), numpoints=1, fontsize='smaller', fancybox=True)
    plt.xlim(300, 1100)
    plt.ylim(0, 1)
    plt.title('Throughput components')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')

    hardware = {}
    system = {}
    m5 = {}
    # Combine components (and individual combination for each detector vendor) by hand.
    for detector in ['combo', 1, 2]:
        core_sb = qevendors[detector].sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
        hardware[detector] = {}
        system[detector] = {}
        m5[detector] = {}
        for f in filters:
            hardware[detector][f] = Bandpass()
            system[detector][f] = Bandpass()
            wavelen = filters[f].wavelen
            hw_sb = core_sb * filters[f].sb
            hardware[detector][f].setBandpass(wavelen, hw_sb)
            system[detector][f].setBandpass(wavelen, hw_sb*atmosphere.sb)
        m5[detector] = calcM5s(hardware[detector], system[detector], atmosphere, title='Vendor %s' %detector)

    # Show what these look like (print m5 limits on throughput curves)
    plt.figure()
    for f in filterlist:
        for det in [1, 2]:
            if det == 1:
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
    for det in [1, 2]:
        for f in filterlist:
            if det == 1:
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
    for det in (1, 2):
        mags[det] = su.calcNatMags(system[det], seds)
    dmags = su.calcDeltaMags(mags[2], mags[1], mmags=True, matchBlue=False)

    gi = su.calcGiColors(mags[1])
    ug = su.calcAnyColor(mags[1], 'u', 'g')

    su.plotDmags(gi, dmags, titletext='Color terms between detectors, as function of vendor1 g-i color')
    su.plotDmagsSingle(gi, dmags, titletext='Color terms between detectors, as function of vendor1 g-i color')
    plt.ylim(-60, 60)
    plt.xlim(-1, 3)
    plt.grid()

    su.plotDmags(ug, dmags, colorname='u-g', xlim=[-1, 3],
                titletext='Color terms between detectors, as function of vendor1 u-g color')
    su.plotDmagsSingle(ug, dmags, colorname='u-g',
                    titletext='Color terms between detectors, as function of vendor1 u-g color')
    plt.ylim(-60, 60)
    plt.xlim(-1, 3)
    plt.grid()

    plt.show()
