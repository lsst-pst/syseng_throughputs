import os
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass, Sed, SignalToNoise, PhotometricParameters, LSSTdefaults
import bandpassUtils as bu

def calcM5s(hardware, system, title='m5'):
    photParams = PhotometricParameters()
    lsstDefaults = LSSTdefaults()
    darksky = Sed()
    darksky.readSED_flambda(os.path.join('../siteProperties', 'darksky.dat'))
    m5 = {}
    for f in system:
        m5[f] = SignalToNoise.calcM5(darksky, system[f], hardware[f], photParams, seeing=lsstDefaults.seeing(f))
    print title
    for f in ('u', 'g' ,'r', 'i', 'z', 'y'):
        print '%s %.2f' %(f, m5[f])
    return m5


if __name__ == '__main__':

    defaultDirs = bu.setDefaultDirs(rootDir = '..')

    addLosses = True

    # Build the separate vendor detectors.
    qevendors = {}
    qevendors[1] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor1'), addLosses)
    qevendors[2] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor2'), addLosses)
    qevendors['Min'] = bu.buildGenericDetector(defaultDirs['detector'], addLosses)
    bu.plotBandpasses(qevendors, title='Vendor Detector Responses')

    # Build each component. Make plot to compare to Chuck's / previous version.
    comparison = Bandpass()

    detector = bu.buildGenericDetector(defaultDirs['detector'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/camera/detThroughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-', label='comparison detector')
    plt.plot(detector.wavelen, detector.sb, 'r-', label='calculated detector')
    plt.title('detector')

    lens1 = bu.buildLens(defaultDirs['lens1'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/camera/lens1Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-', label='comparison lens1')
    plt.plot(lens1.wavelen, lens1.sb, 'r-', label='calculated lens1')
    plt.title('l1')

    lens2 = bu.buildLens(defaultDirs['lens2'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/camera/lens2Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-', label='comparison lens2')
    plt.plot(lens2.wavelen, lens2.sb, 'r-', label='calculated lens2')
    plt.title('l2')

    lens3 = bu.buildLens(defaultDirs['lens3'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/camera/lens3Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-', label='comparison lens3')
    plt.plot(lens3.wavelen, lens3.sb, 'r-', label='calculated lens3')
    plt.title('l3')

    filters = bu.buildFilters(defaultDirs['filters'], addLosses)
    plt.figure()
    for f in filters:
        comparison.readThroughput('../components/camera/'+f+'BandThroughput.dat')
        plt.plot(comparison.wavelen, comparison.sb, 'k-')
        plt.plot(filters[f].wavelen, filters[f].sb, 'r-')
    plt.title('filters')

    mirror1 = bu.buildMirror(defaultDirs['mirror1'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/telescope/m1Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-')
    plt.plot(mirror1.wavelen, mirror1.sb, 'r-')
    plt.title('m1')

    mirror2 = bu.buildMirror(defaultDirs['mirror2'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/telescope/m2Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-')
    plt.plot(mirror2.wavelen, mirror2.sb, 'r-')
    plt.title('m2')

    mirror3 = bu.buildMirror(defaultDirs['mirror3'], addLosses)
    plt.figure()
    comparison.readThroughput('../components/telescope/m3Throughput.dat')
    plt.plot(comparison.wavelen, comparison.sb, 'k-')
    plt.plot(mirror3.wavelen, mirror3.sb, 'r-')
    plt.title('m3')

    atmosphere = bu.buildAtmosphere(defaultDirs['atmosphere'])

    # Plot all components.
    plt.figure()
    plt.plot(detector.wavelen, detector.sb, 'k-', linewidth=2, label='Detector')
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

    # Combine components (and individual combination for each detector vendor) by hand.
    core_sb = detector.sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        wavelen = filters[f].wavelen
        hw_sb = core_sb * filters[f].sb
        hardware[f].setBandpass(wavelen, hw_sb)
        system[f].setBandpass(wavelen, hw_sb*atmosphere.sb)
    calcM5s(hardware, system, 'Min m5')

    # Combine components (and individual combination for each detector vendor) by hand.
    core_sb = qevendors[1].sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        wavelen = filters[f].wavelen
        hw_sb = core_sb * filters[f].sb
        hardware[f].setBandpass(wavelen, hw_sb)
        system[f].setBandpass(wavelen, hw_sb*atmosphere.sb)
    calcM5s(hardware, system, 'Vendor1 m5')

    core_sb = qevendors[2].sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        wavelen = filters[f].wavelen
        hw_sb = core_sb * filters[f].sb
        hardware[f].setBandpass(wavelen, hw_sb)
        system[f].setBandpass(wavelen, hw_sb*atmosphere.sb)
    calcM5s(hardware, system, 'Vendor2 m5')

    # Build the whole system and hardware in one go.
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)


    bu.plotBandpasses(hardware)

    bu.plotBandpasses(system)
    plt.plot(atmosphere.wavelen, atmosphere.sb, 'k:')
    plt.figtext(0.22, 0.75, 'Airmass 1.2')


    plt.show()
