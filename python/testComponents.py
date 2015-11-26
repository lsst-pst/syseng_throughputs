import os
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass
import bandpassUtils as bu

# This script will combine the individual components of the camera and telescope,
# generating all intermediate files as well as final combined system and hardware curves.
#  Note that it does NOT replace the files in $SYSENG_THROUGHPUTS/intermediateFiles,
#  but merely uses those files as a comparison point, making plots that compare the
#  new values to the old ones.
# If you DO you want to replace the intermediate files, look at the ipython notebooks.


if __name__ == '__main__':

    # Get the names of the directories containing each component.
    # (returns a dictionary).
    defaultDirs = bu.setDefaultDirs(rootDir = '..')

    # Add losses into the throughputs?
    addLosses = True

    throughputs = {}

    # Build the separate vendor detectors.
    qevendors = {}
    # Vendor 1
    qevendors[1] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor1'), addLosses)
    # Vendor 2
    qevendors[2] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor2'), addLosses)
    # Generic 'minimum' detector throughputs.
    qevendors['Min'] = bu.buildDetector(defaultDirs['detector'], addLosses)
    throughputs['detector'] = qevendors['Min']
    bu.plotBandpasses(qevendors, title='Combining Vendor Detector Responses')

    # Read the previously generated 'intermediate' file and plot the comparison.
    comparison = Bandpass()
    oldDetector = '../intermediateFiles/components/camera/detThroughput.dat'
    comparison.readThroughput(oldDetector)
    # Plot old and new detectors.
    plotDict = {'New Detector':throughputs['detector'], 'Old Detector':comparison}
    bu.plotBandpasses(plotDict, title='Compare Detector Throughputs')

    # Build and compare the lens throughput curves, for lens1/2/3.
    for lens in ('lens1', 'lens2', 'lens3'):
        throughputs[lens] = bu.buildLens(defaultDirs[lens], addLosses)
        # Read the old intermediate file version.
        comparison.readThroughput('../intermediateFiles/components/camera/%sThroughput.dat'
                                  %(lens))
        plotDict = {'New %s' %(lens): throughputs[lens], 'Old %s' %(lens): comparison}
        bu.plotBandpasses(plotDict, title='Compare %s' %(lens))

    # Build and compare the filter curves.
    filters = bu.buildFilters(defaultDirs['filters'], addLosses)
    oldfilters = {}
    for f in filters:
        oldfilters[f] = Bandpass()
        oldfilters[f].readThroughput('../intermediateFiles/components/camera/'+f+'BandThroughput.dat')
    bu.plotBandpasses(filters)
    bu.plotBandpasses(oldfilters, newfig=False, linestyle=':', title='Compare filters')
    # Put the individual filter bandpasses straight into 'throughputs', for easier use below.
    for f in filters:
        throughputs[f] = filters[f]

    # Build and compare the mirror curves, for mirrors 1/2/3.
    for mirror in ('mirror1', 'mirror2', 'mirror3'):
        throughputs[mirror] = bu.buildMirror(defaultDirs[mirror], addLosses)
        mfile = ('../intermediateFiles/components/telescope/' + mirror[0] + mirror[-1] +
                 'Throughput.dat')
        comparison.readThroughput(mfile)
        plotDict = {'New %s' %(mirror):throughputs[mirror], 'Old %s' %(mirror):comparison}
        bu.plotBandpasses(plotDict, title='Compare %s' %(mirror))

    # Read the atmosphere file.
    throughputs['atmosphere'] = bu.readAtmosphere(defaultDirs['atmosphere'])

    # Plot all components.
    bu.plotBandpasses(throughputs, title='All components')

    # Combine components by hand. Compare to combination returned from bandpassUtils.
    core_sb = (throughputs['detector'].sb * throughputs['lens1'].sb * throughputs['lens2'].sb
            * throughputs['lens3'].sb * throughputs['mirror1'].sb
            * throughputs['mirror2'].sb * throughputs['mirror3'].sb)
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        wavelen = filters[f].wavelen
        hw_sb = core_sb * filters[f].sb
        hardware[f].setBandpass(wavelen, hw_sb)
        system[f].setBandpass(wavelen, hw_sb*throughputs['atmosphere'].sb)
    # Get combination from bandpassUtils.
    bU_hardware, bU_system = bu.buildHardwareAndSystem(defaultDirs)

    bu.plotBandpasses(hardware)
    bu.plotBandpasses(bU_hardware, linestyle=':', newfig=False,
                      title='BandpassUtils components, hardware')

    bu.plotBandpasses(system)
    plt.plot(throughputs['atmosphere'].wavelen, throughputs['atmosphere'].sb, 'k:')
    plt.figtext(0.22, 0.75, 'Airmass 1.2')
    bu.plotBandpasses(bU_system, linestyle=':', newfig=False,
                      title='BandpassUtils components, system')

    plt.show()
