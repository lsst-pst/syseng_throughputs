# Calculate and print to screen the effective wavelengths for the bandpasses
#   Builds the bandpasses from the individual component files, using bandpassUtils
#   Reads the standard atmosphere from siteProperties/pachonModtranAtm_12.dat
#   Reads additional atmospheres from the $LSST_THROUGHPUTS repo, if available.
#
#  Set environment variable SYSENG_THROUGHPUTS_DIR to be the root of syseng_throughputs,
#   (automatically done if package set up via eups).

import os
import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, LSSTdefaults
import bandpassUtils as bu

def calcEffWavelen(hardware, title, throughputDir=None):
    photParams = PhotometricParameters()
    lsstDefaults = LSSTdefaults()
    atmos = {}
    stdAtmoFile = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'siteProperties/pachonModtranAtm_12.dat')
    atmos['std'] = Bandpass()
    atmos['std'].readThroughput(stdAtmoFile)
    multiAtmos = False
    if throughputDir is None:
        throughputDir = os.getenv('THROUGHPUTS_DIR')
    if throughputDir is not None:
        multiAtmos = True
        Xarr = np.arange(1.0, 2.55, 0.1)
        for X in Xarr:
            atmos['%.1f' %X] = Bandpass()
            atmos['%.1f' %X].readThroughput(os.path.join(throughputDir,
                                                         'atmos/atmos_%d.dat' %(int(X*10))))

    atmoskeys = sorted(atmos.keys())
    print title
    print '    %s' %('    '.join(atmoskeys))
    system = Bandpass()
    effsb = {}
    for f in ['u', 'g', 'r', 'i', 'z', 'y']:
        writestring = '%s ' %f
        for k in atmoskeys:
            system.wavelen, system.sb = hardware[f].multiplyThroughputs(atmos[k].wavelen, atmos[k].sb)
            effphi, effsb[k] = system.calcEffWavelen()
            writestring += '%.2f ' %(effsb[k])
        print writestring

if __name__ == '__main__':

    defaultDirs = bu.setDefaultDirs(rootDir = '..')
    # To use a particular vendor's detector response curve, uncomment one of the lines below
    defaultDetector = defaultDirs['detector']
    #defaultDirs['detector'] = os.path.join(defaultDetector, 'vendor1')
    #defaultDirs['detector'] = os.path.join(defaultDetector, 'vendor2')

    addLosses = True

    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Min: Std Atmo')


