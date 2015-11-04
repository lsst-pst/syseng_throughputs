import os
import numpy as np
from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, LSSTdefaults
import bandpassUtils as bu

def calcEffWavelen(hardware, title):
    photParams = PhotometricParameters()
    lsstDefaults = LSSTdefaults()
    atmos = {}
    stdAtmoFile = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'siteProperties/pachonModtranAtm_12.dat')
    atmos['std'] = Bandpass()
    atmos['std'].readThroughput(stdAtmoFile)
    for X in np.arange(1.0, 2.55, 0.1):
        atmos['%.1f' %X] = Bandpass()
        atmos['%.1f' %X].readThroughput(os.path.join(os.getenv('THROUGHPUTS_DIR'), 'atmos/atmos_%d.dat' %(int(X*10))))
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

    addLosses = True

    # Minimum detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Min: Std Atmo')

    # vendor1 detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector/vendor1')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Vendor 1: Std Atmo')

        # Minimum detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector/vendor2')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Vendor 2: Std Atmo')

