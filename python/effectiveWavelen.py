import os
import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters, LSSTdefaults
import bandpassUtils as bu

def calcEffWavelen(hardware, title):
    photParams = PhotometricParameters()
    lsstDefaults = LSSTdefaults()
    atmos = {}
    stdAtmoFile = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'siteProperties/pachonModtranAtm_12.dat')
    atmos['std'] = Bandpass()
    atmos['std'].readThroughput(stdAtmoFile)
    Xarr = np.arange(1.0, 2.55, 0.1)
    for X in Xarr:
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

        effsb_vals = np.zeros(len(Xarr), float)
        for i, X in enumerate(Xarr):
            effsb_vals[i] = effsb['%.1f' %X]
        p = np.polyfit(Xarr, effsb_vals, 1)
        print p
        effsb_eval = np.polyval(p, Xarr) 
        print (effsb_vals - effsb_eval).max(), (effsb_vals-effsb_eval).std()
        plt.plot(Xarr, effsb_vals-effsb_eval, label=f)
    plt.legend()
    plt.show()

if __name__ == '__main__':

    defaultDirs = bu.setDefaultDirs(rootDir = '..')

    addLosses = True

    # Minimum detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Min: Std Atmo')
    exit()
    # vendor1 detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector/vendor1')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Vendor 1: Std Atmo')

        # Minimum detector throughputs
    defaultDirs['detector'] = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'components/camera/detector/vendor2')
    hardware, system = bu.buildHardwareAndSystem(defaultDirs, addLosses)
    calcEffWavelen(hardware, 'Vendor 2: Std Atmo')

