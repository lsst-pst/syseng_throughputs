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

    # Build the separate vendor detectors.
    qevendors = {}
    # Add losses into the throughputs?
    addLosses=True
    qevendors['loss'] = {}
    # Vendor 1
    qevendors['loss'][1] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor1'), addLosses)
    # Vendor 2
    qevendors['loss'][2] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor2'), addLosses)
    # Generic 'minimum' detector throughputs.
    qevendors['loss']['Min'] = bu.buildDetector(defaultDirs['detector'], addLosses)
    bu.plotBandpasses(qevendors['loss'], title='Combining Vendor Detector Responses Losses')


    qevendors['noloss'] = {}
    # Add losses into the throughputs?
    addLosses = False
    # Build the separate vendor detectors.
    # Vendor 1
    qevendors['noloss'][1] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor1'), addLosses)
    # Vendor 2
    qevendors['noloss'][2] = bu.buildVendorDetector(os.path.join(defaultDirs['detector'], 'vendor2'), addLosses)
    # Generic 'minimum' detector throughputs.
    qevendors['noloss']['Min'] = bu.buildDetector(defaultDirs['detector'], addLosses)
    bu.plotBandpasses(qevendors['noloss'], title='Combining Vendor Detector Responses No Losses')


    compare = {}
    #compare['min_loss'] = qevendors['loss']['Min']
    #compare['min_noloss'] = qevendors['noloss']['Min']
    compare['v1_losses'] = qevendors['loss'][1]
    compare['v1_nolosses'] = qevendors['noloss'][1]
    #compare['v2_losses'] = qevendors['loss'][2]
    #compare['v2_nolosses'] = qevendors['noloss'][2]
    bu.plotBandpasses(compare, title='No losses vs losses')

    plt.show()
