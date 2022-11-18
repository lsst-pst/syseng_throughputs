#!/usr/bin/env python

from __future__ import print_function
# Calculate m5 values using the *throughputs* repo files.
import os
import pandas as pd
from rubin_sims.phot_utils import Bandpass, Sed, SignalToNoise
import syseng_throughputs as st

if __name__ == '__main__':

    # Set the max number of columns for pandas output of the m5 values.
    pd.set_option('display.max_columns', 500)

    # Set the directories for each component.
    # Note that this sets the detector to be the 'generic detector' (minimum of each vendor).
    throughputDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    defaultDirs = st.setDefaultDirs()

    # Build the system and hardware throughput curves (without aerosols, with X=1.0).
    print("Calculating m5 using X=1.0 atmosphere")
    atmosphere = st.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
    hardware = {}
    system = {}
    for f in ['u', 'g', 'r', 'i', 'z', 'y']:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        hardware[f].readThroughputList(componentList=['detector.dat', 'filter_'+f+'.dat','lens1.dat',
                                                      'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat'],
                                       rootDir=throughputDir)
        system[f].wavelen, system[f].sb = hardware[f].multiplyThroughputs(atmosphere.wavelen, atmosphere.sb)

    m5 = st.makeM5(hardware, system, X=1.0)
    print(m5)
