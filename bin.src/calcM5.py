#!/usr/bin/env python

# Calculate m5 values using the SYSENG_THROUGHPUTS files.
from __future__ import print_function
import lsst.syseng.throughputs as st

if __name__ == '__main__':

    # Set the directories for each component.
    # Note that this sets the detector to be the 'generic detector' (minimum of each vendor).
    defaultDirs = st.setDefaultDirs()
    # To use a particular vendor, uncomment one of the following lines or edit as necessary.
    # defaultDirs['detector'].replace('joint_mininimum', 'itl')
    # defaultDirs['detector'].replace('joint_minimum', 'e2v')

    # Add losses to each component?
    addLosses = True

    # Build the system and hardware throughput curves (with aerosols, with X=1.2).
    atmosphere = st.readAtmosphere(defaultDirs['atmosphere'], atmosFile='pachonModtranAtm_12_aerosol.dat')
    hardware, system = st.buildHardwareAndSystem(defaultDirs, addLosses, atmosphereOverride=atmosphere)
    m5 = st.makeM5(hardware, system, X=1.2)
    print(m5)

    # At X=1.0
    atmosphere = st.readAtmosphere(defaultDirs['atmosphere'], atmosFile='atmos_10_aerosol.dat')
    hardware, system = st.buildHardwareAndSystem(defaultDirs, addLosses, atmosphereOverride=atmosphere)
    m5 = st.makeM5(hardware, system, X=1.0)
    print(m5)
