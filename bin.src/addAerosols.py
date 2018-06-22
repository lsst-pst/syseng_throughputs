#!/usr/bin/env python

import os, argparse
import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass
import lsst.sims.syseng.throughputs.bandpassUtils as bu

def addAerosol(atmosphere, X, tau0=0.05, alpha=1.0, wavelen0=440.0, plotAtmosphere=True):
    # Calculate the aerosol contribution -- sb with aerosols = sb*exp(-tau)
    tau = tau0 * np.power((wavelen0/atmosphere.wavelen), alpha)
    # Generate new atmosphere bandpass with aerosols.
    atmosphere_aerosol = Bandpass()
    atmosphere_aerosol.setBandpass(wavelen = atmosphere.wavelen,
                                   sb = atmosphere.sb * np.exp(-tau*X))

    if plotAtmosphere:
        # Plot for a check:
        atmodict = {'Original atmosphere':atmosphere,
                    'With aerosols': atmosphere_aerosol}
        bu.plotBandpasses(atmodict)

    return atmosphere_aerosol

if __name__=='__main__':

    parser = argparse.ArgumentParser(description=
                                     'Add a standard aerosol component to an atmosphere without aerosols')
    parser.add_argument('atmosphereFile', type=str, default=None,
                        help='The atmosphere file to add aerosols to.')
    parser.add_argument('airmass', type=float, default=None, help='The airmass of the atmosphere file.')
    args = parser.parse_args()
    atmosDir, atmosFile = os.path.split(args.atmosphereFile)

    defaultDirs = bu.setDefaultDirs()
    atmosphere = bu.readAtmosphere(atmosDir, atmosFile)
    atmosphere_aerosol = addAerosol(atmosphere, args.airmass)

    outfile = atmosphere.bandpassname.replace('.dat', '') + '_aerosol.dat'
    atmosphere_aerosol.writeThroughput(outfile, print_header='Added aerosol component')

    plt.show()
