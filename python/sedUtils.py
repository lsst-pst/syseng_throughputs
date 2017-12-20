from __future__ import print_function
import os
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed


def readPhotSeds(sedDir=None):
    """
    Read all the seds provided by this package, storing them in a
    nested dictionary (sedDict['sedtype']['sedname']).
    """
    if sedDir == None:
        sedDir = os.path.join(getPackageDir('SYSENG_THROUGHPUTS_DIR'), 'seds')
    # Set up the lists for the photometry reference SEDs.
    #  (I'm doing this by hand to keep these lists well-known over time.)
    sedlists = {}
    sedlists['quasar'] = ['quasar.dat',]
    sedlists['stars'] = ['km10_7250.fits_g45', 'km10_6500.fits_g45',
                         'km10_6000.fits_g45', 'km10_5250.fits_g45',
                         'km10_4500.fits_g45', 'm3.0Full.dat']
    sedlists['white_dwarf'] = ['wd_H_100000_80.dat', 'wd_H_15000_80.dat', 'wd_H_50000_80.dat',
                               'wd_H_5500_80.dat', 'wd_He_10000_80.dat', 'wd_He_15000_80.dat', 'wd_He_5500_80.dat']
    sedlists['sn'] = ['sn1a_15.0.dat', 'sn1a_20.0.dat', 'sn1a_10.0.dat']
    sedlists['galaxies'] = ['Sa_template_norm.sed.dat', 'Sdm_template_norm.sed0.dat',
                            'Ell2_template_norm.sed.dat']
    sedlists['photoZ_outliers'] = ['xspec_172.sed.dat', 'xspec_173.sed.dat',
                                   'xspec_175.sed.dat', 'xspec_176.sed.dat',
                                   'xspec_90.sed.dat', 'xspec_91.sed.dat']
    # Let's go read the files.
    sedDict = {}
    # Loop through quasar, stars, galaxies, sn in the sedlist dictionary.
    for objtype in sedlists.keys():
        sedDict[objtype] = {}
        # And for each type of object, loop through and read each SED.
        for s in sedlists[objtype]:
            sedDict[objtype][s] = Sed()
            sedDict[objtype][s].readSED_flambda(os.path.join(sedDir, objtype, s))
    return sedDict

def readAnySeds(inputfileList, sedDir='.'):
    """Read the seds in a list and store them in a format appropriate for use with the
    rest of these routines. sedDir (if set) can be the root directory for the files in the list."""
    # Read the files. We'll store them in a key called 'any' so that this dictionary looks like sedDict from above.
    sedDict = {}
    sedDict['any'] = {}
    for f in inputfileList:
        sedDict['any'][f] = Sed()
        sedDict['any'][f].readSED_flambda(os.path.join(sedDir, filename))
    return sedDict

def makeRedshiftedSeds(sedDict, redshifts):
    """Redshift the quasar, galaxies and SN by the amounts given in redshifts.
    If redshift is a numpy array, all objects except stars are redshifted to the same values.
    If redshift is a dictionary containing numpy arrays in the 'quasar', 'sn' or 'galaxies' keys, then
    those redshifts are applied to those objects only. """
    # Check what kind of redshift object we received.
    if isinstance(redshifts, dict):
        keys = list(redshifts.keys())
        # Loop over all the object types in the redshift dictionary.
        for objtype in keys:
            sedlist = list(sedDict[objtype].keys())
            for sed in sedlist:
                for z in redshifts[objtype]:
                    newsedname = sed + '_Z_%.3f' %(z)
                    sedDict[objtype][newsedname] = redshiftSingleSED(sedDict[objtype][sed], z)
    else:
        # Using the same redshift list for everything (except stars).
        for objtype in sedDict:
            if objtype != 'stars':
                sedlist = list(sedDict[objtype].keys())
                for sed in sedlist:
                    for z in redshifts:
                        newsedname = sed + '_Z_%.3f' %(z)
                        sedDict[objtype][newsedname] = redshiftSingleSED(sedDict[objtype][sed], z)
    return sedDict

def redshiftSingleSED(sed_in, z):
    # Make a copy of the input SED, as otherwise we will overwrite original arrays.
    sed_out = deepcopy(sed_in)
    # Set dimming False just to keep approximate fnu normalization constant (makes it easier to plot seds on same figure)
    sed_out.redshiftSED(z, dimming=False)
    return sed_out


def matchSedsBp(sedDict, bpDict, refFilter=None):
    """
    Match the wavelength ranges for all the Seds and the bandpass dictionary.
    This will speed up later calculation of magnitudes (if you're doing a lot of them),
     but is not strictly necessary.
     """
    if refFilter == None:
        refFilter = bpDict.keys()[0]
    wavelen_match = bpDict[refFilter].wavelen
    # Check all filters in bpDict match in wavelength space (note bandpasses must be regular grid).
    for f in bpDict.keys():
        if np.any(bpDict[f].wavelen != wavelen_match):
            bpDict[f].resampleBandpass(wavelen_min=wavelen_match.min(), wavelen_max=wavelen_match.max(),
                                       wavelen_step = wavelen_match[1] - wavelen_match[0])
    # Check all seds in sedDict match in wavelength space.
    for objtype in sedDict:
        for s in sedDict[objtype]:
            sedDict[objtype][s].resampleSED(wavelen_match=wavelen_match)
            sedDict[objtype][s].flambdaTofnu()
    return sedDict, bpDict


def calcNatMags(bandpassDict, sedDict):
    """
    Calculate (zeropoint-calibrated) natural magnitudes for each SED,
    in all bandpasses. Changes in this magnitude includes only color/wavelength-dependent effects.
    """
    # Create a dictionary to hold all of the magnitude information, for all filters.
    mags = {}
    # Loop over SEDs:
    for objtype in sedDict:
        mags[objtype] = {}
        for s in sedDict[objtype]:
            # Create a dictionary for each object to hold the multiple filter information.
            mags[objtype][s] = {}
            for f in bandpassDict.keys():
                # Calculate the magnitudes.
                mags[objtype][s][f] = sedDict[objtype][s].calcMag(bandpassDict[f])
    return mags

def calcInstMags(bandpassDict, sedDict, photParams):
    """
    Calculate instrumental magnitudes for each SED, in all bandpasses.
    Changes in this magnitude includes gray-scale effects as well as color/wavelength dependent effects.
    """
    mags = {}
    for objtype in sedDict:
        mags[objtype] = {}
        for s in sedDict[objtype]:
            mags[objtype][s] = {}
            for f in bandpassDict.keys():
                mags[objtype][s][f] = sedDict[objtype][s].calcADU(bandpassDict[f], photParams)
                mags[objtype][s][f] = -2.5*numpy.log10(mags[objtype][s][f])
    return mags


def calcDeltaMags(mags1, mags2, mmags=True, matchBlue=False):
    """
    Calculate the difference in magnitudes between two magnitudes, mags calculated as above (mag[objtype][sedname][filter])
    If 'mmags' is True, then returns delta mags in mmags.
    If 'matchBlue' is True, then scales change in magnitude between mags1 and mags2 so that the blue
     star 'km10_7250.fits_g45' has zero magnitude change.
    """
    dmags = {}
    for objtype in mags1:
        if objtype in mags2:
            dmags[objtype] = {}
            for sedname in mags1[objtype]:
                if sedname in mags2[objtype]:
                    dmags[objtype][sedname] = {}
                    for f in mags1[objtype][sedname]: # assume filters are the same between mag dicts
                        dmags[objtype][sedname][f] = mags1[objtype][sedname][f] - mags2[objtype][sedname][f]
                        if mmags:
                            dmags[objtype][sedname][f] *= 1000.0
    # Apply scaling so that bluest star remains constant, if desired.
    if matchBlue:
        # Calculate offset to apply.
        offset = {}
        for f in dmags_filters:
            offset[f] = dmags['stars']['km10_7250.fits_g45'][f]
        # Apply offset.
        for objtype in dmags:
            for sedname in dmags[objtype]:
                for f in dmags[objtype][sedname]:
                    dmags[objtype][sedname][f] = dmags[objtype][sedname][f] - offset[f]
    return dmags

def calcGiColors(mags):
    """Calculate the g-i colors of objects in a set of magnitudes, mags calculated as above."""
    gi = {}
    for objtype in mags:
        gi[objtype] = {}
        for s in mags[objtype]:
            gi[objtype][s] = mags[objtype][s]['g'] - mags[objtype][s]['i']
    return gi

def calcAnyColor(mags, color1, color2):
    """Calculate any color of objects in a set of magnitudes. Specify color1 (e.g. 'g') and color2 (e.g. 'i')."""
    color = {}
    for objtype in mags:
        color[objtype] = {}
        for s in mags[objtype]:
            color[objtype][s] = mags[objtype][s][color1] - mags[objtype][s][color2]
    return color

def printDmags(dmags, filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Print changes in magnitudes to the screen."""
    print('Delta mmag:')
    writestring = "object"
    for f in filterlist:
        writestring += '\t %s ' %(f)
    print(writestring)
    for objtype in dmags:
        print('Object type: ', objtype)
        for s in dmags[objtype]:
            writestring = 'dm %s ' %(s)
            for f in filterlist:
                writestring += ' %f ' %(dmags[objtype][s][f])
            print(writestring)
    return


def plotDmags(color, dmags, colorname = 'g-i', xlim=None,
              newfig=True, titletext=None, savefig=False,
              filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Generate a plot of the change in magnitudes as a function of color. """
    symbs = {'quasar':'o', 'stars':'*', 'white_dwarf':'*', 'sn':'x', 'galaxies':'s', 'photoZ_outliers':'+', 'any':'^'}
    colors = {'quasar':'g', 'stars':'r', 'white_dwarf':'m','sn':'b', 'galaxies':'g', 'photoZ_outliers':'g', 'any':'k'}
    if newfig:
        plt.figure()
    plt.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    for f, i in zip(filterlist, range(1, len(filterlist)+1)):
        plt.subplot(3,2,i)
        for objtype in dmags:
            for s in dmags[objtype]:
                plt.plot(color[objtype][s], dmags[objtype][s][f], color=colors[objtype], marker=symbs[objtype])
        plt.xlabel(colorname)
        plt.ylabel(r'$\Delta$%s (mmag)' %(f))
        if xlim is not None:
            plt.xlim(xlim)
    plt.suptitle(titletext)
    if savefig:
        if titletext != None:
            plt.savefig('%s.%s' %(titletext, figformat), format=figformat)
        else:
            plt.savefig('Dmag.%s' %(figformat), format=figformat)
    return


def plotDmagsSingle(color, dmags, colorname='g-i', plotFilter='u', newfig=True, titletext=None, savefig=False):
    """Generate a plot of the change in magnitudes, in a single filter only. """
    symbs = {'quasar':'o', 'stars':'*', 'white_dwarf':'*', 'sn':'x', 'galaxies':'s', 'photoZ_outliers':'+', 'any':'^'}
    colors = {'quasar':'g', 'stars':'r', 'white_dwarf':'m', 'sn':'b', 'galaxies':'g', 'photoZ_outliers':'g', 'any':'k'}
    if newfig:
        plt.figure()
    f = plotFilter
    for objtype in dmags:
        for s in dmags[objtype]:
            plt.plot(color[objtype][s], dmags[objtype][s][f], color=colors[objtype], marker=symbs[objtype])
    plt.xlabel(colorname)
    plt.ylabel(r'$\Delta$%s (mmag)' %(f))
    plt.title(titletext)
    handles= []
    for objtype in dmags:
        myline = mlines.Line2D([], [], color=colors[objtype], marker=symbs[objtype], linestyle='', label=objtype)
        handles.append(myline)
    plt.legend(handles=handles, loc=(0.9, .2), fontsize='small', numpoints=1, fancybox=True)
    if savefig:
        if titletext != None:
            plt.savefig('%s.%s' %(titletext, figformat), format=figformat)
        else:
            plt.savefig('Dmag.%s' %(figformat), format=figformat)
    return
