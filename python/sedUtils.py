import os
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lsst.sims.photUtils import Sed


def readPhotSeds(sedDir=None):
    """Read all the seds provided by this package, storing them in a
    dictionary and saving lists of the SEDs of each type (so that
    they can be separated later if desired). """
    # The environment variable PHOT_REF_SEDS_DIR is set if the file
    #  'phot_ref.csh' is sourced.
    if sedDir == None:
        sedDir = os.path.join(os.getenv('SYSENG_THROUGHPUTS_DIR'), 'seds')
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
        # And for each type of object, loop through and read each SED.
        for s in sedlists[objtype]:
            sedDict[s] = Sed()
            sedDict[s].readSED_flambda(os.path.join(sedDir, objtype, s))
    # Return the sed dictionary and the dictionary containing the lists of each
    #  type of object.
    return sedDict, sedlists

def readAnySeds(inputfileList, sedDir=None):
    """Read the seds in a list and store them in a format appropriate for use with the
    rest of these routines. sedDir (if set) can be the root directory for the files in the list."""
    # Set up the sedlists dictionary to hold the sed names. We'll store them all keyed under 'any'.
    sedlists = {}
    sedlists['any'] = deepcopy(inputfileList)
    # If the root directory is set, add it to the input file names.
    if sedDir != None:
        ifiles = []
        for i in inputfileList:
            ifiles.append(os.path.join(sedDir, i))
        inputfileList = ifiles
    # Read the files.
    sedDict = {}
    for filename, s in zip(inputfileList, sedlists['any']):
        sedDict[s] = Sed()
        sedDict[s].readSED_flambda(filename)
    # Return the sed dictionary and the dictionary containing the lists of each type of object.
    return sedDict, sedlists


def makeRedshiftedSeds(sedDict, sedlists, redshifts):
    """Redshift the quasar, galaxies and SN by the amounts given in redshifts.
    If redshift is a numpy array, all objects except stars are redshifted to the same values.
    If redshift is a dictionary containing numpy arrays in the 'quasar', 'sn' or 'galaxies' keys, then
    those redshifts are applied to those objects only. """
    # Check what kind of redshift object we received.
    if isinstance(redshifts, dict):
        # Loop over all the object types in the redshift dictionary.
        for objtype in redshifts.keys():
            # Then loop over all redshifts for this type.
            newSedsForList = []
            for z in redshifts[objtype]:
                # And add new SEDs redshifted to this z for this type of object.
                for s in sedlists[objtype]:
                    # Make the new name for the new SED at this redshift.
                    newsedname = s + '_Z_%.3f' %(z)
                    # Add it to the SED dictionary.
                    sedDict[newsedname] = redshiftSingleSED(sedDict[s], z)
                    newSedsForList.append(newsedname)
            # Done redshifting all objects of this type; updated sedlist.
            sedlists[objtype] += newSedsForList
    else:
        # Using the same redshift list for everything (except stars).
        for z in redshifts:
            for objtype in sedlists.keys():
                newSedsForList = []
                if objtype == 'stars':
                    # Skip stars.
                    continue
                for s in sedlists[objtype]:
                    newsedname = s + '_Z_%.3f' %(z)
                    sedDicts[newsedname] = redshiftSingleSED(sedDict[s], z)
                    newSedsForList.append(newsedname)
                sedlists[objtype] += newSedsForList
    return sedDict, sedlists


def redshiftSingleSED(sed_in, z):
    # Make a copy of the input SED, as otherwise we will overwrite original arrays.
    sed_out = deepcopy(sed_in)
    sed_out.redshiftSED(z, dimming=False)
    return sed_out


def matchSedsBp(sedDict, bpDict, refFilter=None):
    """Match the wavelength ranges for all the Seds and the bandpass dictionary.
    This will speed up later calculation of magnitudes (if you're doing a lot of them),
     but is not strictly necessary. """
    if refFilter == None:
        refFilter = bpDict.keys()[0]
    wavelen_match = bpDict[refFilter].wavelen
    # Check all filters in bpDict match in wavelength space (note bandpasses must be regular grid).
    for f in bpDict.keys():
        if numpy.any(bpDict[f].wavelen != wavelen_match):
            bpDict[f].resampleBandpass(wavelen_min=wavelen_match.min(), wavelen_max=wavelen_match.max(),
                                       wavelen_step = wavelen_match[1] - wavelen_match[0])
    # Check all seds in sedDict match in wavelength space.
    for s in sedDict.keys():
        if sedDict[s].needResample(wavelen_match=wavelen_match):
            sedDict[s].resampleSED(wavelen_match=wavelen_match)
    return sedDict, bpDict


def calcNatMags(bandpassDict, sedDict, sedlists):
    """Calculate (zeropoint-calibrated) natural magnitudes for each SED,
    in all bandpasses. Changes in this magnitude includes only color/wavelength-dependent effects. """
    # Create a dictionary to hold all of the magnitude information, for all filters.
    mags = {}
    # Loop over SEDs:
    for o in sedlists.keys():
        for s in sedlists[o]:
            # Create a dictionary for each object to hold the multiple filter information.
            mags[s] = {}
            for f in bandpassDict.keys():
                # Calculate the magnitudes.
                mags[s][f] = sedDict[s].calcMag(bandpassDict[f])
    return mags

def calcInstMags(bandpassDict, sedDict, sedlists, photParams):
    """Calculate instrumental magnitudes for each SED, in all bandpasses.
    Changes in this magnitude includes gray-scale effects as well as color/wavelength dependent effects. """
    mags = {}
    for o in sedlists.keys():
        for s in sedlists[o]:
            mags[s] = {}
            for f in bandpassDict.keys():
                mags[s][f] = sedDict[s].calcADU(bandpassDict[f], photParams)
                mags[s][f] = -2.5*numpy.log10(mags[s][f])
    return mags


def calcDeltaMags(mags1, mags2, mmags=True, matchBlue=False):
    """Calculate the difference in magnitudes between two sets of magnitudes, mags calculated as above.
    If 'mmags' is True, then returns delta mags in mmags. 
    If 'matchBlue' is True, then scales change in magnitude between mags1 and mags2 so that the blue 
     star 'km10_7250.fits_g45' has zero magnitude change."""
    dmags_seds = list(set(mags1) & set(mags2))
    s = mags1.keys()[0]
    dmags_filters = list(set(mags1[s]) & set(mags2[s]))
    dmags= {}
    for s in dmags_seds:
        dmags[s] = {}
        for f in dmags_filters:
            dmags[s][f] = mags1[s][f] - mags2[s][f]
            # Convert to millimags if desired.
            if mmags:
                dmags[s][f] *= 1000.0
    # Apply scaling so that bluest star remains constant, if desired. 
    if matchBlue:
        # Calculate offset to apply.
        offset = {}
        for f in dmags_filters:
            offset[f] = dmags['km10_7250.fits_g45'][f]
        # Apply offset.
        for s in dmags_seds:
            for f in dmags_filters:
                dmags[s][f] = dmags[s][f] - offset[f]
    return dmags

def calcGiColors(mags):
    """Calculate the g-i colors of objects in a set of magnitudes, mags calculated as above."""
    gi = {}
    for s in mags.keys():
        gi[s] = mags[s]['g'] - mags[s]['i']
    return gi

def calcAnyColor(mags, color1, color2):
    """Calculate any color of objects in a set of magnitudes. Specify color1 (e.g. 'g') and color2 (e.g. 'i')."""
    color = {}
    for s in mags.keys():
        color[s] = mags[s][color1] - mags[s][color2]
    return color

def printDmags(sedlists, dmags, filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Print changes in magnitudes to the screen."""
    print 'Delta mmag:'
    writestring = "object"
    for f in filterlist:
        writestring += '\t %s ' %(f)
    print writestring
    for objtype in sedlists.keys():
        print 'Object type: ', objtype
        for s in sedlists[objtype]:
            writestring = 'dm %s ' %(s)
            for f in filterlist:
                writestring += ' %f ' %(dmags[s][f])
            print writestring
    return


def plotDmags(sedlists, color, dmags, colorname = 'g-i', xlim=None,
              newfig=True, titletext=None, savefig=False,
              filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Generate a plot of the change in magnitudes as a function of color. """
    symbs = {'quasar':'o', 'stars':'*', 'white_dwarf':'*', 'sn':'x', 'galaxies':'s', 'photoZ_outliers':'+', 'any':'^'}
    colors = {'quasar':'g', 'stars':'r', 'white_dwarf':'m','sn':'b', 'galaxies':'g', 'photoZ_outliers':'g', 'any':'k'}
    objlist = ['quasar', 'galaxies', 'photoZ_outliers', 'sn', 'stars', 'white_dwarf']
    for objtype in sedlists:
        if objtype not in objlist:
            objlist.append(objtype)
    if newfig:
        plt.figure()
    plt.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    for f, i in zip(filterlist, range(1, len(filterlist)+1)):
        plt.subplot(3,2,i)
        for objtype in objlist:
            for s in sedlists[objtype]:
                plt.plot(color[s], dmags[s][f], color=colors[objtype], marker=symbs[objtype])
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


def plotDmagsSingle(sedlists, color, dmags, colorname='g-i', plotFilter='u', newfig=True, titletext=None, savefig=False):
    """Generate a plot of the change in magnitudes, in a single filter only. """
    symbs = {'quasar':'o', 'stars':'*', 'white_dwarf':'*', 'sn':'x', 'galaxies':'s', 'photoZ_outliers':'+', 'any':'^'}
    colors = {'quasar':'g', 'stars':'r', 'white_dwarf':'m', 'sn':'b', 'galaxies':'g', 'photoZ_outliers':'g', 'any':'k'}
    objlist = ['quasar', 'galaxies', 'photoZ_outliers', 'sn', 'stars', 'white_dwarf']
    for objtype in sedlists:
        if objtype not in objlist:
            objlist.append(objtype)
    if newfig:
        plt.figure()
    f = plotFilter
    for objtype in objlist:
        for s in sedlists[objtype]:
            plt.plot(color[s], dmags[s][f], color=colors[objtype], marker=symbs[objtype])
    plt.xlabel(colorname)
    plt.ylabel(r'$\Delta$%s (mmag)' %(f))
    plt.title(titletext)
    handles= []
    for s in objlist:
        myline = mlines.Line2D([], [], color=colors[s], marker=symbs[s], linestyle='', label=s)
        handles.append(myline)
    plt.legend(handles=handles, loc=(0.9, .2), fontsize='small', numpoints=1, fancybox=True)
    if savefig:
        if titletext != None:
            plt.savefig('%s.%s' %(titletext, figformat), format=figformat)
        else:
            plt.savefig('Dmag.%s' %(figformat), format=figformat)
    return
