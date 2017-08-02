import os, re
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import lsst.utils
from lsst.sims.photUtils import Bandpass

# Input components:

# $SYSENG_THROUGHPUTS_DIR / components
# In each component subdirectory, there is a *_Losses subdirectory containing the loss files
#  (all to be combined together)
# In some component subdirectories, there is a *_Coatings subdirecotry containing the coating
#  information for the throughput curve, all to be combined together
# In each component subdirectory, there is a file (or files in the filters component)
#   containing the throughput response. The name of this file varies. For the glass for
#   lens components, the glass throughput curve must be smoothed by a savitzky_golay function.

belowZeroThreshhold = -1.0e-15

def setDefaultDirs(rootDir=None):
    """
    Returns a dictionary with the default directory locations of each component of the system throughput.
    The default value for each component will mirror the expected values in the syseng_throughputs repository, with
     the defaultDirs['detector'] value pointing to the directory common to all vendors (thus would build a generic detector using the minimum values).
    """
    # Set SYSENG_THROUGHPUTS_DIR to the root dir for syseng_throughputs ('setup syseng_throughputs' will do this automatically)
    defaultDirs = {}
    if rootDir is None:
        rootDir = lsst.utils.getPackageDir('syseng_throughputs')
    defaultDirs['detector'] = os.path.join(rootDir, 'components/camera/detector')
    for lens in ('lens1', 'lens2', 'lens3'):
        defaultDirs[lens] = os.path.join(rootDir, 'components/camera', lens)
    defaultDirs['filters'] = os.path.join(rootDir, 'components/camera/filters')
    for mirror in ('mirror1', 'mirror2', 'mirror3'):
        defaultDirs[mirror] = os.path.join(rootDir, 'components/telescope', mirror)
    defaultDirs['atmosphere'] = os.path.join(rootDir, 'siteProperties')
    return defaultDirs

def _readLosses(componentDir):
    """
    Read and combine the losses in all files from a _Losses subdirectory in componentDir.
    Return a bandpass object with all losses combined.
    """
    lossDir = glob(os.path.join(componentDir, '*_Losses'))
    if len(lossDir) > 1:
        errmsg = 'Expect a single *_Losses subdirectory for component %s.'%(componentDir)
        errmsg += ' Found %s.' %(lossDir)
        raise ValueError(errmsg)
    lossDir = lossDir[0]
    if not os.path.isdir(lossDir):
        errmsg = 'Expect %s to be a subdirectory containing loss files for component %s.' \
          %(lossDir, componentDir)
        raise ValueError(errmsg)
    lossfiles = glob(os.path.join(lossDir, '*.dat'))
    if len(lossfiles) == 0:
        errmsg = 'Expect to find at least one loss file in %s for component %s.'\
          %(lossDir, componentDir)
        errmsg += ' Found no loss files.'
        raise ValueError(errmsg)
    loss = Bandpass()
    loss.readThroughputList(lossfiles)
    return loss

def _readCoatings(componentDir):
    """
    Read and combine the coatings in all files from a _Coatings subdirectory in componentDir.
    Return a bandpass object with all coating throughputs combined.
    """
    coatingDir = glob(os.path.join(componentDir, '*_Coatings'))
    if len(coatingDir) > 1:
        errmsg = 'Expect a single *_Coatings subdirectory for component %s.'%(componentDir)
        errmsg += ' Found %s.' %(coatingDir)
        raise ValueError(errmsg)
    coatingDir = coatingDir[0]
    if not os.path.isdir(coatingDir):
        errmsg = 'Expect %s to be a subdirectory containing coating files for component %s.' \
          %(coatingDir, componentDir)
        raise ValueError(errmsg)
    coatingfiles = glob(os.path.join(coatingDir, '*.dat'))
    if len(coatingfiles) == 0:
        errmsg = 'Expect to find at least one coating file in %s for component %s.'\
          %(coatingDir, componentDir)
        errmsg += ' Found no coating files.'
        raise ValueError(errmsg)
    coatings = Bandpass()
    coatings.readThroughputList(coatingfiles)
    return coatings


def buildVendorDetector(vendorDir, addLosses=True):
    """
    Builds a detector response from the files in vendorDir, by reading the *_QE.dat
      and *_Losses subdirectory for a single version of the detector.
    Returns a Bandpass object.
    If addLosses is True, the QE curve is multiplied by the losses in the *Losses.dat files.
    If addLosses is False, the QE curve does not have any losses included.
    """
    # Read the QE file.
    qefile = glob(os.path.join(vendorDir, '*_QE.dat'))
    if len(qefile) != 1:
        raise ValueError('Expected a single QE file in this directory, found: ', qefile)
    qefile = qefile[0]
    qe = Bandpass()
    qe.readThroughput(qefile)
    if addLosses:
        loss = _readLosses(vendorDir)
        wavelength, sb = qe.multiplyThroughputs(loss.wavelen, loss.sb)
        qe.setBandpass(wavelength, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(qe.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if qe.sb[belowzero] < belowZeroThreshhold:
        raise ValueError('Found values in QE response significantly below zero.')
    # If they are just small errors in interpolation, set to zero.
    qe.sb[belowzero] = 0
    return qe

def buildDetector(detectorDir, addLosses=True):
    """
    Builds a detector response from 'detectorDir', potentially from multiple vendors.
    Returns a bandpass object.
    Behavior depends on contents of 'detectorDir':
      if 'detectorDir' contains a *_QE.dat (and a *_Losses subdirectory, if addLosses=True)
       it will be treated as an individual vendor.
      if 'detectorDir' does not contain these files, but does contain subdirectories which do,
       all of the subdirectories which contain *_QE.dat and _Losses subdirectories will
       be read and assumed to be individual vendors; the resulting response curve will
       be the *MINIMUM* value of the response at each wavelength.
    In both cases, the *QE.dat and *Losses files will be read and combined using
     bandpassUtils.buildVendorDetector.
    The value of addLosses will be passed to buildVendorDetector - if True, losses are
      included in the returned response curve.
    """
    try:
        # Try to treat this as a single vendor.
        qe = buildVendorDetector(detectorDir, addLosses=addLosses)
    except ValueError:
        # But it wasn't a single vendor, so treat it as multiple vendors and look for minimum.
        # Find vendor subdirectories:
        tmp = os.listdir(detectorDir)
        vendorDirs = []
        for t in tmp:
            if os.path.isdir(os.path.join(detectorDir, t)) and not t.endswith('_Losses'):
                vendorDirs.append(os.path.join(detectorDir, t))
        if len(vendorDirs) == 0:
            errmsg = 'Could not find the files required for a single vendor '\
              'and could not find subdirectories for multiple vendors, in component %s.'\
              %(detectorDir)
            raise ValueError(errmsg)
        # There are subdirectories; we should take the minimum value of all QE curves.
        sbAll = []
        for vendorDir in vendorDirs:
            qe = buildVendorDetector(vendorDir, addLosses=addLosses)
            sbAll.append(qe.sb)
        wavelen = qe.wavelen
        sbMin = (np.array(sbAll)).min(axis=0)
        qe.setBandpass(wavelen, sbMin)
    return qe

def buildFilters(filterDir, addLosses=True):
    """
    Build a filter throughput curve from the files in filterDir.
    Assumes there are files [filtername]-bandResponse.dat, together with a 'filterLosses' subdirectory containing loss files.
    Returns a dictionary (keyed by filter name) of the bandpasses for each filter.
    If addLosses is True, the filter throughput curves are multiplied by the loss files.
    """
    # Read the filter files.
    filterfiles = glob(os.path.join(filterDir, '*_band_Response.dat'))
    filters = {}
    for f in filterfiles:
        fname = os.path.split(f)[1].split('_')[0]
        filters[fname] = Bandpass()
        filters[fname].readThroughput(f)
    if addLosses:
        # Read and multiply the losses.
        loss = _readLosses(filterDir)
        for f in filters:
            wavelen, sb = filters[f].multiplyThroughputs(loss.wavelen, loss.sb)
            filters[f].setBandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    for f in filters:
        belowzero = np.where(filters[f].sb < 0)
        # If there are QE values significantly < 0, raise an exception.
        if filters[f].sb[belowzero] < belowZeroThreshhold:
            raise ValueError('Found values in filter response significantly below zero in %s filter' % f)
        # If they are just small errors in interpolation, set to zero.
        filters[f].sb[belowzero] = 0
    return filters

def savitzky_golay(y, window_size=31, order=3, deriv=0, rate=1):
    """
    Method brought from Chuck Claver's makeLens*.ipynb notebook.
    Smoothes the wavelength response
     of the borosilicate glass and returns a smoothed throughput curve.
    """
    # y = throughput for lenses
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except(ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.asmatrix([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * np.math.factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

def buildLens(lensDir, addLosses=True):
    """
    Build the lens throughput curve from the files in lensDir.
    Returns a bandpass object.
    The coatings for the lens are in *_Coatings, the loss files are in *_Losses.
    The borosilicate glass throughput is in l*_Glass.dat; this file is smoothed using the
     savitzsky_golay function.
    The glass response is multiplied by the coatings and (if addLosses is True),
      also the loss curves.
    """
    lens = Bandpass()
    # Read the glass base file.
    glassfile = glob(os.path.join(lensDir, 'l*_Glass.dat'))
    if len(glassfile) != 1:
        raise ValueError('Expected a single glass file in this directory, found: ', glassfile)
    glassfile = glassfile[0]
    glass = Bandpass()
    glass.readThroughput(glassfile)
    # Smooth the glass response.
    smoothSb = savitzky_golay(glass.sb, 31, 3)
    lens = Bandpass()
    lens.setBandpass(glass.wavelen, smoothSb)
    # Read the broad band antireflective (BBAR) coatings files.
    bbars = _readCoatings(lensDir)
    # Multiply the bbars by the glass.
    wavelen, sb = lens.multiplyThroughputs(bbars.wavelen, bbars.sb)
    lens.setBandpass(wavelen, sb)
    # Add losses.
    if addLosses:
        loss = _readLosses(lensDir)
        wavelen, sb = lens.multiplyThroughputs(loss.wavelen, loss.sb)
        lens.setBandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(lens.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if lens.sb[belowzero] < belowZeroThreshhold:
        raise ValueError('Found values in lens throughput significantly below zero.')
    # If they are just small errors in interpolation, set to zero.
    lens.sb[belowzero] = 0
    return lens

def buildMirror(mirrorDir, addLosses=True):
    """
    Build a mirror throughput curve.
    Assumes there are *Losses.dat subdirectory with loss files
       and a m*_Ideal.dat file with the mirror throughput.
    Returns a bandpass object.
    If addLosses is True, the *_Ideal.dat file is multiplied by the *_Losses/*.dat files.
    """
    # Read the mirror reflectance curve.
    mirrorfile = glob(os.path.join(mirrorDir, 'm*Ideal.dat'))
    if len(mirrorfile) != 1:
        raise ValueError('Expected a single mirror file in directory %s, found: ' %mirrorDir, mirrorfile)
    mirrorfile = mirrorfile[0]
    mirror = Bandpass()
    mirror.readThroughput(mirrorfile)
    if addLosses:
        loss = _readLosses(mirrorDir)
        wavelen, sb = mirror.multiplyThroughputs(loss.wavelen, loss.sb)
        mirror.setBandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(mirror.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if mirror.sb[belowzero] < belowZeroThreshhold:
        raise ValueError('Found values in mirror response significantly below zero')
    # If they are just small errors in interpolation, set to zero.
    mirror.sb[belowzero] = 0
    return mirror

def readAtmosphere(atmosDir, atmosFile='pachonModtranAtm_12.dat'):
    """
    Read an atmosphere throughput curve, from the default location 'atmosDir'
     and default filename 'pachonModtranAtm_12.dat'
    Returns a bandpass object.
    """
    atmofile = os.path.join(atmosDir, atmosFile)
    atmo = Bandpass()
    atmo.readThroughput(atmofile)
    # Verify that no values go significantly below zero.
    belowzero = np.where(atmo.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if atmo.sb[belowzero] < belowZeroThreshhold:
        raise ValueError('Found values in atmospheric transmission significantly below zero')
    # If they are just small errors in interpolation, set to zero.
    atmo.sb[belowzero] = 0
    return atmo

def buildHardwareAndSystem(defaultDirs, addLosses=True, atmosphereOverride=None):
    """
    Using directories for each component set by 'defaultDirs',
     builds the system (including atmosphere) and hardware throughput curves for each filter.
     Allows optional override of the default atmosphere by one passed in atmosphereOverride (a bandpass object).
    Returns dictionaries of the hardware and system in bandpass objects, keyed per filtername.

    defaultDirs is a dictionary containing the directories of each throughput component:
      it can be built automatically using the setDefaultDirs function.
    For each component, the relevant method above is called - to build the 'lens1' component, buildLens is called, etc.
    The detector is built using 'buildDetector':
      if the directory defaultDirs['detector'] points to multiple subdirectories
       (one per vendor) then the detector response curve corresponds to the minimum
       of all curves at each wavelength.
      if defaultsDirs['detector'] points to a single vendor directory
        (with a *QE.dat curve present), then this is treated a single vendor directory.
    The value of addLosses is passed to each component as it is built (i.e. losses will be multiplied into each throughput curve).
    """
    # Build each component.
    detector = buildDetector(defaultDirs['detector'], addLosses)
    lens1 = buildLens(defaultDirs['lens1'], addLosses)
    lens2 = buildLens(defaultDirs['lens2'], addLosses)
    lens3 = buildLens(defaultDirs['lens3'], addLosses)
    filters = buildFilters(defaultDirs['filters'], addLosses)
    mirror1 = buildMirror(defaultDirs['mirror1'], addLosses)
    mirror2 = buildMirror(defaultDirs['mirror2'], addLosses)
    mirror3 = buildMirror(defaultDirs['mirror3'], addLosses)
    if atmosphereOverride is None:
        atmosphere = readAtmosphere(defaultDirs['atmosphere'])
    else:
        atmosphere = atmosphereOverride
        if np.all(atmosphere.wavelen != detector.wavelen):
            atmosphere.resampleBandpass(wavelen_min = detector.wavelen_min,
                                        wavelen_max = detector.wavelen_max,
                                        wavelen_step = detector.wavelen_step)
    # Combine the individual components.
    # Note that the process of reading in the files above would have put them onto the same wavelength grid.
    core_sb = (detector.sb * lens1.sb * lens2.sb * lens3.sb
               * mirror1.sb * mirror2.sb * mirror3.sb)
    wavelen = detector.wavelen
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        hw_sb = core_sb * filters[f].sb
        hardware[f].setBandpass(wavelen, hw_sb)
        system[f].setBandpass(wavelen, hw_sb*atmosphere.sb)
    return hardware, system

def plotBandpasses(bandpassDict, title=None, newfig=True, savefig=False, addlegend=True,
                   linestyle='-', linewidth=2):
    """
    Plot the bandpass throughput curves, provided as a dictionary in bandpassDict.

    title = plot title
    newfig = (True/False) start a new figure or use the active matplotlib figure
    savefig = (True/False) save the figure to disk with default name title.png or throughputs.png if title not defined
    addLegend = (True/False) add a legend to the plot
    linestyle = matplotlib linestyle (default '-') for the throughput curve lines
    linewidth = matplotlib linewidth (default 2) for the throughput curve lines.
    """
    # Generate a new figure, if desired.
    if newfig:
        plt.figure()
    # Plot the bandpass curves.  Try to sort by filter name ugrizy if those are the keys.
    names = bandpassDict.keys()
    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    if set(names).issubset(set(filterlist)):
        newnames = []
        for f in filterlist:
            if f in names:
                newnames.append(f)
        names = newnames
    for f in names:
        plt.plot(bandpassDict[f].wavelen, bandpassDict[f].sb, marker="", linestyle=linestyle,
                   linewidth=linewidth, label=f)
    # Only draw the legend if desired (many bandpassDicts plotted together could make the legend unwieldy).
    if addlegend:
        plt.legend(loc='lower right', numpoints=1, fancybox=True, fontsize='smaller')
    # Limit wavelengths to the LSST range.
    plt.xlim(300, 1150)
    plt.ylim(0, 1)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fractional Throughput Response')
    # Only add the grid if it's a new figure (otherwise, it toggles on/off).
    plt.grid(True)
    # Add a plot title.
    if title != None:
        plt.title(title)
    # Save the figure, if desired.
    if savefig:
        figformat = 'png'
        if title is not None:
            plt.savefig('%s.%s' %(title, figformat), format=figformat)
        else:
            plt.savefig('throughputs.%s' %(figformat), format=figformat)
    return
