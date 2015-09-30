import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import lsst.utils
from lsst.sims.photUtils import Bandpass

# Input components:

# $SYSENG_THROUGHPUTS_DIR / components
#   camera /
#      detector /   vendor[1,2] / ven[1,2]QE.dat, ven[1,2]Losses.dat  (multiply for a particular vendor, use min between vendors for 'average')
#      lens[1,2,3] / l[1,2,3]Glass.dat, l[1,2,3]SBBAR.dat, l[1,2,3]S2BBAR.dat  (pass Glass through savitzky_golay and multiply by BBARs)
#              l[1,2,3]Losses   / l[1,2,3]S(1,2)(Condensation,Contamination).dat   (multiply all)
#      filters / [u,g,r,i,z,y]-bandResponse.dat
#              filterLosses/    filterS(1,2)(Condensation,Contamination).dat (multiply all, for all filters)
#   telescope /
#       mirror[1,2,3] / m[1,2,3]Losses.dat, m[1,2,3]ProtAlIdeal.dat   (multiply ideal and loss for each mirror, multiply all mirrors)

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


def buildVendorDetector(vendorDir, addLosses=True):
    """
    Assumes there is a file *QE.dat and *Losses.dat to represent the QE response and the Losses.
    Returns a Bandpass object.
    If addLosses is Frue, the QE curve is multiplied by the losses in the *Losses.dat files.
    If addLosses is False, the QE curve does not have any losses included.
    """
    # Read the QE file.
    qefile = glob(os.path.join(vendorDir, '*QE.dat'))
    if len(qefile) != 1:
        raise ValueError('Expected a single QE file in this directory, found: ', qefile)
    qefile = qefile[0]
    qe = Bandpass()
    qe.readThroughput(qefile)
    if addLosses:
        # Read the Losses file(s). If more than one, multiply.
        lossfiles = glob(os.path.join(vendorDir, '*Losses.dat'))
        loss = Bandpass()
        loss.readThroughputList(lossfiles)
        wavelength, sb = qe.multiplyThroughputs(loss.wavelen, loss.sb)
        qe.setBandpass(wavelength, sb)
    return qe

def buildGenericDetector(detectorDir, addLosses=True):
    """
    Combine detector throughputs from multiple vendors, using the minimum value of QE*losses at each wavelength.
    Returns a bandpass object.
    If 'detectorDir' contains subdirectories, then this method will read the vendor curves in each subdirectory
      (calling buildVendorDetector) for each directory and combine the results so that the resulting response is the MINIMUM
      value at each wavelength.
    If 'detectorDir' is an actual vendor subdirectory containing a *QE.dat throughput curve and no subdirectories, then
      buildVendorDetector will be called directly and the bandpass object returned will corresond to the actual value from that vendor.
    The value of addLosses is passed to buildVendorDetector.
    """
    tmp = os.listdir(detectorDir)
    vendorDirs = []
    for t in tmp:
        if os.path.isdir(os.path.join(detectorDir, t)):
            vendorDirs.append(os.path.join(detectorDir, t))
    if len(vendorDirs) == 0:
        # If there were no subdirectories, this could be the actual vendor.
        # Provide a shortcut so that "buildHardwareAndSystem" can operate on a per-vendor basis.
        qefile = glob(os.path.join(detectorDir, '*QE.dat'))
        if len(qefile) == 1:
            print 'Using %s as the root directory for a single vendor' %(detectorDir)
            qe = buildVendorDetector(detectorDir, addLosses=addLosses)
    else:
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
    Assumes there are files [filtername]-bandResponse.dat, together with a 'filterLosses' subdirectory containing loss files.
    Returns a dictionary (keyed by filter name) of the bandpasses for each filter.
    If addLosses is True, the filter throughput curves are multiplied by the loss files. 
    """
    # Read the filter files.
    filterfiles = glob(os.path.join(filterDir, '*-bandResponse.dat'))
    filters = {}
    for f in filterfiles:
        fname = os.path.split(f)[1].split('-')[0]
        filters[fname] = Bandpass()
        filters[fname].readThroughput(f)
    if addLosses:
        # Read and multiply the losses.
        lossfiles = glob(os.path.join(filterDir, 'filterLosses', '*.dat'))
        loss = Bandpass()
        loss.readThroughputList(lossfiles)
        for f in filters:
            wavelen, sb = filters[f].multiplyThroughputs(loss.wavelen, loss.sb)
            filters[f].setBandpass(wavelen, sb)
    return filters

def savitzky_golay(y, window_size=31, order=3, deriv=0, rate=1):
    """
    Method brought from Chuck Claver's makeLens*.ipynb notebook. Smoothes the wavelength response
    of the borosilicate glass and returns a smoothed throughput curve.
    """
    # y = throughput for lenses
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
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
    Build each individual lens. Assumes there are files l*Glass.dat, l*BBAR.dat, and a subdirectory with loss files called 'l*Losses'.
    Returns a bandpass object.
    The l*Glass throughput file is smoothed using the savitzsky_golay function, and then multiplied by the two l*BBAR coatings.
    If addLosses is True, the throughput curve is then multiplied by the loss curves in l*Losses.
    """
    lens = Bandpass()
    # Read the glass base file.
    glassfile = glob(os.path.join(lensDir, 'l*Glass.dat'))
    if len(glassfile) != 1:
        raise ValueError('Expected a single glass file in this directory, found: ', glassfile)
    glassfile = glassfile[0]
    glass = Bandpass()
    glass.readThroughput(glassfile)
    # Smooth the glass response.
    smoothGlassSb = savitzky_golay(glass.sb, 31, 3)
    lens.setBandpass(glass.wavelen, smoothGlassSb)
    # Read the broad band antireflective (BBAR) coatings files.
    bbarfiles = glob(os.path.join(lensDir, 'l*BBAR.dat'))
    bbars = Bandpass()
    bbars.readThroughputList(bbarfiles)
    # Multiply the bbars by the glass.
    wavelen, sb = lens.multiplyThroughputs(bbars.wavelen, bbars.sb)
    lens.setBandpass(wavelen, sb)
    # Add losses.
    if addLosses:
        lossfiles = glob(os.path.join(lensDir, 'l*Losses', '*.dat'))
        loss = Bandpass()
        loss.readThroughputList(lossfiles)
        wavelen, sb = lens.multiplyThroughputs(loss.wavelen, loss.sb)
        lens.setBandpass(wavelen, sb)
    return lens

def buildMirror(mirrorDir, addLosses=True):
    """
    Build a mirror throughput curve. Assumes there are *Losses.dat file(s) and a *Ideal.dat file.
    Returns a bandpass object.
    If addLosses is True, the *Ideal.dat file is multiplied by the *Losses.dat files.
    """
    # Read the mirror reflectance curve.
    mirrorfile = glob(os.path.join(mirrorDir, 'm*Ideal.dat'))
    if len(mirrorfile) != 1:
        raise ValueError('Expected a single mirror file in directory %s, found: ' %mirrorDir, mirrorfile)
    mirrorfile = mirrorfile[0]
    mirror = Bandpass()
    mirror.readThroughput(mirrorfile)
    if addLosses:
        lossfiles = glob(os.path.join(mirrorDir, 'm*Losses.dat'))
        loss = Bandpass()
        loss.readThroughputList(lossfiles)
        wavelen, sb = mirror.multiplyThroughputs(loss.wavelen, loss.sb)
        mirror.setBandpass(wavelen, sb)
    return mirror

def buildAtmosphere(atmosDir):
    """
    Read the atmosphere throughput curve.
    Returns a bandpass object.
    """
    atmofile = os.path.join(atmosDir, 'pachonModtranAtm_12.dat')
    atmo = Bandpass()
    atmo.readThroughput(atmofile)
    return atmo

def buildHardwareAndSystem(defaultDirs, addLosses=True):
    """
    Go through and build all of the default files into a system (includes atmosphere) and a hardware only throughput.
    Returns dictionaries of the hardware and system (including atmosphere) in bandpass objects, keyed per filtername.

    defaultDirs is a dictionary containing the directories of each throughput component: it can be built automatically using
     the setDefaultDirs function.
    For each component, the relevant method above is called - to build the 'lens1' component, buildLens is called, etc.
    The detector is built using the conservative 'buildGenericDetector': if the directory defaultDirs['detector'] points to
     multiple subdirectories (one per vendor) then the detector response curve corresponds to the minimum of all curves at each wavelength.
     If defaultsDirs['detector'] points to a single vendor directory (with a *QE.dat curve present), then this is treated a single vendor directory.
    The value of addLosses is passed to each component as it is built (i.e. losses will be multiplied into each throughput curve).
    """
    # Build each component.
    detector = buildGenericDetector(defaultDirs['detector'], addLosses)
    lens1 = buildLens(defaultDirs['lens1'], addLosses)
    lens2 = buildLens(defaultDirs['lens2'], addLosses)
    lens3 = buildLens(defaultDirs['lens3'], addLosses)
    filters = buildFilters(defaultDirs['filters'], addLosses)
    mirror1 = buildMirror(defaultDirs['mirror1'], addLosses)
    mirror2 = buildMirror(defaultDirs['mirror2'], addLosses)
    mirror3 = buildMirror(defaultDirs['mirror3'], addLosses)
    atmosphere = buildAtmosphere(defaultDirs['atmosphere'])
    # Combine the individual components.
    # Note that the process of reading in the files above would have put them onto the same wavelength grid.
    core_sb = detector.sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        wavelen = filters[f].wavelen
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
        plt.legend(loc=(0.95, 0.7), numpoints=1, fancybox=True, fontsize='smaller')
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
