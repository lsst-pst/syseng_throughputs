import math
import os
import warnings
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
from rubin_sim.phot_utils import Bandpass

__all__ = [
    "setDefaultDirs",
    "buildVendorDetector",
    "buildDetector",
    "buildFilters",
    "savitzky_golay",
    "buildLens",
    "buildMirror",
    "readAtmosphere",
    "buildHardwareAndSystem",
    "plotBandpasses",
]

"""This module combines the individual base-level components to build a total throughput curve.

This includes the telescope (M1/M2/M3), with their appropriate loss curves, 
the optics (L1/L2/L3) (calculated from the glass + coatings file) with their appropriate loss curves, 
the filters (ugrizy), 
and the QE of the detector, with appropriate loss curve 
(QE can be computed as a minimum between vendors, or each vendor individually). 
The losses for each component can be chosen to be added or not. 
"""

# Input components:

# $SYSENG_THROUGHPUTS_DIR / components
# In each component subdirectory, there is a *_Losses subdirectory containing the loss files
#  (all to be combined together)
# In some component subdirectories, there is a *_Coatings subdirecotry containing the coating
#  information for the throughput curve, all to be combined together
# In each component subdirectory, there is a file (or files in the filters component)
#   containing the throughput response. The name of this file varies. For the glass for
#   lens components, the glass throughput curve must be smoothed by a savitzky_golay function.

belowZeroThreshhold = -0.02
filterlist = ("u", "g", "r", "i", "z", "y")
filtercolors = {"u": "b", "g": "c", "r": "g", "i": "y", "z": "r", "y": "m"}
WAVELEN_MIN = 300
WAVELEN_MAX = 1150
WAVELEN_STEP = 0.1


def findRootDir(rootDir=None):
    """Find the location of the syseng_throughputs data.
    First this looks to see if the rootDir was simply passed a kwarg.
    If not, then it looks to see if SYSENG_THROUGHPUTS_DIR was set as an environment variable.
    If not, then it looks for the location of the installed syseng_throughputs package.
    """
    if rootDir is None:
        rootDir = os.getenv("SYSENG_THROUGHPUTS_DIR")
        if rootDir is None:
            rootDir = os.path.split(pkg_resources.resource_filename("syseng_throughputs", "."))[0]
            # Remove the last syseng_throughputs where the python lives and go to the 'components' level
            rootDir = os.path.split(rootDir)[0]
    return rootDir


def setDefaultDirs(rootDir=None):
    """
    Returns a dictionary with the default directory locations of each component of the system throughput.
    The default value for each component will mirror the expected values in the syseng_throughputs repository,
    with the defaultDirs['detector'] value pointing to the directory common to all vendors
    (thus would build a generic detector using the minimum values).

    Parameters
    ----------
    rootDir : 'str', opt
        Path to top level ('syseng_throughputs' root directory) of throughput data.
        Default None - uses 'SYSENG_THROUGHPUTS_DIR' environment variable or the location of this package.

    Returns
    -------
    defaultDirs : `dict`
        Dictionary with keys = component name, value = default root directory for that component.
    """
    defaultDirs = {}
    rootDir = findRootDir(rootDir)
    defaultDirs["detector"] = os.path.join(rootDir, "components/camera/detector/joint_minimum")
    for lens in ("lens1", "lens2", "lens3"):
        defaultDirs[lens] = os.path.join(rootDir, "components/camera", lens)
    defaultDirs["filters"] = os.path.join(rootDir, "components/camera/filters")
    for mirror in ("mirror1", "mirror2", "mirror3"):
        defaultDirs[mirror] = os.path.join(rootDir, "components/telescope", mirror)
    defaultDirs["atmosphere"] = os.path.join(rootDir, "siteProperties")
    return defaultDirs


def _readLosses(componentDir):
    """
    Read and combine the losses in all files from a _Losses subdirectory in componentDir.

    Parameters
    ----------
    componentDir : `str`
        The directory containing the _Losses subdirectory.

    Returns
    -------
    loss : `Bandpass`
        A rubin_sim.photUtils.Bandpass object with all losses combined.
    """
    lossDir = glob(os.path.join(componentDir, "*_Losses"))
    if len(lossDir) > 1:
        errmsg = "Expect a single *_Losses subdirectory for component %s." % (componentDir)
        errmsg += " Found %s." % (lossDir)
        raise ValueError(errmsg)
    lossDir = lossDir[0]
    if not os.path.isdir(lossDir):
        errmsg = "Expect %s to be a subdirectory containing loss files for component %s." % (
            lossDir,
            componentDir,
        )
        raise ValueError(errmsg)
    lossfiles = glob(os.path.join(lossDir, "*.dat"))
    if len(lossfiles) == 0:
        errmsg = "Expect to find at least one loss file in %s for component %s." % (
            lossDir,
            componentDir,
        )
        errmsg += " Found no loss files."
        raise ValueError(errmsg)
    loss = Bandpass()
    loss.read_throughput_list(
        lossfiles,
        wavelen_min=WAVELEN_MIN,
        wavelen_max=WAVELEN_MAX,
        wavelen_step=WAVELEN_STEP,
    )
    return loss


def _readCoatings(componentDir):
    """
    Read and combine the coatings in all files from a _Coatings subdirectory in componentDir.

    Parameters
    ----------
    componentDir : `str`
        The directory containing the *_Coatings subdirectory.

    Returns
    -------
    coatings : `Bandpass`
        A bandpass object with all coating throughputs combined.
    """
    coatingDir = glob(os.path.join(componentDir, "*_Coatings"))
    if len(coatingDir) > 1:
        errmsg = "Expect a single *_Coatings subdirectory for component %s." % (componentDir)
        errmsg += " Found %s." % (coatingDir)
        raise ValueError(errmsg)
    coatingDir = coatingDir[0]
    if not os.path.isdir(coatingDir):
        errmsg = "Expect %s to be a subdirectory containing coating files for component %s." % (
            coatingDir,
            componentDir,
        )
        raise ValueError(errmsg)
    coatingfiles = glob(os.path.join(coatingDir, "*.dat"))
    if len(coatingfiles) == 0:
        errmsg = "Expect to find at least one coating file in %s for component %s." % (
            coatingDir,
            componentDir,
        )
        errmsg += " Found no coating files."
        raise ValueError(errmsg)
    coatings = Bandpass()
    coatings.read_throughput_list(
        coatingfiles,
        wavelen_min=WAVELEN_MIN,
        wavelen_max=WAVELEN_MAX,
        wavelen_step=WAVELEN_STEP,
    )
    return coatings


def buildVendorDetector(vendorDir, addLosses=True):
    """
    Builds a detector response from the files in vendorDir, by reading the *_QE.dat
    and *_Losses subdirectory for a single version of the detector.

    Parameters
    ----------
    vendorDir : `str`
        The directory from which to read a particular vendor QE curve (_QE.dat)
    addLosses : `bool`, opt
        Flag to determine whether to read and then combine the losses from the *_Losses.dat files with the QE
        Default True.

    Returns
    -------
    qe : `Bandpass`
        A bandpass object with the QE curve, including losses if addLosses=True.
    """
    # Read the QE file.
    qefile = glob(os.path.join(vendorDir, "*_QE.dat"))
    if len(qefile) != 1:
        raise ValueError("Expected a single QE file in this directory, found: ", qefile)
    qefile = qefile[0]
    qe = Bandpass()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        qe.read_throughput(qefile)
    qe.resample_bandpass(wavelen_min=WAVELEN_MIN, wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP)
    if addLosses:
        loss = _readLosses(vendorDir)
        wavelength, sb = qe.multiply_throughputs(loss.wavelen, loss.sb)
        qe.set_bandpass(wavelength, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(qe.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if np.any(qe.sb[belowzero] < belowZeroThreshhold):
        raise ValueError("Found values in QE response significantly below zero.")
    # If they are just small errors in interpolation, set to zero.
    qe.sb[belowzero] = 0
    return qe


def buildDetector(detectorDir, addLosses=True):
    """
    Builds a detector response from 'detectorDir', potentially from multiple vendors.
    Returns a bandpass object.
    Behavior depends on contents of 'detectorDir':
    * if 'detectorDir' contains a *_QE.dat (and a *_Losses subdirectory, if addLosses=True)
    it will be treated as an individual vendor.
    * if 'detectorDir' does not contain these files, but does contain subdirectories which do,
    all of the subdirectories which contain *_QE.dat and _Losses subdirectories will
    be read and assumed to be individual vendors; the resulting response curve will
    be the *MINIMUM* value of the response at each wavelength.

    In both cases, the *QE.dat and *Losses files will be read and combined using
    bandpassUtils.buildVendorDetector.
    The value of addLosses will be passed to buildVendorDetector - if True, losses are
    included in the returned response curve.

    Parameters
    ----------
    detectorDir : `str`
        The directory from which to read either the single vendor or multiple vendor QE curves.
    addLosses : `bool`, opt
        Flag to determine whether to add losses to the final QE curve.
        Default True.

    Returns
    -------
    qe : `Bandpass`
        A bandpass object with the QE curve for the detector (or multiple vendors).
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
            if os.path.isdir(os.path.join(detectorDir, t)) and not t.endswith("_Losses"):
                vendorDirs.append(os.path.join(detectorDir, t))
        if len(vendorDirs) == 0:
            errmsg = (
                "Could not find the files required for a single vendor "
                "and could not find subdirectories for multiple vendors, in component %s." % (detectorDir)
            )
            raise ValueError(errmsg)
        # There are subdirectories; we should take the minimum value of all QE curves.
        sbAll = []
        for vendorDir in vendorDirs:
            qe = buildVendorDetector(vendorDir, addLosses=addLosses)
            sbAll.append(qe.sb)
        wavelen = qe.wavelen
        sbMin = (np.array(sbAll)).min(axis=0)
        qe.set_bandpass(wavelen, sbMin)
    return qe


def buildFilters(filterDir, addLosses=True, shiftFilters=None):
    """
    Build a filter throughput curve from the files in filterDir.
    Assumes there are files [filtername]-band_Response.dat, together with a 'filterLosses' subdirectory
    containing loss files.

    Parameters
    ----------
    filterDir : `str`
        Directory from which to read the [filtername]-band_Response.dat files
    addLosses : `bool`, opt
        Flag to determine whether to add losses to the filter response curves, from the _Losses subdirectory.
        Default True.
    shiftFilters : `float`, opt
        Add a simple offset to the filter wavelength/throughput values.
        The value for shiftFilters is the % of the effective wavelength to shift by.
        Default None, no shift. Typical likely shifts may be 1-3% or so.

    Returns
    -------
    filters : `dict`
        A dictionary {'filter': 'Bandpass`} containing the filter throughputs.
    """
    # Read the filter files.
    filterfiles = glob(os.path.join(filterDir, "*_band_Response.dat"))
    filters = {}
    for f in filterfiles:
        fname = os.path.split(f)[1].split("_")[0]
        filters[fname] = Bandpass()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            filters[fname].read_throughput(f)
        filters[fname].resample_bandpass(
            wavelen_min=WAVELEN_MIN, wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP
        )
    if addLosses:
        # Read and multiply the losses.
        loss = _readLosses(filterDir)
        for f in filters:
            wavelen, sb = filters[f].multiply_throughputs(loss.wavelen, loss.sb)
            filters[f].set_bandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    for f in filters:
        belowzero = np.where(filters[f].sb < 0)
        # If there are QE values significantly < 0, raise an exception.
        if np.any(filters[f].sb[belowzero] < belowZeroThreshhold):
            raise ValueError("Found values in filter response significantly below zero in %s filter" % f)
        # If they are just small errors in interpolation, set to zero.
        filters[f].sb[belowzero] = 0
    if shiftFilters is not None:
        # This is an extremely simple shift. See photoCal / FilterShift.py for more realistic shifts.
        for f in filters:
            effphi, effsb = filters[f].calcEffWavelen()
            shift = shiftFilters / 100.0 * effsb
            filters[f].wavelen = filters[f].wavelen + shift
            filters[f].resample_bandpass()
            filters[f].sbTophi()
    return filters


def savitzky_golay(y, window_size=31, order=3, deriv=0, rate=1):
    """
    Method brought from Chuck Claver's makeLens*.ipynb notebook.
    Smoothes the wavelength response of the borosilicate glass and returns a smoothed throughput curve.
    """
    # y = throughput for lenses
    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except (ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.asmatrix([[k**i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * math.factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1 : half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1 : -1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode="valid")


def buildLens(lensDir, addLosses=True):
    """
    Build the lens throughput curve from the files in lensDir.
    Returns a bandpass object.
    The coatings for the lens are in *_Coatings, the loss files are in *_Losses.
    The borosilicate glass throughput is in l*_Glass.dat;
    the glass throughput values are smoothed using the savitzsky_golay function.
    The glass response is multiplied by the coatings and (if addLosses is True), also the loss curves.

    Parameters
    -----------
    lensDir : `str`
        The directory from which to read the l*_Glass.dat file, along with subdirectories for
        the Coatings and Losses.
    addLosses : `bool`, opt
        Flag to add losses from the _Losses subdirectory.
        Default True.

    Returns
    -------
    lens : `Bandpass`
        A bandpass object with the lens throughput, including coatings and optionally losses.
    """
    # Read the glass base file.
    glassfile = glob(os.path.join(lensDir, "l*_Glass.dat"))
    if len(glassfile) != 1:
        raise ValueError("Expected a single glass file in this directory, found: ", glassfile)
    glassfile = glassfile[0]
    glass = Bandpass()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        glass.read_throughput(glassfile)
    glass.resample_bandpass(wavelen_min=WAVELEN_MIN, wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP)
    # Smooth the glass response.
    smoothSb = savitzky_golay(glass.sb, 31, 3)
    lens = Bandpass()
    lens.set_bandpass(glass.wavelen, smoothSb)
    # Read the broad band antireflective (BBAR) coatings files.
    bbars = _readCoatings(lensDir)
    # Multiply the bbars by the glass.
    wavelen, sb = lens.multiply_throughputs(bbars.wavelen, bbars.sb)
    lens.set_bandpass(wavelen, sb)
    # Add losses.
    if addLosses:
        loss = _readLosses(lensDir)
        wavelen, sb = lens.multiply_throughputs(loss.wavelen, loss.sb)
        lens.set_bandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(lens.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if np.any(lens.sb[belowzero] < belowZeroThreshhold):
        raise ValueError("Found values in lens throughput significantly below zero.")
    # If they are just small errors in interpolation, set to zero.
    lens.sb[belowzero] = 0
    return lens


def buildMirror(mirrorDir, addLosses=True):
    """
    Build a mirror throughput curve.
    Assumes there are *Losses.dat subdirectory with loss files
    and a m*_Ideal.dat file with the mirror throughput.

    Parameters
    ----------
    mirrorDir : `str`
        Path to mirror directory. Must contain a file m*Ideal.dat with mirror throughput.
    addLosses : `bool`, opt
        Flag to add loss contributions from the *_Losses subdirorectires.
        Default True.

    Returns
    -------
    mirror : `Bandpass`
        A bandpass object with the mirror throughput, optionally with added losses.
    """
    # Read the mirror reflectance curve.
    mirrorfile = glob(os.path.join(mirrorDir, "m*.dat"))
    if len(mirrorfile) != 1:
        raise ValueError(
            "Expected a single mirror file in directory %s, found: " % mirrorDir,
            mirrorfile,
        )
    mirrorfile = mirrorfile[0]
    mirror = Bandpass()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mirror.read_throughput(mirrorfile)
    mirror.resample_bandpass(wavelen_min=WAVELEN_MIN, wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP)
    if addLosses:
        loss = _readLosses(mirrorDir)
        wavelen, sb = mirror.multiply_throughputs(loss.wavelen, loss.sb)
        mirror.set_bandpass(wavelen, sb)
    # Verify that no values go significantly below zero.
    belowzero = np.where(mirror.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if np.any(mirror.sb[belowzero] < belowZeroThreshhold):
        raise ValueError("Found values in mirror response significantly below zero")
    # If they are just small errors in interpolation, set to zero.
    mirror.sb[belowzero] = 0
    return mirror


def readAtmosphere(atmosDir, atmosFile="pachonModtranAtm_12_aerosol.dat"):
    """
    Read an atmosphere throughput curve.

    Parameters
    ----------
    atmosDir : `str`
        The directory containing the atmosphere throughput file.
    atmosFile : `str`, opt
        The filename of the atmospheric throughput curve.
        Default pachonModtranAtm_12.dat

    Returns
    -------
    atmos : `Bandpass`
        A bandpass object containing the atmospheric throughput curve.
    """
    atmofile = os.path.join(atmosDir, atmosFile)
    atmo = Bandpass()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        atmo.read_throughput(atmofile)
    atmo.resample_bandpass(wavelen_min=WAVELEN_MIN, wavelen_max=WAVELEN_MAX, wavelen_step=WAVELEN_STEP)
    # Verify that no values go significantly below zero.
    belowzero = np.where(atmo.sb < 0)
    # If there are QE values significantly < 0, raise an exception.
    if np.any(atmo.sb[belowzero] < belowZeroThreshhold):
        raise ValueError("Found values in atmospheric transmission significantly below zero")
    # If they are just small errors in interpolation, set to zero.
    atmo.sb[belowzero] = 0
    return atmo


def buildHardwareAndSystem(defaultDirs, addLosses=True, atmosphereOverride=None, shiftFilters=None):
    """
    Using directories for each component set by 'defaultDirs',
    build the system (including atmosphere) and hardware throughput curves for each filter.

    Parameters
    ----------
    defaultDirs : `dict` of `str`
        Dictionary containing the directories of each throughput component.
        This can be built for the default values, using the setDefaultDirs function.
        There should be 'detector', 'lens1', 'lens2', 'lens3', 'filters', 'mirror1', 'mirror2', 'mirror3',
        and 'atmosphere' keys for this dictionary.
        The 'detector' is built using 'buildDetector' - if this points to a single vendor directory,
        the single vendor QE will be used, but if it is multiple vendors/subdirectories then this is
        the minimum of each QE curve at each wavelength.
    addLosses : `bool`, opt
        Add losses for each component. These correspond to contamination or aging losses.
        Default True.
    atmosphereOverride : `Bandpass`, opt
        A bandpass object containing a user-specified atmosphere throughput file.
        Default None uses the atmosphere file specified in defaultDirs.
    shiftFilters : `float`, opt
        If specified, then the filter throughputs will be shifted by this percent relative to their
        effective wavelengths. Default None (0).

    Returns
    -------
    hardware, system : `dict`, `dict`
        Dictionary of bandpass objects containing the hardware and system throughputs, keyed per filter.
    """
    # Build each component.
    detector = buildDetector(defaultDirs["detector"], addLosses)
    lens1 = buildLens(defaultDirs["lens1"], addLosses)
    lens2 = buildLens(defaultDirs["lens2"], addLosses)
    lens3 = buildLens(defaultDirs["lens3"], addLosses)
    filters = buildFilters(defaultDirs["filters"], addLosses, shiftFilters=shiftFilters)
    mirror1 = buildMirror(defaultDirs["mirror1"], addLosses)
    mirror2 = buildMirror(defaultDirs["mirror2"], addLosses)
    mirror3 = buildMirror(defaultDirs["mirror3"], addLosses)
    # set above was set to default wavelength binning
    if atmosphereOverride is None:
        atmosphere = readAtmosphere(defaultDirs["atmosphere"])
    else:
        atmosphere = atmosphereOverride
        if np.all(atmosphere.wavelen != detector.wavelen):
            atmosphere.resample_bandpass()
    # Combine the individual components.
    core_sb = detector.sb * lens1.sb * lens2.sb * lens3.sb * mirror1.sb * mirror2.sb * mirror3.sb
    wavelen = detector.wavelen
    hardware = {}
    system = {}
    for f in filters:
        hardware[f] = Bandpass()
        system[f] = Bandpass()
        hw_sb = core_sb * filters[f].sb
        hardware[f].set_bandpass(wavelen, hw_sb)
        system[f].set_bandpass(wavelen, hw_sb * atmosphere.sb)
    return hardware, system


def plotBandpasses(
    bandpassDict,
    title=None,
    newfig=True,
    savefig=False,
    addlegend=True,
    linestyle="-",
    linewidth=2,
):
    """
    Plot the bandpass throughput curves.

    Parameters
    -----------
    bandpassDict : `dict` of `Bandpass`
        Dictionary of the bandpass objects, keyed per filter
    title : `str`, opt
        Title for the plot
    newfig : `bool`, opt
        Start a new figure or reuse the active matplotlib figure. Default True.
    savefig : `bool`, opt
        Save the figure to disk with default name title.png or throughputs.png. Default False.
    addlegend : `bool`, opt
        Add a legend to the plot. Default True.
    linestyle : `str`, opt
        The matplotlib linestyle for the throughput curves. Default '-'.
    linewidth : `int`, opt
        The matplotlib linewidth for the throughput curves. Default 2.
    """
    # Generate a new figure, if desired.
    if newfig:
        plt.figure()
    # Plot the bandpass curves.  Try to sort by filter name ugrizy if those are the keys.
    names = bandpassDict.keys()
    if set(names).issubset(set(filterlist)):
        newnames = []
        for f in filterlist:
            if f in names:
                newnames.append(f)
        names = newnames
    # just set other things to k
    for key in names:
        if key not in filtercolors.keys():
            filtercolors[key] = "k"
    for f in names:
        plt.plot(
            bandpassDict[f].wavelen,
            bandpassDict[f].sb,
            marker="",
            linestyle=linestyle,
            linewidth=linewidth,
            color=filtercolors[f],
            label=f,
        )
    # Only draw the legend if desired (many bandpassDicts plotted together could make the legend unwieldy).
    if addlegend:
        plt.legend(loc="lower right", numpoints=1, fancybox=True, fontsize="smaller")
    # Limit wavelengths to the LSST range.
    plt.xlim(300, 1150)
    plt.ylim(0, 1)
    plt.xlabel("Wavelength (nm)", fontsize="x-large")
    plt.ylabel("Fractional Throughput Response", fontsize="x-large")
    # Only add the grid if it's a new figure (otherwise, it toggles on/off).
    plt.grid(True)
    # Add a plot title.
    if title != None:
        plt.title(title, fontsize="x-large")
    # Save the figure, if desired.
    if savefig:
        figformat = "png"
        if title is not None:
            plt.savefig("%s.%s" % (title, figformat), format=figformat)
        else:
            plt.savefig("throughputs.%s" % (figformat), format=figformat)
    return
