# syseng_throughputs #
SysEng-approved LSST throughput curves
The latest m5 depths are available in the notebooks, such as in [notebooks/Overview Paper.ipynb](./notebooks/Overview%20Paper.ipynb).

This repository provides the ultimate source of the throughput curves in the repository [lsst/throughputs](https://github.com/lsst/throughputs).

The [components](./components) directory contains the response curves
for each individual component of the camera and telescope. In each
directory, there is also a `*_Losses` directory that contains the
time-averaged ten-year losses due to contamination or condensation on
the surfaces of the component. In some directories, there is also a
`*_Coatings` directory, which contains information on coatings applied
to the surface, such as the Broad Band Anti-Reflection coatings on the
lenses.

These components curves are maintained and updated by the LSST system
engineering team.

Python utilties to read and combine these various
curves appropriately are maintained in this repository, in the
[syseng_throughputs](.syseng_throughputs) directory. In particular, note the utilities
provided in [bandpassUtils.py](./syseng_throughputs/bandpassUtils.py). At this
time, we expect most users to use the throughputs repository instead
of this repository directly - the curves in the throughputs repository
are constructed from these curves, and can be traced through the git
SHA1 and release tags.

The python code requires [rubin_sim](https://github.com/lsst/rubin_sim) to run.
After installation of rubin_sim, install syseng_throughputs into the same python environment using
```pip install -e .```

# Release 1.9 #

This updates uses the same throughput components are previous versions, but moves from Al-Ag-Al mirror coatings to Ag-Ag-Ag (triple silver, or 3Ag) mirror coatings.
This results in increased throughput in redder bands, at the cost of lower throughput in u band. However, since more survey time is spent in r and i bands (and redder bands in general) than u, the overall impact on survey efficiency is positive.

# Release 1.8 #

This update includes as-measured filter throughput curves and as-measured glass and coating measurements for the lenses.
The mirror reflectivities are also all updated to as-measured curves for Al and Ag, but the mirror coatings are assumed to remain Al-Ag-Al in this tag. 


# Release 1.7 #

  The M2 reflectivity was updated based on witness sample measurements from M2 coating run in July 2019. The PR for this update is https://github.com/lsst-pst/syseng_throughputs/pull/12. The notebooks showing what has changed and the updated m5 calculations are found in the "documentation" subdirectory.

# Release 1.6 #

The mirror reflectivity was updated based on measurements from coating samples from June 2019. The PR for this update is https://github.com/lsst-pst/syseng_throughputs/pull/11. The notebooks showing what has changed and the updated m5 calculations are found in the "documentation" subdirectory.

# Release 1.5 #

This is a minor update for throughputs (the lens2 glass and BBAR coating curves
have been extended in their wavelength information, but the curves themselves
are the same as previously). However it is a major update for documentation
and process information, as reflected in the "documentation" subdirectory.

# Release 1.4 #

The primary update here is in the lens2 response curves. The BBAR coating
has been updated.

Other minor updates include bug fixes in the python code in sedUtils.py,
updating of the jupyter notebooks, and the addition of notebooks evaluating
the effect of the mixed vendor detector focal plane and recreating the
inputs for the LSST Overview Paper.

# Release 1.3: #

The primary update here is in the detector response curves.
The QE response curves here are the result of measurements of multiple
chips provided by each vendor, ITL and E2V. The measurements have been
averaged across multiple CCDs; the default (single) 'generic' curve remains
the minimum QE response at each wavelength between both vendors.
These curves were provided by Steve Ritz in December, 2017.

Other minor updates include additional python code to allow scaling
of the FWHM at different airmasses and wavelengths (according to
details provided in Document-18208 and Document-20160), and a jupyter
notebook which can provide latex-formatted content of Table 2 from the
overview paper.

# Release 1.2: #

This is primarily an update to the python code in the repository, using
corrected and updated readnoise values (which results in corresponding
changes to m5, particularly in the u band).


# As of release 1.1: #

## Camera Components ##
* _Detector_: There are two separate detector response and loss curves,
  corresponding to the expected response (QE response + AR coatings)
  of the CCDs provided by each of the two vendors under
  consideration. For most purposes (including the detector curve reported in
  the throughputs repository), we use a 'generic' detector
  response that is generated by combining both of these throughput
  curves using the *minimum* QE response at each wavelength.
  The response curves from each vendor correpond to a response
  measured in LSST labs, using vendor-provided prototypes. The loss
  curves provided for each vendor represent a simulated effect of
  contamination buildup over time; the loss curves are identical for
  both vendors and are the average expected values over ten
  years. Note that some values in the 'contamination' loss file for
  the detectors are > 1; this is because the contamination is
  primarily a thin film of water, which at some wavelengths can
  enhance the performance of the AR coating on the detector -- this is
  only true for the detector.
* _Lenses_: There are three separate lenses in the camera, each with an
  identical base `*_Glass.dat` curve that represents the fused silica
  throughput of the lens itself. This throughput curve must be smoothed using the
  Savitzy-Golay smoothing function. The fused silica lens transmission curves are
  based on vendor-provided expected transmission curves. The silica base of the len must
  also be combined with the BroadBand AntiReflective (BBAR) coatings
  response in the `*_Coatings` directory. There are two coatings; one
  for each side of the lens. The BBAR coating response is based on vendor-provided
  models, consistent with LSST requested coating requirements. There are small differences between the
  glass components used for each lens; there are also small
  differences in the BBARS, including a difference from one side of
  the lens to the other. In each lens, there are also several files in
  the `*_Losses` directory, representing the time-averaged condensation and
  contamination losses for each surface of each lens. The losses are based on
  models developed by Andy Rasmussen at SLAC. These vary
  depending on the direction the lens is facing and the location of
  the lens in the camera. The final response curves for all lenses are
  similar in shape, however lens3 has a slightly higher overall
  throughput due to slightly lower losses (only by 1-2%).
* _Filters_: For each filter, a goal throughput envelope has been
  provided. This is the goal throughput envelope provided to the
  filter vendors; tolerances on this envelope have also been
  provided. Note that this is not the expected performance for an
  as-manufactured filter, which would likely include some out-of-band throughput leaks
  (within a specified limit), and represents a change compared to
  previously provided throughput curves (which represented one simulation of
  an expected as-provided filter set). In the `*_Losses` directory,
  there are also ten-year-average simulated
  contamination and condensation losses for each surface of the
  filters, based on models developed by Andy Rasmussen.

## Telesope Components ##
* _Mirrors_: Each mirror has a reflectivity curve, which should be
  coupled with the respective losses curve found in the relevant
  `*_Losses` directory. The reflectivity of mirror1 (primary mirror) and mirror3
  (tertiary) is based on using a protected aluminum surface; the
  reflectivity of mirror2 (secondary) is based on using a protected
  silver surface. These mirror reflectivities are based on lab measurements
  of pristine witness samples. The losses represent the ten-year average,
  based on performance degradation measurements from historical telescope performance,
  modified for the expected LSST maintenance schedule.
  Currently mirror cleanings are scheduled yearly, with resurfacing every
  two years.


## Site Properties ##
* _Atmosphere_: The atmosphere throughput is modeled by using MODTRAN to
 produce a 'standard US Atmosphere', which does not include aerosols.
 To better represent the expected atmospheric transmission on site, aerosols
 have been added to the resulting throughput curves, using the python
 script [addAerosols.py](./python/addAerosols.py). The atmospheric
 transmission curves are in the [siteProperties](./siteProperties)
 directory, with an X=1.2 and X=1.0 atmosphere, with and without
 aerosols. To represent 'typical' throughput, the X=1.2, with aerosols
 [atmosphere](./siteProperties/pachonModtranAtm_12_aerosol.dat) curve
 should be used. To represent zenith, optimum throughputs, the X=1.0,
 with aerosols [atmosphere](./siteProperties/atmos_10_aerosol.dat)
 curve should be used.
* _Dark sky_: The expected dark sky, zenith, background spectrum can
   be found in [darksky.dat](./siteProperties/darksky.dat). This is
   used to calculate expected zenith, dark-sky limiting magnitude
   values. The dark sky SED is based on data from UVES and Gemini Near-IR,
   combined with ESO sky data from Ferdinand Patat, modified slightly at
   the red and blue ends to match observed dark sky broadband skybrightness
   values reported by SDSS.
