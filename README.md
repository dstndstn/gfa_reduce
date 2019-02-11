# directories

## etc/
auxiliary files such as master biases, master flats, bad pixel masks, dark images; will not be checked in to git due to file sizes involved

## py/scripts/
miscellaneous scripts not intended to be used as part of production pipeline; includes code that generates calibration files in etc/, such as master biases, master flats, bad pixel masks, dark images

## py/ci_reduce/analysis/
core image analysis such as source detection, centroid measurements, PSF fitting, flux measurements

## py/ci_reduce/imred/
convert raw images to reduced images

## py/ci_reduce/plotting/
plotting utilities, for example overplotting circles at the locations of Gaia stars

## py/ci_reduce/xmatch/
cross-matching utilities, such as for matching to Gaia or other external catalogs

# top-level py/ci_reduce Python files

## py/ci_reduce/common.py
miscellaneous utilities, such as conversions between CI numbering/labeling schemes

## py/ci_reduce/dark_current.py
utilities related to dark current

## py/ci_reduce/exposure.py
class that encapsulates a single CI exposure consisting of multiple single-camera images

## py/ci_reduce/image.py
class that encapsulates a single-camera image (and its metadata) drawn from a single CI exposure

## py/ci_reduce/io.py
input/output utilities

## py/ci_reduce/platescale.py
utility for retrieving desimodel platescale as a function of CI pixel coordinates

# convention for pixel coordinates
unless otherwise stated in a particular portion of the code...

Throughout the codebase, pixel coordinate (0, 0) will be the **center** of the lower left pixel (where "lower left" applies in the context of the IDL/DS9 convention for 2D image display). "x" will refer to the "AXIS1" coordinate (long axis for CI cameras, a.k.a. "width"), and "y" will refer to the "AXIS2" coordinate (short axis for CI cameras, a.k.a. "height").

# environment
* to access auxiliary calibration files, it will be necessary to set the CI_REDUCE_ETC environment variable to the location of the etc/ directory
  * the authoritative copy of these files can be found at /project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/
* to enable Gaia cross-matching, the GAIA_CAT_DIR environment variable must be set to the full path of "chunks-gaia-dr2-astrom"
  * the recommended path at NERSC is /global/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom
* this code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
  * https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC
