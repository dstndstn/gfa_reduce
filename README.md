# directories

## etc/
auxiliary files such as master biases, master flats, bad pixel masks, dark images; will not be checked in to git due to file sizes involved

## py/scripts
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

## py/ci_reduce/ci_wcs.py
utilities related to CI WCS

## py/ci_reduce/common.py
miscellaneous utilities, such as conversions between CI numbering/labeling schemes

## py/ci_reduce/exposure.py
will contain class that encapsulates a single CI exposure consisting of multiple single-camera images

## py/ci_reduce/image.py
will contain class that encapsulates a single-camera image (and its metadata) drawn from a single CI exposure

# environment
* to access auxiliary calibration files, it will be necessary to set $CI\_REDUCE\_ETC to the location of the etc/ directory
* to enable Gaia cross-matching, the GAIA_CAT_DIR environment variable must be set to the full path of "chunks-gaia-dr2-astrom"
* this code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC
