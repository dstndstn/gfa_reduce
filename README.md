# directories

## etc/
auxiliary files such as master biases, master flats, bad pixel masks, dark images ; will not be checked in to git due to file sizes

## py/analysis/
core image analysis such as source detection, centroid measurements, PSF fitting, flux measurements

## py/calib/
code that generates calibration files in etc/, such as master biases, master flats, bad pixel masks, dark images

## py/imred/
convert raw images to reduced images

## py/plotting/
plotting utilities, for example overplotting circles at the locations of Gaia stars

## py/xmatch/
cross-matching utilities, such as for matching to Gaia and other external catalogs

# top-level py/ Python files

## py/ci.py
will contain class that will serve as repository for nominal CI parameters (gain, readnoise, ...)

## py/common.py
miscellaneous utilities, such as conversions between CI numbering/labeling schemes

## py/exposure.py
will contain class that encapsulates a single CI exposure consisting of multiple single-camera images

## py/image.py
will contain class that encapsulates a single-camera image (and its metadata) of a single CI exposure

## py/ci_wcs.py
utilities related to CI WCS

# environment
to access auxiliary calibration files, it will be necessary to set \$CI\_REDUCE\_ETC to the location of the etc/ directory
this code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC
