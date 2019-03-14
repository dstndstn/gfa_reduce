# directories

## py/scripts/
miscellaneous scripts not intended to be used as part of production pipeline; includes code that generates calibration files in `CI_REDUCE_ETC`, such as master biases, master flats, static bad pixel masks

## py/ci_reduce/analysis/
core image analysis such as source detection, centroid measurements, flux measurements

## py/ci_reduce/imred/
utilities for converting raw images to reduced images

## py/ci_reduce/plotting/
plotting utilities

## py/ci_reduce/xmatch/
cross-matching utilities, such as for matching to Gaia or other external catalogs

# top-level py/ci_reduce Python files

## py/ci_reduce/ci_proc.py
primary driver that gets called to run the CI reduction pipeline end to end

## py/ci_reduce/ci_wcs.py
WCS-related utilities

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

# convention for pixel coordinates
unless otherwise stated in a particular portion of the code...

Throughout the codebase, pixel coordinate (0, 0) will be the **center** of the lower left pixel (where "lower left" applies in the context of the IDL/DS9 convention for 2D image display). "x" will refer to the "AXIS1" coordinate (long axis for CI cameras, a.k.a. "width"), and "y" will refer to the "AXIS2" coordinate (short axis for CI cameras, a.k.a. "height").

# CI camera naming convention

Throughout the codebase, the CI camera names CIE, CIN, CIS, and CIW refer to sky directions (not geographic directions).

# environment
* to access auxiliary calibration files, it will be necessary to set the `CI_REDUCE_ETC` environment variable to their directory's location
  * examples of auxiliary calibration files include master biases, master flats, static bad pixel masks, master dark images
  * these files will not be checked in to git due to the file sizes involved
  * the authoritative copy of this directory can be found at `/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc`
* to enable Gaia cross-matching, the `GAIA_CAT_DIR` environment variable must be set to the full path of "chunks-gaia-dr2-astrom"
  * the recommended path at NERSC is `/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom`
  * this directory contains a full-sky set of 12,288 FITS "chunk" files, one per nside = 32 HEALPix pixel, with ring-ordered indexing in equatorial coordinates
* this code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
  * https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC

# example environment configuration script

the following is a bash script that can be used to configure one's environment for running the `ci_reduce` package at NERSC

    source /project/projectdirs/desi/software/desi_environment.sh 18.7
    export PYTHONPATH=${PYTHONPATH}:/global/homes/a/ameisner/ci_reduce/py
    export CI_REDUCE_ETC=/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc
    export GAIA_CAT_DIR=/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom

reading auxiliary files from other places or running a `ci_reduce` checkout located elsewhere on NERSC would require modifications of these example paths
