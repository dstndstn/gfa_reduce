# directories

## py/scripts/
miscellaneous scripts not intended to be used as part of the production pipeline

## py/gfa_reduce/analysis/
core image analysis such as source detection, centroid measurements, flux measurements

## py/gfa_reduce/imred/
utilities for converting raw images to reduced images

## py/gfa_reduce/xmatch/
cross-matching utilities, such as for matching to Gaia or other external catalogs

# top-level py/gfa_reduce Python files

## py/gfa_reduce/gfa_red.py
primary driver that gets called to run the GFA reduction pipeline end to end

## py/gfa_reduce/gfa_wcs.py
WCS-related utilities

## py/gfa_reduce/common.py
miscellaneous utilities

## py/gfa_reduce/dark_current.py
utilities related to dark current

## py/gfa_reduce/exposure.py
class that encapsulates a single GFA exposure consisting of multiple single-camera images

## py/gfa_reduce/image.py
class that encapsulates a single-camera image (and its metadata) drawn from a single GFA exposure

## py/gfa_reduce/io.py
input/output utilities

# convention for pixel coordinates
unless otherwise stated in a particular portion of the code...

Throughout the codebase, pixel coordinate (0, 0) will be the **center** of the lower left pixel (where "lower left" applies in the context of the IDL/DS9 convention for 2D image display). "x" will refer to the "AXIS1" coordinate (long axis for GFA cameras, a.k.a. "width"), and "y" will refer to the "AXIS2" coordinate (short axis for GFA cameras, a.k.a. "height").

# environment
* to access auxiliary calibration files, it will be necessary to set the `GFA_REDUCE_META` environment variable to their directory's location
  * examples of auxiliary calibration files include master biases, master flats, static bad pixel masks
  * these auxiliary files will not be checked in to git due to the file sizes involved
  * the authoritative copy of this directory can be found on NERSC at `/global/cfs/cdirs/desi/users/ameisner/GFA/gfa_reduce_meta`
* to enable Gaia cross-matching, the `GAIA_CAT_DIR` environment variable must be set to the full path of "chunks-gaia-dr2-astrom"
  * the recommended path at NERSC is `/global/cfs/cdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom`
  * this directory contains a full-sky set of 12,288 FITS "chunk" files, one per nside = 32 HEALPix pixel, with ring-ordered indexing in equatorial coordinates
* this code is intended to be run at NERSC using the DESI software environment, which is Python 3 based:
  * https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC

# example environment configuration script

the following is a bash script that can be used to configure one's environment for running the `gfa_reduce` package at NERSC

    source /global/cfs/cdirs/desi/software/desi_environment.sh master
    export PYTHONPATH=${PYTHONPATH}:/global/homes/a/ameisner/gfa_reduce/py
    export GFA_REDUCE_META=/global/cfs/cdirs/desi/users/ameisner/GFA/gfa_reduce_meta
    export GAIA_CAT_DIR=/global/cfs/cdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom

reading auxiliary files from other places or running a `gfa_reduce` checkout located elsewhere on NERSC would require modifications of these example paths

it would be recommended to run `gfa_reduce` using a checkout of this repository's code rather than the exact NERSC location given above, as code development may be taking place within that `/global/homes/a/ameisner/gfa_reduce/py` directory