import astropy.io.fits as fits
import numpy as np
from astropy import wcs
import os
import gfa_reduce.common as common

def nominal_tan_wcs(telra, teldec, extname):
    # Create a new WCS object.  The number of axes must be set
    # from the start

    par = common.ci_misc_params()

    fname = os.path.join(os.environ[par['etc_env_var']],
                         par['headers_dummy_filename'])

    h = fits.getheader(fname, extname=extname)

    h['CRVAL1'] = telra
    h['CRVAL2'] = teldec

    w = wcs.WCS(h)

    return w

def ccd_center_radec(_wcs):
    y_pixel_center = 515.5 # need to double check this relative to convention
    x_pixel_center = 1023.5 # need to double check this relative to convention

    ra, dec = _wcs.all_pix2world(x_pixel_center, y_pixel_center, 0)

    return ra, dec
