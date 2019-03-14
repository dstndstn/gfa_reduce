import astropy.io.fits as fits
import numpy as np
from astropy import wcs
import os
import ci_reduce.common as common

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
