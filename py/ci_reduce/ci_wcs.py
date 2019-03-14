import astropy.io.fits as fits
import numpy as np
from astropy import wcs
import os

def nominal_tan_wcs(telra, teldec, extname):
    # Create a new WCS object.  The number of axes must be set
    # from the start

    # eventually avoid this hardcoding by appropriately
    # updating/using common.py
    fname = os.path.join(os.environ['CI_REDUCE_ETC'],
                         'dummy_with_headers.bigtan.fits.gz')

    h = fits.getheader(fname, extname=extname)

    h['CRVAL1'] = telra
    h['CRVAL2'] = teldec

    w = wcs.WCS(h)

    return w
