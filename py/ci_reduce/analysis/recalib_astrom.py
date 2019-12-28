from .kentools_center import kentools_center, gaia_cat_for_exp
import numpy as np
import astropy.io.fits as fits
# run astrometric calibration given a catalog with the centroids and
# an initial guess (SKYRA, SKYDEC) of the field of view center
# mainly going to be a wrapper for kentools_center

def recalib_astrom(cat, fname_raw):
    # cat should be the catalog for an entire exposure

    extnames = np.unique(cat['camera'])

    try:
        h = fits.getheader(fname_raw, extname='GFA')
    except:
        h = fits.getheader(fname_raw, extname='GUIDER')

    gaia = gaia_cat_for_exp(h['SKYRA'], h['SKYDEC'])
    result = []

    arcmin_max = 6.0
    print('astrometry search using ' + '{:.1f}'.format(arcmin_max) + ' arcminute radius')
    for extname in extnames:
        result.append(kentools_center(cat[cat['camera'] == extname],
                                      h['SKYRA'], h['SKYDEC'],
                                      extname=extname, gaia=gaia,
                                      arcmin_max=arcmin_max))

    return result
