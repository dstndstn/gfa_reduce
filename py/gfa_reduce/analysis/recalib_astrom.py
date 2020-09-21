from .asterisms import pattern_match, gaia_cat_for_exp
import numpy as np
import astropy.io.fits as fits
# run astrometric calibration given a catalog with the centroids and
# an initial guess (SKYRA, SKYDEC) of the field of view center
# mainly going to be a wrapper for asterisms.pattern_match

def recalib_astrom(cat, fname_raw, mjd=None, h=None):
    # cat should be the catalog for an entire exposure

    if cat is None:
        return None
        
    extnames = np.unique(cat['camera'])

    if h is None:
        try:
            h = fits.getheader(fname_raw, extname='GFA')
        except:
            h = fits.getheader(fname_raw, extname='GUIDER')

    gaia = gaia_cat_for_exp(h['SKYRA'], h['SKYDEC'], mjd=mjd)
    # conserve memory
    gaia = gaia[['ra', 'dec']] # only columns needed for pattern matching
    result = []

    arcmin_max = 6.0
    print('astrometry search using ' + '{:.1f}'.format(arcmin_max) + ' arcminute radius')
    for extname in extnames:
        _cat = cat[(cat['camera'] == extname) & (cat['sig_major_pix'] > 1.0) & (cat['min_edge_dist_pix'] > 3)] # 3 is a fairly arbitrary guess
        if len(_cat) < 2:
            _cat = cat[cat['camera'] == extname]
        result.append(pattern_match(_cat,
                                    h['SKYRA'], h['SKYDEC'],
                                    extname=extname, gaia=gaia,
                                    arcmin_max=arcmin_max))

    return result
