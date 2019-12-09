import ci_reduce.common as common
import astropy.io.fits as fits
import healpy
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

# this is intended to mirror how the DESI imaging surveys access
# Gaia, namely through the HEALPix-elized full-sky catalog at:
#     /global/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom
#
# each catalog contains one Nside = 32 HEALPix pixel worth of Gaia souces
# the HEALPix indices are determined using RA/Dec as longitude/latitude
# HEALPix indexing is ring-ordered

nside = 32

def gaia_chunknames(ipix, ps1=False):
    # could add checks to make sure that all ipix values are 
    # sane HEALPix pixel indices
    # RIGHT NOW THIS ASSUMES IPIX IS AN ARRAY !! 
    # should eventually make this also work for scalar ipix

    par = common.ci_misc_params()

    env_var = par['ps1_env_var'] if ps1 else par['gaia_env_var']
    gaia_dir = os.environ[env_var]

    flist = [os.path.join(gaia_dir, 'chunk-' + str(i).zfill(5) + 
                                    '.fits') for i in ipix]
    return flist

def read_gaia_cat(ra, dec, ps1=False):
    # should add checks to make sure that ra and dec have compatible dimensions
    # should also check that this works for both scalar and array ra/dec

    ipix_all = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)

    ipix_u = np.unique(ipix_all)

    flist = gaia_chunknames(ipix_u, ps1=ps1)

    # for the case of PS1, should eventually add checking/handling of cases
    # where a HEALPix pixel not present in the PS1 chunks (dec < -30)
    # is requested
    
    tablist = []
    for f in flist:
        print('READING : ', f)
        tab = fits.getdata(f)
        tablist.append(tab)

    result = np.hstack(tuple(tablist))

    if ps1:
        # dumb stuff about capitalization of column names
        result.dtype.names = tuple([n.lower() for n in result.dtype.names])

    return result

def gaia_xmatch(ra, dec, ps1=False):
    gaia_cat = read_gaia_cat(ra, dec, ps1=ps1)

    assert(len(gaia_cat) > 0)
    assert(type(gaia_cat).__name__ == 'ndarray')

    ci_catalog = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    
    gaia_catalog = SkyCoord(ra=gaia_cat['ra']*u.degree, \
                            dec=gaia_cat['dec']*u.degree)

    idx, ang_sep_deg, _ = ci_catalog.match_to_catalog_sky(gaia_catalog)

    gaia_matches = Table(gaia_cat[idx])

    gaia_matches['ang_sep_deg'] = ang_sep_deg

    return gaia_matches
