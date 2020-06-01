import gfa_reduce.common as common
import astropy.io.fits as fits
import healpy
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

# this is intended to mirror how the DESI imaging surveys access
# Gaia, namely through the HEALPix-elized full-sky catalog at:
#     /global/cfs/cdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom
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

    par = common.gfa_misc_params()

    env_var = par['ps1_env_var'] if ps1 else par['gaia_env_var']
    gaia_dir = os.environ[env_var]

    flist = [os.path.join(gaia_dir, 'chunk-' + str(i).zfill(5) + 
                                    '.fits') for i in ipix]
    return flist

def read_gaia_cat(ra, dec, ps1=False, mjd=None):
    # should add checks to make sure that ra and dec have compatible dimensions
    # should also check that this works for both scalar and array ra/dec

    ipix_all = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)

    ipix_u = np.unique(ipix_all)

    flist = gaia_chunknames(ipix_u, ps1=ps1)

    # for the case of PS1, should eventually add checking/handling of cases
    # where a HEALPix pixel not present in the PS1 chunks (dec < -30)
    # is requested
    
    tablist = []
    for i, f in enumerate(flist):
        if ps1:
            # PS1 is not full-sky ...
            if not os.path.exists(f):
                print('SKIPPING PS1 CHUNK: ', f)
                continue
        
        print('READING : ', f)
        tab = fits.getdata(f)
        _ipix = fits.ColDefs([fits.Column(name=('ps1' if ps1 else 'gaia') + '_heal32_ring', format='I', array=np.ones(len(tab))*ipix_u[i])])
        tab = (fits.BinTableHDU.from_columns(tab.columns + _ipix)).data
        if (not ps1) and (mjd is not None):
            # save a copy of the "original" Gaia (RA, Dec) coords straight from the Gaia DR2 catalog
            radec_orig = fits.ColDefs([fits.Column(name='ra_gaia_orig', format='D', array=tab['ra']), fits.Column(name='dec_gaia_orig', format='D', array=tab['dec'])])

            tab = (fits.BinTableHDU.from_columns(tab.columns + radec_orig)).data

        tablist.append(tab)
        
    if ps1:
        if len(tablist) == 0:
            return None
    
    result = np.hstack(tuple(tablist))

    if ps1:
        # dumb stuff about capitalization of column names
        result.dtype.names = tuple([n.lower() for n in result.dtype.names])
    else:
        # if you have an actual table of Gaia rows, then
        # correct Gaia positions to relevant Mayall epoch
        if mjd is not None:
            assert(mjd >= 58757.0) # 2019Oct01 at 00:00:00
            assert(np.isfinite(mjd))
            assert(isinstance(mjd, float))
            
            print('CORRECTING GAIA POSITIONS FOR PROPER MOTION WHEN POSSIBLE')
            
            # 2015-01-01 00:00:00.000 UTC <-> MJD = 57023
            # a year is 365.25 days according to some definition
            # (may not be the relevant definition but can't be too far off)
            # so mjd_gaia must be (57023 + 365.25/2.0) = 57205.625
            
            mjd_gaia = 57205.625

            # only adjust RA, DEC based on PMRA, PMDEC for rows
            # that have full five-parameter astrometric solutions !!

            # the division by cos(Dec) could cause problems if there were a
            # Gaia source at Dec of exactly +/- 90 deg
            
            ra_corr = result['ra'] + ((mjd - mjd_gaia)/365.25)*result['pmra']/np.cos(result['dec']/(180.0/np.pi))/(3600.0*1000.0)
            dec_corr = result['dec'] + ((mjd - mjd_gaia)/365.25)*result['pmdec']/(3600.0*1000.0)
            
            full_solution = np.isfinite(result['pmra'])
            result['ra'][full_solution] = ra_corr[full_solution]
            result['dec'][full_solution] = dec_corr[full_solution]

            print('adjusted ', np.sum(full_solution), ' of ', len(result),
                  ' Gaia source positions for proper motion')

            assert(np.sum(np.isfinite(result['ra'])) == len(result))
            assert(np.sum(np.isfinite(result['dec'])) == len(result))
            assert(np.sum(np.abs(result['dec']) > 90) == 0)

            # eventually also ensure that ra is bounded between 0 and 360 deg
            
            
    return result

def gaia_xmatch(ra, dec, ps1=False, mjd=None):
    gaia_cat = read_gaia_cat(ra, dec, ps1=ps1, mjd=mjd)

    if ps1 and (gaia_cat is None):
        return None
    
    assert(len(gaia_cat) > 0)
    assert(type(gaia_cat).__name__ == 'ndarray')

    catalog = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    
    gaia_catalog = SkyCoord(ra=gaia_cat['ra']*u.degree, \
                            dec=gaia_cat['dec']*u.degree)

    idx, ang_sep_deg, _ = catalog.match_to_catalog_sky(gaia_catalog)

    gaia_matches = Table(gaia_cat[idx])

    gaia_matches['ang_sep_deg'] = ang_sep_deg

    return gaia_matches
