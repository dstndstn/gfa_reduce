import ci_reduce.xmatch.gaia as gaia_xmatch
import astropy.io.fits as fits
import os
import healpy
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import wcs
from amm_2dhist import amm_2dhist
from scipy.ndimage import gaussian_filter

def downselected_star_sample(cat, n_desi_max):
    assert(len(np.unique(cat['CAMERA'])) == 1)

    _cat = cat[cat['DQ_FLAGS'] == 0]

    sind = np.argsort(-1.0*_cat['APER_SUM_BKGSUB_3'])

    n = len(sind)

    result = _cat[sind[0:(min(n_desi_max, n))]]

    return result
    
def kentools_center(fname_cat, extname='CIC', arcmin_max=2.0):
    assert(os.path.exists(fname_cat))

    # output needs to include, at a minimum:
    #     xshift_best
    #     yshift_best
    #     contrast
    #     extname
    # probably also want
    #     expid retrieved from the catalog table

    fname_reduced = fname_cat.replace('_catalog', '_reduced')

    assert(os.path.exists(fname_reduced))

    h = fits.getheader(fname_reduced, extname=extname)

    expid = h['EXPID']

    cat = fits.getdata(fname_cat)
    racen = h['SKYRA']
    deccen = h['SKYDEC']

    nside = 32
    ipix = healpy.pixelfunc.get_all_neighbours(nside, racen, phi=deccen,
                                               lonlat=True)

    # probably also want to include the neighbors to the neighbors
    ipix_all = []
    for _ipix in ipix:
        ipix_all.append(healpy.pixelfunc.get_all_neighbours(nside,
                                                            _ipix))

    ipix = np.unique(ipix_all)
    
    ra_pixcenters, dec_pixcenters = healpy.pixelfunc.pix2ang(nside, ipix,
                                                              lonlat=True)

    gaia = gaia_xmatch.read_gaia_cat(ra_pixcenters, dec_pixcenters)

    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)
    c = SkyCoord(racen*u.deg, deccen*u.deg)

    dangle = g.separation(c)

    gaia = gaia[dangle.degree < 2]
    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)
    
    cat = cat[cat['CAMERA'] == extname]

    assert(len(cat) > 0)
    
    n_desi_max = 150

    if len(cat) > n_desi_max:
        cat = downselected_star_sample(cat, n_desi_max)

    h = fits.getheader('/project/projectdirs/desi/users/ameisner/CI/ci_reduce_etc/dummy_with_headers-as_built.bigtan.fits.gz', extname=extname)

    assert(h['EXTNAME'] == extname)
    
    astrom = wcs.WCS(h)
    astrom.wcs.crval = [racen, deccen]

    x_gaia_guess, y_gaia_guess = astrom.wcs_world2pix(gaia['ra'], gaia['dec'],
                                                      0)
    dx_all = np.array([], dtype='float64')
    dy_all = np.array([], dtype='float64')

    cat = cat[np.argsort(cat['DEC'])] # not really necessary

    for i in range(len(cat)):
        print(i+1, ' of ', len(cat))
        c = SkyCoord(cat[i]['RA']*u.deg, cat[i]['DEC']*u.deg)
        dangle = c.separation(g)
        w = (np.where(dangle.arcminute < arcmin_max))[0]
        print(len(w), len(g), len(gaia))
        if len(w) == 0:
            continue
        dx = cat[i]['XCENTROID'] - x_gaia_guess[w]
        dy = cat[i]['YCENTROID'] - y_gaia_guess[w]
        dx_all = np.concatenate((dx_all, dx))
        dy_all = np.concatenate((dy_all, dy))

    assert(len(dx_all) == len(dy_all))
    print(np.min(dy_all), np.max(dy_all))

    axlim = max(np.round(arcmin_max*60.0/0.128), 1000.0)

    print(axlim)

    dx = 1.0
    dy = 1.0
    nx = ny = 2*axlim
    counts, x_edges_left, y_edges_left = amm_2dhist(-1.0*axlim, -1.0*axlim,
                                                    nx, ny, dx, dy, 
                                                    dx_all, dy_all)

    counts = counts.astype(float)

    smth = gaussian_filter(counts, 7.5, mode='constant')

    counts_shape = counts.shape
    sidelen = counts_shape[0]
    #plt.imshow(smth, cmap='gray_r')
    #plt.show()
    indmax = np.argmax(smth)

    indx = (indmax // sidelen).astype(int)
    indy = (indmax % sidelen).astype(int)

    print(indmax, indx, indy, counts_shape, x_edges_left[indx],
          y_edges_left[indy])
    #print(x_edges_left[indx], y_edges_left[indy], ' ~~~~~~~~~~~~~~~~')
    
    return smth

def _test(extname='CIC', arcmin_max=3.5):
    fname_cat = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190403/ci-00003697/ci-00003697_catalog.fits'

    
    gaia = kentools_center(fname_cat, extname=extname, arcmin_max=arcmin_max)

    return gaia
