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
from center_contrast import center_contrast
#import fitsio # for testing

def downselected_star_sample(cat, n_desi_max):
    assert(len(np.unique(cat['CAMERA'])) == 1)

    _cat = cat[cat['DQ_FLAGS'] == 0]

    sind = np.argsort(-1.0*_cat['APER_SUM_BKGSUB_3'])

    n = len(sind)

    result = _cat[sind[0:(min(n_desi_max, n))]]

    return result
    
def kentools_center(fname_cat, extname='CIC', arcmin_max=3.5):
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

    # this apparently happened for the CI in some cases...
    if isinstance(racen, str) or isinstance(deccen, str):
        return None
    
    nside = 32
    ipix = healpy.pixelfunc.get_all_neighbours(nside, racen, phi=deccen,
                                               lonlat=True)

    # remove dummy -1 values from ipix

    ipix = ipix[ipix >= 0]
    
    # probably also want to include the neighbors to the neighbors
    ipix_all = []
    for _ipix in ipix:
        ipix_all.append(healpy.pixelfunc.get_all_neighbours(nside,
                                                            _ipix))
    
    ipix = np.unique(ipix_all)

    ipix = ipix[ipix >= 0]
    
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

    # 0.128 value is tailored to the CI -- revisit for GFAs !!
    axlim = max(np.round(arcmin_max*60.0/0.128), 1000.0)

    print(axlim)

    dx = 1.0
    dy = 1.0
    nx = ny = 2*axlim
    counts, x_edges_left, y_edges_left = amm_2dhist(-1.0*axlim, -1.0*axlim,
                                                    nx, ny, dx, dy, 
                                                    dx_all, dy_all)

    counts = counts.astype(float)

    # 7.5 value is tailored to the CI -- revisit for GFAs !!
    fwhm_pix = 7.5
    sigma_pix = fwhm_pix/(2*np.sqrt(2*np.log(2)))
    smth = gaussian_filter(counts, sigma_pix, mode='constant')

    counts_shape = counts.shape
    sidelen = counts_shape[0]
    #plt.imshow(smth, cmap='gray_r')
    #plt.show()
    indmax = np.argmax(smth)

    indx = (indmax // sidelen).astype(int)
    indy = (indmax % sidelen).astype(int)

    #print(indmax, indx, indy, counts_shape, x_edges_left[indx],
    #      y_edges_left[indy])
    #print(x_edges_left[indx], y_edges_left[indy], ' ~~~~~~~~~~~~~~~~')

    xshift_best = x_edges_left[indx] + 0.5*dx
    yshift_best = y_edges_left[indy] + 0.5*dy

# assumes dx = dy = 1 !!
    ycen_grid, xcen_grid = np.meshgrid(x_edges_left[0:(len(x_edges_left)-1)] + 0.5*dx, y_edges_left[0:(len(y_edges_left)-1)] + 0.5*dy)
    
    d = np.sqrt(np.power(xcen_grid - xshift_best, 2) + np.power(ycen_grid - yshift_best, 2))

    #fitsio.write('/global/cscratch1/sd/ameisner/smth.fits', smth)
    #fitsio.write('/global/cscratch1/sd/ameisner/d.fits', d)
    
    r_max = 4*7.5 # roughly 4 asec radius for the CI

    sind = np.argsort(np.ravel(smth))

    val_90 = np.ravel(smth)[sind[int(np.round(0.925*len(sind)))]]

    high = [(d < r_max) & ((smth == np.max(smth)) | (smth > val_90))]

    assert(np.sum(high) > 0)
    
    print(val_90, '^^^^^^^^^^^^^^^^^^^^^^^^^^', np.sum(high))
    
    print(xshift_best, yshift_best)

    wt = np.sum(high*smth)
    xshift_best = np.sum(high*xcen_grid*smth)/wt
    yshift_best = np.sum(high*ycen_grid*smth)/wt

    print(xshift_best, yshift_best)
    contrast = center_contrast(smth)

    print(contrast, ' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

    print(xcen_grid.shape)
    print(counts.shape, smth.shape)

    result = {'xshift_best': xshift_best,
              'yshift_best': yshift_best,
              'contrast': contrast,
              'expid': expid,
              'extname': extname,
              'astr_guess': astrom,
              'fname_cat': fname_cat}

    return result

def _test(extname='CIC', arcmin_max=3.5):
    fname_cat = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190403/ci-00003697/ci-00003697_catalog.fits'

    
    gaia = kentools_center(fname_cat, extname=extname, arcmin_max=arcmin_max)

    return gaia

def _loop(indstart, nproc):
    tab = fits.getdata('/global/homes/a/ameisner/ci/pro/ci_quality_summary.fits')
    # tab = tab[tab['GOOD'] == 1]

    # shuffle
    np.random.seed(seed=99)
    sind = np.argsort(np.random.rand(len(tab)))

    tab = tab[sind]

    results = []

    indend = min(indstart + nproc, len(tab))
    for i in range(indstart, indend):
        t = tab[i]
        fname_cat = t['FNAME'].replace('_psf-a', '_catalog')
        extname = t['EXTNAME']
        print('WORKING ON: ', i, '   ', fname_cat, '   ', extname)
        result = kentools_center(fname_cat, extname=extname, arcmin_max=3.5)
        results.append(result)

    # then write out the results to a pickle file
    # need to construct the output file name first
    outname = 'center_' + str(indstart).zfill(5) + '_' + str(indend-1).zfill(5) + '.pkl'
    assert(not os.path.exists(outname))
    import pickle
    pickle.dump(results, open(outname, 'wb'))
