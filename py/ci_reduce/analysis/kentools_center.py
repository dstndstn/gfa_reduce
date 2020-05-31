import ci_reduce.xmatch.gaia as gaia_xmatch
import astropy.io.fits as fits
import os
import healpy
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import wcs
from .amm_2dhist import amm_2dhist
from scipy.ndimage import gaussian_filter
from .center_contrast import center_contrast
import ci_reduce.common as common
import copy

def downselected_star_sample(cat, n_desi_max):
    assert(len(np.unique(cat['camera'])) == 1)

    try:
        _cat = cat[cat['DQ_FLAGS'] == 0]
    except:
        _cat = cat[cat['dq_flags'] == 0]

    try:
        flux = _cat['APER_SUM_BKGSUB_3']
    except:
        flux = _cat['aper_sum_bkgsub_3']
    
    sind = np.argsort(-1.0*flux)
    
    n = len(sind)

    result = _cat[sind[0:(min(n_desi_max, n))]]

    return result

def gaia_cat_for_exp(racen, deccen, mjd=None):

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

    gaia = gaia_xmatch.read_gaia_cat(ra_pixcenters, dec_pixcenters, mjd=mjd)
    return gaia

def kentools_center(catalog, skyra, skydec, extname='GUIDE0', arcmin_max=3.5,
                    gaia=None):

    # cat needs to have fields xcentroid and ycentroid
    # skyra, skydec are the initial guesses of the 
    # actual center of the field of view; these need to be within
    # arcmin_max of the true FOV center for this routine to succeed
    
    #assert(os.path.exists(fname_cat))

    # output needs to include, at a minimum:
    #     xshift_best
    #     yshift_best
    #     contrast
    #     extname
    # probably also want
    #     expid retrieved from the catalog table

    #fname_reduced = fname_cat.replace('_catalog', '_reduced')

    #assert(os.path.exists(fname_reduced))

    #h = fits.getheader(fname_reduced, extname=extname)

    #expid = h['EXPID']

    #cat = fits.getdata(fname_cat)
    cat = copy.deepcopy(catalog)
 # should just entirely rename racen, deccen to skyr, skydec ...
    racen = skyra
    deccen = skydec

    # this apparently happened for the CI in some cases...
    if isinstance(racen, str) or isinstance(deccen, str):
        return None

    if gaia is None:
        gaia = gaia_cat_for_exp(racen, deccen)

    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)
    c = SkyCoord(racen*u.deg, deccen*u.deg)

    dangle = g.separation(c)

    gaia = gaia[dangle.degree < 2]
    g = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg)
    
    cat = cat[cat['camera'] == extname]

    assert(len(cat) > 0)
    
    n_desi_max = 150

    if len(cat) > n_desi_max:
        cat = downselected_star_sample(cat, n_desi_max)

    par = common.ci_misc_params()
    fname_wcs_templates = os.environ['GFA_REDUCE_ETC'] + '/' + par['headers_dummy_filename']

    assert(os.path.exists(fname_wcs_templates))

    h = fits.getheader(fname_wcs_templates, extname=extname)

    assert(h['EXTNAME'] == extname)
    
    astrom = wcs.WCS(h)
    astrom.wcs.crval = [racen, deccen]

    x_gaia_guess, y_gaia_guess = astrom.wcs_world2pix(gaia['ra'], gaia['dec'],
                                                      0)
    dx_all = np.array([], dtype='float64')
    dy_all = np.array([], dtype='float64')

    cat = cat[np.argsort(cat['dec'])] # not really necessary

    _g = None
    for i in range(len(cat)):
        print(i+1, ' of ', len(cat))
        c = SkyCoord(cat[i]['ra']*u.deg, cat[i]['dec']*u.deg)
        dangle = c.separation(g if _g is None else _g)
        w = (np.where(dangle.arcminute < arcmin_max))[0]
        print(len(w), len(dangle))
        if len(w) == 0:
            continue

        assert(len(dangle) == len(x_gaia_guess))
        assert(len(dangle) == len(y_gaia_guess))
        
        dx = cat[i]['xcentroid'] - x_gaia_guess[w]
        dy = cat[i]['ycentroid'] - y_gaia_guess[w]
        dx_all = np.concatenate((dx_all, dx))
        dy_all = np.concatenate((dy_all, dy))

        if _g is None:
            padfac = 2.0
            keep = (dangle.arcminute < (8.1 + padfac*arcmin_max))
            x_gaia_guess = x_gaia_guess[keep]
            y_gaia_guess = y_gaia_guess[keep]
            _g = g[keep]
        
    assert(len(dx_all) == len(dy_all))
    #print(np.min(dy_all), np.max(dy_all))

    axlim = max(np.round(arcmin_max*60.0/0.214), 1000.0)

    #print(axlim)

    dx = 1.0
    dy = 1.0
    nx = ny = 2*axlim
    counts, x_edges_left, y_edges_left = amm_2dhist(-1.0*axlim, -1.0*axlim,
                                                    nx, ny, dx, dy, 
                                                    dx_all, dy_all)

    counts = counts.astype(float)

    # 7.5 value is tailored to the CI -- revisit for GFAs !!
    fwhm_pix = 4.7
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
    ycen_grid, xcen_grid = np.meshgrid(x_edges_left[0:(len(x_edges_left)-1)], y_edges_left[0:(len(y_edges_left)-1)])
    
    d = np.sqrt(np.power(xcen_grid - xshift_best, 2) + np.power(ycen_grid - yshift_best, 2))

    #fitsio.write('/global/cscratch1/sd/ameisner/smth.fits', smth)
    #fitsio.write('/global/cscratch1/sd/ameisner/d.fits', d)
    
    r_max = 4*4.7 # roughly 4 asec radius for the CI

    sind = np.argsort(np.ravel(smth))

    val_90 = np.ravel(smth)[sind[int(np.round(0.925*len(sind)))]]

    high = [(d < r_max) & ((smth == np.max(smth)) | (smth > val_90))]

    assert(np.sum(high) > 0)
    
    print(xshift_best, yshift_best)

    wt = np.sum(high*smth)
    xshift_best = np.sum(high*xcen_grid*smth)/wt
    yshift_best = np.sum(high*ycen_grid*smth)/wt

    print(xshift_best, yshift_best)
    contrast = center_contrast(smth)


    #print(xcen_grid.shape)
    #print(counts.shape, smth.shape)

    result = {'xshift_best': xshift_best,
              'yshift_best': yshift_best,
              'contrast': contrast,
              'extname': extname,
              'astr_guess': astrom}

    return result

def _test(extname='CIC', arcmin_max=3.5):
    fname_cat = '/project/projectdirs/desi/users/ameisner/CI/reduced/v0001/20190403/ci-00003697/ci-00003697_catalog.fits'

    
    gaia = kentools_center(fname_cat, extname=extname, arcmin_max=arcmin_max)

    return gaia

def _test_gfa():
    fname_cat = '/project/projectdirs/desi/users/ameisner/GFA/reduced/v0000/20191116/00028537/gfa-00028537_catalog.fits'

    result = kentools_center(fname_cat, extname='GUIDE2')

    return result

def __test_gfa():
    fname_cat = '/project/projectdirs/desi/users/ameisner/GFA/reduced/v0000/20191116/00028537/gfa-00028537_catalog.fits'

    fname_raw = '/project/projectdirs/desi/spectro/data/20191116/00028537/gfa-00028537.fits.fz'

    h = fits.getheader(fname_raw, extname='GFA')
    
    cat = fits.getdata(fname_cat)

    result = kentools_center(cat, h['SKYRA'], h['SKYDEC'], extname='GUIDE2')

    return result

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
