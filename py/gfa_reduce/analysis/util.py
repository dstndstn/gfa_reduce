import gfa_reduce.common as common
import numpy as np
import os
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.ndimage.interpolation import shift
import astropy.io.fits as fits
from scipy.optimize import minimize
import gfa_reduce.xmatch.gaia as gaia

def use_for_fwhm_meas(tab, bad_amps=None, snr_thresh=20,
                      no_sig_major_cut=False):
    # return a boolean mask indicating whether each source in a catalog
    # should / should not contribute to its image's reported FWHM
    # measurement

    assert(len(np.unique(tab['camera'])) == 1)
        
    good = ((tab['dq_flags'] == 0) & (tab['min_edge_dist_pix'] > 30) &
            (tab['detmap_peak'] >= snr_thresh))

    if not no_sig_major_cut:
        good = (good & (tab['sig_major_pix'] > 1) & np.isfinite(tab['sig_major_pix']))

        # if bad amps specified, it should be a list 
    # containing the amps thought to be in a state of bad readout
    if (bad_amps is not None) and (len(bad_amps) > 0):
        amp_ok = np.array([(t['amp'] not in bad_amps) for t in tab])
        good = (good & amp_ok)

    return good
    
def has_wrong_dimensions(exp):
    # check meant to catch simulated data
    # or other rare anomalies where GFA images do not have the
    # correct dimensions

    for im in exp.images.values():
        if im is None:
            continue
        sh = im.image.shape
        if (sh[0] != 1032) or (sh[1] != 2048):
            return True
    return False
    
def nominal_pixel_area_sq_asec(extname):
    par = common.gfa_misc_params()

    pixel_area_sq_asec = \
        (par['nominal_mer_cd']*3600.0)*(par['nominal_sag_cd']*3600.0)

    return pixel_area_sq_asec

def nominal_pixel_sidelen_arith():
    # calculate/return the nominal pixel sidelength in arcseconds
    # using the arithmetic mean of the x and y platescales

    par = common.gfa_misc_params()

    return np.mean([par['nominal_mer_cd'], par['nominal_sag_cd']])*3600.0

def nominal_pixel_sidelen_geom(extname):
    # calculate/return the nominal pixel sidelength in arcseconds
    # using the geometric mean of the x and y platescales

    return np.sqrt(nominal_pixel_area_sq_asec(extname))

def gfa_pixel_xmax(pix_center=False, quadrant=None):
    """
    "x" here is in GFA pixel coordinates
    could imagine adding a "binfac" keyword here for use in processing
    steps where I've performed an integer downbinning
    """
    par = common.gfa_misc_params()

    # right edge of rightmost pixel
    xmax = par['width_pix_native'] - 0.5

    if pix_center:
        xmax -= 0.5 # center of rightmost pixel

    if (quadrant == 2) or (quadrant == 3):
        # haven't thought about whether assumption of even width matters here
        xmax -= par['width_pix_native']/2

    return xmax

def gfa_pixel_ymax(pix_center=False, quadrant=None):
    """
    "y" here is in GFA pixel coordinates
    """
    par = common.gfa_misc_params()

    # right edge of rightmost pixel
    ymax = par['height_pix_native'] - 0.5

    if pix_center:
        ymax -= 0.5 # center of rightmost pixel

    if (quadrant == 3) or (quadrant == 4):
        # haven't thought about whether assumption of even width matters here
        ymax -= par['height_pix_native']/2

    return ymax

def gfa_pixel_xmin(pix_center=False, quadrant=None):
    """
    "x" here is in GFA pixel coordinates
    """

    # left edge of leftmost pixel
    xmin = -0.5

    if pix_center:
        xmin += 0.5 # center of leftmost pixel

    if (quadrant == 1) or (quadrant == 4):
        par = common.gfa_misc_params()
        # haven't thought about whether assumption of even width matters here
        xmin += par['width_pix_native']/2

    return xmin

def gfa_pixel_ymin(pix_center=False, quadrant=None):
    """
    "y" here is in GFA pixel coordinates
    """

    # left edge of leftmost pixel
    ymin = -0.5

    if pix_center:
        ymin += 0.5 # center of leftmost pixel

    if (quadrant == 1) or (quadrant == 2):
        par = common.gfa_misc_params()
        # haven't thought about whether assumption of even width matters here
        ymin += par['height_pix_native']/2

    return ymin

def gfa_center_pix_coords():
    # native binning, this is the exact center of the image, 
    # which is at the corner of four pixels because of even sidelengths

    par = common.gfa_misc_params()

    x_pix_center = par['width_pix_native']*0.5 + 0.5
    y_pix_center = par['height_pix_native']*0.5 + 0.5

    return x_pix_center, y_pix_center

def gfa_boundary_pixel_coords(pix_center=True):
    par = common.gfa_misc_params()

    x_top = np.arange(gfa_pixel_xmin(pix_center=pix_center), 
                      gfa_pixel_xmax(pix_center=pix_center) + 1)
    x_left = np.zeros(par['height_pix_native'] + 1*(not pix_center)) + \
                      gfa_pixel_xmin(pix_center=pix_center)
    y_left = np.arange(gfa_pixel_ymin(pix_center=pix_center),
                      gfa_pixel_ymax(pix_center=pix_center) + 1)
    y_bottom = np.zeros(par['width_pix_native'] + 1*(not pix_center)) + \
                        gfa_pixel_ymin(pix_center=pix_center)
    y_top = y_bottom + par['height_pix_native'] - 1 + 1*(not pix_center)
    x_right = x_left + par['width_pix_native'] - 1 + 1*(not pix_center)
    y_right = np.flip(y_left, axis=0)
    x_bottom = np.flip(x_top, axis=0)

    x_bdy = np.concatenate((x_left, x_top, x_right, x_bottom))
    y_bdy = np.concatenate((y_left, y_top, y_right, y_bottom))

    return x_bdy, y_bdy

def gfa_corner_pixel_coords(pix_center=False, wrap=False):
    # LL -> UL -> UR -> LR
    x_pix = [gfa_pixel_xmin(pix_center=pix_center),
             gfa_pixel_xmin(pix_center=pix_center), 
             gfa_pixel_xmax(pix_center=pix_center),
             gfa_pixel_xmax(pix_center=pix_center)]

    y_pix = [gfa_pixel_ymin(pix_center=pix_center),
             gfa_pixel_ymax(pix_center=pix_center),
             gfa_pixel_ymax(pix_center=pix_center),
             gfa_pixel_ymin(pix_center=pix_center)]

    if wrap:
        x_pix.append(x_pix[0])
        y_pix.append(y_pix[0])

    return x_pix, y_pix

# should probably also have something available for the case of upbinning
def gfa_downbinned_shape(binfac):
    # assume integer rebinning until I come across a case where
    # arbitrary rebinning would be valuable

    # assume same rebinning factor in both dimensions for now, until
    # I come across a case where I would want otherwise

    assert((type(binfac).__name__ == 'int') or binfac.is_integer())

    par = common.gfa_misc_params()

    width_native = par['width_pix_native']
    height_native = par['height_pix_native']

    width_downbinned = float(width_native)/float(binfac)
    height_downbinned = float(height_native)/float(binfac)

    assert(width_downbinned.is_integer())
    assert(height_downbinned.is_integer())

    # note Python convention for (height, width)
    return int(height_downbinned), int(width_downbinned)

def min_edge_dist_pix(x, y):
    # minimum distance to any image edge
    # for now inputs are meant to be scalar, not array/list

    min_edge_dist = 10000

    min_edge_dist = min(min_edge_dist, x-gfa_pixel_xmin())
    min_edge_dist = min(min_edge_dist, y-gfa_pixel_ymin())
    min_edge_dist = min(min_edge_dist, gfa_pixel_xmax()-x)
    min_edge_dist = min(min_edge_dist, gfa_pixel_ymax()-y)

    return min_edge_dist

def create_det_ids(catalog, extname, fname_in, add_col=True, cube_index=None):
    # catalog should be an astropy table
    # catalog should pertain to just one **extension**
    # watch out for case where there are no extracted sources in an
    # image

    basename = os.path.basename(fname_in)

    # strip out file extension
    basename = basename.replace('.fz', '')
    basename = basename.replace('.gz', '')
    basename = basename.replace('.fits', '')

    det_ids = [('o' + str(i).zfill(6) + 'e' + extname) for i in range(len(catalog))]

    det_ids = [(basename + det_id) for det_id in det_ids]

    if cube_index is not None:
        det_ids = [(det_id + 'g' + str(cube_index).zfill(5)) for det_id in det_ids]
    
    if add_col:
        catalog['det_id'] = det_ids
    else:
        return det_ids

def slice_indices_for_quadrant(quadrant):

    xmin = int(gfa_pixel_xmin(pix_center=True, quadrant=quadrant))
    xmax = int(gfa_pixel_xmax(pix_center=True, quadrant=quadrant)) + 1
    ymin = int(gfa_pixel_ymin(pix_center=True, quadrant=quadrant))
    ymax = int(gfa_pixel_ymax(pix_center=True, quadrant=quadrant)) + 1

    return xmin, xmax, ymin, ymax

def expid_from_raw_filename(fname):
    # fname should be a single string not a list/array of strings
    
    f = os.path.split(fname)[-1]

    f = f.split('-', 1)[1]

    string = f[0:8]

    # hack for special PlateMaker acquisition image file names
    pos = string.find('.')

    if pos != -1:
        string = string[0:pos]
    
    return int(string)

def average_bintable_metadata(tab):

    result = Table()

    result['EXPTIME'] = [np.mean(tab['EXPTIME'])]

    # eventually do a better job of dealing with missing REQTIME column
    try:
        result['REQTIME'] = [np.mean(tab['REQTIME'])]
    except:
        print('no REQTIME column in guide cube binary table???')

    result['NIGHT'] = [tab['NIGHT'][0]]
    
    columns_to_average = ['MJD-OBS',
                          'GAMBNTT',
                          'GFPGAT',
                          'GFILTERT',
                          'GCOLDTEC',
                          'GHOTTEC',
                          'GCCDTEMP',
                          'GCAMTEMP',
                          'GHUMID2',
                          'GHUMID3']

    for col in columns_to_average:
        if col in tab.columns.names:
            # weighted average...
            result[col] = [np.sum(tab[col]*tab['EXPTIME'])/np.sum(tab['EXPTIME'])]

    return result[0]

def sanity_check_catalog(cat):
    # can build more checks into this as time goes on...

    print('Sanity checking source catalog...')
    assert(np.sum(np.isfinite(cat['xcentroid'])) == len(cat))
    assert(np.sum(np.isfinite(cat['ycentroid'])) == len(cat))

def xy_to_ampname(x, y):
    # x and y should be zero-indexed pixel coordinates within
    # the image area (prescan/overscan stripped out)
    # x and y should be scalar

    if (x <= 1023.5) and (y <= 515.5):
        return 'E'
    if (x > 1023.5) and (y <= 515.5):
        return 'F'
    if (x > 1023.5) and (y > 515.5):
        return 'G'
    if (x <= 1023.5) and (y > 515.5):
        return 'H'

    # should never get here...
    assert(False)

def add_ampname_to_catalog(cat):
    # assumes x and y pixel coordinates conform to expectations of
    # xy_to_ampname function above and that they're stored in columns
    # named 'xcentroid' and 'ycentroid'

    amps = [xy_to_ampname(cat['xcentroid'][i], cat['ycentroid'][i]) for i in range(len(cat))]

    cat['amp'] = amps

def moon_separation(ra_moon, dec_moon, ra, dec):

    assert(len(ra_moon) == len(dec_moon))
    assert(len(ra_moon) == len(ra))
    assert(len(ra_moon) == len(dec))

    c = SkyCoord(ra*u.deg, dec*u.deg)
    m = SkyCoord(ra_moon*u.deg, dec_moon*u.deg)

    dangle = c.separation(m)

    return dangle

def _shift_stamp(stamp, dx, dy):
    order = 4
    mode = 'nearest'
    output = stamp.dtype

    # note ordering of dx and dy ...
    return shift(stamp, [dy, dx], order=order, mode=mode, output=output)

def _stamp_radius_mask(sidelen, return_radius=False,
                       x_centroid=None, y_centroid=None):
    # assume square for now
    # assume odd side length

    assert(sidelen == np.round(sidelen))
    assert(sidelen % 2 == 1)
    
    ybox = np.arange(sidelen*sidelen, dtype=int) // sidelen
    xbox = np.arange(sidelen*sidelen, dtype=int) % sidelen

    if x_centroid is None:
        x_centroid = sidelen // 2

    if y_centroid is None:
        y_centroid = sidelen // 2

    xbox = xbox.astype('float')
    ybox = ybox.astype('float')
        
    xbox -= x_centroid
    ybox -= y_centroid

    dist = np.sqrt(np.power(xbox, 2) + np.power(ybox, 2))

    mask = dist.reshape((sidelen, sidelen)) > (sidelen // 2)

    if not return_radius:
        return mask
    else:
        return mask, dist.reshape(sidelen, sidelen)

def _test_shift():
    import astropy.io.fits as fits
    im = fits.getdata('gaussian.fits')

    result = _shift_stamp(im, 0.5, 0.0)
    
def _resize(arr, fac):
    assert(np.round(fac) == fac)

    fac = int(fac)
    assert(fac >= 1)
    
    return np.repeat(np.repeat(arr, fac, axis=0), fac, axis=1)

def _fiber_fracflux(psf, x_centroid=None, y_centroid=None,
                    fib_rad_pix=None):

    # not really sure if this edge case will ever happen ??
    if (np.sum(psf) <= 0):
        return np.nan
        
    binfac = 11

    if fib_rad_pix is None:
        fib_diam = 107.0/15.0 # GFA pixels

        fib_rad = fib_diam/2.0 # GFA pixels
    else:
        fib_rad = fib_rad_pix

    sidelen = psf.shape[0] # assume square..

    if x_centroid is not None:
        x_centroid *= binfac
    if y_centroid is not None:
        y_centroid *= binfac

    mask, radius = _stamp_radius_mask(sidelen*binfac, return_radius=True,
                                      x_centroid=x_centroid,
                                      y_centroid=y_centroid)
    
    _psf = _resize(psf, binfac)

    _mask = np.logical_not(mask)

    in_fiber = (radius <= fib_rad*binfac)
    
    frac = np.sum(_psf*_mask*in_fiber)/np.sum(_psf*_mask)

    return frac

def _aperture_corr_fac(psf, x_centroid=None, y_centroid=None):

    # would be good to not have this hardcoded...
    rad_asec = 1.5 # corresponds to my _3 aperture fluxes

    asec_per_pix = nominal_pixel_sidelen_arith()
    rad_pix = rad_asec/asec_per_pix

    return _fiber_fracflux(psf, x_centroid=x_centroid,
                           y_centroid=y_centroid,
                           fib_rad_pix=rad_pix)
    
def zenith_zeropoint_photometric_1amp(extname, amp):
    par = common.gfa_misc_params()

    fname = os.path.join(os.environ[par['meta_env_var']], par['zp_filename'])

    # would be better to cache this, but it's of order ~10 kb ..
    tab = fits.getdata(fname)

    good = (tab['EXTNAME'] == extname) & (tab['AMP'] == amp)

    assert(np.sum(good) == 1)

    return tab[good][0]['ZP_ADU_PER_S']
    
def median_zenith_camera_zeropoint(extname):
    # wraps zenith_zeropoint_photometric_1amp
    # eventually want to do a better job of dealing with amp-to-amp
    # zeropoint variation so this is hopefully a temporary hack

    amps = common.valid_amps_list()

    zps = [zenith_zeropoint_photometric_1amp(extname, amp) for amp in amps]

    return np.median(zps)

def zp_photometric_at_airmass(extname, airmass, amp=None):

    # for now don't worry about vectorization

    assert(airmass > 0.99) # allow for some roundoff to < 1

    if amp is None:
        zp_zenith = median_zenith_camera_zeropoint(extname)
    else:
        zp_zenith = zenith_zeropoint_photometric_1amp(extname, amp)

    par = common.gfa_misc_params()

    # account for airmass (k term from DESI-5418-v2)
    
    # "photometric" here means 'in photometric conditions' at this airmass
    zp_photometric = zp_zenith - (airmass - 1)*par['kterm']

    return zp_photometric

def transparency_from_zeropoint(zp_image, airmass, extname):
    # zp_image should be the r band magnitude corresponding to a source
    # with total detected flux of 1 ADU per second in the single-camera image
    # of interest

    if not np.isfinite(airmass):
        return np.nan
    if not np.isfinite(zp_image):
        return np.nan
    
    zp_photometric = zp_photometric_at_airmass(extname, airmass)
    
    transparency = 10**((zp_image - zp_photometric)/2.5)

    return transparency

def _gauss2d_profile(sidelen, xcen, ycen, peak_val, sigma, bg=0):
    
    ybox = np.arange(sidelen*sidelen, dtype=int) // sidelen
    xbox = np.arange(sidelen*sidelen, dtype=int) % sidelen

    xbox = xbox.astype('float')
    ybox = ybox.astype('float')

    xcen = float(xcen)
    ycen = float(ycen)
    sigma = float(sigma)
    
    dx2 = np.power(xbox - xcen, 2)
    dy2 = np.power(ybox - ycen, 2)

    r2 = dx2 + dy2

    prof = np.exp(-1.0*r2/(2*(sigma**2)))
    prof = peak_val*prof/np.max(prof)

    prof += bg

    prof = prof.reshape((sidelen, sidelen))

    return prof

def _moffat2d_profile(sidelen, xcen, ycen, peak_val, fwhm, bg=0, beta=3.5):


    alpha = fwhm/(2.0*np.sqrt(2**(1/beta) - 1))
    # beta = 2.5 is the IRAF default value apparently
    
    ybox = np.arange(sidelen*sidelen, dtype=int) // sidelen
    xbox = np.arange(sidelen*sidelen, dtype=int) % sidelen

    xbox = xbox.astype('float')
    ybox = ybox.astype('float')

    xcen = float(xcen)
    ycen = float(ycen)
    alpha = float(alpha)
    
    dx2 = np.power(xbox - xcen, 2)
    dy2 = np.power(ybox - ycen, 2)

    r2 = dx2 + dy2

    prof = np.power(1 + r2/(alpha**2), -1.0*beta)
    prof = peak_val*prof/np.max(prof)

    prof += bg

    prof = prof.reshape((sidelen, sidelen))
    
    return prof

def _moffat2d_metric(p, xcen, ycen, image):
    # p[0] : fwhm (pixels)
    # p[1] : peak value

    sh = image.shape

    assert(sh[0] == sh[1])

    sidelen = sh[0]

    model = _moffat2d_profile(sidelen, xcen, ycen, p[1], p[0])

    return np.sum(np.power(image-model, 2))
    
def _gauss2d_metric(p, xcen, ycen, image):
    # p[0] : sigma (pixels)
    # p[1] : peak value

    sh = image.shape

    assert(sh[0] == sh[1])

    sidelen = sh[0]

    model = _gauss2d_profile(sidelen, xcen, ycen, p[1], p[0])

    return np.sum(np.power(image-model, 2))

def _fit_gauss2d(xcen, ycen, image):

    # would be good to specify the initial simplex here at some point
    res = minimize(_gauss2d_metric, [6.0, 1.0], args=(xcen, ycen, image), method='Nelder-Mead', options={'maxfev': 200, 'disp': False, 'adaptive': False, 'fatol': 1.0e-5})

    return res

def _fit_moffat2d(xcen, ycen, image):
    res = minimize(_moffat2d_metric, [6.0, 1.0], args=(xcen, ycen, image), method='Nelder-Mead', options={'maxfev': 200, 'disp': False, 'adaptive': False, 'fatol': 1.0e-5})

    return res
    
def _test_gauss2d_fit():
    tab = fits.getdata('/global/cfs/cdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid/20200131/00045485/gfa-00045485_ccds.fits')

    psf = fits.getdata('/global/cfs/cdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid/20200131/00045485/gfa-00045485_psfs.fits')

    res = _fit_gauss2d(tab[0]['XCENTROID_PSF'], tab[0]['YCENTROID_PSF'], psf)

    return res

# maybe this belongs in "io" ...
def load_lst():
    par = common.gfa_misc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['ephem_filename'])

    print('READING EPHEMERIS FILE : ', fname)
    assert(os.path.exists(fname))

    eph = fits.getdata(fname)

    return eph

def interp_lst(mjd, eph=None):

    # for now assume that mjd is a scalar, can deal with vectorization later..

    # LST value returned is in degrees
    
    if (mjd is None) or (mjd == 0) or (np.isnan(mjd)):
        return np.nan
        
    if eph is None:
        eph = load_lst()

    ind_upper = np.searchsorted(eph['MJD'], mjd)

    assert(ind_upper > 0)
    assert(ind_upper != len(eph))

    ind_lower = ind_upper - 1

    mjd_upper = eph['MJD'][ind_upper]
    mjd_lower = eph['MJD'][ind_lower]

    assert(mjd_upper >= mjd)
    assert(mjd_lower <= mjd)

    lst_upper = eph['LST_DEG'][ind_upper]
    lst_lower = eph['LST_DEG'][ind_lower]

    if (lst_lower > lst_upper):
        lst_lower -= 360.0

    lst = ((mjd - mjd_lower)*lst_upper + (mjd_upper - mjd)*lst_lower)/(mjd_upper-mjd_lower)
    
    # bound to be within 0 -> 360

    assert(lst < 360)

    if (lst < 0):
        lst += 360.0

    return lst

def _zenith_distance(ra, dec, lst_deg):

    # output value should be in degrees

    if np.isnan(ra) or np.isnan(dec) or np.isnan(lst_deg):
        return np.nan
    
    # for now assume scalar inputs, can work on vectorization later if desired

    par = common.gfa_misc_params()
    
    kpno_latitude =  par['kpno_lat_deg']
    
    c = SkyCoord(ra*u.deg, dec*u.deg)
    zenith = SkyCoord(lst_deg*u.deg, kpno_latitude*u.deg)

    dangle = c.separation(zenith)

    return dangle.deg

def _get_ha(ra_deg, lst_deg, mountdec):

    # add/implement hours boolean kw arg at some point

    # for now assume scalar inputs, can work on vectorization later if desired

    if mountdec > 90:
        ra_deg -= 180.0

    if ra_deg >= 360:
        ra_deg -= 360.0
    elif ra_deg < 0:
        ra_deg += 360.0
    
    if (np.abs(lst_deg - ra_deg) > 180) and ((lst_deg - ra_deg) > 0):
        lst_deg -= 360.0
    if (np.abs(lst_deg - ra_deg) > 180) and ((lst_deg - ra_deg) < 0):
        lst_deg += 360.0

    ha = lst_deg - ra_deg

    return ha
        
def pm_pi_corr_fiberassign(gfa_targets, mjd):
    # correct fiberassing TARGET_RA, TARGET_DEC to
    # relevant DESI observation epoch based on available parallaxes
    # and proper motions

    assert(mjd is not None)

    # fiberassign files appear to use dummy values of 0 for parallax
    # dra is in true angular mas not RA coordinate mas
    dra_pi_mas, ddec_pi_mas = gaia.parallax_offsets(gfa_targets['TARGET_RA'],
                                                    gfa_targets['TARGET_DEC'],
                                                    mjd, gfa_targets['PARALLAX'])

    def year_to_mjd(yr):
        mjd2k = 51545. # J2000
        return mjd2k + 365.25 * (yr - 2000.)

    #assert(np.sum(gfa_targets['REF_EPOCH'] != 2015.5) == 0)
    #ref_mjd =  57205.625 # 2015.5
    assert(np.all(gfa_targets['REF_EPOCH'] > 2000.))
    ref_mjd = year_to_mjd(gfa_targets['REF_EPOCH'])

    ra_corr = gfa_targets['TARGET_RA'] + np.array(dra_pi_mas)/(np.cos(gfa_targets['TARGET_DEC']/(180.0/np.pi))*3600.0*1000.0)
    dec_corr = gfa_targets['TARGET_DEC'] + np.array(ddec_pi_mas)/(3600.0*1000.0)

    ra_corr += ((mjd - ref_mjd)/365.25)*gfa_targets['PMRA']/(np.cos(gfa_targets['TARGET_DEC']/(180.0/np.pi))*3600.0*1000.0)

    dec_corr += ((mjd - ref_mjd)/365.25)*gfa_targets['PMDEC']/(3600.0*1000.0)

    gfa_targets['TARGET_RA'] = ra_corr
    gfa_targets['TARGET_DEC'] = dec_corr
