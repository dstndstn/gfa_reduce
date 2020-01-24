import ci_reduce.common as common
import ci_reduce.imred.dq_mask as dq_mask
import ci_reduce.analysis.sky as sky
import ci_reduce.analysis.segment as segment
import ci_reduce.analysis.phot as phot
from ci_reduce.ci_wcs import nominal_tan_wcs
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
from astropy.stats import mad_std
from scipy.stats import scoreatpercentile
import ci_reduce.analysis.util as util

class PSF:
    def __init__(self, cube, extname):
        self.psf_image = np.median(cube, 2) # failure for nstars = 1?

        sh = self.psf_image.shape
        bgmask = util._stamp_radius_mask(sh[0])
        self.psf_image -= np.median(self.psf_image[bgmask])

        self.psf_image /= np.max(self.psf_image)
        
        self.extname = extname
        self.nstars = cube.shape[2]
        self.cube = cube # maybe get rid of this eventually to save memory

    def to_hdu(self, primary=False):
         # still need to put metadata into header
         f = (fits.PrimaryHDU if primary else fits.ImageHDU)
         hdu = f(self.psf_image)
         hdu.header['EXTNAME'] = self.extname
         hdu.header['NSTARS'] = self.nstars
         return hdu

class Overscan:
    """Object to encapsulate single-camera worth of overscan and prescan"""

    def __init__(self, image):
        # image should be a 2D numpy array with dimensions
        # 2248 x 1032 in the case of DESI GFA cameras

        par = common.ci_misc_params()

        sh = image.shape
        assert(sh[0] == par['height_with_prescan_overscan'])
        assert(sh[1] == par['width_with_prescan_overscan'])

        amps = common.valid_amps_list()

        self.overscan_cutouts = {}
        self.prescan_cutouts = {}

        for amp in amps:
            bdy = common.overscan_bdy_coords(amp)
            self.overscan_cutouts[amp] = image[bdy['y_l']:bdy['y_u'], bdy['x_l']:bdy['x_u']]
            bdy = common.prescan_bdy_coords(amp)
            self.prescan_cutouts[amp] = image[bdy['y_l']:bdy['y_u'], bdy['x_l']:bdy['x_u']]

        self.n_badpix_overscan = self.count_badpixels()
        self.n_badpix_prescan = self.count_badpixels(prescan=True)

        # still per-amp but summing prescan and overscan counts together
        self.n_badpix = dict([(amp, self.n_badpix_overscan[amp] + self.n_badpix_prescan[amp]) for amp in amps])

        # including all amps and lumping together prescan and overscan
        self.n_badpix_all = np.sum([n for n in self.n_badpix.values()])

        # units are raw ADU
        self.overscan_medians = dict([(amp, np.median(self.overscan_cutouts[amp])) for amp in amps])

        # units are raw ADU
        self.prescan_medians = dict([(amp, np.median(self.prescan_cutouts[amp])) for amp in amps])

    def count_badpixels(self, thresh=10000, prescan=False):
        amps = common.valid_amps_list()
        
        if prescan:
            return dict([(amp, int(np.sum(self.prescan_cutouts[amp] > thresh))) for amp in amps])
        else:
            return dict([(amp, int(np.sum(self.overscan_cutouts[amp] > thresh))) for amp in amps])
    
class CI_image:
    """Single CI image from one CI exposure"""

    def __init__(self, image, header, cube_index=None):
        if cube_index is None:
            self.image = image.astype('float32')
            self.nframe = 1
        elif cube_index == -1:
            self.image = np.mean(image.astype('float32'), axis=0)
            self.nframe = image.shape[0]
        else:
            self.image = image[cube_index, :, :].astype('float32')
            self.nframe = 1

        self.overscan = Overscan(self.image)
        self.remove_overscan()

        par = common.ci_misc_params()
        # important that this saturation threshold be done on truly raw image..
        self.satmask = (self.image > par['sat_thresh_adu']).astype('byte')
        
        self.cube_index = cube_index
        self.header = header
        self.extname = self.header['EXTNAME'].replace(' ', '')
        self.header['CONTRAST'] = 0.0 # typically overwritten with actual value
        self.initialize_wcs()

        # lazily compute bitmask image only as requested
        self.bitmask = None

        # may want to be more clever about this to allow for the possibility
        # of loading in already reduced images in the future
        self.bias_subtracted = False
        self.dark_subtracted = False
        self.flatfielded = False

        self.var_e_sq = None
        # the _adu in ivar_adu is meant to indicate the units are 1/(ADU^2)
        self.ivar_adu = None
        self.sky_mag = None
        self.segmap = None # may want to get rid of this entirely
        self.empirical_bg_sigma = None
        self.sky_level_adu = None
        self.bintable_row = None

        # record CCD temperature used for dark current removal
        self.t_c_for_dark = None

        # should become a boolean indicating whether CCD temperature
        # used for dark removal was a guess (1) or is thought to
        # be correct (0)
        self.t_c_for_dark_is_guess = None
        
        # record exposure time used for dark current removal
        self.time_s_for_dark = None

        self.psf = None
        
    def create_dq_mask(self, dark_image):
        if self.bitmask is not None:
            return

        d = common.mask_bit_dict()
        thresh = scoreatpercentile(dark_image, 99.5)
        if self.bitmask is None:
            self.bitmask = ((dark_image > thresh)*(2**d['HOTDARK'])).astype('byte')
        else:
            self.bitmask += ((dark_image > thresh)*(2**d['HOTDARK'])).astype('byte')
        self.bitmask += self.satmask*(2**d['SATUR'])

        self.bitmask = self.bitmask.astype('byte') # just to make sure...
        del self.satmask

    def update_bitmask_flat(self, flatfield):
        # doing this to avoid having to keep flatfield images in memory

        thresh = 0.6 # very little thought put into this choice..
        d = common.mask_bit_dict()
        if self.bitmask is None:
            self.bitmask = ((flatfield < thresh)*(2**d['FLATBAD'])).astype('byte')
        else:
            self.bitmask += ((flatfield < thresh)*(2**d['FLATBAD'])).astype('byte')
        
    def calc_variance_e_squared(self):
        # at this stage the image ought to have been bias subtracted
        # but not flatfielded or dark subtracted

        assert(self.bias_subtracted)

        assert((not self.dark_subtracted) and (not self.flatfielded))

        ci_extname = self.header['EXTNAME']
        gain = common.ci_camera_gain(ci_extname)

        var_e_sq = (common.ci_camera_readnoise(ci_extname)**2 + \
                    self.image*(self.image >= 0)*gain)
        var_e_sq /= self.nframe

        assert(np.sum(var_e_sq <= 0) == 0)
        self.var_e_sq = var_e_sq

    def calc_variance_adu(self, flatfield):
        assert(self.flatfielded)

        ci_extname = self.header['EXTNAME']
        gain = common.ci_camera_gain(ci_extname)

        variance_adu_sq = self.var_e_sq/(gain**2)

        del self.var_e_sq
        self.var_e_sq = None

        variance_adu_sq *= (flatfield**2)

        assert(np.sum(variance_adu_sq <= 0) == 0)

        # note that I'm not currently taking into account uncertainty due
        # to the uncertainty on the flatfield; not sure if doing so would
        # be desirable eventually

        ivar_adu = 1.0/variance_adu_sq

        assert(self.bitmask is not None)

        # zero out inverse variance of pixels with anything flagged in
        # data quality bitmask
        ivar_adu *= (self.bitmask == 0)

        self.ivar_adu = ivar_adu

    def to_hdu(self, primary=False, flavor=''):
        # convert this image to an HDU

        # currently expect flavor to be one of
        #     REDUCED - reduced image
        #     BITMASK - data quality bitmask
        #     INVVAR - inverse variance image

        # if no flavor is specified then assume ".image" attribute is desired
        # data for this HDU

        f = (fits.PrimaryHDU if primary else fits.ImageHDU)

        if (flavor == '') or (flavor == 'REDUCED'):
            hdu = f(self.image.astype('float32'), header=self.header)
        elif (flavor == 'BITMASK'):
            hdu = f(self.bitmask.astype('int'), header=self.header)
            hdu.header = dq_mask.add_dq_bitmask_header_cards(hdu.header)
        elif (flavor == 'INVVAR'):
            hdu = f(self.ivar_adu.astype('float32'), header=self.header)

        hdu.header['FLAVOR'] = flavor

        ci_extname = self.header['EXTNAME']
        gain = common.ci_camera_gain(ci_extname)
        hdu.header['GAINA'] = (gain, '[e-/ADU] assumed gain')

        hdu.header['BUNIT'] = common.reduced_flavor_to_bunit(flavor)

        ci_num = common.ci_extname_to_ci_number(ci_extname)
        hdu.header['CINUM'] = (ci_num, 'CI# number from DESI-3347')

        return hdu

    def are_pixels_calibrated(self):
        return (self.bias_subtracted and self.dark_subtracted and 
                self.flatfielded)


    def estimate_sky_level(self, careful_sky=False):
        # do something dumb for now, return to this later with something
        # more sophisticated, possibly involving segmentation
        # and/or finding the mode

        assert(self.are_pixels_calibrated())

        if careful_sky:
            if self.segmap is None:
                self.set_segmap()

            self.sky_level_adu = np.median(self.image[self.segmap.array == 0])
        else:
            self.sky_level_adu = np.median(self.image)

        return self.sky_level_adu

    def compute_segmap(self):
        print('Attempting to compute segmentation map for ' + 
              self.header['EXTNAME'])

        segmap = segment.segmentation_map(self.image, self.header['EXTNAME'])
        return segmap

    def set_segmap(self):
        self.segmap = self.compute_segmap()

    def compute_empirical_bg_sigma(self, careful_sky=False):

        if careful_sky:
            if self.segmap is None:
                self.set_segmap()

            # this could go wrong in pathological case that
            # segmap is nonzero for all pixels
            return mad_std(self.image[self.segmap.array == 0])
        else:
            return mad_std(self.image)

    def set_empirical_bg_sigma(self, careful_sky=False):
        if self.empirical_bg_sigma is None:
            self.empirical_bg_sigma = self.compute_empirical_bg_sigma(careful_sky=careful_sky)

    def estimate_sky_mag(self, careful_sky=False):
        # calculate sky brightness in mag per sq asec
        # this is meant to be run on the reduced image in ADU

        assert(self.are_pixels_calibrated())

        sky_adu_per_pixel = self.estimate_sky_level(careful_sky=careful_sky)

        extname = self.header['EXTNAME']

        try:
            acttime = self.header['EXPTIME']
        except:
            print('could not find EXPTIME keyword !!!!')
            acttime = self.header['REQTIME']
        
        sky_mag = sky.adu_to_surface_brightness(sky_adu_per_pixel, 
                                                acttime, extname)

        print(self.header['EXTNAME'] + ' sky mag per square asec AB : ' +  
              '{:.3f}'.format(sky_mag))

        self.sky_mag = sky_mag

        return sky_mag

    def catalog_add_radec(self, catalog):
        # use wcs attribute to convert from pixel to world coordinates
        # be careful about 1-indexed versus 0-indexed convention
        # also be careful about any swapping of x and y pixel coordinates
        
        ra, dec = self.wcs.all_pix2world(catalog['xcentroid'],
                                         catalog['ycentroid'], 0)

        catalog['ra'] = ra
        catalog['dec'] = dec

        return catalog

    def catalog_sources(self):
        print('Attempting to catalog sources in ' + self.header['EXTNAME'] +  
              ' image')

        tab = phot.get_source_list(self.image, self.bitmask, 
                                   self.header['EXTNAME'], self.ivar_adu)

        n_sources = (len(tab) if tab is not None else 0)

        print('Found ' + str(n_sources) + ' sources in ' +
              self.header['EXTNAME'] + ' image')

        if tab is None:
            return tab

        tab = self.catalog_add_radec(tab)

        mjd_obs = self.try_retrieve_meta_keyword('MJD-OBS')
        if mjd_obs is None:
            print('could not find MJD-OBS header keyword !!!')
        else:
            tab['mjd_obs'] = mjd_obs

        util.add_ampname_to_catalog(tab)
        util.sanity_check_catalog(tab)
        return tab

    def initialize_wcs(self):
        telra = self.header['SKYRA']
        teldec = self.header['SKYDEC']
        extname = self.header['EXTNAME']

        print('Attempting to initialize WCS guess for ' + extname)
        self.wcs = nominal_tan_wcs(telra, teldec, extname)

    def remove_overscan(self):
        sh = self.image.shape
        if sh[1] == 2248:
            _image = np.zeros((1032, 2048), dtype=float)
            _image[:, 0:1024] = self.image[:, 50:1074]
            _image[:, 1024:2048] = self.image[:, 1174:2198]
            self.image = _image

    def update_wcs(self, d):
        # d is a dictionary with xshift_best, yshift_best, astr_guess
        # for this EXTNAME

        assert(d['extname'] == self.header['EXTNAME'])

        wcs = d['astr_guess']
        wcs.wcs.crpix = wcs.wcs.crpix + np.array([d['xshift_best'], d['yshift_best']])
        self.wcs = wcs

        # also want to update the header

        new_wcs_header_cards = wcs.to_header()
        new_wcs_header_cards['CONTRAST'] = d['contrast']

        for k,v in new_wcs_header_cards.items():
            if k == 'LATPOLE':
                continue
            self.header[k] = v

        self.header['CD1_1'] = self.header['PC1_1']
        self.header['CD2_1'] = self.header['PC2_1']
        self.header['CD1_2'] = self.header['PC1_2']
        self.header['CD2_2'] = self.header['PC2_2']

        del self.header['PC1_1']
        del self.header['PC2_1']
        del self.header['PC1_2']
        del self.header['PC2_2']

    def try_retrieve_meta_keyword(self, keyword, placeholder=None):
        # examples are MJD-OBS and GCCDTEMP, which
        # are found in different places in the raw data depending
        # on gfa*.fits.fz versus guide*.fits.fz

        # because guider cube metadata has evolved over time, won't
        # always be guaranteed to get e.g., GCCDTEMP at all
        
        # first look in the image header (could be dangerous for EXPTIME
        # in the case of guider cubes)
        
        if keyword in self.header.keys():
            return self.header[keyword]
        elif self.bintable_row is not None:
            return self.bintable_row[keyword]
        else:
            print('could not find ' + keyword + ' !!')
            return placeholder
            
        
    def extract_psf_cutouts(self, _catalog, sidelen=51):
        # sidelen should be an integer...
        assert(np.round(sidelen) == sidelen)

        half = sidelen // 2
        keep = _catalog['used_for_fwhm_meas'].astype('bool') & (_catalog['min_edge_dist_pix'] > (half + 0.5)) & (_catalog['camera'] == self.extname) & (_catalog['aper_sum_bkgsub_3'] > 0)
        
        if np.sum(keep) == 0:
            return None

        n = np.sum(keep)
        #cube = np.zeros((sidelen, sidelen, n))

        catalog = _catalog[keep]

        assert(np.sum(catalog['aper_sum_bkgsub_3'] <= 0) == 0)

        cutouts = []
        for i in range(n):
            ixcentroid = int(np.round(catalog['xcentroid'][i]))
            iycentroid = int(np.round(catalog['ycentroid'][i]))
            cutout = self.image[(iycentroid-half):(iycentroid+half+1),
                                (ixcentroid-half):(ixcentroid+half+1)]

            # background subtract
            bgmask = util._stamp_radius_mask(sidelen)

            bg = np.median(cutout[bgmask])

            cutout -= bg
            
            cutout = cutout/catalog['aper_sum_bkgsub_3'][i]

            # hack to try removing saturated sources
            if np.sum(cutout >= 30000.0) > 1:
                continue
            cutouts.append(cutout)
            
            # background subtraction
            # subpixel shifting
            # outlier rejection

        ncutouts = len(cutouts)
        
        if ncutouts == 0:
            return None

        cube = np.zeros((sidelen, sidelen, ncutouts))
        for i, cutout in enumerate(cutouts):
            cube[:, :, i] = cutout
        
        return cube

    def create_psf(self, catalog, sidelen=51):
        cube = self.extract_psf_cutouts(catalog, sidelen=sidelen)
        print('computing PSF for ' + self.extname)
        if cube is None:
            self.psf = None
        self.psf = PSF(cube, self.extname)        
