import ci_reduce.common as common
import ci_reduce.imred.dq_mask as dq_mask
import ci_reduce.analysis.sky as sky
import ci_reduce.analysis.segment as segment
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
from astropy.stats import mad_std

class CI_image:
    """Single CI image from one CI exposure"""

    def __init__(self, image, header):
        self.image = image.astype('float')
        self.header = header
        self.wcs = wcs.WCS(header)

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
        self.segmap = None
        self.empirical_bg_sigma = None

    def create_dq_mask(self):
        if self.bitmask is not None:
            return

        self.bitmask = dq_mask.dq_bitmask(self.image, self.header['EXTNAME'])

    def calc_variance_e_squared(self):
        # at this stage the image ought to have been bias subtracted
        # but not flatfielded or dark subtracted

        assert(self.bias_subtracted)

        assert((not self.dark_subtracted) and (not self.flatfielded))

        ci_extname = self.header['EXTNAME']
        gain = common.ci_camera_gain(ci_extname)

        var_e_sq = (common.ci_camera_readnoise(ci_extname)**2 + \
                    self.image*(self.image >= 0)*gain)

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
            hdu = f(self.image, header=self.header)
        elif (flavor == 'BITMASK'):
            hdu = f(self.bitmask, header=self.header)
        elif (flavor == 'INVVAR'):
            hdu = f(self.ivar_adu, header=self.header)

        hdu.header['FLAVOR'] = flavor

        return hdu

    def are_pixels_calibrated(self):
        return (self.bias_subtracted and self.dark_subtracted and 
                self.flatfielded)


    def estimate_sky_level(self):
        # do something dumb for now, return to this later with something
        # more sophisticated, possibly involving segmentation
        # and/or finding the mode

        assert(self.are_pixels_calibrated())

        if self.segmap is None:
            self.set_segmap()

        return np.median(self.image[self.segmap.array == 0])

    def compute_segmap(self):
        print('Attempting to compute segmentation map for ' + 
              self.header['EXTNAME'])

        segmap = segment.segmentation_map(self.image, self.header['EXTNAME'])
        return segmap

    def set_segmap(self):
        self.segmap = self.compute_segmap()

    def compute_empirical_bg_sigma(self):
        if self.segmap is None:
            self.set_segmap()

        # this could go wrong in pathological case that
        # segmap is nonzero for all pixels
        self.empirical_bg_sigma = mad_std(image[self.segmap.array == 0])

    def estimate_sky_mag(self):
        # calculate sky brightness in mag per sq asec
        # this is meant to be run on the reduced image in ADU

        assert(self.are_pixels_calibrated())

        sky_adu_per_pixel = self.estimate_sky_level()

        extname = self.header['EXTNAME']
        acttime = self.header['ACTTIME']

        sky_mag = sky.adu_to_surface_brightness(sky_adu_per_pixel, 
                                                acttime, extname)

        print(self.header['EXTNAME'] + ' sky mag per square asec AB : ' +  
              '{:.3f}'.format(sky_mag))

        self.sky_mag = sky_mag

        return sky_mag
