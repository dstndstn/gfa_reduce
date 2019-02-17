import ci_reduce.common as common
import ci_reduce.imred.dq_mask as dq_mask
import numpy as np
from astropy import wcs

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
