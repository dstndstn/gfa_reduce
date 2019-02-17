import ci_reduce.common as common
import ci_reduce.imred.dq_mask as dq_mask
from astropy import wcs

class CI_image:
    """Single CI image from one CI exposure"""

    def __init__(self, image, header):
        self.image = image
        self.header = header
        self.wcs = wcs.WCS(header)
        # lazily compute bitmask image only as requested
        self.bitmask = None

    def create_dq_mask(self):
        if self.bitmask is not None:
            return

        self.bitmask = dq_mask.dq_bitmask(self.image, self.header['EXTNAME'])
