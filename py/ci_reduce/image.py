import astropy.io.fits as fits
import ci_reduce.common as common
from astropy import wcs

class CI_image:
    """Single CI image from one CI exposure"""

    def __init__(self, image, header):
        self.image = image
        self.header = header
        self.wcs = wcs.WCS(header)
