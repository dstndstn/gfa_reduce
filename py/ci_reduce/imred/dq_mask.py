import ci_reduce.common as common
import ci_reduce.imred.load_calibs as load_calibs
import numpy as np

def create_satmask(im, extname):
    # im is just a 2D array of pixels, not a CI_image object

    par = common.ci_misc_params()

    gain = common.ci_camera_gain(extname)

    sat_thresh = par['full_well_electrons']/gain

    satmask = (im >= sat_thresh)

    return satmask

def create_nanmask(im):
    # im is just a 2D array of pixels, not a CI_image object

    return np.logical_not(np.isfinite(im))

def dq_bitmask(im, extname):
    """
    bit definitions
    2^0 = bad pixels based on master flat (inherit from static bad pixel mask)
    2^1 = questionable based on master flat (inherit from static bad pixel mask)
    2^2 = saturated in image
    2^3 = non-finite values in image
    """

    # im is just a 2D array of pixels, not a CI_image object

    mask = load_calibs.read_static_mask_image(extname)

    satmask = create_satmask(im, extname)
    print(satmask.shape, np.sum(satmask != 0), np.min(satmask), np.max(satmask))
    mask = np.bitwise_or(mask, (2**2)*satmask)

    nanmask = create_nanmask(im)
    mask = np.bitwise_or(mask, (2**3)*nanmask)

    return mask

# should add a function for creating the header needed for writing the
# output dq mask -- in particular the self-documentation of the mask bit
# meanings
