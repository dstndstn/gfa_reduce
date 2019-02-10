import ci_reduce.common as common
import numpy as np

def ci_pixel_xmax(pix_center=False):
    """
    "x" here is in CI pixel coordinates
    could imagine adding a "binfac" keyword here for use in processing
    steps where I've performed an integer downbinning
    """
    par = common.ci_misc_params()

    # right edge of rightmost pixel
    xmax = par['width_pix_native'] - 0.5

    if pix_center:
        xmax -= 0.5 # center of rightmost pixel

    return xmax

def ci_pixel_ymax(pix_center=False):
    """
    "y" here is in CI pixel coordinates
    """
    par = common.ci_misc_params()

    # right edge of rightmost pixel
    ymax = par['height_pix_native'] - 0.5

    if pix_center:
        ymax -= 0.5 # center of rightmost pixel

    return ymax

def ci_pixel_xmin(pix_center=False):
    """
    "x" here is in CI pixel coordinates
    """

    # left edge of leftmost pixel
    xmin = -0.5

    if pix_center:
        xmin += 0.5 # center of leftmost pixel

    return xmin

def ci_pixel_ymin(pix_center=False):
    """
    "y" here is in CI pixel coordinates
    """

    # left edge of leftmost pixel
    ymin = -0.5

    if pix_center:
        ymin += 0.5 # center of leftmost pixel

    return ymin

def ci_center_pix_coords():
    # native binning, this is the exact center of the image, 
    # which is at the corner of four pixels because of even sidelengths

    par = common.ci_misc_params()

    x_pix_center = par['width_pix_native']*0.5 + 0.5
    y_pix_center = par['height_pix_native']*0.5 + 0.5

    return x_pix_center, y_pix_center

def ci_boundary_pixel_coords(pix_center=True):
    par = common.ci_misc_params()

    x_top = np.arange(ci_pixel_xmin(pix_center=pix_center), 
                      ci_pixel_xmax(pix_center=pix_center) + 1)
    x_left = np.zeros(par['height_pix_native'] + 1*(not pix_center)) + \
                      ci_pixel_xmin(pix_center=pix_center)
    y_left = np.arange(ci_pixel_ymin(pix_center=pix_center),
                      ci_pixel_ymax(pix_center=pix_center) + 1)
    y_bottom = np.zeros(par['width_pix_native'] + 1*(not pix_center)) + \
                        ci_pixel_ymin(pix_center=pix_center)
    y_top = y_bottom + par['height_pix_native'] - 1 + 1*(not pix_center)
    x_right = x_left + par['width_pix_native'] - 1 + 1*(not pix_center)
    y_right = np.flip(y_left, axis=0)
    x_bottom = np.flip(x_top, axis=0)

    x_bdy = np.concatenate((x_left, x_top, x_right, x_bottom))
    y_bdy = np.concatenate((y_left, y_top, y_right, y_bottom))

    return x_bdy, y_bdy

def ci_corner_pixel_coords(pix_center=False, wrap=False):
    # LL -> UL -> UR -> LR
    x_pix = [ci_pixel_xmin(pix_center=pix_center),
             ci_pixel_xmin(pix_center=pix_center), 
             ci_pixel_xmax(pix_center=pix_center),
             ci_pixel_xmax(pix_center=pix_center)]

    y_pix = [ci_pixel_ymin(pix_center=pix_center),
             ci_pixel_ymax(pix_center=pix_center),
             ci_pixel_ymax(pix_center=pix_center),
             ci_pixel_ymin(pix_center=pix_center)]

    if wrap:
        x_pix.append(x_pix[0])
        y_pix.append(y_pix[0])

    return x_pix, y_pix

# should probably also have something available for the case of upbinning
def ci_downbinned_shape(binfac):
    # assume integer rebinning until I come across a case where
    # arbitrary rebinning would be valuable
    # the native CI image dimensions are very favorable for integer rebinning

    # assume same rebinning factor in both dimensions for now, until
    # I come across a case where I would want otherwise

    assert((type(binfac).__name__ == 'int') or binfac.is_integer())

    par = common.ci_misc_params()

    width_native = par['width_pix_native']
    height_native = par['height_pix_native']

    width_downbinned = float(width_native)/float(binfac)
    height_downbinned = float(height_native)/float(binfac)

    assert(width_downbinned.is_integer())
    assert(height_downbinned.is_integer())

    # note Python convention for (height, width)
    return int(height_downbinned), int(width_downbinned)
