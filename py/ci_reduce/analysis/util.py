import ci_reduce.common as common
import numpy as np
import os

def nominal_pixel_area_sq_asec(extname):
    par = common.ci_misc_params()

    if extname == 'CIC':
        pixel_area_sq_asec = (par['nominal_cen_cd']*3600.0)**2
    else:
        pixel_area_sq_asec = \
            (par['nominal_mer_cd']*3600.0)*(par['nominal_sag_cd']*3600.0)

    return pixel_area_sq_asec

def nominal_pixel_sidelen_arith(extname):
    # calculate/return the nominal pixel sidelength in arcseconds
    # using the arithmetic mean of the x and y platescales

    par = common.ci_misc_params()

    if extname == 'CIC':
        return par['nominal_cen_cd']*3600.0
    else:
        return np.mean([par['nominal_mer_cd'], par['nominal_sag_cd']])*3600.0

def nominal_pixel_sidelen_geom(extname):
    # calculate/return the nominal pixel sidelength in arcseconds
    # using the geometric mean of the x and y platescales

    return np.sqrt(nominal_pixel_area_sq_asec(extname))

def ci_pixel_xmax(pix_center=False, quadrant=None):
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

    if (quadrant == 2) or (quadrant == 3):
        # haven't thought about whether assumption of even width matters here
        xmax -= par['width_pix_native']/2

    return xmax

def ci_pixel_ymax(pix_center=False, quadrant=None):
    """
    "y" here is in CI pixel coordinates
    """
    par = common.ci_misc_params()

    # right edge of rightmost pixel
    ymax = par['height_pix_native'] - 0.5

    if pix_center:
        ymax -= 0.5 # center of rightmost pixel

    if (quadrant == 3) or (quadrant == 4):
        # haven't thought about whether assumption of even width matters here
        ymax -= par['height_pix_native']/2

    return ymax

def ci_pixel_xmin(pix_center=False, quadrant=None):
    """
    "x" here is in CI pixel coordinates
    """

    # left edge of leftmost pixel
    xmin = -0.5

    if pix_center:
        xmin += 0.5 # center of leftmost pixel

    if (quadrant == 1) or (quadrant == 4):
        par = common.ci_misc_params()
        # haven't thought about whether assumption of even width matters here
        xmin += par['width_pix_native']/2

    return xmin

def ci_pixel_ymin(pix_center=False, quadrant=None):
    """
    "y" here is in CI pixel coordinates
    """

    # left edge of leftmost pixel
    ymin = -0.5

    if pix_center:
        ymin += 0.5 # center of leftmost pixel

    if (quadrant == 1) or (quadrant == 2):
        par = common.ci_misc_params()
        # haven't thought about whether assumption of even width matters here
        ymin += par['height_pix_native']/2

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

def min_edge_dist_pix(x, y):
    # minimum distance to any image edge
    # for now inputs are meant to be scalar, not array/list

    min_edge_dist = 10000

    min_edge_dist = min(min_edge_dist, x-ci_pixel_xmin())
    min_edge_dist = min(min_edge_dist, y-ci_pixel_ymin())
    min_edge_dist = min(min_edge_dist, ci_pixel_xmax()-x)
    min_edge_dist = min(min_edge_dist, ci_pixel_ymax()-y)

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

    xmin = int(ci_pixel_xmin(pix_center=True, quadrant=quadrant))
    xmax = int(ci_pixel_xmax(pix_center=True, quadrant=quadrant)) + 1
    ymin = int(ci_pixel_ymin(pix_center=True, quadrant=quadrant))
    ymax = int(ci_pixel_ymax(pix_center=True, quadrant=quadrant)) + 1

    return xmin, xmax, ymin, ymax
