import ci_reduce.common as common
import numpy as np
import ci_reduce.analysis.util as util

def adu_to_surface_brightness(sky_adu_1pixel, acttime, extname):
    """
    convert from ADU (per pixel) to mag per square asec (AB)

    note that this is meant to be applied to an average sky value across
    an entire CI camera; this function does not take into account
    platescale variations within a camera
    """

    if (sky_adu_1pixel <= 0) or (acttime <= 0):
        return np.nan

    par = common.ci_misc_params()

    pixel_area_sq_asec = util.nominal_pixel_area_sq_asec(extname)

    sky_adu_per_sq_asec = sky_adu_1pixel/pixel_area_sq_asec

    sky_adu_per_sec_sq_asec = sky_adu_per_sq_asec/acttime

    sky_e_per_sec_sq_asec = sky_adu_per_sec_sq_asec*common.ci_camera_gain(extname)

    return (par['nominal_zeropoint'] - 2.5*np.log10(sky_e_per_sec_sq_asec))
    
