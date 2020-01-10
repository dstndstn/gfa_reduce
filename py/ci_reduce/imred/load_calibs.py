import ci_reduce.common as common
import astropy.io.fits as fits
import os
import numpy as np

def remove_overscan(image):
    sh = image.shape
    if sh[1] == 2248:
        _image = np.zeros((1032, 2048), dtype=float)
        _image[:, 0:1024] = image[:, 50:1074]
        _image[:, 1024:2048] = image[:, 1174:2198]
    return _image


def read_bias_image(ci_extname):
    assert(common.is_valid_extname(ci_extname))

    par = common.ci_misc_params()
    bias_fname = os.path.join(os.environ[par['etc_env_var']], \
                              par['master_bias_filename'])

    print('Attempting to read master bias : ' + bias_fname + 
          ', extension name : ' + ci_extname)

    assert(os.path.exists(bias_fname))

    bias = fits.getdata(bias_fname, extname=ci_extname)

    bias = remove_overscan(bias)
    return bias

def read_flat_image(ci_extname):
    # at some point should add option to return master flat's
    # inverse variance as well
    assert(common.is_valid_extname(ci_extname))

    par = common.ci_misc_params()
    flat_fname = os.path.join(os.environ[par['etc_env_var']], \
                              par['master_flat_filename'])

    print('Attempting to read master flat : ' + flat_fname + 
          ', extension name : ' + ci_extname)

    assert(os.path.exists(flat_fname))

    flat = fits.getdata(flat_fname, extname=ci_extname)

    flat = remove_overscan(flat)
    return flat

def read_static_mask_image(ci_extname):
    assert(common.is_valid_extname(ci_extname))

    par = common.ci_misc_params()
    mask_fname = os.path.join(os.environ[par['etc_env_var']], \
                              par['static_mask_filename'])

    print('Attempting to read static bad pixel mask : ' + mask_fname + 
          ', extension name : ' + ci_extname)

    assert(os.path.exists(mask_fname))

    mask = fits.getdata(mask_fname, extname=ci_extname)

    mask = remove_overscan(mask)
    return mask
