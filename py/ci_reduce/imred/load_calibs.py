import ci_reduce.common as common
import astropy.io.fits as fits
import os

def read_bias_image(ci_extname):
    assert(common.is_valid_extname(ci_extname))

    par = common.ci_misc_params()
    bias_fname = os.path.join(os.environ[par['etc_env_var']], \
                              par['master_bias_filename'])

    assert(os.path.exists(bias_fname))

    bias = fits.getdata(bias_fname, extname=ci_extname)

    return bias
