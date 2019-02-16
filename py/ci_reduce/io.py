from ci_reduce.image import CI_image
from ci_reduce.exposure import CI_exposure
import ci_reduce.common as common
import astropy.io.fits as fits
import os

# in the context of this file, "image" and "exposure" generally refer to 
# CI_image and CI_exposure objects

def loading_image_extension_message(extname):
    print('Attempting to load image extension : ' + extname)

def load_image_from_hdu(hdu):
    loading_image_extension_message(hdu.header['EXTNAME'])
    return CI_image(hdu.data, hdu.header)

def load_image_from_filename(fname, extname):
    assert(os.path.exists(fname))

    loading_image_extension_message(extname)
    assert(common.is_valid_image_extname(extname))

    data, header = fits.getdata(fname, extname=extname, header=True)
    return CI_image(data, header)

def load_exposure(fname):
    assert(os.path.exists(fname))

    print('Attempting to load exposure : ' + fname)

    par = common.ci_misc_params()

    hdul = fits.open(fname)

    dummy_fz_header = None

    for hdu in hdul:
        if (hdu.header['EXTNAME']).strip() == par['fz_dummy_extname']:
            dummy_fz_header = hdu.header
            hdul.remove(hdu)

    imlist = [load_image_from_hdu(hdu) for hdu in hdul]

    exp = CI_exposure(imlist, dummy_fz_header=dummy_fz_header)

    print('Successfully loaded exposure : ' + fname)
    print('Exposure has ' + str(exp.num_images_populated()) + 
          ' image extensions populated')

    return exp
