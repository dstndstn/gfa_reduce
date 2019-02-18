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
    print('Populated image extension names are : ' + 
          str(exp.populated_extnames()))

    return exp

def reduced_image_fname(outdir, fname_in, flavor, gzip=True):
    assert(os.path.exists(outdir))

    outname = os.path.join(outdir, os.path.basename(fname_in))

    # get rid of any ".fz" or ".gz" present in input filename
    outname = outname.replace('.fz', '')
    outname = outname.replace('.gz', '')

    assert(outname[-5:] == '.fits')

    outname = outname.replace('.fits', 
        common.reduced_image_filename_label(flavor) + '.fits')

    if gzip:
        outname += '.gz'

    assert(not os.path.exists(outname))

    return outname


def write_image_level_outputs(exp, outdir, fname_in, gzip=True):
    # exp is a CI_exposure object
    # outdir is the output directory (string)
    # fname_in is the input filename (string)

    par = common.ci_misc_params()

    for flavor in par['reduced_image_flavors']:
        outname = reduced_image_fname(outdir, fname_in, flavor, gzip=gzip)

        hdulist = exp.to_hdulist(flavor=flavor)

        print('Attempting to write ' + flavor + ' image output to ' + 
              outname)

        hdulist.writeto(outname)

        print('Successfully wrote ' + flavor + ' image output to ' + 
              outname)
