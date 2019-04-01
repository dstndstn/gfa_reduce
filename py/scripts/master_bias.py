import ci_reduce.common as common
import astropy.io.fits as fits
import glob
import numpy as np
import os

# currently I have access to 4 independent bias exposures for each of
# CIE, CIN, CIS, and CIW, and none for CIC
# these available data come from "forDK.tar.gz" samples, where
# CIC was in "simulation mode"

engdir = '/project/projectdirs/desi/spectro/data/20190330'

def get_bias_frame_names():
    flist = glob.glob(engdir + '/*/ci*.fits.fz')
    flist = np.array(flist)

    exptime = np.zeros(len(flist))
    flavor = []
    for i, f in enumerate(flist):
        h = fits.getheader(f, extname='CI')
        exptime[i] = h['EXPTIME']
        flavor.append(h['FLAVOR'])

    flavor = np.array(flavor)
    flist = flist[(exptime == 0) & (flavor == 'zero')]
    
    flist.sort()
    return flist

def read_bias_frames(ci_extname):
    flist = get_bias_frame_names()
    imgs = [fits.getdata(f, extname=ci_extname) for f in flist]

    return imgs

def master_bias_1camera(ci_extname):
    # use common to check that CI extention name is valid

    bias_frames = read_bias_frames(ci_extname)
    bias_frames = np.asarray(bias_frames)
    bias_med = np.median(bias_frames, axis=0)
    return bias_med

def master_bias_header_cards(hdu, ci_extname):
    h = hdu.header
    h['FLAVOR'] = 'ZERO'
    h['EXTNAME'] = ci_extname
    h['EXPTIME'] = 0.0
    # zero idea whatsoever what temperature the forDK.tar.gz bias frames were 
    # acquired at; just putting in 0 here since the dark current will be
    # zeroed out by EXPTIME = 0 anyway

    # CAMTEMP header card name comes from DESI-4000v3, "CI Camera" tab
    # night of March 30 biases did not have CAMTEMP available
    #h['CAMTEMP'] = 0.0
    return h

def write_master_bias(outname=None):

    par = common.ci_misc_params()

    if outname is None:
        outname = par['master_bias_filename']
        outname = os.path.join(os.environ[par['etc_env_var']], outname)

    assert(not os.path.exists(outname))

    ci_extnames = common.valid_image_extname_list()

    hdus = []
    for ci_extname in ci_extnames:
        print('Working on master bias for: ' + ci_extname)
        bias_image = master_bias_1camera(ci_extname)
        bias_image = bias_image.astype('float32')
        if len(hdus) == 0:
            hdu = fits.PrimaryHDU(bias_image)
        else:
            hdu = fits.ImageHDU(bias_image)
        hdu.header = master_bias_header_cards(hdu, ci_extname)

        hdus.append(hdu)

    hdul = fits.HDUList(hdus)

    hdul.writeto(outname)
