import numpy as np
import os

# could add utilities re: rotation relative to CS5

def ci_misc_params():
    """
    repository for various constant values that I don't want to 
    end up hardcoding in various places throughout the code base
    """

    width_pix_native = 2048
    height_pix_native = 1032

    par = {'etc_env_var': 'GFA_REDUCE_ETC',
           'gaia_env_var': 'GAIA_CAT_DIR',
           'width_pix_native': width_pix_native,
           'height_pix_native': height_pix_native,
           'n_cameras': 6,
           'nominal_zeropoint': 26.56,
           'fz_dummy_extname': 'CI',
           'master_bias_filename': 'GFA_master_bias.fits', 
           'master_flat_filename': 'GFA_master_flat.fits',
           'master_dark_filename' : 'GFA_master_dark.fits',
           'static_mask_filename': 'GFA_static_badpixels.fits',
           'nominal_sag_cd': (5.0/3.0)*3.55978e-5, 
           'nominal_mer_cd': (5.0/3.0)*3.26627e-5,
           'nominal_cen_cd': (5.0/3.0)*3.70303e-5,
           'full_well_electrons' : 100000.0,
           'nominal_fwhm_asec' : 1.25,
           'headers_dummy_filename' : 'dummy_with_headers_gaia.zenith.fits.gz',
           'reduced_image_flavors' : ['REDUCED', 'INVVAR', 'BITMASK']}

    return par

def ci_camera_gain_dict():
    # HACK !!
    gains = {'GUIDE0': 3.9,
             'GUIDE2': 3.9,
             'GUIDE3': 3.9, 
             'GUIDE5': 3.9,
             'GUIDE7': 3.9,
             'GUIDE8': 3.9}
    return gains

def ci_camera_gain(extname):
    assert(is_valid_extname(extname))

    gains = ci_camera_gain_dict()

    return gains[extname]

def ci_camera_readnoise_dict():
    # units are electrons per pixel
    # HACK !!
    readnoise_electrons = {'GUIDE0': 20.0,
                           'GUIDE2': 20.0,
                           'GUIDE3': 20.0,
                           'GUIDE5': 20.0,
                           'GUIDE7': 20.0,
                           'GUIDE8': 20.0}

    return readnoise_electrons

def ci_camera_readnoise(extname):
    assert(is_valid_extname(extname))

    r = ci_camera_readnoise_dict()
    return r[extname]

def valid_image_extnum_list(fz=True):
    """ 
    When fzipped, valid image extension numbers are 1-5 inclusive.
    When fzipped, there is also a dummy zeroth extension.
    When not fzipped, the valid image extensions are 0-4 inclusive.
    When not fzipped, there is no dummy zeroth extension.

    astropy.io.fits has weird (buggy?) behavior for fzipped files:
        - fits.getheader with ext=0 returns the dummy ext=0 header.
        - fits.getdata with ext=0 and header=True returns the ext=1 header and
          ext=1 image data.

    I want image extension numbers that are valid for both getheader and
    getdata, so for the fzipped case, I will assert that the valid image
    extension numbers are 1-5 inclusive and that ext=0 is invalid.

    hdul.info()
    Filename: FORDK_NOBIN-00000001.fits.fz
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  CI            1 PrimaryHDU       15   ()
      1  CIN           1 CompImageHDU     43   (3072, 2048)   int16
      2  CIW           1 CompImageHDU     43   (3072, 2048)   int16
      3  CIC           1 CompImageHDU     43   (2048, 3072)   int16
      4  CIE           1 CompImageHDU     43   (3072, 2048)   int16
      5  CIS           1 CompImageHDU     43   (3072, 2048)   int16

    hdul.info()
    Filename: FORDK_NOBIN-00000001.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0 CIN            1 PrimaryHDU      47   (3072, 2048)    int16
      1 CIW            1 ImageHDU        42   (3072, 2048)    int16
      2 CIC            1 ImageHDU        42   (2048, 3072)    int16
      3 CIE            1 ImageHDU        42   (3072, 2048)    int16
      4 CIS            1 ImageHDU        42   (3072, 2048)    int16

    Note that the swapped dimensions of CIC above are a bug in the sample
    data, and that the number of header cards will presumably be changing
    relative to the preliminary sample data.
    """

    return (np.arange(1, 6) if fz else np.arange(0, 5))

def valid_extname_list(fz=True):
    """ returned list is in order of CI# in DESI-3347 """

    extnames = ['GUIDE0', 'GUIDE2', 'GUIDE3', 'GUIDE5', 'GUIDE7', 'GUIDE8']

    par = ci_misc_params()
    if fz:
        extnames.append(par['fz_dummy_extname'])

    return extnames

def valid_image_extname_list():
    return valid_extname_list(fz=False)

def valid_ci_number_list():
    """ from DESI-3347 page 2 CI# labels """

    return [0, 2, 3, 5, 7, 8]

def ci_extname_to_extnum_dict(fz=True):
    """
    return a dictionary with the correspondence between CI 
    extension name and CI extension number
    for now I'm including the dummy zeroth extension in the fzipped case    
    """

    par = ci_misc_params()

    if fz:
        d = {par['fz_dummy_extname']: 0, 'CIN': 1, 'CIW': 2, 'CIC': 3, 
                 'CIE': 4, 'CIS': 5}
    else:
        d = {'CIN': 0, 'CIW': 1, 'CIC': 2, 'CIE': 3, 'CIS': 4}

    return d

def ci_extnum_to_extname_dict(fz=True):
    """ 
    make this a wrapper of ci_extname_to_extnum_dict to
    to avoid having the same information typed out explicitly twice, 
    which would be asking for additional bugs/inconsistencies
    """
    d_reverse = ci_extname_to_extnum_dict(fz=fz)

    d = dict((v,k) for k,v in d_reverse.items())

    return d

def ci_extname_to_extnum(extname, fz=True):
    """ convert CI image extension name to extension number"""

    assert(is_valid_extname(extname))

    d = ci_extname_to_extnum_dict(fz=fz)

    return d[extname]

def ci_extname_to_ci_number_dict():
    return dict(zip(valid_image_extname_list(), valid_ci_number_list()))

def ci_number_to_ci_extname_dict():
    return dict(zip(valid_ci_number_list(), valid_image_extname_list()))

def ci_extname_to_ci_number(extname):
    """ convert CI image extension name to integer CI# from DESI-3347"""

    # trim whitespace in case the input somehow has any
    extname = extname.replace(' ', '')

    assert(is_valid_extname(extname))

    d = ci_extname_to_ci_number_dict()

    return d[extname]

def ci_extnum_to_extname(extnum, fz=True):

    assert(is_valid_extnum(extnum, fz=fz))

    d = ci_extnum_to_extname_dict(fz=fz)

    return d[extnum]

def ci_number_to_ci_extname(num):
    """
    CI number refers to CI# as labeled on DESI-3347 page 4
    """
    assert(is_valid_ci_number(num))

    d = ci_number_to_ci_extname_dict()

    return d[num]

def ci_number_to_extnum(num, fz=True):
    # convert ci number to ci extname
    # then use extname to get extnum
    extname = ci_number_to_ci_extname(num)

    d = ci_extname_to_extnum_dict(fz=fz)
    return d[extname]

def is_valid_extnum(extnum, fz=True):
    return (extnum in valid_image_extnum_list(fz=fz))

def is_valid_image_extname(extname):
    return (extname in valid_image_extname_list())

def is_valid_extname(extname, fz=True):
    """
    this is different from the above because when fzipped there's
    a dummy zeroth extension with extname 'CI'
    """

    return (extname in valid_extname_list(fz=fz))

def is_valid_ci_number(num):
    return (num in valid_ci_number_list())

def valid_flavor_list():
    """
    this includes at least 'BIAS', 'LIGHT' based on forDK.tar.gz samples
    capitalization inconsistent in forDK.tar.gz samples
    need to keep an eye out for additional valid flavors to add
    """

    # not sure how to deal with reduced image flavors that I've invented:
    #     REDUCED, INVVAR, BITMASK
    valid_flavors = ['BIAS', 'LIGHT', 'FLAT', 'MASK']

    return valid_flavors

def is_valid_flavor(flavor):
    """ need to find out what the valid flavors (observation types) will be"""
    # be case insensitive for now

    return (flavor.upper() in valid_flavor_list())

def mask_bit_dict():
    # keys are strings with shorthand name for each bit, values
    # are the corresponding powers of 2

    d = {'FLATBAD' : 0, 
         'FLATQ'   : 1, 
         'SATUR'   : 2,
         'NAN'     : 3}

    return d

def mask_bit_from_bitname(bitname):

    d = mask_bit_dict()

    assert(bitname in d.keys())

    return d[bitname]

def mask_bit_description_dict():

    d = {'FLATBAD' : 'bad pixel based on master flat', 
         'FLATQ'   : 'questionable pixel based on master flat', 
         'SATUR'   : 'saturated pixel',
         'NAN'     : 'non-finite pixel value'}

    return d

def mask_bit_description(bitname):
    # bitname is the shorthand name for 

    d = mask_bit_description_dict()

    assert(bitname in d.keys())

    return d[bitname]

def reduced_image_filename_label(flavor):
    par = ci_misc_params()

    assert(flavor in par['reduced_image_flavors'])

    label = '_' + flavor.lower()

    return label

def reduced_flavor_to_bunit(flavor):
    par = ci_misc_params()

    assert(flavor in par['reduced_image_flavors'])

    d = {'REDUCED' : 'ADU', 'INVVAR' : '1/ADU^2', 
         'BITMASK' : 'dimensionless'}

    return d[flavor]

def expid_from_filename(fname):
    # assumes ci-????????.fits.fz pattern
    # this is meant to be run on an input fname that is a single string
    # as opposed to an array/list of strings

    assert(isinstance(fname, str))

    f = os.path.basename(fname)

    expid = int(f[3:11])

    return expid

def get_median_bias_adu(extname):
    # using CI_master_bias.fits with md5sum:
    #     5c48bef5923647ac9646d9462af88032
    # see py/scripts/median_bias_per_camera.py
    # result in in ADU
    bias_med_dict = {'CIE' : 991.0, 
                     'CIN' : 1013.0, 
                     'CIC' : 1020.0, 
                     'CIS' : 984.0, 
                     'CIW' : 990.0}

    assert(extname in bias_med_dict.keys())
    return bias_med_dict[extname]

def get_median_dark_current(extname):
    # return values are in ADU/s at a nominal temperature of 7.5 Celsius
    # the nominal temperature is just a guess since the dark frames
    # that my master dark are based on don't have CCDTEMP available

    # using CI_master_bias.fits with md5sum:
    #     dab062d6f06cc09429efc0f600f44baa
    # see py/scripts/median_dark_current_per_camera.py

    dark_med_dict = {'CIE' : 0.128333, 
                     'CIN' : 0.161667,
                     'CIC' : 0.143333,
                     'CIS' : 0.153333,
                     'CIW' : 0.143333}

    assert(extname in dark_med_dict.keys())
    return dark_med_dict[extname]
