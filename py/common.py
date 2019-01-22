import numpy as np

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

def valid_image_extname_list():
    """ any reason to set order a specific way here? """

    # order of CI# in DESI-3347
    return ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']

def valid_ci_number_list():
    """ from DESI-3347 page 2 CI# labels """

    return np.arange(1, 6)

def ci_extname_to_extnum(extname):
    

def ci_extname_to_ci_number(extname):


def ci_extnum_to_extname(extnum):


def ci_number_to_ci_extname():


def ci_number_to_extnum():


def ci_number_to_extname():


def ci_boundary_coords(pix_center=False):


def ci_corner_coords(pix_center=False):


def check_valid_extnum():


def check_valid_extname():


def check_valid_ci_number(num):
    return (num in valid_ci_number_list())

def valid_flavor_list():
    """
    this includes at least LIGHT, BIAS 
    capitalization inconsistent in sample data
    """

def check_valid_flavor(flavor):
    """ need to find out what the valid flavors (observation types) will be"""
