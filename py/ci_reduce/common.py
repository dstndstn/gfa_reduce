import numpy as np

# could add utilities re: rotation relative to CS5
# add retrieval of nominal gain and readnoise values for each camerax

def ci_misc_params():
    """
    repository for various constant values that I don't want to 
    end up hardcoding in various places throughout the code base
    """

    width_pix = 3072
    height_pix = 2048

    par = {'width_pix': width_pix,
           'height_pix': height_pix,
           'center_pix_coord_x': 0.5*width_pix + 0.5,
           'center_pix_coord_y': 0.5*height_pix + 0.5,
           'n_cameras': 5}

    return par

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

def valid_image_extname_list(fz=True):
    """ returned list is in order of CI# in DESI-3347 """

    extnames = ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']

    if fz:
        extnames.append('CI')

    return extnames

def valid_ci_number_list():
    """ from DESI-3347 page 2 CI# labels """

    return np.arange(1, 6)

def ci_extname_to_extnum_dict(fz=True):
    """
    return a dictionary with the correspondence between CI 
    extension name and CI extension number
    for now I'm including the dummy zeroth extension in the fzipped case    
    """

    if fz:
        d = {'CI': 0, 'CIN': 1, 'CIW': 2, 'CIC': 3, 'CIE': 4, 'CIS': 5}
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

    return (extname in valid_image_extname_list(fz=fz))

def is_valid_ci_number(num):
    return (num in valid_ci_number_list())

def valid_flavor_list():
    """
    this includes at least 'BIAS', 'LIGHT' based on forDK.tar.gz samples
    capitalization inconsistent in forDK.tar.gz samples
    need to keep an eye out for additional valid flavors to add
    """

    valid_flavors = ['BIAS', 'LIGHT']

    return valid_flavors

def is_valid_flavor(flavor):
    """ need to find out what the valid flavors (observation types) will be"""
    # be case insensitive for now

    return (flavor.upper() in valid_flavor_list())

def ci_pixel_xmax(pix_center=False):
    """
    "x" here is in CI pixel coordinates
    could imagine adding a "binfac" keyword here for use in processing
    steps where I've performed an integer downbinning
    """
    par = ci_misc_params()

    # right edge of rightmost pixel
    xmax = par['width_pix'] - 0.5

    if pix_center:
        xmax -= 0.5 # center of rightmost pixel

    return xmax

def ci_pixel_ymax(pix_center=False):
    """
    "y" here is in CI pixel coordinates
    """
    par = ci_misc_params()

    # right edge of rightmost pixel
    ymax = par['height_pix'] - 0.5

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

def ci_boundary_pixel_coords(pix_center=True):
    par = ci_misc_params()

    x_top = np.arange(ci_pixel_xmin(pix_center=pix_center), 
                      ci_pixel_xmax(pix_center=pix_center) + 1)
    x_left = np.zeros(par['height_pix'] + 1*(not pix_center)) + \
                      ci_pixel_xmin(pix_center=pix_center)
    y_left = np.arange(ci_pixel_ymin(pix_center=pix_center),
                      ci_pixel_ymax(pix_center=pix_center) + 1)
    y_bottom = np.zeros(par['width_pix'] + 1*(not pix_center)) + \
                        ci_pixel_ymin(pix_center=pix_center)
    y_top = y_bottom + par['height_pix'] - 1 + 1*(not pix_center)
    x_right = x_left + par['width_pix'] - 1 + 1*(not pix_center)
    y_right = np.flip(y_left, axis=0)
    x_bottom = np.flip(x_top, axis=0)

    print(np.min(x_left), np.max(x_left), len(x_left))
    print(np.min(x_top), np.max(x_top), len(x_top))
    print(np.min(y_left), np.max(y_left), len(y_left))
    print(np.min(y_bottom), np.max(y_bottom), len(y_bottom))

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
