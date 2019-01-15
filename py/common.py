import numpy as np

def valid_image_extnum_list():
    """ 
    valid image extension numbers are 1-5 inclusive
    I guess there is also a dummy zeroth extension...
    """

    return np.arange(1, 6)

def valid_image_extname_list():
    """ any reason to set order a specific way here? """

    # order of CI# in DESI-3347
    return ['CIE', 'CIN', 'CIC', 'CIS', 'CIE']

def valid_ci_number_list():
    return np.arange(1, 6)

def ci_extname_to_extnum(extname):
    

def ci_extname_to_ci_number(extname):


def ci_extnum_to_extname(extnum):


def ci_number_to_ci_extname():


def ci_number_to_extnum():


def ci_number_to_extname():


def ci_boundary_coords():


def ci_corner_coords():


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
