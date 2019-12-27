import ci_reduce.common as common
import numpy as np
import matplotlib.pyplot as plt
import ci_reduce.imred.load_calibs as load_calibs

# return dark current rate in e-/pix/sec as a function of temperature
# based on DESI-3358 slide 9, or else my own dark current measurements once 
# I get access to the relevant engineering data 

# once I get access to the engineering dark frames, would also be good
# to construct *images* of the dark current rate

# assume all cameras have same dark current rate, since I don't currently
# have sufficient data to study any camera-to-camera dark current differences

def desi_3358_measurements():
    # not aware of uncertainties being available for either of these quantities
    t_celsius = np.array([0.0, 5.0, 10.0, 20.0])
    e_per_pix_per_sec = np.array([0.08, 0.16, 0.38, 0.62])

    return t_celsius, e_per_pix_per_sec

def fit_dark_doubling_rate():
    t_celsius, e_per_pix_per_sec = desi_3358_measurements()
    y = np.log(e_per_pix_per_sec)
    coeff = np.polyfit(t_celsius, y, 1)

    print('Best-fit dark current doubling rate (deg C) = ', \
          np.log(2.0)/coeff[0])
    print('Best-fit 0 Celsius dark current (e-/pix/sec) = ', \
          np.e**coeff[1])

def dark_current_rate(t_celsius):
    """
    t_celsius - temperature in deg celsius
    I = I(0 Celsius)*2^(T_celsius/dT), dT = doubling rate in deg C
    I(0 Celsius) and hence output will have units of e-/sec/pix
    parameters I(0 Celsius) and dT determined from fit_dark_doubling_rate()

    should work for both scalar and array t_celsius inputs
    """

    # not sure if there's any value in storing these parameters somewhere
    # more useful, such as allowing them to be retrieved via common.py

    I0 = 0.0957
    dT = 6.774

    I = I0*np.power(2, t_celsius/dT)

    return I

def total_dark_current_electrons(acttime, t_celsius):
    """
    accttime - actual exposure time in seconds
    t_celsius - temperature in deg celsius

    output is in e-/pix
    """
    dark_e_per_pix_per_sec = dark_current_rate(t_celsius)
    dark_e_per_pix = dark_e_per_pix_per_sec*acttime

    return dark_e_per_pix

def total_dark_current_adu(ci_extname, acttime, t_celsius):
    """
    accttime - actual exposure time in seconds
    t_celsius - temperature in deg celsius

    output is in ADU/pix
    """

    assert(common.is_valid_image_extname(ci_extname))

    assert(acttime >= 0)

    gain = common.ci_camera_gain(ci_extname)
    dark_e_per_pix = total_dark_current_electrons(acttime, t_celsius)
    dark_adu_per_pix = dark_e_per_pix/gain

    return dark_adu_per_pix


def get_linear_coeff(extname):

  assert(extname in common.valid_extname_list())

  # https://github.com/desihub/desicmx/blob/master/analysis/gfa/GFA-Dark-Calibration.ipynb
  # 91cb09761ccaa0e57b8a4bf0673bacb6  GFA-Dark-Calibration.ipynb
  d = {'GUIDE0' : 0.223,
       'FOCUS1' : 0.224,
       'GUIDE2' : 0.237,
       'GUIDE3' : 0.222,
       'FOCUS4' : 0.233,
       'GUIDE5' : 0.237,
       'FOCUS6' : 0.209,
       'GUIDE7' : 0.228,
       'GUIDE8' : 0.195,
       'FOCUS9' : 0.225}

  return d[extname]


def dark_rescaling_factor(t_master, t_image, extname):
    f = get_linear_coeff(extname)

    # linear rescaling factors returned by get_linear_coeff are
    # referenced to 11.0 C (would be good to extract this special number
    # to somewhere else)
    fac = (1 + f*(t_image - 11.0))/(1 + f*(t_master - 11.0))

    assert(fac > 0)

    print(fac, '$$$$$$$$$$$$$$$$$$$$$$$$$$')
    return fac

def total_dark_image_adu(extname, exptime, t_celsius):
    # return a predicted image of the total dark current in
    # a CI image by scaling the master dark image to account 
    # for the exposure time and temperature
    # return value will be in ADU !!!

    dark_image, hdark = load_calibs.read_dark_image(extname)

    dark_image *= dark_rescaling_factor(hdark['GCCDTEMP'], t_celsius, extname)

    dark_image *= exptime

    return dark_image
