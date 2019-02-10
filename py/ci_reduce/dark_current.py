import ci_reduce.common as common
import numpy as np
import matplotlib.pyplot as plt

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

    gain = common.ci_camera_gain(ci_extname)
    dark_e_per_pix = total_dark_current_electrons(acttime, t_celsius)
    dark_adu_per_pix = dark_e_per_pix/gain

    return dark_adu_per_pix
