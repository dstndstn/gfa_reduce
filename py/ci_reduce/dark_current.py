import ci_reduce.common as common
import numpy as np
import matplotlib.pyplot as plt
import ci_reduce.imred.load_calibs as load_calibs
import os
import astropy.io.fits as fits
from scipy.optimize import minimize

# return dark current rate in e-/pix/sec as a function of temperature
# based on DESI-3358 slide 9, or else my own dark current measurements once 
# I get access to the relevant engineering data 

# once I get access to the engineering dark frames, would also be good
# to construct *images* of the dark current rate

# assume all cameras have same dark current rate, since I don't currently
# have sufficient data to study any camera-to-camera dark current differences

class DarkCurrentInfo:
    def __init__(self, fname_master_dark, extname, header=None):

        self.fname_master_dark = fname_master_dark
        self.extname = extname
        
        if header is None:
            self.header = fits.getheader(fname_master_dark, extname=extname)
        else:
            self.header = header

        assert(self.header['EXTNAME'].replace(' ', '') == self.extname)

def _objective_function(p, im, dark):

    # im and dark should already be made 1-dimensional before being input !!
    
    # p should have one element
    # p should generally be a factor relatively close to 1 like 1.1 or 0.9
    _sub = im - dark*p[0]

    n_all = _sub.size

    # important that this not be a divisor of the # of rows or cols
    fac = 10

    n_subsample = int(np.floor(float(n_all)/float(fac)))

    ind_subsample = fac*np.arange(n_subsample, dtype=int)

    sub = _sub[ind_subsample]

    sind = np.argsort(sub)

    # think more later on about whether floor is really what i want here
    val_l = sub[sind[int(np.floor(0.16*sub.size))]]
    val_u = sub[sind[int(np.floor(0.84*sub.size))]]

    spread = val_u - val_l

    return spread
        
def fit_dark_scaling_1amp(im, dark_guess_scaled, amp, extname):
    is_focus = extname.find('FOCUS') != -1

    if amp == 'E':
        xmax = 900 if is_focus else 1024
        cutout = im[0:516, 0:xmax]
        dark_cutout = dark_guess_scaled[0:516, 0:xmax]
    elif amp == 'F':
        xmin = 1170 if is_focus else 1024
        cutout = im[0:516, xmin:2048]
        dark_cutout = dark_guess_scaled[0:516, xmin:2048]
    elif amp == 'G':
        xmin = 1170 if is_focus else 1024
        cutout = im[516:1032, xmin:2048]
        dark_cutout = dark_guess_scaled[516:1032, xmin:2048]
    elif amp == 'H':
        xmax = 900 if is_focus else 1024
        cutout = im[516:1032, 0:xmax]
        dark_cutout = dark_guess_scaled[516:1032, 0:xmax]

    _cutout = np.ravel(cutout)
    _dark_cutout = np.ravel(dark_cutout)

    initial_simplex = np.zeros((2, 1))
    initial_simplex[0] = 0.99
    initial_simplex[1] = 1.01
    
    res = minimize(_objective_function, [1.0], args=(_cutout, _dark_cutout), method='Nelder-Mead', options={'maxfev': 200, 'disp': False, 'initial_simplex': initial_simplex, 'adaptive': False, 'fatol': 1.0e-5})

    return res

def fit_dark_scaling(im, dark_guess_scaled, extname):
    # im needs to be bias-subtracted !!!!!!!!!

    # im should have overscan stripped out
    # dark should have overscan stripped out
  
    # dark should be scaled according to exptime and best guess of
    # T dependence so that the additional scaling fit out by this routine
    # is a perturbation

    sh_im = im.shape
    assert((sh_im[0] == 1032) and (sh_im[1] == 2048))

    sh_dark = dark_guess_scaled.shape
    assert((sh_dark[0] == 1032) and (sh_dark[1] == 2048))

    amps = common.valid_amps_list()

    rescale_factors = np.zeros(4)
    ncalls = np.zeros(4, dtype=int)

    for i, amp in enumerate(amps):
        res = fit_dark_scaling_1amp(im, dark_guess_scaled, amp, extname)
        x = res.x
        print(x[0], ' !!!!!!!!!!!!!!!!!!! ', extname, '  ', res.nfev)
        rescale_factors[i] = x[0]
        ncalls[i] = res.nfev
        #rescale_factors[i] = rescale_factor
        #ncalls[i] = _ncalls

    rescale_factor = np.median(rescale_factors)
    print(rescale_factor, ' !!! ')
    return rescale_factor
    
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

def read_dark_image(ci_extname, exptime, t_celsius):
    assert(common.is_valid_extname(ci_extname))

    par = common.ci_misc_params()

    # try getting a master dark with an exactly matching integration time
    dark_fname = choose_master_dark(exptime, ci_extname, t_celsius)

    # if no master dark has an exactly matching integration time
    # then just go back to some 'standard' 5 s master dark
    # REVISIT THIS LATER TO DO BETTER
    if dark_fname is None:
        print('could not find a master dark with ORIGTIME matching EXPTIME')
        dark_fname = os.path.join(os.environ[par['etc_env_var']], \
                                  par['master_dark_filename'])

    print('Attempting to read master dark : ' + dark_fname + 
          ', extension name : ' + ci_extname)

    assert(os.path.exists(dark_fname))

    dark, hdark = fits.getdata(dark_fname, extname=ci_extname, header=True)

    dark = load_calibs.remove_overscan(dark)
    return dark, hdark, dark_fname

def total_dark_image_adu(extname, exptime, t_celsius, im):

    # im is the bias-subtracted, overscan/prescan removed version of the
    # image being reduced
    
    # return a predicted image of the total dark current in
    # a GFA image by scaling the master dark image to account 
    # for the exposure time and temperature
    # return value will be in ADU !!!

    dark_image, hdark, dark_fname = read_dark_image(extname, exptime, t_celsius)

    dark_image *= dark_rescaling_factor(hdark['GCCDTEMP'], t_celsius, extname)

    dark_image *= exptime

    rescale_factor = fit_dark_scaling(im, dark_image, extname)
    
    return dark_image*rescale_factor

def choose_master_dark(exptime, extname, gccdtemp):

    par = common.ci_misc_params()
    
    # eventually could cache the index of master darks...
    fname_index = os.path.join(os.environ[par['etc_env_var']], par['dark_index_filename'])
    
    assert(os.path.exists(fname_index))
    str = fits.getdata(fname_index)

    # cases of potential bad readout should already be removed, but just in case
    str = str[str['READWARN'] == 0]

    str = str[(str['ORIGTIME'] == exptime) & (str['EXTNAME'] == extname)]

    # this is the case where EXPTIME does not have any available
    # master darks in the library of master darks
    if len(str) == 0:
        return None

    indmin = np.argmin(np.abs(str['GCCDTEMP'] - gccdtemp))

    return str[indmin]['FNAME_FULL'].replace(' ', '')
