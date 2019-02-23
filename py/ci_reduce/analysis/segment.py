import ci_reduce.common as common
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources

def segmentation_map(image):
    # in this context image means a 2D numpy array rather than a CI_image
    # object

    threshold = detect_threshold(image, snr=2.)
    
    # fwhm = 9.0 pixels
    sigma = 9.0*gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=9, y_size=9)
    kernel.normalize()

    segm = detect_sources(image, threshold, npixels=5, filter_kernel=kernel)

    return segm
