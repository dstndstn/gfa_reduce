import ci_reduce.analysis.util as util
import numpy as np
from astropy.stats import mad_std
from astropy.table import Table, hstack

# add quadrant kw args to allow computation of each statistic not only
# on each full CCD but also on each quadrant separately

# for now compute these stats on the *raw* image, could consider
# switching to the reduced image, or doing both in the future

def slice_indices_for_quadrant(quadrant):

    xmin = int(util.ci_pixel_xmin(pix_center=True, quadrant=quadrant))
    xmax = int(util.ci_pixel_xmax(pix_center=True, quadrant=quadrant)) + 1
    ymin = int(util.ci_pixel_ymin(pix_center=True, quadrant=quadrant))
    ymax = int(util.ci_pixel_ymax(pix_center=True, quadrant=quadrant)) + 1

    return xmin, xmax, ymin, ymax

def compute_image_median(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return np.median(image[ymin:ymax, xmin:xmax])

def compute_image_mean(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return np.mean(image[ymin:ymax, xmin:xmax])

def n_non_finite(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return np.sum(np.logical_not(np.isfinite(image[ymin:ymax, xmin:xmax])))

def compute_image_max(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return np.max(image[ymin:ymax, xmin:xmax])

def compute_image_min(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return np.min(image[ymin:ymax, xmin:xmax])

def compute_mad_std(image, quadrant=None):

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return mad_std(image[ymin:ymax, xmin:xmax])

def compute_std(image, quadrant=None):
    # non-robust standard deviation

    xmin, xmax, ymin, ymax = slice_indices_for_quadrant(quadrant)

    return mad_std(image[ymin:ymax, xmin:xmax])

# do a meta-analysis of the per-object FWHM values, 
# presumably restricting to relative high s/n sources
# if no good sources are available, then need to return some dummy
# value 
def overall_image_fwhm(tab):
    print('stub')

# do a meta-analysis of the per-object ellipticity values, 
# presumably restricting to relative high s/n sources
# if no good sources are available, then need to return some dummy
# value 
def overall_image_ellipticity(tab):
    print('stub')

def compute_all_stats(image, skip_quadrants=False, extname=None):
    sectors = [None]
    if not skip_quadrants:
        sectors += [1, 2, 3, 4]

    t = None
    for q in sectors:
        med = [compute_image_median(image, quadrant=q)]
        avg = [compute_image_mean(image, quadrant=q)]
        nbad = [n_non_finite(image, quadrant=q)]
        _max = [compute_image_max(image, quadrant=q)]
        _min = [compute_image_min(image, quadrant=q)]
        sig_robust = [compute_mad_std(image, quadrant=q)]
        sig = [compute_std(image, quadrant=q)]

        colnames = ['median', 'mean', 'n_non_finite', 'max', 'min', 
                    'sig_robust', 'sig']

        if q is not None:
             colnames = [(colname + '_q' + str(q)) for colname in colnames]

        if extname is not None:
            colnames = [(colname + '_' + extname) for colname in colnames]

        row = Table([med, avg, nbad, _max, _min, sig_robust, sig], names=tuple(colnames))

        if t is None:
            t = row
        else:
            # hstack the tables
            t = hstack([t, row])

    return t
