from astropy.stats import mad_std
from scipy.ndimage.filters import gaussian_filter
import numpy as np
from scipy.ndimage.measurements import label, find_objects
from astropy.table import Table

def slices_to_table(slices, detsn, extname):
    xcen = []
    ycen = []
    xmin = []
    xmax = []
    ymin = []
    ymax = []
    detsn_value = []

    for slc in slices:
        ycen.append( (float(slc[0].start) + float(slc[0].stop))/2.0)
        xcen.append( (float(slc[1].start) + float(slc[1].stop))/2.0)
        xmin.append(slc[1].start)
        xmax.append(slc[1].stop)
        ymin.append(slc[0].start)
        ymax.append(slc[0].stop)
        detsn_value.append(np.max(detsn[slc]))

    xcen = np.array(xcen)
    ycen = np.array(ycen)
    xmin = np.array(xmin)
    xmax = np.array(xmax)
    ymin = np.array(ymin)
    ymax = np.array(ymax)
    extname = [extname]*len(xcen)
    detsn_value = np.array(detsn_value)

    tab = Table([extname, xcen, ycen, xmin, xmax, ymin, ymax, detsn_value], names=('extname', 'xcen', 'ycen', 'xmin', 'xmax', 'ymin', 'ymax', 'detsn_value'))

    return tab

def detection_map(im, fwhm):
    psf_sigma = fwhm/2.355

    smth = gaussian_filter(im, psf_sigma)

    sig_smth = mad_std(smth)

    detsn = (smth-np.median(smth))/sig_smth

    return detsn

def detect_sources(detsn, thresh):
    peaks = (detsn > thresh)
    blobs, _ = label(peaks)
    slices = find_objects(blobs)
    return slices

def get_nominal_fwhm_pix(extname):
    nominal_fwhm_pix = (9.3766944 if (extname == 'CIC') else 10.192287)

    return nominal_fwhm_pix

def get_source_list(image, bitmask, extname, thresh=5):
    filler_value = np.median(image)

    # should do something like djs_maskinterp instead
    image[bitmask != 0] = filler_value
 
    nominal_fwhm_pix = get_nominal_fwhm_pix(extname)

    detsn = detection_map(image, nominal_fwhm_pix)

    slices = detect_sources(detsn, thresh)

    tab = slices_to_table(slices, detsn, extname)

    tab.rename_column('xcen', 'xcentroid')
    tab.rename_column('ycen', 'ycentroid')

    return tab
