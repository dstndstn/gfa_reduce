import ci_reduce.analysis.util as util
from astropy.stats import mad_std
from scipy.ndimage.filters import gaussian_filter
import numpy as np
from scipy.ndimage.measurements import label, find_objects
from astropy.table import Table
from photutils import centroid_com, centroid_2dg
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus

def slices_to_table(slices, detsn, extname):
    nslc = len(slices)
    xcen = np.zeros(nslc)
    ycen = np.zeros(nslc)
    xmin = np.zeros(nslc, dtype='int')
    xmax = np.zeros(nslc, dtype='int')
    ymin = np.zeros(nslc, dtype='int')
    ymax = np.zeros(nslc, dtype='int')
    detmap_peak = np.zeros(nslc)

    for i, slc in enumerate(slices):
        ycen[i] = (float(slc[0].start) + float(slc[0].stop))/2.0
        xcen[i] = (float(slc[1].start) + float(slc[1].stop))/2.0
        xmin[i] = slc[1].start
        xmax[i] = slc[1].stop
        ymin[i] = slc[0].start
        ymax[i] = slc[0].stop
        detmap_peak[i] = np.max(detsn[slc])

    tab = Table([xcen, ycen, xmin, xmax, ymin, ymax, detmap_peak], names=('xcen_init', 'ycen_init', 'detmap_xmin', 'detmap_xmax', 'detmap_ymin', 'detmap_ymax', 'detmap_peak'))

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

def do_aper_phot(data, catalog, extname):
    # catalog should be the catalog with refined centroids 
    # for **one CI camera**

    print('Attempting to do aperture photometry')
    positions = list(zip(catalog['xcentroid'], catalog['ycentroid']))
    apertures = CircularAperture(positions, r=get_nominal_fwhm_pix(extname)/2.0)
    annulus_apertures = CircularAnnulus(positions, r_in=30.0, r_out=40.0)
    annulus_masks = annulus_apertures.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)
    phot = aperture_photometry(data, apertures)
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * apertures.area()
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']

    # still need to actually append some of this aperture photometry
    # info to the catalog in the form of new columns

    catalog['aper_sum_bkgsub'] = phot['aper_sum_bkgsub']
    catalog['aper_bkg'] = phot['aper_bkg']

def get_nominal_fwhm_pix(extname):
    # this is a nominal FWHM for use as an initial guess
    # when creating the detection map, when nothing is yet known
    # about the FWHM of the particular data being analyzed
    # could revisit this and/or clean it up later on
    nominal_fwhm_pix = (9.3766944 if (extname == 'CIC') else 10.192287)

    return nominal_fwhm_pix

def refine_centroids(tab, image, bitmask):
    # input table tab gets augmented with additional columns

    print('Attempting to refine initial centroids')
    # could scale this based on the typical size of the slices
    boxsize = 11

    half = int(np.floor(boxsize/2))

    med = np.median(image)

    nobj = len(tab)
    xcentroid = np.zeros(nobj)
    ycentroid = np.zeros(nobj)

    no_centroid_refinement = np.zeros(nobj, dtype=bool)
    centroid_refinement_fail = np.zeros(nobj, dtype=bool)

    for i in range(nobj):
        ix_guess = int(round(tab[i]['xcen_init']))
        iy_guess = int(round(tab[i]['ycen_init']))
        min_edge_dist = util.min_edge_dist_pix(ix_guess, iy_guess)

        # if the centering box extends outside of image, don't
        # try to do any refinement of the centroid
        if min_edge_dist < half:
            xcentroid[i] = tab[i]['xcen_init']
            ycentroid[i] = tab[i]['ycen_init']
            no_centroid_refinement[i] = True
            continue

        _xcentroid, _ycentroid = centroid_2dg(image[(iy_guess-half):(iy_guess+half+1), (ix_guess-half):(ix_guess+half+1)] - med)
        if (not np.isfinite(_xcentroid)) or (not np.isfinite(_ycentroid)):
            xcentroid[i] = tab[i]['xcen_init']
            ycentroid[i] = tab[i]['ycen_init']
            centroid_refinement_fail[i] = True
        else:
            xcentroid[i] = _xcentroid + ix_guess - half
            ycentroid[i] = _ycentroid + iy_guess - half

    tab['no_centroid_refinement'] = no_centroid_refinement.astype(int)
    tab['centroid_refinement_fail'] = centroid_refinement_fail.astype(int)
    tab['xcentroid'] = xcentroid
    tab['ycentroid'] = ycentroid

def add_metadata_columns(tab, bitmask):
    # input table tab gets modified

    min_edge_dist = [util.min_edge_dist_pix(c[0], c[1]) for c in zip(tab['xcentroid'], tab['ycentroid'])]
    tab['min_edge_dist_pix'] = min_edge_dist

    # there's probably a scipy function for this but w/e

    xmin = util.ci_pixel_xmin(pix_center=True)
    xmax = util.ci_pixel_xmax(pix_center=True)
    ymin = util.ci_pixel_ymin(pix_center=True)
    ymax = util.ci_pixel_ymax(pix_center=True)

    ixs = [int(min(max(np.round(t['xcentroid']), xmin), xmax)) for t in tab]
    iys = [int(min(max(np.round(t['ycentroid']), ymin), ymax)) for t in tab]

    tab['dq_flags'] = bitmask[iys, ixs]

def get_source_list(image, bitmask, extname, thresh=5):
    filler_value = np.median(image)

    # should do something like djs_maskinterp instead
    image[bitmask != 0] = filler_value
 
    nominal_fwhm_pix = get_nominal_fwhm_pix(extname)

    detsn = detection_map(image, nominal_fwhm_pix)

    slices = detect_sources(detsn, thresh)

    tab = slices_to_table(slices, detsn, extname)

    refine_centroids(tab, image, bitmask)

    do_aper_phot(image, tab, extname)

    add_metadata_columns(tab, bitmask)

    return tab
