import ci_reduce.analysis.util as util
import numpy as np

# do a meta-analysis of the per-object FWHM values, 
# presumably restricting to relative high s/n sources
# if no good sources are available, then need to return some dummy
# value 
def overall_image_fwhm(tab, bad_amps=None, snr_thresh=20):
    # remove bad/dummy FWHM values
    # remove sources that are near image boundaries
    # could remove sources that aren't isolated, but might be fine
    # to leave this as an optimization for the future

    assert(len(np.unique(tab['camera'])) == 1)
        
    good = ((tab['sig_major_pix'] > 1) & np.isfinite(tab['sig_major_pix']) & 
            (tab['dq_flags'] == 0) & (tab['min_edge_dist_pix'] > 30) &
            (tab['detmap_peak'] >= snr_thresh))

        # if bad amps specified, it should be a list 
    # containing the amps thought to be in a state of bad readout
    if (bad_amps is not None) and (len(bad_amps) > 0):
        amp_ok = np.array([(t['amp'] not in bad_amps) for t in tab])
        good = (good & amp_ok)

    # check for case of not enough sources
    if np.sum(good) < 2:
        return -1, -1, -1, -1, 0

    # just take geometric mean of semi-major / semi-minor sigmas ??
    fwhm_major_pix = np.median(tab[good]['sig_major_pix'])*2.355
    fwhm_minor_pix = np.median(tab[good]['sig_minor_pix'])*2.355

    fwhm_pix = np.sqrt(fwhm_major_pix*fwhm_minor_pix)

    # this is a temporary HACK -- should really take into account 
    # the PA and directionality of the platescale when doing this

    asec_per_pix = util.nominal_pixel_sidelen_arith(tab[0]['camera'])

    fwhm_asec= asec_per_pix*fwhm_pix

    return fwhm_major_pix, fwhm_minor_pix, fwhm_pix, fwhm_asec, int(np.sum(good))
