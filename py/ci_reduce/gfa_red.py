#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
from datetime import datetime
import ci_reduce.analysis.util as util
import ci_reduce.common as common
import ci_reduce.analysis.recalib_astrom as wcs
import time

def _proc(fname_in, outdir=None, careful_sky=False, no_cataloging=False,
          no_gaia_xmatch=False, no_ps1_xmatch=False,
          cube_index=None, skip_image_outputs=False,
          realtime=False, no_dark_rescaling=False, 
          dont_write_invvar=False, make_psf_models=False):

    print('Starting GFA reduction pipeline at: ' + str(datetime.utcnow()) + 
          ' UTC')

    t0 = time.time()

    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')
    
    write_outputs = (outdir is not None)

    assert(os.path.exists(fname_in))

    gitrev = io.retrieve_git_rev()

    if write_outputs:
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # fail if ANY of expected outputs already exist
        io.check_image_level_outputs_exist(outdir, fname_in, gzip=True,
                                           cube_index=cube_index)

    exp = io.load_exposure(fname_in, cube_index=cube_index, realtime=realtime)

    # check for simulated data
    if (exp is None) or util.has_wrong_dimensions(exp):
        # point is to not crash, for sake of real time reductions
        print('EXITING: exposure may be a simulation or contain only focus camera images?!')
        return
    
    print('Attempting to compute basic statistics of raw pixel data')
    imstats = io.gather_pixel_stats(exp)

    # create data quality bitmasks
    # exp.create_all_bitmasks() # revisit this later

    # go from "raw" images to "reduced" images
    exp.calibrate_pixels(do_dark_rescaling=(not no_dark_rescaling))

    # calculate sky brightness in mag per sq asec
    exp.estimate_all_sky_mags(careful_sky=careful_sky)
    exp.estimate_all_sky_sigmas(careful_sky=careful_sky)

    par = common.ci_misc_params()

    if not no_cataloging:
        catalogs = exp.all_source_catalogs()

        for extname, cat in catalogs.items():
            if cat is not None:
                util.create_det_ids(cat, extname, fname_in,
                                    cube_index=cube_index)

        # reformat the output catalogs into a single merged astropy Table
        catalog = io.combine_per_camera_catalogs(catalogs)

        # run astrometric recalibration
        print('Attempting astrometric recalibration relative to Gaia DR2')
        astr = wcs.recalib_astrom(catalog, fname_in)
        exp.update_wcs(astr)
        exp.recompute_catalog_radec(catalog)

        if make_psf_models:
            exp.compute_psfs(catalog)
        
        if (not no_ps1_xmatch) and (par['ps1_env_var'] in os.environ):
            # probably should look into dec < -30 handling more at some point
            print('Attempting to perform PS1 cross-matching...')
            ps1 = io.write_ps1_matches(catalog, outdir, fname_in,
                                       cube_index=cube_index)
        else:
            ps1 = None
            
        if (not no_gaia_xmatch) and (par['gaia_env_var'] in os.environ):
            print('Attempting to identify Gaia cross-matches')
            catalog = io.append_gaia_crossmatches(catalog)
        
    # try to write image-level outputs if outdir is specified

    if write_outputs:
        if not skip_image_outputs:
            print('Attempting to write image-level outputs to directory : ' + 
                  outdir)
            # could add command line arg for turning off gzip compression
            io.write_image_level_outputs(exp, outdir, fname_in, gzip=True,
                                         cube_index=cube_index,
                                         dont_write_invvar=dont_write_invvar)

        # make this work correctly in the case that --no_cataloging is set
        io.write_ccds_table(imstats, catalog, exp, outdir, fname_in,
                            cube_index=cube_index, ps1=ps1)
        
        if not no_cataloging:
            io.write_exposure_source_catalog(catalog, outdir, fname_in, exp,
                                             cube_index=cube_index)
            if make_psf_models:
                io.write_psf_models(exp, outdir, fname_in,
                                    cube_index=cube_index)

    print('Successfully finished reducing ' + fname_in)

    dt = time.time() - t0
    print('GFA reduction pipeline took ' + '{:.2f}'.format(dt) + ' seconds')
    print('GFA reduction pipeline completed at: ' + str(datetime.utcnow()) + 
          ' UTC')

if __name__ == "__main__":
    descr = 'run full gfa_reduce pipeline on a GFA exposure'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--careful_sky', default=False, action='store_true',
        help='use image segmentation when deriving sky quantities')

    parser.add_argument('--no_cataloging', default=False, action='store_true', 
                        help='reduce image without cataloging sources')

    parser.add_argument('--no_gaia_xmatch', default=False, action='store_true',
                        help='skip Gaia cross-match')

    parser.add_argument('--no_ps1_xmatch', default=False, action='store_true',
                        help='skip PS1 cross-match')

    parser.add_argument('--cube_index', default=None, type=int,
                        help='guide cube index')

    parser.add_argument('--skip_image_outputs', default=False,
                        action='store_true',
                        help='skip writing of full-frame image outputs')

    parser.add_argument('--realtime', default=False,
                        action='store_true',
                        help='avoid crashing on partially written raw files')

    parser.add_argument('--no_dark_rescaling', default=False,
                        action='store_true',
                        help='skip empirical rescaling of dark current')

    parser.add_argument('--dont_write_invvar', default=False, 
                        action='store_true',
                        help="don't write out invvar maps")

    parser.add_argument('--make_psf_models', default=False,
                        action='store_true',
                        help="generate per-camera PSF models")
    
    args = parser.parse_args()

    fname_in = args.fname_in[0]

    _proc(fname_in, outdir=args.outdir, careful_sky=args.careful_sky,
          no_cataloging=args.no_cataloging, no_gaia_xmatch=args.no_gaia_xmatch,
          no_ps1_xmatch=args.no_ps1_xmatch, cube_index=args.cube_index,
          skip_image_outputs=args.skip_image_outputs, realtime=args.realtime,
          no_dark_rescaling=args.no_dark_rescaling, 
          dont_write_invvar=args.dont_write_invvar,
          make_psf_models=args.make_psf_models)
