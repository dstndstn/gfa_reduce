#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
from datetime import datetime
import ci_reduce.analysis.util as util

if __name__ == "__main__":
    descr = 'run full ci_reduce pipeline on a CI exposure'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default='', type=str,
                        help='directory to write outputs in')

    parser.add_argument('--careful_sky', default=False, action='store_true',
        help='use image segmentation when deriving sky quantities')

    parser.add_argument('--no_cataloging', default=False, action='store_true', 
        help='reduce image without cataloging sources')

    parser.add_argument('--no_gaia_xmatch', default=False, action='store_true',
        help='skip Gaia cross-match')

    args = parser.parse_args()

    print('Starting CI reduction pipeline at: ' + str(datetime.utcnow()) + 
          ' UTC')

    fname_in = args.fname_in[0]

    write_outputs = (len(args.outdir) > 0)

    assert(os.path.exists(fname_in))

    gitrev = io.retrieve_git_rev()

    if write_outputs:
        outdir = args.outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # fail if ANY of expected outputs already exist
        io.check_image_level_outputs_exist(outdir, fname_in, gzip=True)

    exp = io.load_exposure(fname_in)

    # create data quality bitmasks
    exp.create_all_bitmasks()

    # go from "raw" images to "reduced" images
    exp.calibrate_pixels()

    # calculate sky brightness in mag per sq asec
    exp.estimate_all_sky_mags(careful_sky=args.careful_sky)
    exp.estimate_all_sky_sigmas(careful_sky=args.careful_sky)

    if not args.no_cataloging:
        catalogs = exp.all_source_catalogs()

        for extname, cat in catalogs.items():
            util.create_det_ids(cat, extname, fname_in)

        # reformat the output catalogs into a single merged astropy Table
        catalog = io.combine_per_camera_catalogs(catalogs)
        print('Attempting to identify Gaia cross-matches')
        if not args.no_gaia_xmatch:
            catalog = io.append_gaia_crossmatches(catalog)

    # try to write image-level outputs if outdir is specified
    if write_outputs:
        print('Attempting to write image-level outputs to directory : ' + 
              outdir)

    if write_outputs:
        # could add command line arg for turning off gzip compression
        io.write_image_level_outputs(exp, outdir, fname_in, gzip=True)
        if not args.no_cataloging:
            io.write_exposure_source_catalog(catalog, outdir, fname_in)

    print('Successfully finished reducing ' + fname_in)

    print('CI reduction pipeline completed at: ' + str(datetime.utcnow()) + 
          ' UTC')
