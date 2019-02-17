#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io

if __name__ == "__main__":
    descr = 'run full ci_reduce pipeline on a CI exposure'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default='', type=str,
                        help='directory to write outputs in')

    args = parser.parse_args()

    fname_in = args.fname_in[0]

    assert(os.path.exists(fname_in))

    exp = io.load_exposure(fname_in)

    # create data quality bitmasks
    exp.create_all_bitmasks()

    # go from "raw" images to "reduced" images
    exp.calibrate_pixels()

    # try to write image-level outputs if outdir is specified
    if len(args.outdir) > 0:
        outdir = args.outdir
        print('Attempting to write outputs to : ' + 
              outdir)

        if not os.path.exists(outdir):
            os.mkdir(outdir)

    print('Succesfully finished reducing ' + fname_in)
