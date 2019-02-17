#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io

if __name__ == "__main__":
    descr = 'run full ci_reduce pipeline on a CI exposure'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    args = parser.parse_args()

    fname_in = args.fname_in[0]

    assert(os.path.exists(fname_in))

    exp = io.load_exposure(fname_in)

    # go from "raw" images to "reduced" images
    exp.calibrate_pixels()

    # create data quality bitmasks
    exp.create_all_bitmasks()

    print('Succesfully finished reducing ' + fname_in)
