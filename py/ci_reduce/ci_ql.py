#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
import glob

def ql_stats_1exp(fname_in):
    exp = io.load_exposure(fname_in)
    imstats = io.gather_pixel_stats(exp)

    colnames = imstats.colnames
    colnames.insert(0, colnames.pop())

    imstats = imstats[colnames]
    print(imstats)


if __name__=="__main__":
    # works on simulated data
    #fname_in = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/sims/dci-01402.fits'

    # try a real frame
    #fname_in = '/project/projectdirs/desi/spectro/data/20190330/00002930/ci-00002930.fits.fz'

    flist = glob.glob('/project/projectdirs/desi/spectro/data/20190330/*/ci*.fits.fz')

    for f in flist:
        ql_stats_1exp(f)

    #fname = '/project/projectdirs/desi/spectro/data/20190403/00003746/ci-00003746.fits.fz'
    #ql_stats_1exp(fname)
