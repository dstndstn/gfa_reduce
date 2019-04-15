#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
import glob
import time

def print_imstats_1exp(imstats, fname_in, verbose=False):

    if verbose:
        for row in imstats:
            print(row['camera'], fname_in)
            for c in row.colnames:
                print(' '*5, '{:16}'.format(c), ' : ', row[c])

    cols = ['camera', 'median', 'mean', 'max', 'min', 'sig_robust']

    print(imstats[cols])

def ql_stats_1exp(fname_in, verbose=False):

    for i in range(5):
        exp = io.load_exposure(fname_in, realtime=True)
        if exp is not None:
            break
        else:
            time.sleep(2.0)

    imstats = io.gather_pixel_stats(exp)

    colnames = imstats.colnames
    colnames.insert(0, colnames.pop())

    imstats = imstats[colnames]
    print_imstats_1exp(imstats, fname_in, verbose=verbose)


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
