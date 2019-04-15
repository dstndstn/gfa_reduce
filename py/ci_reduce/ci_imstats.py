#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
import glob
import time
import numpy as np

def print_imstats_1exp(imstats, fname_in, verbose=False):

    if verbose:
        for row in imstats:
            print(row['camera'], fname_in)
            for c in row.colnames:
                print(' '*5, '{:16}'.format(c), ' : ', row[c])

    cols = ['camera', 'median', 'max', 'min', 'sig_robust']

    # save horizntal space on printouts
    imstats['median'] = np.round(imstats['median']).astype('int')
    imstats['min'] = np.round(imstats['min']).astype('int')
    imstats['max'] = np.round(imstats['max']).astype('int')
    imstats['sig_robust'] = np.round(imstats['sig_robust']).astype('int')

    print(imstats[cols])

def imstats_1exp(fname_in, verbose=False):

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
        imstats_1exp(f)

    #fname = '/project/projectdirs/desi/spectro/data/20190403/00003746/ci-00003746.fits.fz'
    imstats_1exp(fname)
