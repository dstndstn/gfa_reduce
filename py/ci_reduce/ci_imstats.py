#!/usr/bin/env python

import argparse
import os
import ci_reduce.io as io
import glob
import time
import numpy as np
import copy

def print_imstats_1exp(imstats, fname_in, verbose=False):

    if verbose:
        for row in imstats:
            print(row['camera'], fname_in)
            for c in row.colnames:
                print(' '*5, '{:16}'.format(c), ' : ', row[c])

    cols = ['camera', 'median', 'max', 'min', 'sigma']

    _imstats = copy.deepcopy(imstats)

    # save horizntal space on printouts
    _imstats['median'] = np.round(_imstats['median']).astype('int')
    _imstats['min'] = np.round(_imstats['min']).astype('int')
    _imstats['max'] = np.round(_imstats['max']).astype('int')
    _imstats['sigma'] = np.round(_imstats['sig_robust']).astype('int')

    print(_imstats[cols])
    print('*sigma column is a robust standard deviation measurement')
    print('**all pixel values quoted are in ADU')

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
