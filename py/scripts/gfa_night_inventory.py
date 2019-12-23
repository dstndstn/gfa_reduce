#!/usr/bin/env python3

import argparse
import astropy.io.fits as fits
import glob

def _printouts(night):
    basedir = '/project/projectdirs/desi/spectro/data'
    nightdir = basedir + '/' + night
    flist = glob.glob(nightdir + '/*/gfa*.fits.fz')

    flist.sort()

    for f in flist:
        h = fits.getheader(f, extname='GFA')
        print(f.replace(basedir + '/', '') + '   ', h['FLAVOR'].ljust(10, ' '), '{:.1f}'.format(h['EXPTIME']).ljust(10, ' '), h['PROGRAM'])

if __name__ == "__main__":
    descr = 'print information about GFA prescan/overscan bad pixels'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('night', type=str, nargs=1,
                        help='NIGHT string')

    args = parser.parse_args()

    night = args.night[0]

    _printouts(night)
