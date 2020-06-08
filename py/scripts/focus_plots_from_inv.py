from focus_scan import focus_plots
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

fname_inv = '/global/homes/a/ameisner/gfa/pro/focus_scan_inventory-thru_20200315.fits'

def plot_one_row(ind, inv=None, outdir=None, basedir=None,
                 dont_plot_centroid=False):

    if inv is None:
        inv = fits.getdata(fname_inv)

    row = inv[ind]

    night = row['NIGHT']
    expids = row['EXPIDS'][row['EXPIDS'] != -1]

    if basedir is None:
        basedir = '/global/cfs/cdirs/desi/users/ameisner/GFA/reduced/v0009'

    if outdir is None:
        outdir = '/global/cscratch1/sd/ameisner/focus_plots'

    print('Working on index : ', ind)
    focus_plots(night, expids, basedir=basedir, outdir=outdir,
                no_popups=True, dont_plot_centroid=dont_plot_centroid)

def _loop(indstart=0, nproc=None):

    inv = fits.getdata(fname_inv)

    if nproc is None:
        nproc = len(inv)

    indend = min(indstart + nproc, len(inv))
    
    for ind in np.arange(indstart, indend):
        plot_one_row(ind, inv=inv)
        plt.cla()
