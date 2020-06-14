#!/usr/bin/env python

import astropy.io.fits as fits
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

def color_frame(sidelen):
    x = [0, 0, sidelen-1, sidelen-1, 0]
    y = [0, sidelen-1, sidelen-1, 0, 0]

    plt.plot(x, y, c='orange', linewidth=2)

def focus_plots(night, expids,
                basedir='/n/home/datasystems/users/ameisner/reduced/focus',
                outdir='/n/home/desiobserver/focus_scan', no_popups=False,
                dont_plot_centroid=False, n_stars_min=-1):

    plt.figure(1, figsize=(12.0*(len(expids)/7.0), 9))
    extnames = ['GUIDE0', 'GUIDE2', 'GUIDE3', 'GUIDE5', 'GUIDE7', 'GUIDE8']

    focus_z = []
    fwhm_pix = []

    # PSF stamps plot
    plt.subplots_adjust(hspace=0.01, wspace=0.01)
    n_stamps_plotted = 0
    for i, expid in enumerate(expids):
        fname = basedir + '/' + night + '/' + str(expid).zfill(8) + '/gfa-' + str(expid).zfill(8) + '_psfs.fits'
        print(fname)
        fname_ccds = fname.replace('_psfs.fits', '_ccds.fits')
        if not os.path.exists(fname):
            continue
        ccds = fits.getdata(fname_ccds)

        if np.sum(np.isfinite(ccds['PSF_FWHM_PIX'])) != 0:
            fwhm_pix.append(np.median(ccds['PSF_FWHM_PIX'][np.isfinite(ccds['PSF_FWHM_PIX'])]))
            focus_z.append(float(ccds[0]['FOCUS'].split(',')[2]))
            
        hdul = fits.open(fname)
        extnames_present = [hdu.header['EXTNAME'] for hdu in hdul]
        for j, extname in enumerate(extnames):
            if extname not in extnames_present:
                continue
            print(i, j)
            plt.subplot(6, len(expids), len(expids)*j + i +  1)
            plt.xticks([])
            plt.yticks([])
            im, h = fits.getdata(fname, extname=extname, header=True)
            plt.imshow(im, interpolation='nearest', origin='lower', cmap='gray_r', vmin=0.01)
            n_stamps_plotted += 1
            n_stars = h['NSTARS']
            if n_stars < n_stars_min:
                print('expid = ', h['EXPID'], ' ; extname = ', h['EXTNAME'], ' has too few contributing sources')
                color_frame(im.shape[0])
                
            plt.text(5, 44, str(expid) + '; ' + extname, color='r', fontsize=9)
            plt.text(10, 3.5, 'z = ' + str(int(float(ccds[0]['FOCUS'].split(',')[2]))), color='r')
            
            if np.isfinite(ccds[j]['XCENTROID_PSF']) and np.isfinite(ccds[j]['YCENTROID_PSF']) and (not dont_plot_centroid):
                plt.scatter([ccds[j]['XCENTROID_PSF']], [ccds[j]['YCENTROID_PSF']], marker='.', c='r')

    
    expid_min = int(np.min(expids))

    print(focus_z)
    print(fwhm_pix)

    if n_stamps_plotted > 0:
        plt.savefig(os.path.join(outdir, 'stamps_focus_scan-' + str(expid_min).zfill(8)+'.png'), bbox_inches='tight')
    else:
        print('WARNING : NO PSF MODELS WERE CREATED IN ANY IMAGES !!!')

    # doesn't make sense to fit a parabola to < 3 data points...
    if len(focus_z) < 3:
        print('NOT ENOUGH EXPOSURES AVAILABLE TO FIT A PARABOLA')
        return
    
    plt.figure(200)
    
    asec_per_pix = 0.205

    focus_z = np.array(focus_z)
    fwhm_asec = np.array(fwhm_pix)*asec_per_pix
    plt.scatter(focus_z, fwhm_asec)
    plt.xlabel('focus z (micron)')
    plt.ylabel('FWHM (asec)')

    coeff = np.polyfit(focus_z, fwhm_asec, 2)

    xsamp = np.arange(np.min(focus_z), np.max(focus_z))
    ysamp = coeff[0]*(np.power(xsamp, 2)) + coeff[1]*xsamp + coeff[2]

    plt.title('focus scan starting with EXPID = ' + str(expid_min))

    
    plt.plot(xsamp, ysamp)

    zmin = -coeff[1]/(2*coeff[0])

    min_fwhm_fit_asec = coeff[0]*(zmin**2) + coeff[1]*zmin + coeff[2]
    
    yrange = [np.min(fwhm_asec), np.max(fwhm_asec)]
    plt.text(focus_z[2], yrange[0] + 0.8*(yrange[1]-yrange[0]), 'best FWHM (meas) : ' + '{:.2f}'.format(np.min(fwhm_asec)))
    plt.text(focus_z[2], yrange[0] + 0.7*(yrange[1]-yrange[0]), 'best FWHM (fit) : ' + '{:.2f}'.format(min_fwhm_fit_asec))
    plt.text(focus_z[2], yrange[0] + 0.9*(yrange[1]-yrange[0]), 'best focus : ' + str(int(np.round(zmin))))
    
    plt.savefig(os.path.join(outdir, 'fit_focus_scan-' + str(expid_min).zfill(8) + '.png'), bbox_inches='tight')
    if not no_popups:
        plt.show()
    
def _test():
    night = '20200131'
    expids = 45446 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid', outdir='.')

def _test_missing_cam():
    night = '20200131'
    expids = 45485 + np.arange(7)

    focus_plots(night, expids, basedir='/project/projectdirs/desi/users/ameisner/GFA/run/psf_flux_weighted_centroid')

if __name__ == "__main__":
    descr = 'GFA focus sequence plots/analysis'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('first_expid', type=int, nargs=1)

    parser.add_argument('night', type=str, nargs=1)
    
    parser.add_argument('--basedir', default='/n/home/datasystems/users/ameisner/reduced/focus',
                        type=str, help='base directory for GFA reductions')

    parser.add_argument('--outdir', default='/n/home/desiobserver/focus_scan', 
                        type=str, help='output directory for plot PNGs')

    parser.add_argument('--no_popups', default=False, action='store_true',
                        help='write PNGs without popping up plot windows')

    parser.add_argument('--dont_plot_centroid', default=False, action='store_true',
                        help='do not overplot centroid location on each postage stamp')

    args = parser.parse_args()

    expids = args.first_expid + np.arange(7, dtype=int)

    print(expids)
    print(args.night[0])
    print(args.basedir)
    
    outdir = args.outdir if os.path.exists(args.outdir) else '.'
    focus_plots(args.night[0], expids, basedir=args.basedir, outdir=outdir, no_popups=args.no_popups,
                dont_plot_centroid=args.dont_plot_centroid)
