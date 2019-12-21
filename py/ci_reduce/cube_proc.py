import argparse
from ci_reduce.gfa_red import _proc
import astropy.io.fits as fits
import os
import numpy as np

if __name__ == "__main__":
    descr = 'run gfa_reduce pipeline on all GFA guide cube slices'

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--indstart', default=0, type=int,
                        help='starting frame index')

    parser.add_argument('--nproc', default=None, type=int,
                        help='number of frames to process')

    parser.add_argument('--write_image_outputs', default=False,
                        action='store_true',
                        help='write image-level reduction outputs')
    
    args = parser.parse_args()
    
    fname_in = args.fname_in[0]

    assert(os.path.exists(fname_in))
    
    h = fits.getheader(fname_in, extname='GUIDER')

    nframes = h['FRAMES']
    assert(args.indstart < nframes)

    nproc = args.nproc if args.nproc is not None else nframes
    
    indend = min(nframes, args.indstart + nproc)

    skip_image_outputs = not args.write_image_outputs
    for i in np.arange(args.indstart, indend):
        # note the hardcoding of certain arguments, especially 
        # skipping write-out of image-level outputs (saves disk space...)
        print('WORKING ON FRAME ' + str(i+1) + ' OF ' + str(nframes))
        _proc(fname_in, outdir=args.outdir,
              careful_sky=False, no_cataloging=False,
              no_gaia_xmatch=False, no_ps1_xmatch=False, cube_index=i,
              skip_image_outputs=skip_image_outputs,
              realtime=False)
