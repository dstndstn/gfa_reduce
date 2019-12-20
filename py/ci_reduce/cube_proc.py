import argparse
from ci_reduce.gfa_red import _proc
import astropy.io.fits as fits
import os

if __name__ == "__main__":
    descr = 'run gfa_reduce pipeline on all GFA guide cube slices'

    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('fname_in', type=str, nargs=1)

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')

    args = parser.parse_args()
    
    fname_in = args.fname_in[0]

    assert(os.path.exists(fname_in))
    
    h = fits.getheader(fname_in, extname='GUIDER')

    nframes = h['FRAMES']
    
    for i in range(nframes):
        # note the hardcoding of certain arguments, especially 
        # skipping write-out of image-level outputs (saves disk space...)
        print('WORKING ON FRAME ' + str(i+1) + ' OF ' + str(nframes))
        _proc(fname_in, outdir=args.outdir,
              careful_sky=False, no_cataloging=False,
              no_gaia_xmatch=False, cube_index=i, skip_image_outputs=True,
              realtime=False)
