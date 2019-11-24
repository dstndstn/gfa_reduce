#!/usr/bin/env python

import sys, os, time
import multiprocessing as mp
import argparse
import glob
from ci_reduce.gfa_red import _proc
from ci_reduce.common import expid_from_filename
import numpy as np
import astropy.io.fits as fits
import json

# using Stephen Bailey's "multirunner" template as the basis for this script
# https://raw.githubusercontent.com/sbailey/multirunner/master/multirunner.py

parser = argparse.ArgumentParser(usage = "{prog} [options]")
parser.add_argument("--night", type=str,  help="NIGHT string")
parser.add_argument("-n", "--numworkers", type=int,  default=1, help="number of workers")
parser.add_argument("-w", "--waittime", type=int, default=5, help="wait time between directory checks")
parser.add_argument("-e", "--expid_min", type=int, default=-1, help="start with this EXPID value")
parser.add_argument("--out_basedir", type=str, default='/n/home/datasystems/users/ameisner/reduced/realtime', help="base output directory for GFA reductions")
args = parser.parse_args()

def check_flavor_json(gfa_image_fname):
    gfa_json_fname = gfa_image_fname.replace('gfa-', 'request-')
    gfa_json_fname = gfa_json_fname.replace('.fits.fz', '.json')

    print(gfa_json_fname)
    assert(os.path.exists(gfa_json_fname))
 
    with open(gfa_json_fname) as json_file:
        data = json.load(json_file)

    return data['FLAVOR']

def is_flavor_science(gfa_image_fname):
    return check_flavor_json(gfa_image_fname).lower() == 'science'

indir = '/exposures/desi/' + args.night
assert(os.path.exists(indir))
print(args.out_basedir)
assert(os.path.exists(args.out_basedir))
night_basedir_out = os.path.join(args.out_basedir, args.night)
if not os.path.exists(night_basedir_out):
    os.mkdir(night_basedir_out)

#- Create communication queue to pass files to workers
q = mp.Queue()

#- Function to run for each worker.
#- Listens on Queue q for filenames to process.
def run(workerid, q):
    global night_basedir_out
    print('Worker {} ready to go'.format(workerid))
    while True:
        filename = q.get(block=True)
        print('Worker {} processing {}'.format(workerid, filename))
        sys.stdout.flush()
        #- Do something with that filename
        h = fits.getheader(filename, extname='GFA')
        if h['FLAVOR'].lower() == 'science':
            outdir = os.path.join(night_basedir_out, str(expid_from_filename(filename)).zfill(8))
            _proc(filename, outdir=outdir, realtime=True) # realtime HARDCODED to true
        else:
            print('New GFA file ' + filename + ' is not flavor=science ; skipping')
        print('Worker {} done with {}'.format(workerid, filename))
        sys.stdout.flush()

#- Start workers
for i in range(args.numworkers):
    p = mp.Process(target=run, args=(i, q))
    p.start()

#- Track what files have already been added to queue.
#- TODO: Upon startup, this could compare against files in output dir
#- and only load input files haven't already been processed.
exp_outdirs = glob.glob(night_basedir_out + '/????????')
known_files = set([indir + '/' + os.path.split(d)[-1] + '/gfa-' + os.path.split(d)[-1] + '.fits.fz' for d in exp_outdirs])

print('Number of known files = ', len(known_files))

#- Periodically check for any new files that may have appeared and add them
#- to the queue for a worker to process.

glob_pattern = os.path.join(indir, '????????/gfa*.fits.fz')
while(True):
    flist = glob.glob(glob_pattern)
    flist.sort()
    flist = np.array(flist)
    expids = np.array([expid_from_filename(f) for f in flist])
    flist = flist[expids >= args.expid_min]
    for filename in flist:
        if filename not in known_files:
            if is_flavor_science(filename):
                print('Server putting {} in the queue'.format(filename))
                sys.stdout.flush()
                q.put(filename)
            else:
                print('skipping ' + filename + ' ; NOT flavor=science')
            known_files.add(filename)

    time.sleep(args.waittime)
