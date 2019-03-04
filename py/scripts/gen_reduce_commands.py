#!/usr/bin/env python

import argparse
import glob
import os

def ci_reduce_1command(fname_in, outdir, backgrounded=False):
    basename = os.path.basename(fname_in)

    pos = basename.find('.fits')

    assert(pos != -1)

    expname = basename[0:pos]
    this_outdir = os.path.join(outdir, expname)

    cmd = 'python -u /global/homes/a/ameisner/ci_reduce/py/ci_reduce/ci_proc.py ' + fname_in + ' --outdir ' + this_outdir + ' &> ' + expname + '.log'

    if backgrounded:
        cmd += ' &'

    return cmd

def ci_reduce_commands(data_dir, outdir, backgrounded=False):
    flist = glob.glob(data_dir + '/dci-?????.fits')

    cmds = [ci_reduce_1command(f, outdir, backgrounded=backgrounded) for f in flist]

    return cmds

if __name__ == "__main__":
    descr = 'create Python commands to run ci_reduce pipeline'
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('data_dir', type=str, nargs=1)

    default_outdir = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/results'

    parser.add_argument('--outdir', default=default_outdir, type=str,
                        help='directory to write outputs in')

    parser.add_argument('--backgrounded', default=False, action='store_true', 
                        help='add ampersand to the end of each command')

    args = parser.parse_args()

    data_dir = args.data_dir[0]
    outdir = args.outdir
    backgrounded = args.backgrounded

    cmds = ci_reduce_commands(data_dir, outdir, backgrounded=backgrounded)
    for cmd in cmds:
        print(cmd)
