#!/usr/bin/env python

import sys
import glob
import os

def ci_reduce_1command(fname_in, outdir):
    basename = os.path.basename(fname_in)

    pos = basename.find('.fits')

    assert(pos != -1)

    expname = basename[0:pos]
    this_outdir = os.path.join(outdir, expname)

    cmd = 'python -u /global/homes/a/ameisner/ci_reduce/py/ci_reduce/ci_proc.py ' + fname_in + ' --outdir ' + this_outdir + ' &> ' + expname + '.log &'
    return cmd

def ci_reduce_commands(data_dir, outdir):
    flist = glob.glob(data_dir + '/dci-?????.fits')

    cmds = [ci_reduce_1command(f, outdir) for f in flist]

    return cmds

if __name__ == "__main__":
    # add kw arg for backgrounding ampersand

    data_dir = sys.argv[1]

    outdir = '/project/projectdirs/desi/users/ameisner/CI/ci_data_challenge/results'
    cmds = ci_reduce_commands(data_dir, outdir)
    for cmd in cmds:
        print(cmd)
