"""
    %prog run
"""

import sys
from sys import stderr
import os
import numpy
import esutil as eu
from esutil.misc import wlog
from shapesim import shapesim

from shapesim.shapesim import get_npair_nsplit,get_default_fs, \
        read_config, get_output_url

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-f','--full',action='store_true',
                  help="if true, check number of rows as well")


def check_is2n(run, is2n, flist, full=False):
    fs=get_default_fs()
    c = read_config(run)

    npair,nsplit = get_npair_nsplit(c, is2n)
    ntot=npair*2

    for isplit in xrange(nsplit):
        f=get_output_url(run, 0, is2n, itrial=isplit, fs=fs)

        print >>stderr,'checking:',f
        if f not in flist:
            print 'missing:',f
            continue

        if full:
            t=eu.io.read(f)

            if t.size != ntot:
                print "expected %d, got %d" % (npair,t.size)

def get_flist(run):
    import glob
    fs=get_default_fs()
    f=get_output_url(run, 0, 0, itrial=0, fs=fs)
    d=os.path.dirname(f)

    flist=glob.glob(d+'/*.rec')
    return flist

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]


    c = shapesim.read_config(run)

    ns2n=len(c['s2n_vals'])

    flist=get_flist(run)

    for is2n,s2n in enumerate(c['s2n_vals']):
        check_is2n(run, is2n, flist, full=options.full)

main()
