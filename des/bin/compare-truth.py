#!/usr/bin/env python
"""
    %prog [options] run mock_catalog

Plot comparisons of shears

mock_catalog should be e.g. 'desmocks-3.02'
"""
import os
import sys
from sys import stderr
import deswl
import lensing
import des
import esutil as eu
import numpy

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-n','--nperbin',default=None,
                  help="nperbin when binning by S/N")
parser.add_option('--nperbin-sub',default=None,
                  help="nperbin within a given S/N bin")
parser.add_option('-s','--show',action="store_true",
                  help="show plots")

def set_nperbin(run, nperbin, nperbin_sub):
    if nperbin is None:
        if run[0:2] == 'me':
            nperbin = 300000
        else:
            nperbin = 3000000
    
    if nperbin_sub is None:
        if run[0:2] == 'me':
            nperbin_sub = 10000
        else:
            nperbin_sub = 100000
    return nperbin, nperbin_sub
def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) != 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    mock_catalog=args[1]
    nperbin,nperbin_sub=set_nperbin(run,options.nperbin,options.nperbin_sub)
    c=des.compare_truth.Comparator(run,mock_catalog) 

    c.plot_sheardiff_bys2n(nperbin, nperbin_sub, show=options.show)

main()
