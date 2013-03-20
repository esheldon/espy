#!/usr/bin/env python
"""
    %prog [options] run mock_catalog type

Plot comparisons of shears

mock_catalog should be e.g. 'desmocks-3.02'
type should be one of
    mean_vs_true - plot mean shear vs true, and linear fits
    diff_vs - plot the shear diff as a function of various variables

The above also write out html files.
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
                  help="nperbin when binning")
parser.add_option('--nperbin-sub',default=None,
                  help="nperbin within a given S/N bin for mean shear plot")
parser.add_option('-s','--show',action="store_true",
                  help="show plots")

def set_nperbin(run, type, nperbin, nperbin_sub):
    if nperbin is None:
        if run[0:2] == 'me':
            if type == 'mean_vs_true':
                nperbin = 300000
            else:
                nperbin = 100000
        else:
            if type == 'mean_vs_true':
                nperbin = 3000000
            else:
                nperbin = 1000000
    
    if nperbin_sub is None:
        if run[0:2] == 'me':
            nperbin_sub = 10000
        else:
            nperbin_sub = 100000

    return nperbin, nperbin_sub

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) != 3:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    mock_catalog=args[1]
    types=args[2].split(',')

    c=des.compare_truth.Comparator(run,mock_catalog) 

    for type in types:
        nperbin,nperbin_sub=set_nperbin(run,type,options.nperbin,
                                        options.nperbin_sub)
        if type=='mean_vs_true':
            c.plot_meanshear_vs_trueshear_bys2n(nperbin, nperbin_sub, 
                                                show=options.show)
        if type=='diff_vs':
            c.plot_sheardiff_vs(nperbin, show=options.show)

main()
