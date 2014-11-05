#!/usr/bin/env python
"""
plot dsig for the input run and binning scheme

Note jack also implies corr
"""
from __future__ import print_function
import os,sys
import des
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run", help="run identifier")
parser.add_argument("bin_scheme", help="bin scheme, e.g. bin-lgt05-zwide")
parser.add_argument("-t","--type", 
                    help="type to plot",
                    choices=['bin', 'corr', 'jack'],
                    default='bin')
   
def main():
    args = parser.parse_args()

    if args.type=='bin':
        data=des.files_common.read_binned(args.run,
                                          args.bin_scheme)
    des.plotting.plot_dsig(data)

main()
