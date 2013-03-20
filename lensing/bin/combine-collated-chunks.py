"""
    %prog [options] run

Description:
    Combine the chunks into a single collated file
"""
from __future__ import print_function

import os
import sys
import lensing

import fitsio

import esutil as eu

import numpy
from numpy import where, zeros, ones

from optparse import OptionParser

parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run=args[0]

    outfile = lensing.files.collated_file(sample=run)
    c=lensing.files.cascade_config(run)
    nsplit=c['lens_config']['nsplit']

    print("processing",nsplit,"lens splits")

    # chunks are not the same size, so we need an initial pass
    ntot=0
    ii=0
    for ipass in [1,2]:
        if ipass==1:
            print('pass 1 to get sizes')
        else:
            print('now copying')
        for i in xrange(nsplit):
            t=lensing.files.collated_read(sample=run,lens_split=i,
                                          verbose=False)
            if ipass == 1:
                print(t.size)
                ntot += t.size
            else:
                if i == 0:
                    print('-'*70)
                    print('ntot:',ntot)
                    descr=eu.numpy_util.descr_to_native(t.dtype.descr)
                    data = numpy.zeros(ntot, dtype=descr)
                    print(data.dtype.descr)

                data[ii:ii+t.size] = t
                ii += t.size

    eu.io.write(outfile, data, clobber=True, verbose=True)
main()
