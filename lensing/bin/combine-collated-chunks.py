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

    print(outfile)
    with eu.hdfs.HDFSFile(outfile,verbose=True) as hdfs_obj:
        print("opening:",hdfs_obj.localfile)
        with fitsio.FITS(hdfs_obj.localfile,'rw',clobber=True) as fobj:

            for i in xrange(nsplit):
                t=lensing.files.collated_read(sample=run,
                                              lens_split=i,
                                              verbose=False)
                if i==0:
                    fobj.write(t)
                else:
                    fobj[-1].append(t)

        hdfs_obj.put(clobber=True)
main()
