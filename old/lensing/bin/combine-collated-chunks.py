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
parser.add_option('--fs',default='nfs',
                  help="which file system")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run=args[0]

    outfile = lensing.files.collated_file(sample=run,
                                          fs=options.fs)
    c=lensing.files.cascade_config(run)
    nsplit=c['lens_config']['nsplit']

    print("processing",nsplit,"lens splits")

    print(outfile)

    if options.fs=='hdfs':
        hdfs_file=eu.hdfs.HDFSFile(outfile,verbose=True)
        local_file=hdfs_file.localfile
    else:
        local_file=outfile

    print("opening:",local_file)
    with fitsio.FITS(local_file,'rw',clobber=True) as fobj:

        for i in xrange(nsplit):
            t=lensing.files.collated_read(sample=run,
                                          lens_split=i,
                                          verbose=False,
                                          fs=options.fs)
            if i==0:
                fobj.write(t)
            else:
                fobj[-1].append(t)

    if options.fs=='hdfs':
        hdfs_file.put(clobber=True)
        hdfs_file.cleanup()
main()
