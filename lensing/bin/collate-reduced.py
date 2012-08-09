"""
    %prog [options] run

Description:

    Collate the reduced lensout with the original catalog.

"""

from __future__ import print_function

import sys
import lensing
import esutil as eu
from optparse import OptionParser


parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run = args[0]

    outfile = lensing.files.sample_file(type='collated', sample=run, fs='hdfs')
    print("will collate to file:",outfile)

    conf = lensing.files.cascade_config(run)
    lsample = conf['lens_config']['sample']
    cat = lensing.files.read_original_catalog(type='lens',sample=lsample)
    reduced = lensing.files.reduced_read(sample=run)

    # trim down to the ones we used
    print("Extracting by zindex")
    cat = cat[ reduced['zindex'] ]

    # create a new struct with the lensing outputs at the end
    print('collating')

    # for randoms, we use the input lcat so we don't want zindex
    # in the added tags
    add_dt = [d for d in reduced.dtype.descr if d[0] != 'zindex']
    data = eu.numpy_util.add_fields(cat, add_dt)
    eu.numpy_util.copy_fields(reduced, data)

    eu.io.write(outfile, data, verbose=True, clobber=True)


main()
