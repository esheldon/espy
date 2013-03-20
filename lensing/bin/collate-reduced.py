"""
    %prog [options] run

Description:

    Collate the reduced lensout with the original catalog.

"""

from __future__ import print_function

import sys
import lensing
import esutil as eu
from esutil.hdfs import HDFSFile
import fitsio
from fitsio import FITS
import recfile
from recfile import Recfile

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-s",dest="lens_split",default=None,
                  help="which lens split to collate. %default")


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run = args[0]
    # this can be None for now
    lens_split = options.lens_split


    conf = lensing.files.cascade_config(run)
    lsample = conf['lens_config']['sample']
    nbin = conf['lens_config']['nbin']

    if lens_split is not None:
        lens_split=int(lens_split)
        # in this case zindex must be in the catalog, so we can match
        cat = lensing.lcat.read_original(sample=lsample, lens_split=lens_split)
        if 'zindex' not in cat.dtype.names:
            raise ValueError("when collating splits, zindex must be in "
                             "'original' catalog")
        reduced_file = lensing.files.sample_file(type='src-reduced',
                                                 sample=run,
                                                 lens_split=lens_split,
                                                 fs='hdfs')
        outfile = lensing.files.collated_file(sample=run,lens_split=lens_split)
    else:
        reduced_file = lensing.files.reduced_file(sample=run)
        outfile = lensing.files.collated_file(sample=run)

    print("will collate to file:",outfile)

    dt=lensing.files.lensout_dtype(nbin)

    # we have to stage the reduced file to local disk first
    with HDFSFile(reduced_file,verbose=True) as hdfs_red, \
            HDFSFile(outfile,verbose=True) as hdfs_coll:

        hdfs_red.stage()

        # now open the local files
        print("opening local files")
        with Recfile(hdfs_red.localfile,mode='r',delim=' ',dtype=dt) as robj, \
                FITS(hdfs_coll.localfile,mode='rw',clobber=True) as fits:

            reduced=robj.read()

            if lens_split is not None:
                rind, cind = eu.numpy_util.match(reduced['zindex'],cat['zindex'])
                if rind.size != reduced.size:
                    raise ValueError("all zindex did not match")
                catsub = cat[cind]
                reduced = reduced[rind]  # make sure they line up
            else:
                catsub = cat[ reduced['zindex'] ]

            add_dt = [d for d in reduced.dtype.descr if d[0] != 'zindex']

            data = eu.numpy_util.add_fields(catsub, add_dt)
            eu.numpy_util.copy_fields(reduced, data)

            fits.write(data)

        # this must be outside the FITS context so the fits file gets closed
        # and is in a consistent state
        hdfs_coll.put(clobber=True)



main()
