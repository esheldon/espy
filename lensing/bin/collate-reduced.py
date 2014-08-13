"""
    %prog [options] run

Description:

    Collate the reduced lensout with the original catalog.

"""

from __future__ import print_function

import sys
import numpy
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
parser.add_option("--fs",default='nfs',
                  help="file system. %default")

def add_im3_dt(dt):
    dt += [('pT','f8'),('pe1','f8'),('pe2','f8'),('pe','f8')]
    return dt

def add_ellip(data):
    data['pe']  =  9999.0
    data['pe1'] = -9999.0
    data['pe2'] = -9999.0
    #w,=numpy.where(data['totpairs'] > 10)
    #if w.size > 0:
    #    ninv  = 1.0/(data['totpairs'][w]-1)
    w,=numpy.where( (data['totpairs'] > 10) & (data['weight'] > 0.0) )
    if w.size > 0:
        print("kept:",w.size,"for ellip")
        # weight is also wsum
        winv  = 1.0/data['weight'][w]
        T     = data['x2sum'][w] + data['y2sum'][w]
        x2my2 = data['x2sum'][w] - data['y2sum'][w]
        xy    = data['xysum'][w]

        T     *= winv
        x2my2 *= winv
        xy    *= winv

        data['pT'][w] = T

        w2,=numpy.where(T > 0.0)
        if w2.size > 0:

            pe1 = x2my2[w2]/T[w2]
            pe2 = 2.0*xy[w2]/T[w2]

            data['pe1'][w[w2]] = pe1
            data['pe2'][w[w2]] = pe2
            data['pe'][w[w2]] = numpy.sqrt(pe1**2 + pe2**2)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run = args[0]
    # this can be None for now
    lens_split = options.lens_split
    if lens_split is not None:
        lens_split=int(lens_split)


    conf = lensing.files.cascade_config(run)
    lsample = conf['lens_config']['sample']
    nbin = conf['lens_config']['nbin']

    nsplit=conf['lens_config']['nsplit']
    cat = lensing.lcat.read_original(sample=lsample,
                                     lens_split=lens_split,
                                     fs=options.fs)

    # in this case zindex must be in the catalog, so we can match
    if nsplit > 1 and 'zindex' not in cat.dtype.names:
        raise ValueError("when collating nsplit > 1, zindex must be in "
                         "'original' catalog")
    reduced_file = lensing.files.sample_file(type='src-reduced',
                                             sample=run,
                                             lens_split=lens_split,
                                             fs=options.fs)
    outfile = lensing.files.collated_file(sample=run,
                                          lens_split=lens_split,
                                          fs=options.fs)

    print("will collate to file:",outfile)

    if 'im3' in run:
        dt=lensing.files.lensout_im3shape_dtype(nbin)
    else:
        dt=lensing.files.lensout_dtype(nbin)

    # we have to stage the reduced file to local disk first
    if options.fs=='hdfs':
        hdfs_red=HDFSFile(reduced_file,verbose=True)
        hdfs_coll=HDFSFile(outfile,verbose=True)

        hdfs_red.stage()

        red_local=hdfs_red.localfile
        coll_local=hdfs_coll.localfile
    else:
        red_local=reduced_file
        coll_local=outfile

    # now open the local files
    print("opening local files")
    print(red_local)
    print(coll_local)
    with Recfile(red_local,mode='r',delim=' ',dtype=dt) as robj, \
            FITS(coll_local,mode='rw',clobber=True) as fits:

        print("reading")
        reduced=robj.read()

        print("matching")
        if nsplit > 1:
            rind, cind = eu.numpy_util.match(reduced['zindex'],cat['zindex'])
            if rind.size != reduced.size:
                raise ValueError("all zindex did not match")
            catsub = cat[cind]
            reduced = reduced[rind]  # make sure they line up
        else:
            catsub = cat[ reduced['zindex'] ]

        add_dt = [d for d in reduced.dtype.descr if d[0] != 'zindex']
        if 'im3' in run:
            add_dt = add_im3_dt(add_dt)

        data = eu.numpy_util.add_fields(catsub, add_dt)
        eu.numpy_util.copy_fields(reduced, data)

        if 'im3' in run:
            add_ellip(data)

        print("writing")
        fits.write(data)

    if options.fs=='hdfs':
        # this must be outside the FITS context so the fits file gets closed
        # and is in a consistent state
        hdfs_coll.put(clobber=True)

        hdfs_coll.cleanup()
        hdfs_red.cleanup()



main()
