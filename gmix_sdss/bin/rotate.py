"""
    %prog [options] gmix_run

Description:

    rotate the detrended ellipticities
"""

from __future__ import print_function
import sys
import numpy
import lensing
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    gmix_run=args[0]
    conf=gmix_sdss.files.read_config(gmix_run)

    cols=gmix_sdss.collate.open_columns(gmix_run)

    data=cols.read_columns(['run','camcol','field','g_dt'], verbose=True)
    
    grot = numpy.zeros( (data.size, 2) )
    grot[:,:] = -9999.
    angles = numpy.zeros(data.size)
    angles[:] = -9999.

    # we only detrended in a selected range, and the rest
    # are -9999
    print("selecting")
    w,=numpy.where(  (data['g_dt'][:,0] >= -1)
                   & (data['g_dt'][:,1] >= -1) )
    
    rotator = lensing.rotation.SDSSRotator("eq")

    print("getting rotations")
    tg1, tg2, tangle = rotator.rotate(data['run'][w],
                                      data['camcol'][w],
                                      data['field'][w],
                                      conf['filter'],
                                      data['g_dt'][w,0],
                                      data['g_dt'][w,1],
                                      getrot=True)

    grot[w,0] = tg1
    grot[w,1] = tg2
    angles[w] = tangle

    colname="g_dt_eq"
    print("writing:",colname)
    cols.write_column(colname, grot, create=True)

    colname="rot_eq"
    print("writing:",colname)
    cols.write_column(colname, angles, create=True)


main()
