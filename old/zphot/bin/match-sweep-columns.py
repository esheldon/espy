""" 

    %prog gal_sweep_run

    Match CAS photozs to the reduced sweeps outputs by ra,dec.  These are
    identified by the run number, e.g. gal02.  

"""
import sys
from sys import stdout

import zphot
import sdssgal

import numpy
import esutil

from optparse import OptionParser
parser=OptionParser(__doc__)


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        return

    gal_procrun = args[0]

    zp = zphot.cas.CasZphot(origin='uchicago', type='dr7pofz')

    zpcols = zp.open_columns()
    galcols = sdssgal.open_columns(gal_procrun)

    # match ra/dec between these two.  We will write a new column into the sdss
    # gals directory that is the index into the zphot

    # depth=10 is good for close matching
    depth=10
    radius = 2.0/3600.0 # degrees

    h=esutil.htm.HTM(depth)

    stdout.write("Reading zphot ra,dec\n")
    zra = zpcols['ra'][:]
    zdec = zpcols['dec'][:]

    stdout.write("Reading sdss gal ra,dec\n")
    gra = galcols['ra'][:]
    gdec = galcols['dec'][:]

    crap="""
    n=1000000
    print 'just doing top',n,'as a test'
    stdout.write("Reading zphot ra,dec\n")
    zra = zpcols['ra'][0:n]
    zdec = zpcols['dec'][0:n]

    stdout.write("Reading sdss gal ra,dec\n")
    n1=10000000
    gra = galcols['ra'][n1:n1+n]
    gdec = galcols['dec'][n1:n1+n]
    """

    stdout.write("Doing match\n")
    mz, mg, dist = h.match(zra, zdec, gra, gdec, radius)

    stdout.write("Found {n} matches\n".format(n=mz.size))
    u=numpy.unique1d(mz)
    stdout.write("Found {nuniq} unique matches\n".format(nuniq=u.size))

    ngal = gra.size
    del zra
    del zdec
    del gra
    del gdec

    zid = numpy.zeros(ngal,dtype='i4')
    zid[:] = -1

    # numpy complains because mz is 64-bit
    zid[mg] = mz.astype('i4')

    stdout.write("Writing zid column\n")
    galcols.write_column('zid', zid, create=True, type='col')
    # this is crazy slow over nfs
    stdout.write("Creating index on zid\n")
    galcols['zid'].create_index()


if __name__=='__main__':
    main()

