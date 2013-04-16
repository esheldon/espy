"""
    %prog [options] gmix_run zphot_run

E.g.
    %prog gmix-r02 12

Note for epsf, the correction for <g> vs s2n is applied
"""
from __future__ import print_function
import sys
import numpy
import gmix_sdss
import zphot
import esutil as eu

from optparse import OptionParser
parser=OptionParser(__doc__)


def zphot_match(gmix_run, zphot_run):
    '''
    add a match index into photoz column databases

    The column name will be match_zphot{zphot_run}
    '''

    scols = gmix_sdss.collate.open_columns(gmix_run)
    print("  #rows:",scols['photoid'].size)
    pzcols = zphot.weighting.open_pofz_columns(zphot_run)
    print("  #rows:",pzcols['photoid'].size)

    print("reading num from zphot")
    num = pzcols['num'][:]

    print("reading photoid from sources")
    s_photoid = scols['photoid'][:]
    print("reading photoid from zphot")
    pz_photoid = pzcols['photoid'][:]

    print("matching")
    ms, mpz =  eu.numpy_util.match(s_photoid, pz_photoid)
    print("  matched:",ms.size)


    print("now determining which zphot are recoverable")
    w_recover, = numpy.where(num[mpz] > 0)
    print("  found:",w_recover.size)

    ms = ms[w_recover]
    mpz = mpz[w_recover]

    matches = numpy.empty(s_photoid.size, dtype='i4')
    matches[:] = -1

    # must explicitly convert to i4
    matches[ms] = numpy.array(mpz, dtype='i4')

    match_column = 'match_zphot%s' % zphot_run
    print("Adding zphot match column:",match_column)
    scols.write_column(match_column, matches, create=True)
    print("Creating index")
    scols[match_column].create_index()

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    zphot_run=args[1]

    zphot_match(gmix_run, zphot_run)
 
main()
