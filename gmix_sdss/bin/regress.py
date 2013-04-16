"""
    %prog [options] gmix_run regress_type camcol

type should be
    s2n epsf

Note for epsf, the correction for <g> vs s2n is applied
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--nbin',default=40, help="number of bins")


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    regress_type=args[1]
    camcol=int(args[2])

    nbin=int(options.nbin)

    if regress_type=='s2n':
        reg=gmix_sdss.regress.S2NRegressor(gmix_run=gmix_run,
                                           camcol=camcol,
                                           nbin=nbin)

        reg.doplot()
    elif regress_type=='epsf':
        reg=gmix_sdss.regress.EPSFRegressor(gmix_run=gmix_run,
                                            camcol=camcol)
        reg.doplot()
    else:
        raise ValueError("only support type s2n,epsf for now")
main()
