"""
    %prog [options] gmix_run regress_type

type should be
    s2n epsf
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--start',default=None,
                  help="start processing at this run,camcol csv")

parser.add_option('-r','--run',default=None,
                  help="Do runwise regressions")
parser.add_option('-c','--camcol',default=None,
                  help="Do run/camcol regressions")
parser.add_option('--s2n',default='s2n',
                  help="Field for s/n")
parser.add_option('--objtype',default=None, help="obj type")
parser.add_option('--nbin',default=40, help="number of bins")


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    regress_type=args[1]

    run=None
    if options.run is not None:
        run=int(options.run)
    camcol=None
    if options.camcol is not None:
        camcol=int(options.camcol)
    nbin=int(options.nbin)

    if regress_type=='s2n':
        if run is not None:
            reg=gmix_sdss.regress.RunS2NRegressor(gmix_run=gmix_run,
                                                  run=run,
                                                  camcol=camcol,
                                                  nbin=nbin,
                                                  s2n_field=options.s2n,
                                                  objtype=options.objtype)
        else:
            reg=gmix_sdss.regress.S2NRegressor(gmix_run=gmix_run,
                                               camcol=camcol,
                                               nbin=nbin,
                                               s2n_field=options.s2n,
                                               objtype=options.objtype)

        reg.doplot()
    elif regress_type=='epsf':
        reg=gmix_sdss.regress.EPSFRegressor(gmix_run=gmix_run,
                                           camcol=camcol,
                                           nbin=nbin,
                                           objtype=options.objtype)
        reg.doplot()
    else:
        raise ValueError("only support type s2n,epsf for now")
main()
