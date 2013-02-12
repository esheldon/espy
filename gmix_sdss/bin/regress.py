"""
    %prog [options] gmix_run type nperbin
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--start',default=None,
                  help="start processing at this run,camcol csv")
parser.add_option('--index',action='store_true',
                  help="create the indexes")

parser.add_option('-r','--run',default=None,
                  help="Do runwise regressions")
parser.add_option('-c','--camcol',default=None,
                  help="Do run/camcol regressions")


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    type=args[1]
    nperbin=int(args[2])

    run=None
    if options.run is not None:
        run=int(options.run)
    camcol=None
    if options.camcol is not None:
        camcol=int(options.camcol)

    if type=='s2n':
        if run is not None:
            reg=gmix_sdss.regress.RunS2NRegressor(gmix_run=gmix_run,
                                                  run=run,
                                                  camcol=camcol,
                                                  nperbin=nperbin)
        else:
            reg=gmix_sdss.regress.S2NRegressor(gmix_run=gmix_run,
                                               camcol=camcol,
                                               nperbin=nperbin)

        reg.doplot()
main()
