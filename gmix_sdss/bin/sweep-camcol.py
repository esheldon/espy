"""
    %prog [options] gmix_run run camcol
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    run=int(args[1])
    camcol=int(args[2])

    gmix_sdss.select.sweep_camcol(gmix_run=gmix_run,
                                  run=run,
                                  camcol=camcol)

main()
