"""
    %prog [options] gmix_run sdss_run camcol field
"""
import sys, os
import gmix_sdss

from optparse import OptionParser

parser=OptionParser(__doc__)

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 4:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]
    run=int(args[1])
    camcol = int(args[2])
    field=int(args[3])
    pipe=gmix_sdss.pipe.GMixField(gmix_run, run, camcol, field)

    pipe.go()

main()
