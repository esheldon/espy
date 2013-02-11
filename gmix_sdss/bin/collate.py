"""
    %prog [options] gmix_run
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]

    cm=gmix_sdss.collate.ColumnsMaker(gmix_run=gmix_run)
    cm.make_columns()

main()
