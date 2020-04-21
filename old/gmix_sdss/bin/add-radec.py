"""
    %prog [options] gmix_run

Description:

    Add ra,dec to the columns
"""

from __future__ import print_function
import sys
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    gmix_run = args[0]

    gmix_sdss.collate.add_radec(gmix_run)

main()
