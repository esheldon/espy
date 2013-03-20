"""
    %prog [options] gmix_run
"""
import sys, os
import gmix_sdss

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--start',default=None,
                  help="start processing at this run,camcol csv")
parser.add_option('--index',action='store_true',
                  help="create the indexes")


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    gmix_run=args[0]


    start=options.start
    if start is not None:
        run,camcol=start.split(',')
        run=int(run)
        camcol=int(camcol)
        start={'run':run,'camcol':camcol}
        print 'start:\n',start

    if options.index:
        cols=gmix_sdss.collate.open_columns(gmix_run=gmix_run)
        for icol in gmix_sdss.collate.index_cols:
            print 'index:',icol
            cols[icol].create_index(force=True)
    else:
        cm=gmix_sdss.collate.ColumnsMaker(gmix_run=gmix_run,
                                          start=start)
        cm.make_columns()

main()
