"""
usage:
    %prog [options] run

Run this on tutti or it will take days for SE stuff!
"""
import des
import sys
from sys import stdout


from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-s","--small",
                  default=False,action='store_true',
                  help="Don't write interp_pfs_shapelets or shapelets_prepsf")
parser.add_option("--fits",
                  default=False,action='store_true',
                  help="Convert to fits")

parser.add_option("--index",
                  default=False,action='store_true',
                  help="Create the indexes")

parser.add_option("--html",
                  default=False,action='store_true',
                  help="Create the html")

parser.add_option("-n","--njob", default=None, help="Number of jobs for collate")
parser.add_option("-j","--job", default=None, help="This job number")

options, args = parser.parse_args(sys.argv[1:])


# "all" probably here
if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run=args[0]

if run[0:2] == 'se':
    c = des.collate.SEColumnCollator(run,small=options.small,
                                     njob=options.njob,job=options.job)
elif run[0:2] == 'me':
    c = des.collate.MEColumnCollator(run)
else:
    raise ValueError("Expected run 'me*' or 'se*'")

if options.fits:
    c.convert2fits()
elif options.index:
    c.create_indexes()
elif options.html:
    c.write_html()
else:
    c.collate()
