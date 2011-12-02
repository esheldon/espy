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

parser.add_option("-s","--split",
                  default=None,
                  help="Process the split for SE.  Should be 1 or 2")
parser.add_option("--fits",
                  default=False,action='store_true',
                  help="Convert to fits")

parser.add_option("--index",
                  default=False,action='store_true',
                  help="Create the indexes")

parser.add_option("--html",
                  default=False,action='store_true',
                  help="Create the html")

options, args = parser.parse_args(sys.argv[1:])


# "all" probably here
if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run=args[0]
split=options.split
if split is not None:
    split=int(split)

if run[0:2] == 'se':
    c = des.collate.SEColumnCollator(run,split=split)
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
