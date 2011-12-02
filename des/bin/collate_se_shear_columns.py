"""
usage:
    %prog [options] run

Run this on tutti or it will take days!
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

options, args = parser.parse_args(sys.argv[1:])


# "all" probably here
if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run=args[0]
split=options.split
if split is not None:
    split=int(split)
fits=options.fits

if run[0:2] == 'se':
    c = des.collate.SEColumnCollator(run,split=split)
else:
    c = des.collate.MEColumnCollator(run)
if fits:
    c.convert2fits()
else:
    c.collate()
