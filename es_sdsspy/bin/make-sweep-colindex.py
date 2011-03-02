"""
    %prog [options] type
"""
import sys
import sdsspy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-p","--primary",
                  action="store_true",
                  default=False,
                  help="Only do primary objects")

parser.add_option("-d","--tempdir",
                  default=None,
                  help="Use this tempdir, e.g. /dev/shm")
parser.add_option("-c","--columns",
                  default=None,
                  help="Index these columns, CSV list")

options, args = parser.parse_args(sys.argv[1:])
if len(args) < 1:
    parser.print_help()
    sys.exit(45)

type = args[0]
tempdir = options.tempdir
primary = options.primary
columns=options.columns
if columns is not None:
    columns = columns.split(',')

#make_indices = options.make_indices
sel = sdsspy.sweeps.ColumnSelector(type,primary=primary)
sel.make_indices(tempdir=tempdir,columns=columns)

