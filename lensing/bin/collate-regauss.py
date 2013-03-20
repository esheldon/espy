"""
    %prog procrun type

Description:

    Collate the outputs from regauss into a columns database.
    procrun should be, for example
        %prog 04 gal
"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-i","--add-indices",action="store_true",
                  help="create the indices")
parser.add_option("-r","--rotate",action="store_true",
                  help="Add the rotated e1/e2 in equatorial coords")
parser.add_option("--fs", help="file system for regauss outputs")
parser.add_option("-d", "--coldir", default=None, help="Use this columns dir")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(45)

procrun = args[0]
type = args[1]

c = lensing.regauss.Collator(procrun, sweeptype=type, fs=options.fs,coldir=options.coldir)

if options.add_indices:
    c.create_indices()
elif options.rotate:
    c.add_rotated_e1e2()
else:
    c.collate_as_columns_byband()
