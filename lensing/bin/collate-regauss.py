"""
    %prog procrun

Description:

    Collate the outputs from regauss into a columns database.
    procrun should be, for example
        %prog 04
"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-i","--add-indices",action="store_true",
                  help="create the indices")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)

procrun = args[0]
add_indices = options.add_indices

c = lensing.regauss.Collator(procrun)

if add_indices:
    c.create_indices()
else:
    c.collate_as_columns_byband()
