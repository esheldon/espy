"""
    %prog procrun

Description:

    Collate the outputs from regauss into a columns database.
    procrun should be, for example, 02"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)

procrun = args[0]

c = lensing.regauss.Collator(procrun)
c.collate_as_columns_byband()
c.create_indexes()
