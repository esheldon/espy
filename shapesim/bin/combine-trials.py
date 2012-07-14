"""
    %prog run i1 i2

i1 is usually is2, i2 can be is2n or ie
"""

import sys
import shapesim
from optparse import OptionParser

parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 3:
    parser.print_help()
    sys.exit(1)

run=args[0]
i1 = int(args[1])
i2 = int(args[2])
shapesim.shapesim.combine_trials(run, i1, i2)
