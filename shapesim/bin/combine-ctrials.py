"""
    %prog run is2n

"""

import sys
import shapesim
from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('--allow-missing',action='store_true',
                  help="allow missing splits")
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

run=args[0]
is2n=int(args[1])

shapesim.shapesim.combine_ctrials(run, is2n,
                                  allow_missing=options.allow_missing)

