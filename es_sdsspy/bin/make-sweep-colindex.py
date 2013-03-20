"""
    %prog [options] type

type should be 'gal','primgal','star','primstar'
"""
import sys
import es_sdsspy

from optparse import OptionParser
parser=OptionParser(__doc__)


options, args = parser.parse_args(sys.argv[1:])
if len(args) < 1:
    parser.print_help()
    sys.exit(45)

type = args[0]

collator = es_sdsspy.sweeps_collate.Collator(type)
collator.create_indices()


