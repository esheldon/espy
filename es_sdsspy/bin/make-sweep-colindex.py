"""
    %prog [options] type
"""
import sys
import sdsspy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-p","--primary",
                  default=False,
                  action="store_true",
                  help="Restrict to primary objects")

options, args = parser.parse_args(sys.argv[1:])
if len(args) < 1:
    parser.print_help()
    sys.exit(45)

type = args[0]
primary = options.primary

collator = es_sdsspy.sweeps_collate.Collator(type, primary=primary)
collator.create_indices()


