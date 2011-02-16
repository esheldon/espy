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

options, args = parser.parse_args(sys.argv[1:])
if len(args) < 1:
    parser.print_help()
    sys.exit(45)

type = args[0]
primary = options.primary

#make_indices = options.make_indices
sel = sdsspy.sweeps.Selector(type,primary=primary)
sel.process_all()

