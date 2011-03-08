"""
    %prog proctype procrun sweeptype
"""
import sys
import es_sdsspy

from optparse import OptionParser
parser=OptionParser(__doc__)


options, args = parser.parse_args(sys.argv[1:])
if len(args) < 3:
    parser.print_help()
    sys.exit(45)

proctype=args[0]
procrun=args[1]
sweeptype=args[2]

p = es_sdsspy.sweeps.Proc(proctype,procrun,sweeptype)
# no need to crash
p.verify(nohalt=True)
