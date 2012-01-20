"""
    %prog [options] run njob

Description
    Generate wq scripts for parallel collation
"""
import des
import sys
from sys import stdout


from optparse import OptionParser
parser=OptionParser(__doc__)

options, args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

run=args[0]
njob=int(args[1])

if run[0:2] == 'se':
    c = des.collate.SECollateWQJob(run,njob)
elif run[0:2] == 'me':
    raise ValueError("ME not yet implemented")
else:
    raise ValueError("Expected run 'me*' or 'se*'")

c.write()
