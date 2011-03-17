"""
    %prog run

Description:
    Create config and pbs files for the input run, one for each lens in the
    split.

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run = args[0]
lensing.config.write_config(run)
lensing.pbslens.write_pbs(run)
