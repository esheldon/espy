"""
    %prog [options] run

Description:
    Create config files for the input run, one for each source in the
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
