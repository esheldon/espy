"""
    %prog [options] run

Description:

    Create config, script, and condor files for the input run, one for each
    split of the sources.

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

lensing.condor.write_submit_scripts(run)
lensing.scripts.write_scripts(run)
lensing.config.write_configs(run)
