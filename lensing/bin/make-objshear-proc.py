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
parser.add_option("-t",dest="types",default="config,script,condor",
                  help="types to make.  Default is config,script,condor")
parser.add_option("--manager",dest="manager",default="modules",
                  help="How to manage software.  Default is %default")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run = args[0]
types=options.types.split(',')

if 'config' in types:
    lensing.objshear_config.write_configs(run)
if 'script' in types:
    lensing.scripts.write_scripts(run, manage=options.manager)
if 'condor' in types:
    lensing.condor.write_submit_scripts(run)
