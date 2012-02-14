"""
    %prog [options] run

Description:

    Create config, wq submit and reduce script

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="types",default="config,script,condor",
                  help="types to make.  Default is config,wq")
parser.add_option("-g",dest="groups",default="new,new2",
                  help="machine groups to use.  Default is %default")
parser.add_option("-p",dest="priority",default="med",
                  help="priority use.  Default is %default")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run = args[0]
types=options.types.split(',')

if 'config' in types:
    lensing.objshear_config.write_configs(run)
if 'wq' in types:
    wql=lensing.wqsubmit.WQLens(run,options.groups,options.priority)
    wql.write_scripts()
