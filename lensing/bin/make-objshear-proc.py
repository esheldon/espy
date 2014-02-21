"""
    %prog [options] run version

Description:

    Create config and shear/reduce wq scripts

    version is the code version, e.g. gsens, im3

    types is by default 
        config,shear

    types can be
    ------------
    config:
        Write the config file for sobjshear
    shear:
        Write the wq script for shear for each lens/src split
    src_reduce:
        Write the set of wq scripts for reducing across sources
    collate:
        Write scripts to collate each of the lens splits

    # these are optional
    lens_concat:
        A script to concatenate all the lens splits together
        into a final reduced file.  useful if not collating
        in the lens splits

    all_reduce:
        Write the wq script to reduce all at once and produce
        the reduced file

If you don't have any lens splits you can just use the all_reduce script.

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="types",default="config,shear,src_reduce,collate",
                  help="types to make.  Default is %default")
parser.add_option("-g",dest="groups",default=None,
                  help="machine groups to use.  Default is %default")
parser.add_option("-n",dest="notgroups",default=None,
                  help="machine notgroups to use.  Default is %default")
parser.add_option("-p",dest="priority",default=None,
                  help="priority use.  Default is %default")
parser.add_option("--fs",default='nfs',
                  help="which file system.  Default is '%default'")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run = args[0]
version=args[1]
types=options.types.split(',')


if 'config' in types:
    lensing.objshear_config.write_configs(run)

wql=lensing.wqsubmit.WQLens(run, version,
                            groups=options.groups,
                            notgroups=options.notgroups, 
                            priority=options.priority,
                            fs=options.fs)


if 'shear'in types:
    wql.write_shear_scripts()


if 'src_reduce' in types:
    wql.write_src_reduce_scripts()

if 'collate' in types:
    wql.write_lens_collate_scripts()


if 'lens_concat' in types:
    wql.write_lens_concat_script()



if 'all_reduce' in types:
    wql.write_all_reduce_script()
