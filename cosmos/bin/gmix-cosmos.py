"""
    %prog [options] config_file output_file
"""
import sys, os
import cosmos
import yaml

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--obj-range',default=None,
                  help="csv inclusive range, e.g. 35,66")
parser.add_option('--make-plots',action='store_true',
                  help="show some diagnostic plots")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    config_file=args[0]
    output_file=args[1]

    conf=yaml.load(open(config_file))

    obj_range=options.obj_range
    if obj_range is not None:
        obj_range=[int(x) for x in obj_range.split(',')]

    models=conf['models']
    del conf['models']

    conf['obj_range']=obj_range
    conf['make_plots']=options.make_plots

    cosmos.gmix_cosmos.fit_cosmos(models,
                                  output_file,
                                  **conf)

main()
