"""
    %prog [options] output_file
"""
import sys, os
import cosmos

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--obj-range',None,
                  help="csv inclusive range, e.g. 35,66")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    output_file=args[0]
    obj_range=options.obj_range
    if obj_range is not None:
        obj_range=[int(x) for x in obj_range.split(',')]

    models=['exp','dev']
    cosmos.gmix_cosmos.fit_cosmos(models,
                                  output_file,
                                  obj_range=obj_range)

main()
