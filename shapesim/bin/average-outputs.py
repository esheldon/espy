"""
    %prog run

Description

    Average the trials

"""

import sys
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('--skip1',default=None,
                  help="elements in index 1 to skip")
parser.add_option('--skip2',default=None,
                  help="elements in index 2 to skip")
parser.add_option('--nocum',action='store_true',
                  help="don't make the cumulative sums")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]

skip1=options.skip1
if skip1 is None:
    skip1=[]
else:
    skip1 = [int(v) for v in skip1.split(',')]
skip2=options.skip2
if skip2 is None:
    skip2=[]
else:
    skip2 = [int(v) for v in skip2.split(',')]

if options.nocum:
    docum=False
else:
    docum=True
shapesim.shapesim.make_averaged_outputs(run, docum=docum,
                                        skip1=skip1, skip2=skip2)
