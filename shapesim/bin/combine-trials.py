"""
    %prog run [i1 i2]

i1 is usually is2, i2 can be is2n or ie

Either send just the run or both i1 and i2
"""

import sys
import shapesim
from optparse import OptionParser

parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])


run=args[0]

if len(args) == 3:
    i1 = int(args[1])
    i2 = int(args[2])
    shapesim.shapesim.combine_trials(run, i1, i2)
elif len(args) == 1:
    c = shapesim.read_config(run)
    cs = shapesim.read_config(c['sim'])
    runtype = c['runtype']

    n1 = cs['nums2']
    if runtype == 'byellip':
        n2 = cs['nume']
    else:
        n2 = shapesim.get_nums2n(c)


    for i1 in xrange(n1):
        for i2 in xrange(n2):
            shapesim.shapesim.combine_trials(run, i1, i2)

else:
    if len(args) < 3:
        parser.print_help()
        sys.exit(1)
