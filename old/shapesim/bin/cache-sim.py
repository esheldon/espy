"""
    %prog simname is2 ie [num]

Send num if not a ring test
"""

import sys
import os
from esutil.misc import wlog
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("--start-index",default=None,
                  help="index to start in series")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)

    simname =args[0]
    is2 = int(args[1])
    ie = int(args[2])

    start_index = options.start_index
    if start_index is not None:
        start_index=int(start_index)


    ss=shapesim.ShapeSim(simname)
    orient=ss.get('orient','rand')
    if orient == 'ring':
        num = ss['nring']
        print >>sys.stderr,'making',num,'in ring'
    elif orient == 'ring-rand' or orient=='rand':
        if len(args) < 4:
            parser.print_help()
            sys.exit(45)
        num=int(args[3])

    if start_index is not None:
        trials=range(start_index, num)
    else:
        trials=range(num)

    for i in trials:
        print >>sys.stderr,'-'*70
        print >>sys.stderr,'%d/%d' % (i+1,num)

        if orient == 'ring-rand':
            ss.write_ring_trial(is2, ie)
        else:
            if orient == 'ring':
                itheta=i
            else:
                itheta=None
            ss.write_trial(is2, ie, itheta=itheta)

main()
