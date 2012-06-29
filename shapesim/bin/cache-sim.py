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


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)

    simname =args[0]
    is2 = int(args[1])
    ie = int(args[2])

    ss=shapesim.ShapeSim(simname)
    orient=ss.get('orient','rand')
    if orient == 'ring':
        num = ss['nring']
        print >>sys.stderr,'making',num,'in ring'
    else:
        if len(args) < 4:
            parser.print_help()
            sys.exit(45)
        num=int(args[3])


    for i in xrange(num):
        print >>sys.stderr,'-'*70
        print >>sys.stderr,'%d/%d' % (i+1,num)

        if orient == 'ring':
            itheta=i
        else:
            itheta=None
        ss.write_trial(is2, ie, itheta=itheta)

main()
