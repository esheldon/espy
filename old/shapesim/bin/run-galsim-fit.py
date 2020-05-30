"""
    %prog run profile ellip resolution

"""

import sys
import os
import shapesim

from optparse import OptionParser
parser=OptionParser(__doc__)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(1)

    run=args[0]
    profile=args[1]
    ellip=args[2]
    res=args[3]

    g=shapesim.gmix_fit_sim.GMixGalSim(run,profile,ellip,res)
    g.run()

main()
