"""
    %prog run is2 ie
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

    run=args[0]
    is2 = int(args[1])
    ie = int(args[2])

    if run[0:5] == 'deswl':
        sim=shapesim.deswl_sim.DESWLSim(run)
    elif run[0:4] == 'gmix':
        sim=shapesim.gmix_em_sim.GMixEMSim(run)
    else:
        raise ValueError("Don't know about run '%s'" % run)

    sim.process_trials(is2, ie)


main()
