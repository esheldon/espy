"""
    %prog run objmodel psfmodel

Description:

    Run some sims with the label "run" and the given models"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 3:
    parser.print_help()
    sys.exit(45)


run=args[0]
objmodel=args[1]
psfmodel=args[2]

if len(args) > 3:
    s2=args[3]
    rs = lensing.regauss_sim.RegaussSimulatorHirez(run, s2, objmodel, psfmodel)
    rs.run_many_ellip()
else:
    lensing.regauss_sim.run_many_s2(run, objmodel, psfmodel)
