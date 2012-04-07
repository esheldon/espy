"""
    %prog run [s2]

Description:

    Run some sims with the label "run" and the given models"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-v','--verbose',action='store_true')

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]

if len(args) > 1:
    s2=float(args[1])
    rs = lensing.regauss_sim.RegaussSimulatorRescontrol(run, s2, 
                                                        verbose=options.verbose)
    rs.run_many_ellip()
else:
    lensing.regauss_sim.run_many_s2(run, verbose=options.verbose)
