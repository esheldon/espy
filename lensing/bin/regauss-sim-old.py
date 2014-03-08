"""
    %prog run

Description:

    Run some sims with the label "run".

    Note the ellip option only has effect if s2 is also sent
"""

import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-v','--verbose',action='store_true')
parser.add_option('--s2',default=None,help="use the specified s2")
parser.add_option('-e','--ellip',default=None,help="use the specified ellipticity")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]
s2=options.s2
if s2 is not None:
    s2=float(s2)
ellip=options.ellip
if ellip is not None:
    ellip=float(ellip)

if s2 is not None:
    rs = lensing.regauss_sim.RegaussSimulatorRescontrol(run, s2, 
                                                        verbose=options.verbose)
    if ellip is not None:
        rs.run_ellip(ellip)
    else:
        rs.run_many_ellip()
else:
    lensing.regauss_sim.run_many_s2(run, verbose=options.verbose)
