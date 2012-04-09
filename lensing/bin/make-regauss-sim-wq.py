"""
    %prog run
"""


import sys
import lensing

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('--by-s2',action='store_true')

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)


run=args[0]

if options.by_s2:
    lensing.regauss_sim.create_sim_wq_bys2(run)
else:
    lensing.regauss_sim.create_sim_wq_byellip(run)
