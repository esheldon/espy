"""
    %prog pzrun

Description:
    Sum up the p(z) for all recoverable objects in the given run.

"""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

pzrun = args[0]

wo = zphot.weighting.WeightedOutputs()
wo.make_pofz_hist(pzrun)
wo.plot_pofz_hist(pzrun)
