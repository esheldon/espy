"""
    %prog pzrun

Description:

    Create the official release files."""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)

pzrun = args[0]

wo=zphot.weighting.WeightedOutputs()
wo.make_release(pzrun)
