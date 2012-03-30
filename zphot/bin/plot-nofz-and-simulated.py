"""
    %prog pzrun

Make a plot with the simulated N(z) overplotted
"""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    pzrun = args[0]
    
    zphot.weighting.combine_nofz_and_simulated(pzrun)

main()
