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

    zmin=0.0
    binsize=0.031429
    pzstruct = zphot.weighting.read_pofz(pzrun, 'rand', with_rmag=True)
    zphot.weighting.plot_6rand_pofz(pzstruct, zmin, binsize, seed=25)

main()
