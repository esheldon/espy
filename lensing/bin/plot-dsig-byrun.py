"""
    %prog [options] run type nbin

Description:

    Bin lenses.  Type is the type of binning, e.g. n200 or mz,
    run is the lensing run, e.g. 07, and nbin is the number of bins

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-p",dest="dops",action="store_true",default=False,
                  help="Write a hardcopy encapsulated postscript file.  Default %default")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 3:
    parser.print_help()
    sys.exit(1)

run = args[0]
type = args[1]
nbin = int(args[2])

b = lensing.binning.instantiate_binner(type, nbin)
b.plot_dsig_byrun_1var(run, dops=options.dops)
