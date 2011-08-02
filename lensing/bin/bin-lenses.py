"""
    %prog [options] type run nbin

Description:

    Bin lenses.  Type is the type of binning, e.g. n200 or mz,
    run is the lensing run, e.g. 07, and nbin is the number of bins

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
#parser.add_option("-t",dest="types",default="config,script,condor",
#                  help="types to make.  Default is config,script,condor")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 3:
    parser.print_help()
    sys.exit(1)

type = args[0]
run = args[1]
nbin = int(args[2])

lensing.binning.bin_lenses_byrun(type, run, nbin)
