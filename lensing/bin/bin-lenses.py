"""
    %prog [options] run type nbin

Description:

    Bin lenses.  
    
    run is, e.g., rm05gmix01
    
    type is the type of binning, e.g. lambda, 
    and nbin is the number of bins
"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])


if len(args) < 3:
    parser.print_help()
    sys.exit(1)

run = args[0]
type = args[1]
nbin = int(args[2])

lensing.binning.bin_lenses_byrun(run, type, nbin)
