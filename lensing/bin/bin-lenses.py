"""
    %prog [options] run type nbin
    or
    %prog [options] run bintype

Description:

    %prog rm05gmix01 lambda 10 
    %prog rm05gmix01 lambda04-z01-e01
"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 2:
    parser.print_help()
    sys.exit(1)

run = args[0]
type = args[1]

if len(args) > 2:
    nbin = int(args[2])
    lensing.binning.bin_lenses_byrun(run, type, nbin=nbin)
else:
    print 'using full bin specification'
    lensing.binning.bin_lenses_byrun(run, type)
