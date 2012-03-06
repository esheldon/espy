"""
    %prog [options] run bintype nbin

Description:

    Bin lenses.  Type is the type of binning, e.g. n200 or mz,
    run is the lensing run, e.g. 07, and nbin is the number of bins

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-s",dest="show",action="store_true",default=False,
                  help="Show plot on screen.  Default %default")
parser.add_option("-t",dest="type",default='corrected',
                  help="Should be binned, corrected, jackknife.  Default %default")
parser.add_option("-o","--osig",action="store_true",
                  help="Make plots of osig")
parser.add_option("--compare-osig",action="store_true",
                  help="Make plots of dsig compared to osig")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 3:
    parser.print_help()
    sys.exit(1)

run = args[0]
bintype = args[1]
nbin = int(args[2])

b = lensing.binning.instantiate_binner(bintype, nbin)
if options.compare_osig:
    b.plot_dsig_osig_byrun(run, options.type, show=options.show, range4var=[0.5,100.0])
elif options.osig:
    b.plot_osig_byrun_1var(run, options.type, show=options.show)
else:
    b.plot_dsig_byrun_1var(run, options.type, show=options.show)
