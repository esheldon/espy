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
parser.add_option("-s",'--show',action="store_true",default=False,
                  help="Show plot on screen.  Default %default")
parser.add_option("-t",'--type',default='corrected',
                  help="Should be binned, corrected, corrected-ssh, jackknife.  Default %default")
parser.add_option("-o","--osig",action="store_true",
                  help="Make plots of osig")
parser.add_option("--compare-osig",action="store_true",
                  help="Make plots of dsig compared to osig")
parser.add_option("--run2",default=None, help="Overplot a second run")

parser.add_option("-l",'--linear',action="store_true",default=False,
                  help="only a linear plot")

parser.add_option('-y','--yrange',default=None,
                  help="yrange for plot")
parser.add_option('-x','--xrange',default=None,
                  help="xrange for plot")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 3:
    parser.print_help()
    sys.exit(1)

run = args[0]
bintype = args[1]
nbin = int(args[2])

yrange=options.yrange
if yrange is not None:
    yrange=[float(y) for y in yrange.split(',')]
xrng=options.xrange
if xrng is not None:
    xrng=[float(x) for x in xrng.split(',')]



b = lensing.binning.instantiate_binner(bintype, nbin)
if options.compare_osig:
    b.plot_dsig_osig_byrun(run, options.type, show=options.show, range4var=[0.5,100.0],
                           linear=options.linear, yrange=yrange, xrange=xrng)
elif options.osig:
    b.plot_osig_byrun_1var(run, options.type, show=options.show)
else:
    if options.run2 is not None:
        b.plot_dsig_2runs(run, options.run2, options.type, show=options.show)
    else:
        b.plot_dsig_byrun_1var(run, options.type, show=options.show)
