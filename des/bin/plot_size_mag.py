"""
    %prog [options] run expname

Description

    Make size-mag plots for the given run.  If the collated columns database
    exists it will be used, otherwise data are read from the flat files.

"""
import des
import sys
from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-t","--type",default='png',
                  help="plot type: eps,png,screen default %default")
parser.add_option("-c","--ccd",default=None,
                  help="only do a plot for the specified exposure,ccd")
parser.add_option("--nocombine",action="store_true",
                  help="Don't combine ccds in an exposure")
parser.add_option("-a","--all",action="store_true",
                  help="Combine data from all exposures")
parser.add_option("--fs",default='nfs',
                  help="file system, 'nfs' or 'hdfs'")


def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    run = args[0]
    expname = args[0]
    combine=True
    if options.nocombine:
        combine=False

    if options.all:
        s = des.select.SizeMagSelector(run)
        s.plot()
    else:
        s = des.checksg.SizeMagPlotter(run, ptype=options.type, fs=options.fs)
        if options.ccd is not None:
            s.plot_ccd(expname,int(options.ccd))
        else:
            s.plot_exposure(expname,combine=combine)

main()
