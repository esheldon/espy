"""
    %prog [options] run

Description
    Make size-mag plots for the given run.
"""
import des
import sys
from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-t","--type",default='png',
                  help="plot type: eps or png, default %default")
parser.add_option("-e","--expname",default=None,
                  help="only do a plot for the specified exposure")
parser.add_option("-c","--combine",action="store_true",
                  help="Combine ccds in an exposure")
parser.add_option("-a","--all",action="store_true",
                  help="Combine data from all exposures")

def main():
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    run = args[0]
    if options.all:
        s = des.select.SizeMagSelector(run)
        s.plot()
    else:
        des.checksg.plot_sg_all_exposures(run, 
                                          combine=options.combine, 
                                          ptype=options.type,
                                          expname=options.expname)

main()
