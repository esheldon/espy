"""
    %prog [options] procrun sweeptype

Description:

    Plot the detrend of mean ellipticity vs the R.  You first generate the
    detrending polynomial using rg-diagnostics.py
    
    procrun should be a sweeps regauss run or princeton
"""

import sys
import lensing
import numpy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-f","--filter", default=None, help="the band for rg runs")
parser.add_option("--rmag-min", default=18.0, help="Min mag in r, default %default")
parser.add_option("--rmag-max", default=21.8, help="Max mag in r, default %default")
parser.add_option("-s","--show", action='store_true', help="show plots on screen")


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    procrun = args[0]
    sweeptype = args[1]

    filter = options.filter
    rmag_min=options.rmag_min
    rmag_max=options.rmag_max

    rmag_min=float(rmag_min)
    rmag_max=float(rmag_max)


    t = lensing.regauss_test.Tester(procrun, sweeptype, filter)
    t.detrend(rmag_min, rmag_max, show=options.show)

main()
