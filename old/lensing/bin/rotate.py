"""
    %prog [options] procrun sweeptype

Description:

    rotate ellipticities
"""

import sys
import lensing
import numpy

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option("-d","--detrend", action="store_true", 
                  help="rotate the detrend ellipticities.  "
                  "In this case filters defaults to r,i and "
                  "you can send rmag_max")
parser.add_option("--rmag-max", default=21.8, 
                  help="Max mag in r for detrended shapes, default %default")
parser.add_option("-s","--system", default='eq', help="system")
parser.add_option("-f","--filters", default='u,g,r,i,z', 
                  help="filters to rotate")


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    procrun = args[0]
    sweeptype = args[1]
    filters = options.filters.split(',')
    rmag_max=float(options.rmag_max)

    if options.detrend:
        filters = ['r','i']

    c=lensing.regauss.Collator(procrun,sweeptype)
    c.add_rotated_e1e2(filters=filters, system=options.system, 
                       detrend=options.detrend, rmag_max=rmag_max)

main()
