
"""
    %prog [options] lensrun randrun

Description:

    Correct the delta sigma for 

        - subtract the signal from randoms
        - 1/ssh  the shear polarizability correction
        - clustering of sources with lenses
    Note the subtraction of randoms should happen before
    the other two"""

import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-t",dest="bintype",default=None,
                  help="The type of binning, default %default")
parser.add_option("-n",dest="nbin",default=None,
                  help="The number of bins, default %default")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    lensrun=args[0]
    randrun=args[1]

    bintype=options.bintype
    nbin=int(options.nbin)

    if bintype is None or nbin is None:
        raise ValueError("currently demand some kind of binning")

    b = instantiate_binner(bintype, nbin)

    # read collated lens catalog and select lenses with the
    # bin criteria.  Then match the randoms redshift histogram
    data = lensing.files.sample_read('collated', lensrun)
    rand = lensing.files.sample_read('collated', randrun)

    for binnum in xrange(nbin):
        ii = b.select_bin(data, binnum)

