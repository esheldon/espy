"""
    %prog sample

Description:

    Create mean inverse critical density as a function of lens redshift
    for regauss objects with criteria for the input sample

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

sample = args[0]

c=lensing.scat.instantiate_sample(sample)
c.add_scinv()

