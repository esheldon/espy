"""
    %prog train_sample

Description:
    Convert the matched training sets to inputs usable for the weighting
    method.
"""
import sys
import zphot
from optparse import OptionParser

parser=OptionParser(__doc__)

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

train_sample = args[0]
wt = zphot.weighting.WeightedTraining(train_sample)
wt.convert_training()
