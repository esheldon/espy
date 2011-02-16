"""
    %prog make-objshear-scat.py sample

Description:
    Create an input source catalog for objshear.  The sample
    id implies a config file in ${ESPY_DIR}/lensing/config. 
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
wt = zphot.weighting.WeightedTraining(train_sample)
wt.convert_training()
