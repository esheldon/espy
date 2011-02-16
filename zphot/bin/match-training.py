"""
    %prog train_sample
"""
import sys
import zphot
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("--no-photo-cuts",default=None,action='store_true',
                  help="process the new 'other' files")

options,args = parser.parse_args(sys.argv[1:])
no_photo_cuts = options.no_photo_cuts

if len(args) < 1:
    parser.print_help()
    sys.exit(1)

train_sample = args[0]
zt = zphot.training.Training(train_sample,no_photo_cuts=no_photo_cuts)
zt.match()

zt.plot_all_radec()
