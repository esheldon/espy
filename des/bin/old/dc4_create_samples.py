"""
    %prog [options] serun objclass 
"""
import des
import sys
from sys import stdout

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-d",
        dest="delim",
        default=None,
        help="delimiter for file")
parser.add_option("-f",
        dest="ftype",
        default='rec',
        help="file type ['rec','fits']")


options, args = parser.parse_args(sys.argv[1:])
if len(args) < 2:
    parser.print_help()
    sys.exit(45)

serun=args[0]
objclass=args[1]
delim=options.delim
ftype=options.ftype


indir='~/data/DES/wlbnl/%s/collated' % serun
outdir='~/data/DES/wlbnl/%s/collated' % serun
des.util.dc4_create_shear_samples(serun, objclass, ftype=ftype, delim=delim, 
                                  indir=indir, outdir=outdir)
