"""
    %prog [options] type run

Description:
    Create yaml file(s).  Type can be [weights pofz]
"""
import sys
import zphot
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-n","--nchunk",default=None,
                  help="Number of chunks.  Must be present for pofz")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 2:
    parser.print_help()
    sys.exit(1)

type = args[0]
run = args[1]

if type == 'weights':
    zphot.weighting.create_weights_wq(run)
elif type == 'pofz':
    nchunk = options.nchunk
    if nchunk is None:
        raise ValueError("You must send nchunk for type=pofz")
    nchunk = int(nchunk)
    zphot.weighting.create_pofz_wq(run,nchunk)
else:
    raise ValueError("type must be 'weights' or 'pofz'")
