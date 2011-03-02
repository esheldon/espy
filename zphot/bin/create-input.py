"""
    %prog [options] sample

Description:

    Create input samples for photoz codes.  You need to create the config file
    first in ${ESPY_DIR}/zphot/config/zinput-{sample}.json.  

    For weighting you will want run to make the overall file, and run again to
    split into chunks, e.g.  -n 50

"""
import sys
import zphot

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("-n","--nchunk",default=None,
                  help="Split into N chunks.")

options,args = parser.parse_args(sys.argv[1:])

if len(args) < 1:
    parser.print_help()
    sys.exit(45)

sample = args[0]
nchunk = options.nchunk
if nchunk is not None:
    nchunk = int(nchunk)


zcs = zphot.select.ColumnSelector(sample, nchunk=nchunk)

zcs.select()
zcs.write()

