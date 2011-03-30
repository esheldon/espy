"""
    %prog [options] run

Description:
    Create config and pbs files for the input run, one for each lens in the
    split.

"""
import sys
import lensing
from optparse import OptionParser

parser=OptionParser(__doc__)
parser.add_option("-n","--nthreads", help="max allowed openmp threads")

options,args = parser.parse_args(sys.argv[1:])


if len(args) < 1:
    parser.print_help()
    sys.exit(1)

run = args[0]
nthreads = int(options.nthreads)

lensing.config.write_config(run)
lensing.pbslens.write_pbs(run, nthreads=nthreads)
