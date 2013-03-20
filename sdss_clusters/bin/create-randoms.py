"""
    %prog [options] name version nsamples nper

name: e.g. rm for redmapper
version: e.g. dr8-v3
"""
import sys
import sdss_clusters

from optparse import OptionParser
parser=OptionParser(__doc__)
options, args = parser.parse_args(sys.argv[1:])

if len(args) < 4:
    parser.print_help()
    sys.exit(45)

name=args[0]
version=args[1]
nsamples=int(args[2])
nper=int(args[3])

r = sdss_clusters.random.Random(name, version, nsamples, nper)
r.write_randoms()
