"""
    %prog [options] name version

name: e.g. rm for redmapper
version: e.g. dr8-v2
"""
import sys
import sdss_clusters

from optparse import OptionParser
parser=OptionParser(__doc__)
options, args = parser.parse_args(sys.argv[1:])
if len(args) < 2:
    parser.print_help()
    sys.exit(45)

name=args[0]
version=args[1]

mcs = sdss_clusters.select.Selector(name, version)
mcs.select()
mcs.write_columns()
