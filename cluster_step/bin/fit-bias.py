"""
    %prog [options]
"""

import sys
import os
import cluster_step
from cluster_step.fitters import BiasFitter

from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-p','--psfnums',default=None,
                  help='restrict to these PSFs, comma separated')


parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-f','--field',default='s2n_w',
                  help="field for S/N, default %default")

parser.add_option('--s2n',default=20, help=("threshold in s/n"))

parser.add_option('--s2',default=0.5,
                  help='restrict s2 less than this value, default %d')


parser.add_option('--show',action='store_true',
                  help="show the plot on the screen")

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if options.run is None:
        parser.print_help()
        sys.exit(1)

    doshow  = options.show

    bf=BiasFitter(options.run,
                  psfnums=options.psfnums,
                  objtype=options.type,
                  s2n_field=options.field,
                  s2n_min=float(options.s2n),
                  s2_max=float(options.s2))

    if doshow:
        bf.show()

main()
