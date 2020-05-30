"""
    %prog run is2 ie/is2n

Check for the existence of outputs
"""

import sys
import os
import shapesim

from sh import hadoop, awk

from optparse import OptionParser
parser=OptionParser(__doc__)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 3:
        parser.print_help()
        sys.exit(45)


    run=args[0]
    is2 = int(args[1])
    ie_or_is2n = int(args[2])

    conf=shapesim.read_config(run)
    simconf = shapesim.read_config(conf['sim'])

    pattern=shapesim.get_output_url(run, is2, ie_or_is2n, itrial='*', fs='hdfs')

    flist=awk(hadoop('fs','-ls',pattern), '{print $8}').split()

    nring = simconf['nring']
    for i in xrange(nring):
        f=shapesim.get_output_url(run, is2, ie_or_is2n, itrial=i, fs='hdfs')
        f=f.replace('hdfs://','')
        if f not in flist:
            print f

main()
