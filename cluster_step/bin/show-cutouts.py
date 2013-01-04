"""
    %prog [options] index1 index2 index3 ...

shearnums 1-8
psfnums 1-6
"""

import sys
import os
from numpy import zeros, sqrt

import cluster_step

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-p','--psfnum',default=None,
                  help='the psf number, required')
parser.add_option('-s','--shnum',default=None,
                  help='The shear number, required')
parser.add_option('-c','--ccd',default=None,
                  help='the ccd number, required')



def main():
    options,args = parser.parse_args(sys.argv[1:])

    if (len(args)==0
            or options.run is None 
            or options.psfnum is None
            or options.shnum is None
            or options.ccd is None):
        parser.print_help()
        sys.exit(1)

    run=options.run
    psfnum=int(options.psfnum)
    shnum=int(options.shnum)
    ccd=int(options.ccd)

    indices=[int(a) for a in args]

    pipe=cluster_step.pipe.Pipe(run=run,
                                psfnum=psfnum,
                                shnum=shnum,
                                ccd=ccd)
    pipe.show_many_cutouts(indices)

main()
