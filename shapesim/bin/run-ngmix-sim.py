"""
    %prog run s2n npairs output_file
"""

import sys
from sys import stderr
import shapesim.ngmix_sim
import fitsio

from optparse import OptionParser
parser=OptionParser(__doc__)

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    s2n=float(args[1])
    npairs=int(args[2])
    output_file=args[3]

    sim=shapesim.ngmix_sim.NGMixSim(run, s2n, npairs)

    sim.run_sim()

    data=sim.get_data()

    print >>stderr,'writing:',output_file
    with fitsio.FITS(output_file,'rw',clobber=True) as fobj:
        fobj.write(data)



main()
