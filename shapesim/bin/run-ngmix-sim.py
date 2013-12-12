"""
    %prog run s2n npairs output_file
"""

import sys
import os
from sys import stderr
import shapesim.ngmix_sim
import fitsio

from optparse import OptionParser
parser=OptionParser(__doc__)

def get_checkpoint_data(output_file):
    """
    Read in checkpoint data if it exists
    """
    checkpoint_file=output_file.replace('.fits','-checkpoint.fits')
    data=None

    if os.path.exists(checkpoint_file):
        print >>stderr,'reading checkpoint data:',checkpoint_file
        data=fitsio.read(checkpoint_file)

    return checkpoint_file, data

def cleanup_checkpoint(checkpoint_file):
    """
    if we get this far, we have succeeded in writing the data. We can remove
    the checkpoint file
    """
    if os.path.exists(checkpoint_file):
        print >>stderr,'removing checkpoint file',checkpoint_file
        os.remove(checkpoint_file)


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 4:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    s2n=float(args[1])
    npairs=int(args[2])
    output_file=args[3]

    checkpoint_file, checkpoint_data=get_checkpoint_data(output_file)
    sim=shapesim.ngmix_sim.NGMixSim(run, s2n, npairs,
                                    checkpoint_file=checkpoint_file,
                                    checkpoint_data=checkpoint_data)

    sim.run_sim()

    data=sim.get_data()

    print >>stderr,'writing:',output_file
    with fitsio.FITS(output_file,'rw',clobber=True) as fobj:
        fobj.write(data)

    cleanup_checkpoint(checkpoint_file)

main()
