"""
    %prog [options]
"""

import sys
import os
from numpy import where
import cluster_step
from cluster_step import files
from esutil.numpy_util import strmatch, combine_arrlist

from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')

parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-f','--field',default='s2n_w',
                  help="field for S/N, default %default")

parser.add_option('--s2n',default=20, help=("threshold in s/n"))

parser.add_option('--s2',default=1.0,
                  help='restrict s2 less than this value, default %d')


def read_single(run,psfnum,shnum,ccd):
    fname=files.get_output_path(run=run, 
                                psfnum=psfnum, 
                                shnum=shnum, 
                                ccd=ccd, 
                                ftype='shear')
    data0=None
    if os.path.exists(fname):
        data0=files.read_fits_output(run=run, 
                                     psfnum=psfnum, 
                                     shnum=shnum, 
                                     ccd=ccd, 
                                     ftype='shear',
                                     verbose=False)
    else:
        print 'missing:',fname
    return data0

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if options.run is None:
        parser.print_help()
        sys.exit(1)

    s2_max=float(options.s2)
    s2n_min=float(options.s2n)
    s2n_field=options.field
    objtype=options.type

    for psfnum in files.PSFNUMS:
        for shnum in files.SHNUMS:
            for ccd in files.CCDS:
                data0=read_single(run,psfnum,shnum,ccd)
                if data0 is not None:

                    logic= (data0['s2'] < s2_max)
                    logic = logic & (data0[s2n_field] > s2n_min)
                    if objtype is not None:
                        logic=logic & strmatch(data0['model'],objtype)

                    w,=where(logic)
                    frac=float(w.size)/data0.size
                    print '%s/%s %.2f' % (w.size,data0.size,frac)
                    if w.size > 0:
                        pass

main()
