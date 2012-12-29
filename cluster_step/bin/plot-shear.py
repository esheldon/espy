"""
    %prog [options] run shearnum

shearnums 1-8
psfnums 1-6
"""

import sys
import os

import cluster_step
from cluster_step import files, stats

import esutil as eu
from esutil.numpy_util import aprint
from biggles import FramedArray, Points, \
        SymmetricErrorBarsX, SymmetricErrorBarsY

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--psfnums',default=None,
                  help='restrict to these PSFs, comma separated')
parser.add_option('-f','--field',default='s2n_w',
                  help="bin by this field, default s2n_w")
parser.add_option('-n','--nperbin',default=1000,
                  help="bin by this field, default s2n_w")
parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")



def write_eps(arr):
    raise RuntimeError("implement")
    #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

def doplot(bindata, bin_field, show=False):
    arr=FramedArray(2,1)

    arr.uniform_limits=1
    arr.xlog=True
    arr.xrange=[0.5*bindata[bin_field].min(), 1.5*bindata[bin_field].max()]
    arr.xlabel=bin_field

    xdata=bindata[bin_field]
    xerr=bindata[bin_field+'_err']

    xerrpts1 = SymmetricErrorBarsX(xdata, bindata['g1'], xerr)
    xerrpts2 = SymmetricErrorBarsX(xdata, bindata['g2'], xerr)

    g1pts = Points(xdata, bindata['g1'])
    g1errpts = SymmetricErrorBarsY(xdata, bindata['g1'], bindata['g1_err'])
    g2pts = Points(xdata, bindata['g2'])
    g2errpts = SymmetricErrorBarsY(xdata, bindata['g2'], bindata['g2_err'])

    arr[0,0].add( xerrpts1, g1pts, g1errpts )
    arr[1,0].add( xerrpts2, g2pts, g2errpts )

    if show:
        arr.show()

    return arr

    
def get_psfnums(psfnum):
    if psfnum is None:
        psfnums=range(6)
    else:
        psfnums=[int(p) for p in psfnum.split(',')]

    return psfnums

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    shnum=int(args[1])
    nperbin=int(options.nperbin)
    objtype=options.type
    if objtype:
        print 'selecting type:',objtype

    psfnums=get_psfnums(options.psfnums)
    bin_field=options.field

    data=files.read_output_set(run, psfnums, shnum, objtype=objtype)
    bindata=stats.bin_data(data, bin_field, nperbin)

    aprint(bindata, header=True, page=False, fancy=True)

    doplot(bindata, bin_field, show=options.show)


main()
