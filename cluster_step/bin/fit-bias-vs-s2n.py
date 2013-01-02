"""
    %prog [options]
"""

import sys
import os
from numpy import zeros

import cluster_step
from cluster_step.fitters import BiasFitter
from esutil.numpy_util import aprint

from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY, \
        PlotLabel, Table, PlotKey


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

parser.add_option('--s2',default=0.5,
                  help='restrict s2 less than this value, default %d')


def get_labels(fitter):
    psflab=fitter.get_psf_label()
    typelab=fitter.get_objtype_label()
    sizelab=fitter.get_size_label()

    return [psflab,typelab,sizelab]
    

def doplot(fitters, st, s2n_field):
    tab=Table(2,1)

    color1='blue'
    color2='red'

    m1pts=Points(st['s2n'],st['m1'],type='filled circle',color=color1)
    m1errpts=SymmetricErrorBarsY(st['s2n'],st['m1'],st['m1_err'],color=color1)
    m2pts=Points(st['s2n'],st['m2'],type='filled circle',color=color2)
    m2errpts=SymmetricErrorBarsY(st['s2n'],st['m2'],st['m2_err'],color=color2)

    c1pts=Points(st['s2n'],st['c1'],type='filled circle',color=color1)
    c1errpts=SymmetricErrorBarsY(st['s2n'],st['c1'],st['c1_err'],color=color1)
    c2pts=Points(st['s2n'],st['c2'],type='filled circle',color=color2)
    c2errpts=SymmetricErrorBarsY(st['s2n'],st['c2'],st['c2_err'],color=color2)

    m1pts.label=r'$\gamma_1$'
    m2pts.label=r'$\gamma_2$'

    key=PlotKey(0.9,0.9,[m1pts,m2pts], halign='right')

    mplt=FramedPlot()
    cplt=FramedPlot()

    mplt.xlabel=s2n_field+' threshold'
    mplt.ylabel='m'
    cplt.xlabel=s2n_field+' threshold'
    cplt.ylabel='c'

    zvals=0*st['s2n'].copy()
    zplt=Curve(st['s2n'], zvals)

    labels=get_labels(fitters[0])

    mplt.add( zplt, m1pts, m1errpts, m2pts, m2errpts, key )
    mplt.add(*labels)

    cplt.add( zplt, c1pts, c1errpts, c2pts, c2errpts )

    tab[0,0] = mplt
    tab[1,0] = cplt

    tab.show()

def get_stats(fitters):
    dt=[('s2n','f8'),
        ('m1','f8'),
        ('m1_err','f8'),
        ('c1','f8'),
        ('c1_err','f8'),
        ('m2','f8'),
        ('m2_err','f8'),
        ('c2','f8'),
        ('c2_err','f8')]
    st=zeros(len(fitters), dtype=dt)

    for i,bf in enumerate(fitters):
        st['s2n'][i] = bf.s2n_min
        st['m1'][i] = bf.g1fit.pars[0]
        st['m1_err'][i] = bf.g1fit.perr[0]
        st['c1'][i] = bf.g1fit.pars[1]
        st['c1_err'][i] = bf.g1fit.perr[1]

        st['m2'][i] = bf.g2fit.pars[0]
        st['m2_err'][i] = bf.g2fit.perr[0]
        st['c2'][i] = bf.g2fit.pars[1]
        st['c2_err'][i] = bf.g2fit.perr[1]

    return st

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if options.run is None or len(args) < 1:
        parser.print_help()
        sys.exit(1)

    s2n_vals=[float(s2n) for s2n in args]

    fitters=[]
    for i,s2n in enumerate(s2n_vals):
        bf=BiasFitter(options.run,
                      psfnums=options.psfnums,
                      objtype=options.type,
                      s2n_field=options.field,
                      s2n_min=s2n,
                      s2_max=float(options.s2))
        fitters.append(bf)

    st = get_stats(fitters)
    aprint(st,fancy=True)
    doplot(fitters, st, options.field)
    

main()
