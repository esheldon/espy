"""
    %prog [options] [s2n1 s2n2 s2n3] .... 

Either

    - Send -n nbin, which will result in logarithmic binning between the min
    s/n and 1000.
    - Send the bin edges as args
        Bins will be [s2n1,s2n2], [s2n2,s2n3], ...
"""

import sys
import os
from numpy import zeros, logspace, log10, array

import cluster_step
from cluster_step.fitters import BiasFitter
from esutil.numpy_util import aprint

from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY, \
        PlotLabel, Table, PlotKey, FillBetween


from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-p','--psfnums',default=None,
                  help='restrict to these PSFs, comma separated')

parser.add_option('-n','--nbin',default=None,
                  help='number of bins for logarithmic binning in s/n')
parser.add_option('--s2n',default='10,500',
                  help='s/n range when log binning, default %default')


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
    runlab=fitter.get_run_label()

    return [psflab,typelab,sizelab, runlab]
    

def get_symmetric_range(data1, err1, data2, err2):
    minval1=(data1-err1).min()
    maxval1=(data1+err1).max()
    minval2=(data2-err2).min()
    maxval2=(data2+err2).max()

    minval=min(minval1, minval2)
    maxval=max(maxval1, maxval2)
    
    rng=max(abs(minval), abs(maxval))
    return array([-1.1*rng,1.1*rng])

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

    key=PlotKey(0.9,0.2,[m1pts,m2pts], halign='right')

    xrng=array( [0.5*st['s2n'].min(), 1.5*st['s2n'].max()] )
    cyrng=get_symmetric_range(st['c1'],st['c1_err'],st['c2'],st['c2_err'])
    myrng=get_symmetric_range(st['m1'],st['m1_err'],st['m2'],st['m2_err'])

    mplt=FramedPlot()
    cplt=FramedPlot()
    mplt.xlog=True
    cplt.xlog=True

    mplt.xrange=xrng
    cplt.xrange=xrng
    #mplt.yrange=myrng
    #cplt.yrange=cyrng
    mplt.yrange=[-0.15,0.15]
    cplt.yrange=[-0.01,0.01]

    mplt.xlabel=s2n_field
    mplt.ylabel='m'
    cplt.xlabel=s2n_field
    cplt.ylabel='c'

    zplt=Curve(xrng, [0,0])
    mallow=Curve(xrng, [-0.004, 0.004])
    callow=Curve(xrng, [-0.0004, 0.0004])

    mallow=FillBetween(xrng, [0.004,0.004], 
                       xrng, [-0.004,-0.004],
                       color='grey80')
    callow=FillBetween(xrng, [0.0004,0.0004], 
                       xrng, [-0.0004,-0.0004],
                       color='grey80')



    labels=get_labels(fitters[0])

    mplt.add( mallow, zplt, m1pts, m1errpts, m2pts, m2errpts, key )
    mplt.add(*labels)

    cplt.add( callow, zplt, c1pts, c1errpts, c2pts, c2errpts )

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
        #st['s2n'][i] = bf.s2n_min
        st['s2n'][i] = bf.avg['s2n'].mean()
        st['m1'][i] = bf.g1fit.pars[0]
        st['m1_err'][i] = bf.g1fit.perr[0]
        st['c1'][i] = bf.g1fit.pars[1]
        st['c1_err'][i] = bf.g1fit.perr[1]

        st['m2'][i] = bf.g2fit.pars[0]
        st['m2_err'][i] = bf.g2fit.perr[0]
        st['c2'][i] = bf.g2fit.pars[1]
        st['c2_err'][i] = bf.g2fit.perr[1]

    return st

def get_s2n_ranges(options, args):
    if options.nbin:
        nbin=int(options.nbin)
        s2n_range=options.s2n.split(',')
        s2n_range=[float(s) for s in s2n_range]
        log10min=log10(s2n_range[0])
        log10max=log10(s2n_range[1])

        s2n_vals=logspace(log10min, log10max, nbin+1)
    else:
        s2n_vals=[float(s2n) for s2n in args]

    s2n_minvals=[]
    s2n_maxvals=[]
    for i in xrange(len(s2n_vals)-1):
        s2n_minvals.append(s2n_vals[i])
        s2n_maxvals.append(s2n_vals[i+1])

    return s2n_minvals, s2n_maxvals

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if (options.run is None 
            or (len(args) < 1 and options.nbin is None)):
        parser.print_help()
        sys.exit(1)

    s2n_minvals, s2n_maxvals=get_s2n_ranges(options,args)

    fitters=[]
    for i in xrange(len(s2n_minvals)):
        s2n_range=[ s2n_minvals[i], s2n_maxvals[i] ]
        bf=BiasFitter(options.run,s2n_range,
                      psfnums=options.psfnums,
                      objtype=options.type,
                      s2n_field=options.field,
                      s2_max=float(options.s2))
        #bf.show()
        fitters.append(bf)

    st = get_stats(fitters)
    aprint(st,fancy=True)
    doplot(fitters, st, options.field)
    

main()
