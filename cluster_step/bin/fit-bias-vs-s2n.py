"""
    %prog [options]

Either

    - Send -n nbin, which will result in logarithmic binning between the min
    s/n and 1000.
    - Send the bin edges as args
        Bins will be [s2n1,s2n2], [s2n2,s2n3], ...
"""

import sys
import os
from numpy import zeros, logspace, log10, array, where

import cluster_step
from cluster_step.fitters import BiasFitter
from cluster_step import files
from esutil.numpy_util import aprint

from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY, \
        PlotLabel, Table, PlotKey, FillBetween


from optparse import OptionParser

parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')

parser.add_option('-p','--psfnums',default='1,2,3,4,5,6',
                  help='restrict to these PSFs, comma separated')

parser.add_option('-n','--nbin',default=10,
                  help=('number of bins for logarithmic binning in s/n '
                        ', default %default'))

parser.add_option('--s2n',default='10,200', help="s/n range, %default")
parser.add_option('--Ts2n',default='2,200', help="Ts2n range, %default")
parser.add_option('--sratio',default='1.0,10.0',
                  help='sratio range, %default')
parser.add_option('--Tmean',default='4,20',
                  help='Tmean range, %default')
parser.add_option('--mag',default='0,100',
                  help='mag range, %default')


parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")

parser.add_option('-f','--field',default='s2n_w',
                  help="field for S/N, default %default")

"""
        # starting point for labels
        self.lab1_loc=[1.-0.075, 0.1]
        self.lab1_halign='right'
        self.lab1_yshift=0.075

        self.lab2_loc=[0.075,0.075]
        self.lab2_halign='left'
        self.lab2_yshift=+0.075

        self.lab3_loc=[1-0.075,1-0.075]
        self.lab3_halign='right'
        self.lab3_yshift=-0.075

def get_s2n_label(s2n_range):
    labs=r'$%.2f < %s < %.2f$' % (s2n_range[0],
                                  s2n_field,
                                  s2n_range[1])
    yshift=1*self.lab2_yshift
    lab=PlotLabel(self.lab2_loc[0],self.lab2_loc[1]+yshift,
                  labs,
                  halign=self.lab2_halign)
    return lab

"""

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



    #labels=get_labels(fitters[0])

    mplt.add( mallow, zplt, m1pts, m1errpts, m2pts, m2errpts, key )
    #mplt.add(*labels)

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
    nbin=int(options.nbin)
    s2n_range=options.s2n.split(',')
    s2n_range=[float(s) for s in s2n_range]
    log10min=log10(s2n_range[0])
    log10max=log10(s2n_range[1])

    s2n_vals=logspace(log10min, log10max, nbin+1)

    s2n_minvals=[]
    s2n_maxvals=[]
    for i in xrange(len(s2n_vals)-1):
        s2n_minvals.append(s2n_vals[i])
        s2n_maxvals.append(s2n_vals[i+1])

    return s2n_minvals, s2n_maxvals

def main():
    options,args = parser.parse_args(sys.argv[1:])

    if (options.run is None or options.nbin is None):
        parser.print_help()
        sys.exit(1)

    run=options.run
    s2n_range=[float(s) for s in options.s2n.split(',')]
    sratio_range=[float(s) for s in options.sratio.split(',')]
    Ts2n_range=[float(s) for s in options.Ts2n.split(',')]
    Tmean_range=[float(s) for s in options.Tmean.split(',')]
    mag_range=[float(s) for s in options.mag.split(',')]
    objtype=options.type
    psfnums=[float(s) for s in options.psfnums.split(',')]

    s2n_field=options.field

    reader=files.Reader(run=run, 
                        objtype=objtype,
                        psfnums=psfnums,
                        s2n_range=s2n_range,
                        sratio_range=sratio_range,
                        Ts2n_range=Ts2n_range,
                        Tmean_range=Tmean_range,
                        mag_range=mag_range,
                        progress=True)

    data=reader.get_data()

    s2n_minvals, s2n_maxvals=get_s2n_ranges(options,args)

    fitters=[]
    for i in xrange(len(s2n_minvals)):
        w,=where((data[s2n_field] > s2n_minvals[i])
                 &
                 (data[s2n_field] < s2n_maxvals[i]))

        print s2n_minvals[i],s2n_maxvals[i],w.size

        bf=BiasFitter(data[w], run, s2n_field=options.field)
        fitters.append(bf)

    st = get_stats(fitters)
    aprint(st,fancy=True)
    doplot(fitters, st, s2n_field)
    

main()
