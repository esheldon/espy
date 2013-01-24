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

parser.add_option('--sratio',default='1.0,1.e6',
                  help='sratio range, %default')

parser.add_option('-f','--field',default='s2n_w',
                  help="field for S/N, default %default")

parser.add_option('-s','--set',default='use1',
                  help="selection set")

def get_labels(run, psfnums, setname, nobj, sratio_range):


    run_mess='run: %s' % run
    set_mess='set: %s' % setname
    nobj_mess='nobj: %s' % nobj

    if sratio_range[1] > 1000:
        srat_mess = r'$\sigma_{gal}/\sigma_{psf} > %.2f$' % sratio_range[0]
    else:
        srat_mess=r'$%.2f < \sigma_{gal}/\sigma_{psf} < %.2f$' % tuple(sratio_range)

    pstr=[str(s) for s in psfnums]
    pstr=','.join(pstr)
    pmess="p: %s" % pstr

    halign='left'
    x=0.05
    y=0.9
    inc=-0.075

    runlab=PlotLabel(x, y, run_mess, halign=halign)
    y += inc
    pnumlab=PlotLabel(x, y, pmess, halign=halign)
    y += inc
    setlab=PlotLabel(x, y, set_mess, halign=halign)
    y += inc
    sratlab=PlotLabel(x, y, srat_mess, halign=halign)
    y += inc
    nobjlab=PlotLabel(x, y, nobj_mess, halign=halign)
    
    return (runlab,pnumlab,setlab,sratlab,nobjlab)

def get_symmetric_range(data1, err1, data2, err2):
    minval1=(data1-err1).min()
    maxval1=(data1+err1).max()
    minval2=(data2-err2).min()
    maxval2=(data2+err2).max()

    minval=min(minval1, minval2)
    maxval=max(maxval1, maxval2)
    
    rng=max(abs(minval), abs(maxval))
    return array([-1.1*rng,1.1*rng])

def doplot(run, st, s2n_field, setname, s2n_range, sratio_range, psfnums, nobj):
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

    xrng=array( [0.5*s2n_range[0], 1.5*s2n_range[1]])
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



    labels=get_labels(run,
                      psfnums,
                      setname,
                      nobj,
                      sratio_range)

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
        res1=bf.g1fit.get_result()
        res2=bf.g2fit.get_result()
        if res1['perr'] is not None:
            #st['s2n'][i] = bf.s2n_min
            st['s2n'][i] = bf.avg['s2n'].mean()

            st['m1'][i] = res1['pars'][0]
            st['m1_err'][i] = res1['perr'][0]
            st['c1'][i] = res1['pars'][1]
            st['c1_err'][i] = res1['perr'][1]

            st['m2'][i] = res2['pars'][0]
            st['m2_err'][i] = res2['perr'][0]
            st['c2'][i] = res2['pars'][1]
            st['c2_err'][i] = res2['perr'][1]

    w,=where(st['s2n'] != 0)
    st=st[w]
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
    psfnums=[int(s) for s in options.psfnums.split(',')]

    s2n_field=options.field

    reader=files.Reader(run=run, 
                        psfnums=psfnums,
                        s2n_range=s2n_range,
                        sratio_range=sratio_range,
                        setname=options.set,
                        progress=True)

    data=reader.get_data()

    s2n_minvals, s2n_maxvals=get_s2n_ranges(options,args)

    fitters=[]
    for i in xrange(len(s2n_minvals)):
        import mcmc
        w,=where((data[s2n_field] > s2n_minvals[i])
                 &
                 (data[s2n_field] < s2n_maxvals[i]))

        print s2n_minvals[i],s2n_maxvals[i],w.size

        bf=BiasFitter(data[w], s2n_field=options.field)
        fitters.append(bf)
        print bf.g1fit
        print bf.g2fit

        print 'g1 guess:',bf.g1fit.guess
        print 'g2 guess:',bf.g2fit.guess

        print bf.g1fit.x,bf.g1fit.y
        print bf.g2fit.x,bf.g2fit.y
        stop

        plt1=mcmc.plot_results(bf.g1fit.trials)
        plt2=mcmc.plot_results(bf.g2fit.trials)
        stop



    st = get_stats(fitters)
    aprint(st,fancy=True)
    doplot(run, st, s2n_field, options.set, s2n_range,sratio_range, psfnums, data.size)
    

main()
