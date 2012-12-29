"""
    %prog [options] run
"""

import sys
import os

from numpy import sqrt, linspace

import cluster_step
from cluster_step import files, stats, prior

import esutil as eu
from esutil.numpy_util import aprint
from esutil.stat import histogram

from biggles import FramedPlot, FramedArray, Table, Points, PlotKey, \
        SymmetricErrorBarsX, SymmetricErrorBarsY, Histogram, Curve

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--psfnums',default=None,
                  help='restrict to these PSFs, comma separated')
parser.add_option('--sh',default=None,
                  help='restrict to these shears, comma separated')

parser.add_option('-b','--binsize',default=0.02,
                  help="bin size, default %default")
parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")



def write_eps(arr):
    raise RuntimeError("implement")
    #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

def do_histograms(data, binsize):

    err2=data['gcov'][:,0,0] + data['gcov'][:,1,1]
    wts = 1.0/(stats.SHAPE_NOISE**2 + err2)
    h1=histogram(data['g'][:,0], binsize=binsize, min=-1., max=1., weights=wts)
    h2=histogram(data['g'][:,1], binsize=binsize, min=-1., max=1., weights=wts)

    gtot = sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)

    h=histogram(gtot, binsize=binsize, min=0., max=1., more=True, weights=wts)

    return h1,h2,h

def do_fit(h):
    hvals=h['hist'].astype('f8')
    res=prior.fit_gprior_exp(h['center'], hvals)
    return res

def do_fit_new(data):
    gvals=sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)
    mg=prior.fit_gprior_exp_gmix(gvals)
    return mg



def get_prior_vals(h, fitres):
    gp=prior.GPriorExp(fitres['A'], fitres['a'], fitres['g0'],
                       gmax=fitres['gmax'])
    xvals=linspace(h['low'][0], h['high'][-1], 1000)
    yvals = gp.prior_gabs(xvals)

    return xvals,yvals

def get_prior_vals_new(h, mg):
    xvals=linspace(h['low'][0], h['high'][-1], 1000)
    binsize=h['low'][1]-h['low'][0]
    yvals = h['hist'].sum()*binsize*mg.eval(xvals)
    return xvals,yvals

def oplot_gaussians(plt, h, mg):
    import pcolors

    ngauss=mg.gmm.n_states
    colors=pcolors.rainbow(ngauss, 'hex')

    xvals=linspace(h['low'][0], h['high'][-1], 1000)
    binsize=h['low'][1]-h['low'][0]
    hsum=h['hist'].sum()
    for i in xrange(ngauss):
        yvals = hsum*binsize*mg.evalone(xvals,i)

        line=Curve(xvals, yvals, color=colors[i])

        plt.add(line)
 

def doplot(h1, h2, h, binsize, fitres, show=False):
    tab=Table(2,1)

    #arr.uniform_limits=1
    #arr.xrange=[-1,1]

    pltboth=FramedPlot()
    pltboth.xlabel=r'$g$'

    hplt1=Histogram(h1['hist'], x0=h1['low'][0], binsize=binsize,color='red')
    hplt2=Histogram(h2['hist'], x0=h2['low'][0], binsize=binsize,color='blue')

    hplt1.label=r'$g_1$'
    hplt2.label=r'$g_2$'
    keyboth=PlotKey(0.9,0.9,[hplt1,hplt2],halign='right')

    pltboth.add(hplt1, hplt2, keyboth)
    #arr[0,0].add(hplt1, hplt2, keyboth)
    tab[0,0]=pltboth
    

    plt=FramedPlot()
    plt.xlabel=r'$|g|$'

    hplt=Histogram(h['hist'], x0=h['low'][0], binsize=binsize)

    
    xfit,yfit=get_prior_vals(h, fitres)
    line=Curve(xfit, yfit, color='blue')
    plt.add(line)

    #oplot_gaussians(plt, h, fitres)


    #arr[1,0].add(hplt,key)
    plt.add(hplt)

    tab[1,0]=plt
    
    if show:
        #arr.show()
        #plt.show()
        tab.show()

    #return plt
    #return arr
    return tab


def main():
    options,args = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help()
        sys.exit(45)

    run=args[0]
    binsize=float(options.binsize)
    objtype=options.type
    if objtype:
        print 'selecting type:',objtype

    data=files.read_output_set(run, options.psfnums, options.sh, 
                               objtype=objtype,
                               s2n_min=100,
                               s2_max=1./4.,
                               gsens_min=0.95,
                               gerr_max=0.05)


    h1,h2,h=do_histograms(data, binsize)
    fitres=do_fit(h)
    #fitres=do_fit(data)
    doplot(h1,h2,h, binsize, fitres, show=options.show)


main()
