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

    more=True
    h1=histogram(data['g'][:,0], binsize=binsize, min=-1., max=1., more=more)
    h2=histogram(data['g'][:,1], binsize=binsize, min=-1., max=1., more=more)

    gtot = sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)

    h=histogram(gtot, binsize=binsize, min=0., max=1., more=more)
    #h=histogram(gtot, binsize=binsize, more=more)

    return h1,h2,h

def do_fit(h):
    hvals=h['hist'].astype('f8')
    res=prior.fit_gprior_exp(h['center'], hvals)
    return res


def get_prior_vals(h, fitres):
    gp=prior.GPriorExp(fitres['A'], fitres['a'], fitres['g0'], fitres['gmax'])
    xvals=linspace(h['low'][0], h['high'][-1], 1000)
    yvals = gp.prior1d(xvals)

    return xvals,yvals, gp


def doplot(h1, h2, h, binsize, fitres, show=False, title=''):
    tab=Table(2,1)
    tab.title=title

    xfit,yfit,gprior = get_prior_vals(h, fitres)

    nrand=100000
    #nrand=h['hist'].sum()

    g1rand,g2rand=gprior.sample2d(nrand)
    grand=gprior.sample1d(nrand)

    hrand=histogram(grand, binsize=binsize, min=0., max=1., more=True)
    h1rand=histogram(g1rand, binsize=binsize, min=-1., max=1., more=True)

    fbinsize=xfit[1]-xfit[0]
    hrand['hist'] = hrand['hist']*float(yfit.sum())/hrand['hist'].sum()*fbinsize/binsize
    h1rand['hist'] = h1rand['hist']*float(h1['hist'].sum())/h1rand['hist'].sum()

    #arr.uniform_limits=1
    #arr.xrange=[-1,1]

    pltboth=FramedPlot()
    pltboth.xlabel=r'$g$'

    hplt1=Histogram(h1['hist'], x0=h1['low'][0], binsize=binsize,color='red')
    hplt2=Histogram(h2['hist'], x0=h2['low'][0], binsize=binsize,color='blue')
    hpltrand=Histogram(hrand['hist'], x0=hrand['low'][0], binsize=binsize,
                       color='magenta')
    hplt1rand=Histogram(h1rand['hist'], x0=h1rand['low'][0], binsize=binsize,
                       color='magenta')

    hplt1.label=r'$g_1$'
    hplt2.label=r'$g_2$'
    hplt1rand.label='rand'
    hpltrand.label='rand'
    keyboth=PlotKey(0.9,0.9,[hplt1,hplt2,hplt1rand],halign='right')

    pltboth.add(hplt1, hplt2, hplt1rand, keyboth)
    #arr[0,0].add(hplt1, hplt2, keyboth)
    tab[0,0]=pltboth
    

    plt=FramedPlot()
    plt.xlabel=r'$|g|$'

    hplt=Histogram(h['hist'], x0=h['low'][0], binsize=binsize)
    hplt.label='|g|'

    
    line=Curve(xfit, yfit, color='blue')
    line.label='model'

    key=PlotKey(0.9,0.9,[hplt,line,hpltrand],halign='right')
    plt.add(line, hplt, hpltrand, key)


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
                               gerr_max=0.05,
                               subtract_mean=True,
                               progress=True)


    h1,h2,h=do_histograms(data, binsize)
    fitres=do_fit(h)

    title=run
    if options.psfnums is not None:
        title="%s-p%s" % (title,options.psfnums)
    if options.sh is not None:
        title="%s-s%s" % (title,options.sh)

    doplot(h1,h2,h, binsize, fitres, show=options.show, title=title)


main()
