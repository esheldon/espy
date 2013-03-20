"""
    %prog [options] run objtype
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
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")


class FitRunner(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        if len(args) < 2:
            parser.print_help()
            sys.exit(45)

        self.run=args[0]
        self.objtype=args[1]
        if self.objtype not in ['gexp','gdev']:
            raise ValueError("send objtype gexp or gdev")

        self.binsize=float(options.binsize)
        self.show=options.show


        self.title=self.run
        if options.psfnums is not None:
            self.title="%s-p%s" % (self.title,options.psfnums)
        if options.sh is not None:
            self.title="%s-s%s" % (self.title,options.sh)

        self.options=options

    def go(self):
        print 'selecting type:',self.objtype
        self.data=files.read_output_set(self.run, 
                                        self.options.psfnums, 
                                        self.options.sh, 
                                        objtype=self.objtype,
                                        s2n_min=100,
                                        s2_max=1./4.,
                                        gsens_min=0.95,
                                        gerr_max=0.05,
                                        subtract_mean=True,
                                        progress=True)


        self.do_histograms()
        self.do_fit()

        self.doplot()

    def write_eps(self, arr):
        raise RuntimeError("implement")
        #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

    def do_histograms(self):

        more=True
        data=self.data
        binsize=self.binsize

        h1=histogram(data['g'][:,0], binsize=binsize, min=-1., max=1., more=more)
        h2=histogram(data['g'][:,1], binsize=binsize, min=-1., max=1., more=more)

        gtot = sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)

        h=histogram(gtot, binsize=binsize, min=0., max=1., more=more)
        #h=histogram(gtot, binsize=binsize, more=more)

        self.h=h
        self.h1=h1
        self.h2=h2

    def do_fit(self):
        hvals=self.h['hist'].astype('f8')
        if self.objtype=='gdev':
            res=prior.fit_gprior_dev(self.h['center'], hvals)
        else:
            res=prior.fit_gprior_exp(self.h['center'], hvals)

        self.fitres=res


    def get_prior_vals(self):
        fitres=self.fitres
        h=self.h

        if self.objtype=='gdev':
            pars=fitres['pars']
            gp=prior.GPriorDev(pars)
        else:
            pars=[fitres['A'], fitres['a'], fitres['g0'], fitres['gmax']]
            gp=prior.GPriorExp(pars)

        xvals=linspace(h['low'][0], h['high'][-1], 1000)
        yvals = gp.prior1d(xvals)

        return xvals,yvals, gp


    def doplot(self):
        tab=Table(2,1)
        tab.title=self.title

        xfit,yfit,gprior = self.get_prior_vals()

        nrand=100000
        binsize=self.binsize
        h=self.h
        h1=self.h1
        h2=self.h2


        g1rand,g2rand=gprior.sample2d(nrand)
        grand=gprior.sample1d(nrand)

        hrand=histogram(grand, binsize=binsize, min=0., max=1., more=True)
        h1rand=histogram(g1rand, binsize=binsize, min=-1., max=1., more=True)

        fbinsize=xfit[1]-xfit[0]
        hrand['hist'] = hrand['hist']*float(yfit.sum())/hrand['hist'].sum()*fbinsize/binsize
        h1rand['hist'] = h1rand['hist']*float(h1['hist'].sum())/h1rand['hist'].sum()


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
        
        if self.show:
            tab.show()

        return tab


def main():
    runner=FitRunner()
    runner.go()

main()
