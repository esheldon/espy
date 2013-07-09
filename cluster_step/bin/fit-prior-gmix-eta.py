"""
    %prog [options]

"""

import sys
import os

import numpy
from numpy import sqrt, linspace, zeros, where

import lensing

import cluster_step
from cluster_step import files, stats, prior

import esutil as eu
from esutil.numpy_util import aprint
from esutil.stat import histogram

import biggles
from biggles import FramedPlot, FramedArray, Table, Points, PlotKey, \
        SymmetricErrorBarsX, SymmetricErrorBarsY, Histogram, Curve,\
        ErrorBarsX

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")

       
class FitRunner(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])
        self.show=options.show
        self.options=options
        self.otype = 'nosplit-gmix-eta'

        self.ellip_name='eta'
        self.min_eta_2d = None
        self.max_eta_2d = None
        self.max_eta=None

        self.ellip_name='g'
        self.min_g_2d = -1.0
        self.max_g_2d = 1.0
        self.max_g=1.0

        self.ngauss=3
        #self.ngauss=4

    def get_dtype(self, npars):
        npars=self.ngauss
        dt=[('minmag', 'f8'),
            ('maxmag', 'f8'),
            ('weight','f8',npars),
            ('sigma','f8',npars)]
        return dt

    def get_struct(self,nbin):
        npars=4
        dt=self.get_dtype(npars)
        st=zeros(nbin, dtype=dt)
        return st

    def get_data(self):
        data0=files.read_prior_original()

        
        eta1,eta2 = lensing.util.g1g2_to_eta1eta2(data0['g'][:,0], data0['g'][:,1])

        data=eu.numpy_util.add_fields(data0, [('eta','f8',2)])
        data['eta'][:,0] = eta1
        data['eta'][:,1] = eta2

        self.data=data

    def go(self):

        self.get_data()

        maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
                 21.25,21.50,21.75,22.0,22.25,22.5,22.75,
                 23.0,23.25,23.5,23.75,24.0,24.25,24.5]

        nbin=len(maglims)-1

        st=self.get_struct(nbin)
        for i in xrange(nbin):
            minmag=maglims[i]
            maxmag=maglims[i+1]

            w,=where(  (self.data['mag'] > minmag)
                     & (self.data['mag'] < maxmag))

            gprior=self.do_fit(self.data[w])

            self.doplot(gprior, minmag, maxmag)
            print '\nstats:'
            print gprior.means_.shape
            print gprior.means_
            print gprior.covars_.shape
            print gprior.covars_
            print gprior.weights_.shape
            print gprior.weights_

            st['minmag'][i] = minmag
            st['maxmag'][i] = maxmag
            st['weight'][i] = gprior.weights_
            st['sigma'][i] = sqrt( gprior.covars_[:,0] ) # round

            if False and self.show:
                key=raw_input('hit a key (q to quit): ')
                if key=='q':
                    stop

        #self.plot_fits(st)
        self.write_data(st)

    def plot_fits(self,st):
        import biggles
        biggles.configure( 'default', 'fontsize_min', 2)
        parnames=['A','a','g0','gmax']

        npars=len(parnames)
        tab=Table(npars,1)

        magmiddle=(st['minmag'] + st['maxmag'])/2
        for i in xrange(npars):

            yval=st['pars'][:,i]
            yerr=st['perr'][:,i]

            ymean=yval.mean()
            ystd=yval.std()
            yrange=[ymean-3.5*ystd, ymean+3.5*ystd]
            
            pts=Points(magmiddle, yval, type='filled circle')
            yerrpts=SymmetricErrorBarsY(magmiddle,yval,yerr)
            xerrpts=ErrorBarsX(yval, st['minmag'], st['maxmag'])

            
            plt=FramedPlot()
            plt.yrange=yrange
            plt.add(pts,xerrpts, yerrpts)
            plt.xlabel='mag'
            plt.ylabel=parnames[i]

            tab[i,0] = plt

        if self.show:
            tab.show()

        d=files.get_prior_dir()
        d=os.path.join(d, 'plots')
        epsfile='pofe-pars-%s.eps' % self.otype
        epsfile=os.path.join(d,epsfile)
        eu.ostools.makedirs_fromfile(epsfile)
        print epsfile
        tab.write_eps(epsfile)
        os.system('converter -d 100 %s' % epsfile)


    def write_data(self, st):
        import fitsio
        outfile=files.get_gprior_path(type=self.otype)
        print 'writing:',outfile
        with fitsio.FITS(outfile, mode='rw', clobber=True) as fobj:
            fobj.write(st)

    def write_eps(self, arr):
        raise RuntimeError("implement")
        #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

    def do_histograms(self, minmag, maxmag, ellip_name):

        data=self.data
        w,=where((data['mag'] > minmag) & (data['mag'] < maxmag))
        more=True
        data=self.data

        g1=data[ellip_name][w,0]
        g2=data[ellip_name][w,1]

        gtot = sqrt(g1**2 + g2**2)

        sigma=gtot.std()
        binsize=0.1*sigma

        if ellip_name=='eta':
            mine_2d = self.min_eta_2d
            maxe_2d = self.max_eta_2d
            maxe = self.max_eta
            self.binsize_eta=binsize
        else:
            mine_2d = self.min_g_2d
            maxe_2d = self.max_g_2d
            maxe = self.max_g
            self.binsize_g=binsize

        h1=histogram(g1, binsize=binsize, min=mine_2d, max=maxe_2d, more=more)
        h2=histogram(g2, binsize=binsize, min=mine_2d, max=maxe_2d, more=more)


        h=histogram(gtot, binsize=binsize, min=0., max=maxe, more=more)

        return h1, h2, h

    def do_fit(self, data):
        res=prior.fit_gprior_gmix_em(data['eta'][:,0],
                                     data['eta'][:,1],
                                     self.ngauss)
        return res


    def get_prior_vals(self, fitres, h):

        pars=[fitres['A'], fitres['a'], fitres['g0'], fitres['gmax']]
        gp=prior.GPriorExp(pars)

        xvals=linspace(h['low'][0], h['high'][-1], 1000)
        yvals = gp.prior1d(xvals)

        return xvals,yvals, gp


    def doplot(self, gprior, minmag, maxmag):
        from lensing.util import eta1eta2_to_g1g2
        tab=Table(2,2)
        tab.title='%s %.2f %.2f ' % (self.otype, minmag, maxmag)

        h1_g,h2_g,h_g = self.do_histograms(minmag, maxmag, 'g')
        h1_eta,h2_eta,h_eta = self.do_histograms(minmag, maxmag, 'eta')

        nrand=1000000
        binsize_eta=self.binsize_eta
        binsize_g=self.binsize_g

        rbinsize_eta=binsize_eta*0.2
        rbinsize_g=binsize_g*0.2

        gr = gprior.sample(nrand)
        eta1_rand=gr[:,0]
        eta2_rand=gr[:,1]
        eta_rand = numpy.sqrt( eta1_rand**2 + eta2_rand**2 )

        g1_rand,g2_rand = eta1eta2_to_g1g2(eta1_rand, eta2_rand)
        g_rand = numpy.sqrt(g1_rand**2 + g2_rand**2)


        hrand_eta=histogram(eta_rand, binsize=rbinsize_eta,
                            min=0,max=self.max_eta,
                            more=True)
        h1rand_eta=histogram(eta1_rand, binsize=rbinsize_eta,
                             min=self.min_eta_2d, max=self.max_eta_2d,
                             more=True)


        hrand_g=histogram(g_rand, binsize=rbinsize_g,
                          min=0,max=self.max_g,
                          more=True)
        h1rand_g=histogram(g1_rand, binsize=rbinsize_g,
                           min=self.min_g_2d, max=self.max_g_2d,
                           more=True)


        # eta 2d plots
        bratio_eta = self.binsize_eta/rbinsize_eta
        hrand_eta['hist'] = hrand_eta['hist']*bratio_eta*float(h_eta['hist'].sum())/nrand
        h1rand_eta['hist'] = h1rand_eta['hist']*bratio_eta*float(h1_eta['hist'].sum())/h1rand_eta['hist'].sum()

        pltboth_eta=FramedPlot()
        pltboth_eta.xlabel=r'$\eta$'

        hplt1_eta=Histogram(h1_eta['hist'], x0=h1_eta['low'][0], binsize=binsize_eta,color='darkgreen')
        hplt2_eta=Histogram(h2_eta['hist'], x0=h2_eta['low'][0], binsize=binsize_eta,color='blue')
        hpltrand_eta=Histogram(hrand_eta['hist'], x0=hrand_eta['low'][0], binsize=rbinsize_eta,
                           color='red')
        hplt1rand_eta=Histogram(h1rand_eta['hist'], x0=h1rand_eta['low'][0], binsize=rbinsize_eta,
                           color='red')

        hplt1_eta.label=r'$\eta_1$'
        hplt2_eta.label=r'$\eta_2$'
        hplt1rand_eta.label='rand'
        hpltrand_eta.label='rand'
        keyboth_eta=PlotKey(0.9,0.9,[hplt1_eta,hplt2_eta,hplt1rand_eta],halign='right')

        pltboth_eta.add(hplt1_eta, hplt2_eta, hplt1rand_eta, keyboth_eta)

        tab[0,0]=pltboth_eta

        plt1d_eta=FramedPlot()
        plt1d_eta.xlabel=r'$|\eta|$'

        hplt_eta=Histogram(h_eta['hist'], x0=h_eta['low'][0], binsize=binsize_eta)
        hplt_eta.label=r'$|\eta|$'

        key_eta=PlotKey(0.9,0.9,[hplt_eta,hpltrand_eta],halign='right')
        plt1d_eta.add(hplt_eta, hpltrand_eta, key_eta)

        tab[1,0]=plt1d_eta

        # g plots
        
        bratio_g = self.binsize_g/rbinsize_g
        hrand_g['hist'] = hrand_g['hist']*bratio_g*float(h_g['hist'].sum())/nrand
        h1rand_g['hist'] = h1rand_g['hist']*bratio_g*float(h1_g['hist'].sum())/h1rand_g['hist'].sum()


        pltboth_g=FramedPlot()
        pltboth_g.xlabel=r'$g$'

        hplt1_g=Histogram(h1_g['hist'], x0=h1_g['low'][0], binsize=binsize_g,color='darkgreen')
        hplt2_g=Histogram(h2_g['hist'], x0=h2_g['low'][0], binsize=binsize_g,color='blue')
        hpltrand_g=Histogram(hrand_g['hist'], x0=hrand_g['low'][0], binsize=rbinsize_g,
                           color='red')
        hplt1rand_g=Histogram(h1rand_g['hist'], x0=h1rand_g['low'][0], binsize=rbinsize_g,
                           color='red')

        hplt1_g.label=r'$g_1$'
        hplt2_g.label=r'$g_2$'
        hplt1rand_g.label='rand'
        hpltrand_g.label='rand'
        keyboth_g=PlotKey(0.9,0.9,[hplt1_g,hplt2_g,hplt1rand_g],halign='right')

        pltboth_g.add(hplt1_g, hplt2_g, hplt1rand_g, keyboth_g)

        tab[0,1]=pltboth_g

        plt1d_g=FramedPlot()
        plt1d_g.xlabel=r'$|g|$'

        hplt_g=Histogram(h_g['hist'], x0=h_g['low'][0], binsize=binsize_g)
        hplt_g.label='|g|'

        key_g=PlotKey(0.9,0.9,[hplt_g,hpltrand_g],halign='right')
        plt1d_g.add(hplt_g, hpltrand_g, key_g)

        tab[1,1]=plt1d_g
        
        if self.show:
            tab.show()

        d=files.get_prior_dir()
        d=os.path.join(d, 'plots')
        epsfile='pofe-pofeta-%.2f-%.2f-%s.eps' % (minmag,maxmag,self.otype)
        epsfile=os.path.join(d,epsfile)
        eu.ostools.makedirs_fromfile(epsfile)
        print epsfile
        tab.write_eps(epsfile)
        os.system('converter -d 100 %s' % epsfile)

        return tab


def main():
    runner=FitRunner()
    runner.go()

main()
