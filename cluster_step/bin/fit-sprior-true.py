"""
    %prog [options]
"""

import sys
import os

from numpy import sqrt, linspace, zeros, where

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

#parser.add_option('-b','--binsize',default=0.5,
#                  help="bin size, default %default")
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")
parser.add_option('-t','--type',default=None,
                  help="object type")

        
class FitRunner(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.objtype=options.type
        if self.objtype is None:
            raise ValueError("send -t objtype")
        if self.objtype not in ['gexp','gdev']:
            raise ValueError("send objtype gexp or gdev")

        #self.binsize=float(options.binsize)
        self.show=options.show

        self.options=options

    def get_dtype(self, npars):
        dt=[('minmag', 'f8'),
            ('maxmag', 'f8'),
            ('pars','f8',npars),
            ('perr','f8',npars),
            ('pcov','f8',(npars,npars))]
        return dt

    def get_struct(self,nbin):
        npars=3
        dt=self.get_dtype(npars)
        st=zeros(nbin, dtype=dt)
        return st

    def load_data(self):
        data=files.read_prior_original()
        logic = (data['scale'] < 25)
        if self.objtype=='gdev':
            w,=where((data['n'] > 3) & logic)
        else:
            w,=where((data['n'] < 2) & logic)
        self.data=data[w]
        self.set_sigma()

    def set_sigma(self):
        import sersic
        if self.objtype=='gdev':
            fac = sersic.re2T_factor(4)
        else:
            fac = sersic.re2T_factor(1)

        # T = fac*re^2
        # sigma = sqrt(T/2) = sqrt(fac/2)*re
        self.sigma = sqrt(fac/2)*self.data['scale']

        #self.sigma = fac*self.data['scale']**2

    def go(self):
        print 'fitting:',self.objtype
        self.load_data()

        # it is important that a bin starts at 23 since
        # that is where samples changed from SDSS to HST
        # old, fewer
        #maglims=[18.0,20.0,20.5,21.0,21.5,22.0,22.5,23.0,
        #         23.5,24.0,24.5]
        if self.objtype=='gexp':
            self.binfac=0.15
            maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
                     21.25,21.50,21.75,22.0,22.25,22.5,22.75,
                     23.0,23.25,23.5,23.75,24.0,24.25,24.5]
        else:
            self.binfac=0.1
            maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
                     21.25,21.50,21.75,22.0,22.25,22.5,22.75,
                     23.0,23.5,24.0,24.5]



        nbin=len(maglims)-1

        st=self.get_struct(nbin)
        for i in xrange(nbin):
            minmag=maglims[i]
            maxmag=maglims[i+1]
            h,mean,sigma = self.do_histogram(minmag, maxmag)
            fitter=self.do_fit(h,mean,sigma)
            fitres=fitter.get_result()

            print fitter
            self.doplot(fitter, h, minmag, maxmag)

            st['minmag'][i] = minmag
            st['maxmag'][i] = maxmag
            st['pars'][i] = fitres['pars']
            st['perr'][i] = fitres['perr']
            st['pcov'][i] = fitres['pcov']

        self.plot_fits(st)
        self.write_data(st)

    def plot_fits(self,st):
        biggles.configure( 'default', 'fontsize_min', 2)
        parnames=['A','mean','sigma']

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

        d=files.get_prior_dir()
        d=os.path.join(d, 'plots')
        epsfile='pofs-pars-%s.eps' % self.objtype
        epsfile=os.path.join(d,epsfile)
        eu.ostools.makedirs_fromfile(epsfile)
        print epsfile
        tab.write_eps(epsfile)
        os.system('converter -d 100 %s' % epsfile)

        if self.show:
            tab.show()

    def write_data(self, st):
        import fitsio
        outfile=files.get_sprior_path(type=self.objtype)
        print 'writing:',outfile
        with fitsio.FITS(outfile, mode='rw', clobber=True) as fobj:
            fobj.write(st)

    def write_eps(self, arr):
        raise RuntimeError("implement")
        #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

    def do_histogram(self, minmag, maxmag):

        data=self.data
        w,=where((data['mag'] > minmag) & (data['mag'] < maxmag))
        more=True
        data=self.data

        sigma_vals=self.sigma[w]

        mean=sigma_vals.mean()
        sigma=sigma_vals.std()

        for i in xrange(3):
            w,=where(sigma_vals < (mean+4.*sigma))

            sigma_vals=sigma_vals[w]

            mean=sigma_vals.mean()
            sigma=sigma_vals.std()

        binsize=sigma*self.binfac
        self.binsize=binsize

        h=histogram(sigma_vals, binsize=binsize, more=more)

        if False:
            hp=biggles.Histogram(h['hist'], x0=h['low'][0], binsize=binsize)
            plt=biggles.FramedPlot()
            plt.add(hp)
            plt.show()

        return h, mean, sigma

    def do_fit(self, h, mean, sigma):
        from mcmc import LogNormalFitter
        nwalkers=100
        burnin=400
        nstep=100

        hvals=h['hist'].astype('f8')
        herr=sqrt(hvals)
        herr=herr.clip(1.0, herr.max())

        ntot=hvals.sum()
        guess=[
            ntot*self.binsize,
            mean, # arbitrary
            sigma # arbitrary
        ]
        print guess
        #width=[sqrt(ntot)*self.binsize, 100., 100.]
        width=[1.0, 100., 100.]
        #width=None
        nlf=LogNormalFitter(h['center'], hvals, guess, nwalkers, burnin, nstep,
                            yerr=herr, width=width)
        return nlf

    def doplot(self, fitres, h, minmag, maxmag):
        tab=biggles.Table(2,1)

        plt=FramedPlot()
        plt.title='%s %.2f %.2f ' % (self.objtype, minmag, maxmag)
        plt.xlabel=r'$\sigma$'

        sprior=fitres.get_model()

        nrand=100000
        binsize=self.binsize

        hplt=Histogram(h['hist'], x0=h['low'][0], binsize=binsize)
        hplt.label='data'

        srand=sprior.sample(nrand)
        hrand=histogram(srand, binsize=binsize, min=h['low'][0], max=h['high'][-1], more=True)
        hrand['hist'] = hrand['hist']*float(h['hist'].sum())/nrand

        hpltrand=Histogram(hrand['hist'], x0=hrand['low'][0], binsize=binsize,
                           color='blue')
        hpltrand.label='rand'


        key=PlotKey(0.9,0.9,[hplt,hpltrand],halign='right')

        plt.add(hplt, hpltrand, key)

        
        tplt=fitres.plot_trials(show=False,fontsize_min=0.88888888)
        
        tab[0,0] = plt
        tab[1,0] = tplt

        if self.show:
            tab.show()

        d=files.get_prior_dir()
        d=os.path.join(d, 'plots')
        epsfile='pofs-%.2f-%.2f-%s.eps' % (minmag,maxmag,self.objtype)
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
