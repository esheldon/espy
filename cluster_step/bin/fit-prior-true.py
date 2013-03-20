"""
    %prog [options]

always use the "exp" p(e) shape
"""

import sys
import os

from numpy import sqrt, linspace, zeros, where

import cluster_step
from cluster_step import files, stats, prior

import esutil as eu
from esutil.numpy_util import aprint
from esutil.stat import histogram

from biggles import FramedPlot, FramedArray, Table, Points, PlotKey, \
        SymmetricErrorBarsX, SymmetricErrorBarsY, Histogram, Curve,\
        ErrorBarsX

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")
parser.add_option('-t','--type',default=None,
                  help="object type")
parser.add_option('--evals',action='store_true',
                  help="assume the values are e instead of g")

        
class FitRunner(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.objtype=options.type
        if self.objtype not in ['gexp','gdev']:
            raise ValueError("send objtype gexp or gdev")

        self.show=options.show

        self.options=options
        self.evals=options.evals


    def get_dtype(self, npars):
        dt=[('minmag', 'f8'),
            ('maxmag', 'f8'),
            ('pars','f8',npars),
            ('perr','f8',npars),
            ('pcov','f8',(npars,npars))]
        return dt

    def get_struct(self,nbin):
        npars=4
        dt=self.get_dtype(npars)
        st=zeros(nbin, dtype=dt)
        return st

    def get_data(self):
        data=files.read_prior_original()
        if self.objtype=='gdev':
            w,=where(data['n'] > 3)
        else:
            w,=where(data['n'] < 2)
        self.data=data[w]

    def go(self):
        print 'fitting:',self.objtype
        self.get_data()



        # it is important that a bin starts at 23 since
        # that is where samples changed from SDSS to HST
        #maglims=[0.0,20.0,20.5,21.0,21.5,22.0,22.5,23.0,
        #         23.5,24.0,30.0]
        #maglims=[18.0,20.0,20.5,21.0,21.5,22.0,22.5,23.0,
        #         23.5,24.0,24.5]
        #maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
        #         21.25,21.50,21.75,22.0,22.25,22.5,22.75,
        #         23.0,23.25,23.5,23.75,24.0,24.25,24.5]
        #maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
        #         21.25,21.50,21.75,22.0,22.25,22.5,22.75,
        #         23.0,23.5,24.5]


        if self.objtype=='gexp':
            maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
                     21.25,21.50,21.75,22.0,22.25,22.5,22.75,
                     23.0,23.25,23.5,23.75,24.0,24.25,24.5]
        else:
            maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
                     21.25,21.50,21.75,22.0,22.25,22.5,22.75,
                     23.0,23.5,24.0,24.5]
            #maglims=[18.0,19.0,20.0,20.25,20.5,20.75,21.0,
            #         21.25,21.50,21.75,22.0,22.25,22.5,22.75,
            #         23.0,23.5,24.5]



        nbin=len(maglims)-1

        st=self.get_struct(nbin)
        for i in xrange(nbin):
            minmag=maglims[i]
            maxmag=maglims[i+1]
            h1,h2,h = self.do_histograms(minmag, maxmag)
            fitres=self.do_fit(h)
            self.doplot(fitres, h1, h2, h, minmag, maxmag)

            st['minmag'][i] = minmag
            st['maxmag'][i] = maxmag
            st['pars'][i] = fitres['pars']
            st['perr'][i] = fitres['perr']
            st['pcov'][i] = fitres['pcov']

        self.plot_fits(st)
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
        epsfile='pofe-pars-%s.eps' % self.objtype
        epsfile=os.path.join(d,epsfile)
        eu.ostools.makedirs_fromfile(epsfile)
        print epsfile
        tab.write_eps(epsfile)
        os.system('converter -d 100 %s' % epsfile)


    def write_data(self, st):
        import fitsio
        outfile=files.get_gprior_path(type=self.objtype)
        print 'writing:',outfile
        with fitsio.FITS(outfile, mode='rw', clobber=True) as fobj:
            fobj.write(st)

    def write_eps(self, arr):
        raise RuntimeError("implement")
        #epsfile=files.get_output_path(ftype='sizemag', ext='eps', **self)

    def do_histograms(self, minmag, maxmag):

        data=self.data
        w,=where((data['mag'] > minmag) & (data['mag'] < maxmag))
        more=True
        data=self.data

        if self.evals:
            import lensing
            # assum g is actually e
            e1=data['g'][w,0]
            e2=data['g'][w,1]
            g1=zeros(e1.size,dtype='f8')
            g2=zeros(e1.size,dtype='f8')
            for i in xrange(g1.size):
                g1[i],g2[i] =lensing.util.e1e2_to_g1g2(e1[i],e2[i])
        else:
            g1=data['g'][w,0]
            g2=data['g'][w,1]

        gtot = sqrt(g1**2 + g2**2)

        sigma=gtot.std()
        binsize=0.2*sigma
        self.binsize=binsize

        h1=histogram(g1, binsize=binsize, min=-1., max=1., more=more)
        h2=histogram(g2, binsize=binsize, min=-1., max=1., more=more)


        #h=histogram(gtot, binsize=binsize, min=0., max=1., more=more)
        h=histogram(gtot, binsize=binsize, more=more)

        if False:
            import biggles
            hp=biggles.Histogram(h['hist'], x0=h['low'][0], binsize=binsize)
            plt=biggles.FramedPlot()
            plt.add(hp)
            plt.show()
        return h1, h2, h

    def do_fit(self, h):
        hvals=h['hist'].astype('f8')

        yerr=sqrt(hvals)
        yerr=yerr.clip(1.0, yerr.max())
        ivar=1./yerr**2

        res=prior.fit_gprior_exp_mcmc(h['center'], hvals, ivar, gmax=0.7, Awidth=0.1, gmax_min=0.5)

        return res


    def get_prior(self, fitres):
        pars=fitres['pars']
        gp=prior.GPriorExp(pars)

        return gp

    def get_prior_vals(self, fitres, h):

        pars=[fitres['A'], fitres['a'], fitres['g0'], fitres['gmax']]
        gp=prior.GPriorExp(pars)

        xvals=linspace(h['low'][0], h['high'][-1], 1000)
        yvals = gp.prior1d(xvals)

        return xvals,yvals, gp


    def doplot(self, fitres, h1, h2, h, minmag, maxmag):
        tab=Table(2,1)
        tab.title='%s %.2f %.2f ' % (self.objtype, minmag, maxmag)

        #xfit,yfit,gprior = self.get_prior_vals(fitres, h)
        gprior=self.get_prior(fitres)

        nrand=100000
        binsize=self.binsize


        g1rand,g2rand=gprior.sample2d(nrand)
        grand=gprior.sample1d(nrand)

        #hrand=histogram(grand, binsize=binsize, min=0., max=1., more=True)
        hrand=histogram(grand, binsize=binsize, min=h['low'][0], max=h['high'][-1], more=True)
        h1rand=histogram(g1rand, binsize=binsize, min=-1., max=1., more=True)

        #fbinsize=xfit[1]-xfit[0]
        #hrand['hist'] = hrand['hist']*float(yfit.sum())/hrand['hist'].sum()*fbinsize/binsize
        hrand['hist'] = hrand['hist']*float(h['hist'].sum())/nrand
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

        
        #line=Curve(xfit, yfit, color='blue')
        #line.label='model'

        #key=PlotKey(0.9,0.9,[hplt,line,hpltrand],halign='right')
        #plt.add(line, hplt, hpltrand, key)
        key=PlotKey(0.9,0.9,[hplt,hpltrand],halign='right')
        plt.add(hplt, hpltrand, key)


        tab[1,0]=plt
        
        if self.show:
            tab.show()

        d=files.get_prior_dir()
        d=os.path.join(d, 'plots')
        epsfile='pofe-%.2f-%.2f-%s.eps' % (minmag,maxmag,self.objtype)
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