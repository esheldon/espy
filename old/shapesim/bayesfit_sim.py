"""
MCMC with bayesian approach
"""
import os
import pprint
import numpy
from numpy import sqrt, cos, sin, exp, log, log10, pi, zeros, ones, empty, \
        random, where, array, linspace, diag, median
from numpy import tanh, arctanh
from numpy.random import randn
from numpy.random import random as randu
from . import shapesim
import lensing
from lensing.shear import Shear
import fimage
from fimage.convolved import NoisyConvolvedImage

import gmix_image
from gmix_image import print_pars, GMix, gmix2pars
from gmix_image.gmix_mcmc import MixMCSimple
from gmix_image.priors import GPriorOld, CenPrior

import images
import esutil as eu
from esutil.random import srandu
from esutil.misc import wlog



import math
import time

from sys import stderr

import scipy.stats

LOWVAL=-9999.9e9

class BayesFitSim(shapesim.BaseSim):
    def __init__(self, run, extra=None):
        """
        use config files 

        can over-ride any values with extra= for quick
        tests. 

        """
        super(BayesFitSim,self).__init__(run)
        if 'verbose' not in self:
            self['verbose'] = False

        if extra:
            for k,v in extra.iteritems():
                print k,v
                self[k] = v

        self.gprior = GPriorOld(A=self.simc['A'],
                             B=self.simc['B'],
                             C=self.simc['C'],
                             D=self.simc['D'])

    def process_trials_by_s2n(self, is2, is2n):
        """
        ring test

        Generates a random total ellipticity from the 
        assumed prior

        Run this many times to sample the prior distribution
        """

        s2n = shapesim.get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']

        nellip = self.get_nellip(is2n)

        nring = self.simc['nring']

        ntot = nring*nellip
        out = zeros(ntot, dtype=self.out_dtype())

        ii = 0
        for i in xrange(nring):
            itheta=i

            dolog=False
            if i==0:
                dolog=True
            st = self.process_trial_by_s2n(is2, is2n, itheta, dolog=dolog)
            out[ii:ii+nellip] = st
            ii += nellip

        shapesim.write_output(self['run'], is2, is2n, out, fs=self.fs)
        return out


    def get_fitmodels(self):
        fitmodels=self['fitmodel']
        if not isinstance(fitmodels,list):
            fitmodels=[fitmodels]
        return fitmodels

    def process_trial_by_s2n(self, is2, is2n, itheta,
                             dowrite=False, 
                             dolog=False):
        """
        Process a singe element in the ring, with nellip
        possible noise realizations
        """

        s2 = linspace(self.simc['mins2'],
                      self.simc['maxs2'], 
                      self.simc['nums2'])[is2]
        s2n_psf = self['s2n_psf']
        s2n = shapesim.get_s2n(self, is2n)
        theta = shapesim.get_theta(self.simc, itheta=itheta)

        s2n_fac = self['s2n_fac']
        s2n_method = self['s2n_method']
        s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

        nellip=self.get_nellip(is2n)

        # these are generated on the same series every itheta for a given run
        # seed so that each itheta gets the same ellip values; otherwise no
        # ring
        gvals = self.get_gvals(is2, is2n, nellip)


        fitmodels=self.get_fitmodels()
        if fitmodels[0] not in ['gexp','gdev','gauss']:
            ngauss=self['ngauss_obj']
            npars=2*ngauss+4
        else:
            npars=6

        out = zeros(nellip, dtype=self.out_dtype(npars))
        out['s2'] = s2
        self['s2']=s2

        if dolog:
            wlog('s2n:',s2n,'s2n_psf:',s2n_psf,'s2n_method:',s2n_method)


        for i,g in enumerate(gvals):

            ellip=lensing.util.g2e(g)

            ci_nonoise = self.shapesim.get_trial(s2, ellip, theta)

            retrim = self['retrim']
            if retrim:
                if 'retrim_fluxfrac' not in self:
                    raise ValueError("you must set fluxfrac for a retrim")
                retrim_fluxfrac=self['retrim_fluxfrac']
                olddims=ci_nonoise.image.shape
                ci_nonoise.trim(fluxfrac=retrim_fluxfrac)
                if self['verbose']:
                    wlog("re-trimming with fluxfrac: %.12g" % retrim_fluxfrac)
                    wlog("old dims:",str(olddims),"new dims:",str(ci_nonoise.image.shape))

            e1true=ci_nonoise['e1true']
            e2true=ci_nonoise['e2true']
            g1true,g2true=lensing.util.e1e2_to_g1g2(e1true,e2true)

            #wlog("g1true:",g1true,"g2true:",g2true)
            out['gtrue'][i,0] = g1true
            out['gtrue'][i,1] = g2true
            out['shear_true'][i,0] = ci_nonoise['shear1']
            out['shear_true'][i,1] = ci_nonoise['shear2']
            
            if self['verbose']:
                stderr.write('-'*70 + '\n')
            # we always write this, although slower when not verbose
            if (nellip > 1) and (( (i+1) % 10) == 0 or i== 0):
                stderr.write("  %s/%s ellip done\n" % ((i+1),nellip))
            #if i==20:
            #    stop
            ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                     s2n_method=s2n_method,
                                     fluxfrac=s2ncalc_fluxfrac)
            if self['verbose']:
                wlog("s2n_admom:",ci['s2n_admom'],"s2n_uw:",ci['s2n_uw'],
                     "s2n_matched:",ci['s2n_matched'])

            out['s2n_admom'][i] = ci['s2n_admom']
            out['s2n_matched'][i] = ci['s2n_matched']
            out['s2n_uw'][i] = ci['s2n_uw']
            out['s2n_admom_psf'][i] = ci['s2n_admom_psf']
            out['s2n_matched_psf'][i] = ci['s2n_matched_psf']
            out['s2n_uw_psf'][i] = ci['s2n_uw_psf']

            res=self._run_models(ci, fitmodels)

            out['model'][i] = res['model']
            out['pars'][i,:] = res['pars']
            out['pcov'][i,:,:] = res['pcov']

            if 'gcov' in res:
                out['g'][i,:] = res['g']
                out['gsens'][i,:] = res['gsens']
                out['gcov'][i,:,:] = res['gcov']
            else:
                out['e'][i,:] = res['e']
                out['ecov'][i,:,:] = res['ecov']
                out['emed'][i,:] = res['emed']

            if 'Ts2n' in res:
                for tn in ['Tmean','Terr','Ts2n']:
                    out[tn][i] = res[tn]

            out['s2n_meas_w'][i] = res['s2n_w']
            out['loglike'][i] = res['loglike']
            out['chi2per'][i] = res['chi2per']
            out['dof'][i] = res['dof']
            out['fit_prob'][i] = res['fit_prob']
            if 'arate' in res:
                out['arate'][i] = res['arate']
            
            if 'when_prior' in self and self['when_prior'] == 'after':
                if self['eta']:
                    print 'doing eta'
                    out['eta0'][i,:] = res['eta0']
                    out['etacov0'][i,:] = res['etacov0']
                else:
                    out['g0'][i,:] = res['g0']
                    out['gcov0'][i,:] = res['gcov0']

        if dowrite:
            shapesim.write_output(self['run'], is2, is2n, out, itrial=itheta,
                         fs=self.fs)
        return out

    def _run_models(self, ci, fitmodels):
        # fit models, keep the one that most looks like random error
        probrand=-9999.
        for fitmodel in fitmodels:
            self._run_fitter(ci, fitmodel)
            res0 = self.fitter.get_result()
            if len(fitmodels) > 1:
                print '  model:',fitmodel,'probrand:',res0['fit_prob']
            if res0['fit_prob'] > probrand:
                res=res0
                probrand=res0['fit_prob']

        if len(fitmodels) > 1:
            print '    best model:',res['model']
        return res

    def _gmix_fit_psf(self, ci):
        import admom
        counts=ci.psf.sum()
        psfres = admom.admom(ci.psf,
                             ci['cen_uw'][0],
                             ci['cen_uw'][1], 
                             sigsky=ci['skysig_psf'],
                             guess=2.,
                             nsub=1)

        ngauss=self['ngauss_psf']
        npars=2*ngauss+4

        if self['ngauss_psf']==1:
            psf=[{'p':1,
                  'row':psfres['row'],
                  'col':psfres['col'],
                  'irr':psfres['Irr'],
                  'irc':psfres['Irc'],
                  'icc':psfres['Icc']}]
        elif ngauss==2:

            Texamp=array([12.6,3.8])
            pexamp=array([0.30, 0.70])

            Tfrac=Texamp/Texamp.sum()
            pfrac=pexamp/pexamp.sum()

            prior=zeros(npars)
            width=zeros(npars) + 100

            Tpsf=psfres['Irr']+psfres['Icc']

            prior[0]=psfres['row']
            prior[1]=psfres['col']
            prior[2]=psfres['e1']
            prior[3]=psfres['e2']
            prior[4:4+2] = Tpsf*Tfrac
            prior[6:6+2] = counts*pfrac

            # randomize
            prior[0] += 0.01*srandu()
            prior[1] += 0.01*srandu()

            e1start=prior[2]
            e2start=prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)

            prior[4:npars] = prior[4:npars]*(1+0.05*srandu(2*ngauss))

            gm = gmix_image.GMixFitCoellip(ci.psf, ci['skysig_psf'],
                                           prior,width,
                                           Tpositive=True)
            psf=gm.get_gmix()
        elif ngauss==3:
            # these are good for guessing, but the final answer is
            # often a bit off from here
            Texamp=array([0.46,5.95,2.52])
            pexamp=array([0.1,0.7,0.22])

            Tfrac=Texamp/Texamp.sum()
            pfrac=pexamp/pexamp.sum()


            prior=zeros(npars)
            width=zeros(npars) + 100

            Tpsf=psfres['Irr']+psfres['Icc']

            prior[0]=psfres['row']
            prior[1]=psfres['col']
            prior[2]=psfres['e1']
            prior[3]=psfres['e2']

            prior[4:4+3] = Tpsf*Tfrac
            prior[7:7+3] = counts*pfrac

            # randomize
            prior[0] += 0.01*srandu()
            prior[1] += 0.01*srandu()
            e1start=prior[2]
            e2start=prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)


            prior[4:npars] = prior[4:npars]*(1+0.05*srandu(2*ngauss))

            gm = gmix_image.GMixFitCoellip(ci.psf, ci['skysig_psf'],
                                           prior,width,
                                           Tpositive=True)
            psf=gm.get_gmix()

        else:
            raise ValueError("bad ngauss_psf: %s" % ngauss)

        return psf


    def _gmix_fit_convolved_simple(self, ci, psf_gmix, Tguess, fitmodel):
        """
        Use the fitter to get a starting guess
        """

        counts=ci.image.sum()

        npars=6
        prior=zeros(npars)
        width=zeros(npars) + 1000

        width[0] = 0.1
        width[1] = 0.1
        width[5]=0.1
        if self['Tprior']:
            # this is intended to be lognormal but I
            # haven't implemented that yet in the fitter
            width[4]=Tguess*self['Twidthfrac']

        ntry=10
        for i in xrange(ntry):
                             
            prior[0]=ci['cen_uw'][0]*(1.+ .1*srandu())
            prior[1]=ci['cen_uw'][1]*(1.+ .1*srandu())
            prior[2],prior[3] = randomize_e1e2(None,None)
            prior[4] = Tguess*(1. + .10*srandu())
            prior[5] = counts*(1. + .01*srandu())

            gm = gmix_image.GMixFitCoellip(ci.image, 
                                           ci['skysig'],
                                           prior,width,
                                           psf=psf_gmix,
                                           model=fitmodel)

            if 0==gm.get_flags():
                pars=gm.get_pars()
                perr=gm.get_perr()
                return pars,perr

        wlog("could not find a fit after %s tries, returning None" % ntry)
        return None,None

    def _gmix_fit_convolved_ngauss(self, ci, psf_gmix, Tguess, ngauss, fitmodel, onedelta=True):
        """
        Use the fitter to get a starting guess
        """
        counts=ci.image.sum()

        npars=2*ngauss+4
        prior=zeros(npars)
        width=zeros(npars) + 10000.

        width[0] = 0.1
        width[1] = 0.1
        
        pfrac=array([0.26,0.55,0.18])

        if onedelta:
            Tfrac=array([0.77, 0.23, 7.7e-08])
            width[4+ngauss-1]= 7.7e-10
        else:
            Tfrac=array([0.77, 0.22, .01])

        Tvals = Tguess*Tfrac
        pvals = counts*pfrac

        ntry=10
        for i in xrange(ntry):
                             
            prior[0]=ci['cen_uw'][0]*(1.+ .1*srandu())
            prior[1]=ci['cen_uw'][1]*(1.+ .1*srandu())
            prior[2],prior[3] = randomize_e1e2(None,None)
            prior[4:4+ngauss] = Tvals*(1. + .10*srandu(ngauss))
            prior[4+ngauss:4+2*ngauss] = pvals*(1. + .10*srandu(ngauss))

            gm = gmix_image.GMixFitCoellip(ci.image, 
                                           ci['skysig'],
                                           prior,width,
                                           psf=psf_gmix,
                                           model=fitmodel)

            if 0==gm.get_flags():
                pars=gm.get_pars()
                perr=gm.get_perr()
                return pars,perr

        wlog("could not find a fit after %s tries, returning simple guess" % ntry)
        return prior,None


    def _run_fitter(self, ci, fitmodel):
        from gmix_image.gmix_em import GMixEMPSF

        psf_ivar=1./ci['skysig_psf']**2
        gmpsf=GMixEMPSF(ci.psf, psf_ivar, self['ngauss_psf'],
                        cen=ci['cen_psf'],
                        maxiter=self['em_maxiter'],tol=self['em_tol'])

        psf_gmix=gmpsf.get_gmix()

        Tguess = ci['Ttrue']*(1. + 0.1*srandu())
        ivar=1./ci['skysig']**2

        """
        e1true=ci['e1true']
        e2true=ci['e2true']
        g1true,g2true=lensing.util.e1e2_to_g1g2(e1true,e2true)
        print 'Ttrue:',ci['Ttrue']
        print 'gtrue:',g1true,g2true
        """

        cenprior=CenPrior(ci['cen'], [0.1]*2)

        raise RuntimeError("convert to use new mixmc")
        self.fitter=MixMCSimple(ci.image, ivar, 
                                psf_gmix, self.gprior, fitmodel,
                                cen=ci['cen'],
                                nwalkers=self['nwalkers'],
                                nstep=self['nstep'], 
                                burnin=self['burnin'],
                                mca_a=self['mca_a'],
                                iter=self.get('iter',False),
                                draw_gprior=self['draw_gprior'])

    def _run_fitter_old(self, ci, fitmodel):
        """
        cheat on psf, T and cen for now
        """

        if 'Ttrue' in ci:
            T=ci['Ttrue']
        else:
            cov=ci['covtrue']
            T = cov[0] + cov[2]
        cen = ci['cen']

        psf_gmix=self._gmix_fit_psf(ci)
        if 'mixmc' in self['run'] and fitmodel=='coellip':
            start_pars,start_err=self._gmix_fit_convolved_ngauss(ci,
                                                                 psf_gmix,
                                                                 T,
                                                                 self['ngauss_obj'],
                                                                 fitmodel,
                                                                 onedelta=self['onedelta'])
        else:
            start_pars,start_err=self._gmix_fit_convolved_simple(ci,psf_gmix,T,fitmodel)

        burnin = self['burnin']
        if start_pars is None:
            # we might have a very bad guess, increase burnin
            burnin=self['burnin']*4

        if 'fixcen' in self and self['fixcen']:
            raise ValueError("fixcen no longer supported")

        elif 'mixmc' in self['run'] and fitmodel=='coellip':
            #print_pars(start_pars,front='pars:')
            #print_pars(start_err,front='perr:')

            cenprior=CenPrior(ci['cen'], [0.1]*2)
            doiter=self.get('iter',False)
            self.fitter=EmceeNGaussFitter(ci.image,
                                          1./ci['skysig']**2,
                                          psf_gmix,
                                          start_pars, # cenprior will over-ride
                                          cenprior,
                                          self.gprior,
                                          nwalkers=self['nwalkers'],
                                          nstep=self['nstep'], 
                                          burnin=burnin,
                                          mca_a=self['mca_a'],
                                          onedelta=self['onedelta'],
                                          iter=doiter)


        else:
            cenprior=CenPrior(ci['cen'], [0.1]*2)
            if self['Tprior']:
                # need to guess this width for real data
                Twidth=T*self['Twidthfrac']
                Tsend = eu.random.LogNormal(T, Twidth)
            else:
                Tsend = T

            self.fitter=MixMCOld(ci.image,
                              1./ci['skysig']**2,
                              psf_gmix,
                              cenprior,
                              Tsend,
                              self.gprior,
                              fitmodel,
                              nwalkers=self['nwalkers'],
                              nstep=self['nstep'], 
                              burnin=burnin,
                              mca_a=self['mca_a'],
                              iter=self.get('iter',False),
                              start_pars=start_pars) # Tprior/cenprior over-ride

            """
            self.fitter=EmceeFitter(ci.image,
                                    1./ci['skysig']**2,
                                    psf_gmix,
                                    cenprior,
                                    Tsend,
                                    self.gprior,
                                    fitmodel,
                                    nwalkers=self['nwalkers'],
                                    nstep=self['nstep'], 
                                    burnin=burnin,
                                    logT=logT,
                                    mca_a=self['mca_a'],
                                    eta=eta,
                                    when_prior=self['when_prior'],
                                    iter=self.get('iter',False),
                                    start_pars=start_pars) # Tprior/cenprior over-ride
            """


    def get_nellip(self, is2n):
        s2n = shapesim.get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']
        nellip = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        if nellip < self['min_gcount']:
            nellip=self['min_gcount']
        return nellip

    def get_gvals(self, is2, is2n, nellip):
        if self['seed'] == None:
            raise ValueError("can't use null seed for bayesfit")

        seed=self['seed']
        allseed= seed*10000 + is2n*100 + is2

        print 'seed,is2n,is2,allseed:',seed,is2n,is2,allseed
        # always use same seed for a given is2/is2n and config seed so we use
        # the same g values at given theta in ring

        numpy.random.seed(allseed)
        gvals = self.gprior.sample1d(nellip)

        # now random seed
        numpy.random.seed(None)

        return gvals


    def out_dtype(self, npars):

        if 'mca' in self['run']:
            dt=[('model','S20'),
                ('s2n_admom','f8'),
                ('s2n_matched','f8'),
                ('s2n_uw','f8'),
                ('s2n_admom_psf','f8'),
                ('s2n_matched_psf','f8'),
                ('s2n_uw_psf','f8'),
                ('s2','f8'),
                ('shear_true','f8',2),
                ('gtrue','f8',2),
                ('e','f8',2),
                ('ecov','f8',(2,2)),
                ('pars','f8',npars),
                ('pcov','f8',(npars,npars)),
                ('emed','f8',2),
                ('s2n_meas_w','f8'),  # weighted s/n based on most likely point
                ('loglike','f8'),     # loglike of fit
                ('chi2per','f8'),     # chi^2/degree of freedom
                ('dof','i4'),         # degrees of freedom
                ('fit_prob','f8'),    # probability of the fit happening randomly
                ('arate','f8'),
               ]

        else:
            dt=[('model','S20'),
                ('s2n_admom','f8'),
                ('s2n_matched','f8'),
                ('s2n_uw','f8'),
                ('s2n_admom_psf','f8'),
                ('s2n_matched_psf','f8'),
                ('s2n_uw_psf','f8'),
                ('s2','f8'),
                ('shear_true','f8',2),
                ('gtrue','f8',2),
                ('g','f8',2),
                ('gsens','f8',2),
                ('gcov','f8',(2,2)),
                ('pars','f8',npars),
                ('pcov','f8',(npars,npars)),
                ('Tmean','f8','f8'),
                ('Terr','f8','f8'),
                ('Ts2n','f8','f8'),
                ('s2n_meas_w','f8'),  # weighted s/n based on most likely point
                ('loglike','f8'),     # loglike of fit
                ('chi2per','f8'),     # chi^2/degree of freedom
                ('dof','i4'),         # degrees of freedom
                ('fit_prob','f8')     # probability of the fit happening randomly
               ]
            if 'bayesfit' not in self['run']:
                dt += [('arate','f8')]

            when_prior=self.get('when_prior','during')
            if when_prior == 'after':
                if self['eta']:
                    dt += [('eta0','f8',2),
                           ('etacov0','f8',(2,2))]
                else:
                    dt += [('g0','f8',2),
                           ('gcov0','f8',(2,2))]
        return dt




class EmceeNGaussFitter:
    def __init__(self, 
                 image, 
                 ivar, 
                 psf,
                 start_pars,
                 cenprior,
                 gprior,
                 onedelta=True, # one is a delta function
                 nwalkers=20,
                 nstep=200, 
                 burnin=400,
                 mca_a=2.0,
                 iter=True):  # tprior,cenprior take precedence
        """
        mcmc sampling of posterior.

        parameters
        ----------
        image:
            sky subtracted image as a numpy array
        ivar:
            1/(Error per pixel)**2
        psf:
            The psf gaussian mixture
        cenprior:
            The center prior object.
        T:
            Starting value for ixx+iyy of main component
            or a LogNormal object
        gprior:
            The prior on the g1,g2 surface.
        nstep:
            Number of steps in MCMC chain.
        burnin:
            Number of burn in steps.
        """
        
        self.make_plots=False

        # cen1,cen2,e1,e2,Ti,pi
        self.npars=len(start_pars)
        self.start_pars=array(start_pars) # makes a copy
        self.ngauss=(len(start_pars)-4)/2

        self.image=image
        self.ivar=float(ivar)
        self.cenprior=cenprior
        self.delta_prior=None

        self.mca_a=mca_a

        self.nwalkers=nwalkers
        self.nstep=nstep
        self.burnin=burnin
        self.gprior=gprior
        self.iter=iter

        self.onedelta=onedelta
        if self.onedelta:
            self._set_onedelta()

        self._set_psf(psf)

        self.tpars=zeros(self.npars,dtype='f8')

        self._go()


    def get_result(self):
        return self._result

    def get_maxlike_model(self):
        """
        Get the model representing the maximum likelihood point in the chain
        Is this useful?
        """
        w=self.lnprobs.argmax()
        pars = self.trials[w,:].copy()
        e1,e2,ok=g1g2_to_e1e2(pars[2],pars[3])
        if not ok:
            raise ValueError("bad e1,e2")

        pars[2],pars[3]=e1,e2
        gmix=self._get_convolved_gmix(pars)
        model=gmix_image.gmix2image(gmix, self.image.shape)
        return model


    def _go(self):
        import emcee

        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._calc_lnprob,
                                        a=self.mca_a)
        
        guess=self._get_guess()

        if self.iter:
            pos, prob, state = sampler.run_mcmc(guess, self.burnin)
            sampler.reset()
            while True:
                pos, prob, state = sampler.run_mcmc(pos, self.nstep)
                try:
                    acor=sampler.acor
                    tau = (sampler.acor/self.nstep).max()
                    if tau > 0.1:
                        wlog("tau",tau,"greater than 0.1")
                    else:
                        break
                except:
                    # something went wrong with acor, run some more
                    print 'error'
                    pass

        else:
            pos, prob, state = sampler.run_mcmc(guess, self.burnin)
            sampler.reset()
            pos, prob, state = sampler.run_mcmc(pos, self.nstep)

        self.trials  = sampler.flatchain

        lnprobs = sampler.lnprobability.reshape(self.nwalkers*self.nstep)
        self.lnprobs = lnprobs - lnprobs.max()

        self._emcee_sampler=sampler

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if self.make_plots:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop


    def _calc_lnprob(self, pars):
        """
        pars are [cen1,cen2,g1,g2,Ti,pi]
        """

        g1,g2=pars[2],pars[3]

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if not self._check_pvals_Tvals(pars):
            return LOWVAL

        self.tpars[0] = pars[0]
        self.tpars[1] = pars[1]
        self.tpars[2] = e1
        self.tpars[3] = e2
        self.tpars[4:] = pars[4:]

        gmix=self._get_convolved_gmix(self.tpars)
        T=gmix.get_T()
        if T < 0:
            return LOWVAL

        logprob = self._get_loglike_c(gmix)

        gp = self._get_lngprior(g1,g2)
        logprob += gp

        cp = self.cenprior.lnprob(pars[0:2])
        logprob += cp

        if self.onedelta:
            val=pars[self.delta_index]
            dprior = self.delta_prior.lnprob(val)
            logprob += dprior

        return logprob

    def _check_pvals_Tvals(self,pars):
        Tvals=pars[4:4+self.ngauss]
        pvals=pars[4+self.ngauss:]

        w,=where(Tvals <= 0)
        if w.size > 0:
            return False
 
        w,=where(pvals <= 0)
        if w.size > 0:
            return False

        return True
 
    def _get_loglike_c(self, gmix):
        loglike,s2n,flags=\
            gmix_image.render._render.loglike(self.image, 
                                              gmix,
                                              self.ivar)

        if flags != 0:
            return LOWVAL
        return loglike

    def _get_convolved_gmix(self,pars):
        gmix0=gmix_image.GMixCoellip(pars)
        gmix=gmix0.convolve(self.psf_gmix)
        return gmix


    def _get_lngprior(self, g1, g2):
        g=sqrt(g1**2 + g2**2)
        gp = self.gprior.prior2d_gabs_scalar(g)
        if gp > 0:
            gp = log(gp)
        else:
            gp=LOWVAL
        return gp

    def _calc_result(self):
        """
        We marginalize over all parameters but g1,g2, which
        are index 0 and 1 in the pars array
        """
        import mcmc

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g1vals=self.trials[:,2]
        g2vals=self.trials[:,3]

        prior = self.gprior(g1vals,g2vals)
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        g0=None
        gcov0=None

        # prior is already in the distribution of points.  This is simpler for
        # most things but for sensitivity we need a factor of (1/P)dP/de

        pars,pcov = mcmc.extract_stats(self.trials)

        #print_pars(pars,front='pars:')
        #print_pars(sqrt(diag(pcov)),front='perr:')

        g[:] = pars[2:4]
        gcov[:,:] = pcov[2:4, 2:4]

        g1diff = g[0]-g1vals
        g2diff = g[1]-g2vals

        w,=where(prior > 0)
        if w.size == 0:
            raise ValueError("no prior values > 0!")

        gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
        gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()


        # generate p1*T1 + p2*T2 + ...
        Tsums=zeros(self.trials.shape[0])
        psums=Tsums.copy()
        for i in xrange(self.ngauss):
            Tvals = self.trials[:,4+i]
            pvals = self.trials[:,4+self.ngauss+i]

            Tsums += Tvals*pvals
            psums += pvals
 
        self.Ttots = Tsums/psums
        self.ptots = psums

        Tmean = self.Ttots.mean()
        Terr = self.Ttots.std()
        Ts2n = Tmean/Terr

        #print '%g +/- %g   s/n: %g' % (Tmean,Terr,Ts2n)

        arates = self._emcee_sampler.acceptance_fraction
        arate = arates.mean()

        # weighted s/n based on the most likely point
        s2n,loglike,chi2per,dof,prob=self._calculate_maxlike_stats()
        self._result={'model':self.model,
                      'g':g,
                      'gcov':gcov,
                      'gsens':gsens,
                      'g0':g0,
                      'gcov0':gcov0,
                      'pars':pars,
                      'pcov':pcov,
                      'Tmean':Tmean,
                      'Terr':Terr,
                      'Ts2n':Ts2n,
                      'arate':arate,
                      's2n_w':s2n,
                      'loglike':loglike,
                      'chi2per':chi2per,
                      'dof':dof,
                      'fit_prob':prob}

    def _calculate_maxlike_stats(self):
        """
        Stats Based on the most likely point
        """

        w=self.lnprobs.argmax()
        pars = self.trials[w,:].copy()
        e1,e2,ok=g1g2_to_e1e2(pars[2],pars[3])
        if not ok:
            raise ValueError("bad e1,e2")
        pars[2],pars[3]=e1,e2

        gmix=self._get_convolved_gmix(pars)

        loglike,s2n,flags=\
            gmix_image.render._render.loglike(self.image, 
                                              gmix,
                                              self.ivar)
        chi2=loglike/(-0.5)
        dof=self.image.size-pars.size
        chi2per = chi2/dof

        prob = scipy.stats.chisqprob(chi2, dof)
        return s2n, loglike, chi2per, dof, prob


    def _get_guess(self):
        nwalkers=self.nwalkers
        npars=self.npars
        guess=zeros( (nwalkers,npars) )

        guess[:,0]=self.cenprior.cen[0] + 0.01*srandu(nwalkers)
        guess[:,1]=self.cenprior.cen[1] + 0.01*srandu(nwalkers)

        sh=Shear(e1=self.start_pars[2],e2=self.start_pars[3])
        for i in xrange(nwalkers):
            g1s,g2s=randomize_e1e2(sh.g1,sh.g2,width=0.05)
            guess[i,2] = g1s
            guess[i,3] = g2s

        for i in xrange(4,npars):
            guess[:,i] = self.start_pars[i]*(1. + 0.1*srandu(nwalkers))

        return guess

    def _set_onedelta(self):
        # place delta function on smallest from start_pars
        Tvals=self.start_pars[4:4+self.ngauss]
        argmin=Tvals.argmin()

        self.delta_index=4+argmin

        delta_mn=Tvals[argmin]

        if delta_mn < 0:
            delta_mn=1.e-7
            self.start_pars[self.delta_index] = delta_mn

        delta_width = delta_mn*0.1
        self.delta_prior = eu.random.LogNormal(delta_mn,delta_width)


    def _set_psf(self, psf):
        if psf is not None:
            self.psf_gmix = GMix(psf)
            self.psf_pars = gmix2pars(self.psf_gmix)
        else:
            self.psf_gmix = None
            self.psf_pars = None


    def _doplots(self):
        import mcmc
        import biggles

        biggles.configure("default","fontsize_min",1.2)

        plot_logT=False
        tab=biggles.Table(6,2)


        cen1vals=self.trials[:,0]
        cen2vals=self.trials[:,1]
        Tvals=self.Ttots
        g1vals=self.trials[:,2]
        g2vals=self.trials[:,3]
        g1lab=r'$g_1$'
        g2lab=r'$g_2$'

        ampvals=self.ptots

        ind=numpy.arange(g1vals.size)

        burn_cen=biggles.FramedPlot()
        cen1p=biggles.Curve(ind, cen1vals, color='blue')
        cen2p=biggles.Curve(ind, cen2vals, color='red')
        cen1p.label=r'$x_1$'
        cen2p.label=r'$x_2$'
        burn_cen.add(cen1p)
        burn_cen.add(cen2p)
        key=biggles.PlotKey(0.9,0.9,[cen1p,cen2p],halign='right')
        burn_cen.add(key)
        burn_cen.ylabel='cen'

        burn_g1=biggles.FramedPlot()
        burn_g1.add(biggles.Curve(ind, g1vals))
        burn_g1.ylabel=r'$g_1$'

        burn_g2=biggles.FramedPlot()
        burn_g2.add(biggles.Curve(ind, g2vals))
        burn_g2.ylabel=r'$g_2$'

        burn_T=biggles.FramedPlot()
        burn_T.add(biggles.Curve(ind, Tvals))
        burn_T.ylabel='T'

        burn_amp=biggles.FramedPlot()
        burn_amp.add(biggles.Curve(ind, ampvals))
        burn_amp.ylabel='Amplitide'



        likep = biggles.FramedPlot()
        likep.add( biggles.Curve(ind, self.lnprobs) )
        likep.ylabel='ln( prob )'


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))

        res=self.get_result()
        print 's2n weighted maxlike:',res['s2n_w']
        print 'chi^2/dof: %.3f/%i = %f' % (res['chi2per']*res['dof'],res['dof'],res['chi2per'])
        print 'prob:',res['fit_prob']
        print 'acceptance rate:',res['arate']
        print 'T:  %.16g +/- %.16g' % (Tvals.mean(), Tvals.std())

        print_pars(self._result['pars'])
        print_pars(sqrt(diag(self._result['pcov'])))
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'median g1:  %.16g ' % median(g1vals)
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]

        if self._result['g0'] is not None:
            g0=self._result['g0']
            err0=sqrt(diag(self._result['gcov0']))
            print 'g1_0: %.16g +/- %.16g' % (g0[0],err0[0])
            print 'g2_0: %.16g +/- %.16g' % (g0[1],err0[1])

        cenw = cen1vals.std()
        cen_bsize=cenw*0.2
        hplt_cen0 = eu.plotting.bhist(cen1vals,binsize=cen_bsize,
                                      color='blue',
                                      show=False)
        hplt_cen = eu.plotting.bhist(cen2vals,binsize=cen_bsize,
                                     color='red',
                                     show=False, plt=hplt_cen0)
        hplt_cen.add(key)

        bsize1=g1vals.std()*0.2 #errs[0]*0.2
        bsize2=g2vals.std()*0.2 # errs[1]*0.2
        hplt_g1 = eu.plotting.bhist(g1vals,binsize=bsize1,
                                  show=False)
        hplt_g2 = eu.plotting.bhist(g2vals,binsize=bsize2,
                                  show=False)

        Tsdev = Tvals.std()
        Tbsize=Tsdev*0.2
        #hplt_T = eu.plotting.bhist(Tvals,binsize=Tbsize,
        #                          show=False)

        if plot_logT:
            Tvals2plot=log10(Tvals)
            hTlab=r'$log_{10}T$'
        else:
            Tvals2plot=Tvals
            hTlab='T'
        Tsdev = Tvals2plot.std()
        Tbsize=Tsdev*0.2
        hplt_T = eu.plotting.bhist(Tvals2plot,binsize=Tbsize,
                                   show=False)



        amp_sdev = ampvals.std()
        amp_bsize=amp_sdev*0.2
        hplt_amp = eu.plotting.bhist(ampvals,binsize=amp_bsize,
                                     show=False)



        hplt_cen.xlabel='center'
        hplt_g1.xlabel=g1lab
        hplt_g2.xlabel=g2lab
        hplt_T.xlabel=hTlab
        hplt_amp.xlabel='Amplitude'

        tab[0,0] = burn_cen
        tab[1,0] = burn_g1
        tab[2,0] = burn_g2
        tab[3,0] = burn_T
        tab[4,0] = burn_amp

        tab[0,1] = hplt_cen
        tab[1,1] = hplt_g1
        tab[2,1] = hplt_g2
        tab[3,1] = hplt_T
        tab[4,1] = hplt_amp
        tab[5,0] = likep
        tab.show()

        if False: 
            nx = ny = 40
            levels=8
            h2d = eu.stat.histogram2d(Tvals, g1vals, nx=nx, ny=ny,more=True)
            images.view(h2d['hist'], type='cont',
                        xdr=[h2d['xcenter'][0], h2d['xcenter'][-1]],
                        ydr=[h2d['ycenter'][0], h2d['ycenter'][-1]],
                        xlabel='T', ylabel='g1', levels=levels)
            h2d = eu.stat.histogram2d(Tvals, g2vals, nx=nx, ny=ny,more=True)
            images.view(h2d['hist'], type='cont',
                        xdr=[h2d['xcenter'][0], h2d['xcenter'][-1]],
                        ydr=[h2d['ycenter'][0], h2d['ycenter'][-1]],
                        xlabel='T', ylabel='g2', levels=levels)



class EmceeFitter:
    def __init__(self, 
                 image, 
                 ivar, 
                 psf,
                 cenprior,
                 T,
                 gprior,
                 model,
                 nwalkers=10,
                 nstep=100, 
                 burnin=400,
                 logT=False,
                 eta=False,
                 mca_a=2.0,
                 when_prior='during', 
                 iter=False,
                 start_pars=None):  # tprior,cenprior take precedence
        """
        mcmc sampling of posterior.

        parameters
        ----------
        image:
            sky subtracted image as a numpy array
        ivar:
            1/(Error per pixel)**2
        psf:
            The psf gaussian mixture
        cenprior:
            The center prior object.
        T:
            Starting value for ixx+iyy of main component
            or a LogNormal object
        gprior:
            The prior on the g1,g2 surface.
        nstep:
            Number of steps in MCMC chain.
        burnin:
            Number of burn in steps.
        when_prior:
            'during' or 'after'
        """
        
        self.make_plots=False

        # cen1,cen2,e1,e2,T,p
        self.npars=6

        self.image=image
        self.ivar=float(ivar)
        self.T=T
        self.cenprior=cenprior

        self.logT=logT
        self.eta=eta
        self.mca_a=mca_a

        self.model=model

        self.nwalkers=nwalkers
        self.nstep=nstep
        self.burnin=burnin
        self.gprior=gprior
        self.iter=iter
        self.when_prior=when_prior

        self.start_pars=start_pars

        self._set_psf(psf)

        if isinstance(T, eu.random.LogNormal):
            self.T_is_prior=True
        else:
            self.T_is_prior=False

        self.tpars=zeros(6,dtype='f8')

        self._go()

    def get_result(self):
        return self._result

    def get_maxlike_model(self):
        """
        Get the model representing the maximum likelihood point in the chain
        Is this useful?
        """
        w=self.lnprobs.argmax()
        pars = self.trials[w,:].copy()
        if self.logT:
            pars[4] = 10.0**(pars[4])
        if self.eta:
            e1,e2=eta1eta2_to_e1e2(pars[2],pars[3])
            if not ok:
                raise ValueError("bad e1,e2")
        else:
            e1,e2,ok=g1g2_to_e1e2(pars[2],pars[3])
            if not ok:
                raise ValueError("bad e1,e2")

        pars[2],pars[3]=e1,e2
        gmix=self._get_convolved_gmix(pars)
        model=gmix_image.gmix2image(gmix, self.image.shape)
        return model


    def _go(self):
        import emcee

        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._calc_lnprob,
                                        a=self.mca_a)
        
        guess=self._get_guess()

        if self.iter:
            pos, prob, state = sampler.run_mcmc(guess, self.burnin)
            sampler.reset()
            while True:
                pos, prob, state = sampler.run_mcmc(pos, self.nstep)
                try:
                    acor=sampler.acor
                    tau = (sampler.acor/self.burnin).max()
                    if tau > 0.1:
                        wlog("tau",tau,"greater than 0.1")
                    else:
                        break
                except:
                    # something went wrong with acor, run some more
                    pass

        else:
            pos, prob, state = sampler.run_mcmc(guess, self.burnin)
            sampler.reset()
            pos, prob, state = sampler.run_mcmc(pos, self.nstep)

        self.trials  = sampler.flatchain

        lnprobs = sampler.lnprobability.reshape(self.nwalkers*self.nstep)
        self.lnprobs = lnprobs - lnprobs.max()

        self._emcee_sampler=sampler

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if self.make_plots:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop


    def _calc_lnprob(self, pars):
        """
        pars are [g1,g2,T]

        if using logT we convert to linear here
        """
        # hard priors first

        if self.logT:
            #T=exp(pars[4])
            T=10.0**(pars[4])
        else:
            T=pars[4]

        if self.eta:
            """
            try:
                s=Shear(eta1=pars[2],eta2=pars[3])
                e1,e2=s.e1,s.e2
                ok=True
            except:
                ok=False
            """
            g1,g2,ok = eta1eta2_to_g1g2(pars[2],pars[3])
            if not ok:
                return LOWVAL
        else:
            g1,g2=pars[2],pars[3]

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0] = pars[0]
        self.tpars[1] = pars[1]
        self.tpars[2] = e1
        self.tpars[3] = e2
        self.tpars[4] = T
        self.tpars[5:] = pars[5:]

        logprob = self._get_loglike_c(self.tpars)

        if self.when_prior=='during':
            gp = self._get_lngprior(g1,g2)
            logprob += gp

        cp = self.cenprior.lnprob(pars[0:2])
        logprob += cp

        if self.T_is_prior:
            Tp = self.T.lnprob(T)
            logprob += Tp

        return logprob

    def _get_convolved_gmix(self,pars):
        """
        This should have T linear
        """
        if self.model=='gexp':
            gmix0=gmix_image.GMixExp(pars)
        elif self.model=='gdev':
            gmix0=gmix_image.GMixDev(pars)
        elif self.model=='gauss':
            gmix0=gmix_image.GMixCoellip(pars)
        else:
            raise ValueError("bad model: '%s'" % self.model)
        gmix=gmix0.convolve(self.psf_gmix)
        return gmix

 
    def _get_loglike_c(self, pars):
        gmix=self._get_convolved_gmix(pars)

        loglike,s2n,flags=\
            gmix_image.render._render.loglike(self.image, 
                                              gmix,
                                              self.ivar)

        if flags != 0:
            return LOWVAL
        return loglike

    def _get_lngprior(self, g1, g2):
        g=sqrt(g1**2 + g2**2)
        gp = self.gprior.prior2d_gabs_scalar(g)
        if gp > 0:
            gp = log(gp)
        else:
            gp=LOWVAL
        return gp

    def _calc_result(self):
        """
        if when_prior='after', We apply the shape prior here

        We marginalize over all parameters but g1,g2, which
        are index 0 and 1 in the pars array
        """
        import mcmc

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        if self.eta:
            eta1vals,eta2vals=self.trials[:,2], self.trials[:,3]
            g1vals,g2vals=self._get_gvals_from_eta(eta1vals,eta2vals)
        else:
            g1vals=self.trials[:,2]
            g2vals=self.trials[:,3]

        prior = self.gprior(g1vals,g2vals)
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        g0=None
        gcov0=None
        if self.when_prior=='after':
            # this will be eta not g if using that parametrization
            g0,gcov0 = mcmc.extract_stats(self.trials[:,2:4])

            pars, pcov = mcmc.extract_stats(self.trials,weights=prior)
            # we need to multiply each by the prior
            """
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum


            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum
            """

            g[:] = pars[2:4]
            gcov[:,:] = pcov[2:4, 2:4]

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals
            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de

            pars,pcov = mcmc.extract_stats(self.trials)
            #g, gcov = mcmc.extract_stats(self.trials[:,2:2+2])

            g[:] = pars[2:4]
            gcov[:,:] = pcov[2:4, 2:4]

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
            gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()

 
        arates = self._emcee_sampler.acceptance_fraction
        arate = arates.mean()
        #print 'acceptance rate:',w.size/float(self.trials.size)

        # weighted s/n based on the most likely point
        s2n,loglike,chi2per,dof,prob=self._calculate_maxlike_stats()
        if self.eta:
            g0name='eta'
        else:
            g0name='g'

        Tmean=pars[4]
        Terr=sqrt(pcov[4,4])
        Ts2n=pars[4]/sqrt(pcov[4,4])
        #print 'T s/n:',Ts2n
        self._result={'model':self.model,
                      'g':g,
                      'gcov':gcov,
                      'gsens':gsens,
                      g0name+'0':g0,
                      g0name+'cov0':gcov0,
                      'pars':pars,
                      'pcov':pcov,
                      'Tmean':Tmean,
                      'Terr':Terr,
                      'Ts2n':Ts2n,
                      'arate':arate,
                      's2n_w':s2n,
                      'loglike':loglike,
                      'chi2per':chi2per,
                      'dof':dof,
                      'fit_prob':prob}
        #wlog("arate:",self._result['arate'])

    def _calculate_maxlike_stats(self):
        """
        Stats Based on the most likely point
        """

        w=self.lnprobs.argmax()
        pars = self.trials[w,:].copy()
        if self.logT:
            pars[4] = 10.0**(pars[4])
        if self.eta:
            e1,e2,ok=eta1eta2_to_e1e2(pars[2],pars[3])
            if not ok:
                raise ValueError("bad e1,e2")
        else:
            e1,e2,ok=g1g2_to_e1e2(pars[2],pars[3])
            if not ok:
                raise ValueError("bad e1,e2")
        pars[2],pars[3]=e1,e2

        gmix=self._get_convolved_gmix(pars)

        loglike,s2n,flags=\
            gmix_image.render._render.loglike(self.image, 
                                              gmix,
                                              self.ivar)
        chi2=loglike/(-0.5)
        dof=self.image.size-pars.size
        chi2per = chi2/dof

        prob = scipy.stats.chisqprob(chi2, dof)
        return s2n, loglike, chi2per, dof, prob


    def _get_gvals_from_eta(self, eta1, eta2):
        g1vals=zeros(eta1.size)
        g2vals=zeros(eta2.size)
        for i in xrange(eta1.size):
            g1,g2,ok=eta1eta2_to_g1g2(eta1[i], eta2[i])
            if not ok:
                raise ValueError("bad g1,g2")
            g1vals[i]=g1
            g2vals[i]=g2
        return g1vals,g2vals

    def _get_guess(self):
        guess=zeros( (self.nwalkers,self.npars) )

        guess[:,0]=self.cenprior.cen[0] + 0.01*(randu(self.nwalkers)-0.5)
        guess[:,1]=self.cenprior.cen[1] + 0.01*(randu(self.nwalkers)-0.5)

        if self.start_pars is not None:
            sh=Shear(e1=self.start_pars[2],e2=self.start_pars[3])
            for i in xrange(self.nwalkers):
                e1s,e2s=randomize_e1e2(sh.g1,sh.g2,width=0.05)
                guess[i,2] = e1s
                guess[i,3] = e2s
        else:
            # (0,0) with some scatter
            guess[:,2]=0.1*(randu(self.nwalkers)-0.5)
            guess[:,3]=0.1*(randu(self.nwalkers)-0.5)

        # guess for T is self.T with scatter
        if self.T_is_prior:
            T=self.T.mean
            guess[:,4] = T + T*0.1*(randu(self.nwalkers)-0.5)
        else:
            if self.logT:
                if self.T < 0.01:
                    raise ValueError("guess for T must be > 0.01")
                #logT=log(self.T)
                logT=log10(self.T)
                # rand range +/- 0.005
                guess[:,4] = logT + .01*(randu(self.nwalkers)-0.5)
            else:
                T=self.T
                guess[:,4] = T + T*0.1*(randu(self.nwalkers)-0.5)

        # first guess at amp is the total flux
        if self.start_pars is not None:
            pguess=self.start_pars[5]
            guess[:,5] = pguess*(1+0.1*(randu(self.nwalkers)-0.5))
        else:
            imtot=self.image.sum()
            guess[:,5] = imtot + imtot*0.1*(randu(self.nwalkers)-0.5)

        return guess

    def _set_psf(self, psf):
        if psf is not None:
            self.psf_gmix = GMix(psf)
            self.psf_pars = gmix2pars(self.psf_gmix)
        else:
            self.psf_gmix = None
            self.psf_pars = None


    def _doplots(self):

        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(6,2)

        cen1vals=self.trials[:,0]
        cen2vals=self.trials[:,1]
        Tvals=self.trials[:,4]
        if self.eta:
            eta1vals=self.trials[:,2]
            eta2vals=self.trials[:,3]
            g1vals,g2vals=self._get_gvals_from_eta(eta1vals,eta2vals)

            g1lab=r'$\eta_1$'
            g2lab=r'$\eta_2$'
        else:
            g1vals=self.trials[:,2]
            g2vals=self.trials[:,3]
            g1lab=r'$g_1$'
            g2lab=r'$g_2$'

        ampvals=self.trials[:,5]

        ind=numpy.arange(g1vals.size)

        burn_cen=biggles.FramedPlot()
        cen1p=biggles.Curve(ind, cen1vals, color='blue')
        cen2p=biggles.Curve(ind, cen2vals, color='red')
        cen1p.label=r'$x_1$'
        cen2p.label=r'$x_2$'
        burn_cen.add(cen1p)
        burn_cen.add(cen2p)
        key=biggles.PlotKey(0.9,0.9,[cen1p,cen2p],halign='right')
        burn_cen.add(key)
        burn_cen.ylabel='cen'

        burn_g1=biggles.FramedPlot()
        burn_g1.add(biggles.Curve(ind, g1vals))
        burn_g1.ylabel=r'$g_1$'

        burn_g2=biggles.FramedPlot()
        burn_g2.add(biggles.Curve(ind, g2vals))
        burn_g2.ylabel=r'$g_2$'

        burn_T=biggles.FramedPlot()
        burn_T.add(biggles.Curve(ind, Tvals))
        burn_T.ylabel='T'

        burn_amp=biggles.FramedPlot()
        burn_amp.add(biggles.Curve(ind, ampvals))
        burn_amp.ylabel='Amplitide'



        likep = biggles.FramedPlot()
        likep.add( biggles.Curve(ind, self.lnprobs) )
        likep.ylabel='ln( prob )'


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))

        res=self.get_result()
        print 's2n weighted maxlike:',res['s2n_w']
        print 'chi^2/dof: %.3f/%i = %f' % (res['chi2per']*res['dof'],res['dof'],res['chi2per'])
        print 'prob:',res['fit_prob']
        print 'acceptance rate:',res['arate']
        if self.logT:
            tmp=10.0**(Tvals)
            print 'T:  %.16g +/- %.16g' % (tmp.mean(), tmp.std())
        else:
            print 'T:  %.16g +/- %.16g' % (Tvals.mean(), Tvals.std())

        print_pars(self._result['pars'])
        print_pars(sqrt(diag(self._result['pcov'])))
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'median g1:  %.16g ' % median(g1vals)
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]

        if self.eta:
            print 'eta1: %.16g +/- %.16g' % (eta1vals.mean(),eta1vals.std())
            print 'eta2: %.16g +/- %.16g' % (eta2vals.mean(),eta2vals.std())
        elif self._result['g0'] is not None:
            g0=self._result['g0']
            err0=sqrt(diag(self._result['gcov0']))
            print 'g1_0: %.16g +/- %.16g' % (g0[0],err0[0])
            print 'g2_0: %.16g +/- %.16g' % (g0[1],err0[1])

        cenw = cen1vals.std()
        cen_bsize=cenw*0.2
        hplt_cen0 = eu.plotting.bhist(cen1vals,binsize=cen_bsize,
                                      color='blue',
                                      show=False)
        hplt_cen = eu.plotting.bhist(cen2vals,binsize=cen_bsize,
                                     color='red',
                                     show=False, plt=hplt_cen0)
        hplt_cen.add(key)

        if self.eta:
            bsize1=eta1vals.std()*0.2 #errs[0]*0.2
            bsize2=eta2vals.std()*0.2 # errs[1]*0.2
            hplt_g1 = eu.plotting.bhist(eta1vals,binsize=bsize1,
                                      show=False)
            hplt_g2 = eu.plotting.bhist(eta2vals,binsize=bsize2,
                                      show=False)

        else:
            bsize1=g1vals.std()*0.2 #errs[0]*0.2
            bsize2=g2vals.std()*0.2 # errs[1]*0.2
            hplt_g1 = eu.plotting.bhist(g1vals,binsize=bsize1,
                                      show=False)
            hplt_g2 = eu.plotting.bhist(g2vals,binsize=bsize2,
                                      show=False)

        Tsdev = Tvals.std()
        Tbsize=Tsdev*0.2
        #hplt_T = eu.plotting.bhist(Tvals,binsize=Tbsize,
        #                          show=False)

        if self.logT:
            Tsdev = Tvals.std()
            Tbsize=Tsdev*0.2
            hplt_T = eu.plotting.bhist(Tvals,binsize=Tbsize,
                                       show=False)
        else:
            logTvals=log10(Tvals)
            Tsdev = logTvals.std()
            Tbsize=Tsdev*0.2
            hplt_T = eu.plotting.bhist(logTvals,binsize=Tbsize,
                                       show=False)



        amp_sdev = ampvals.std()
        amp_bsize=amp_sdev*0.2
        hplt_amp = eu.plotting.bhist(ampvals,binsize=amp_bsize,
                                     show=False)



        hplt_cen.xlabel='center'
        hplt_g1.xlabel=g1lab
        hplt_g2.xlabel=g2lab
        hplt_T.xlabel=r'$log_{10}T$'
        hplt_amp.xlabel='Amplitude'

        tab[0,0] = burn_cen
        tab[1,0] = burn_g1
        tab[2,0] = burn_g2
        tab[3,0] = burn_T
        tab[4,0] = burn_amp

        tab[0,1] = hplt_cen
        tab[1,1] = hplt_g1
        tab[2,1] = hplt_g2
        tab[3,1] = hplt_T
        tab[4,1] = hplt_amp
        tab[5,0] = likep
        tab.show()

        if False: 
            nx = ny = 40
            levels=8
            h2d = eu.stat.histogram2d(Tvals, g1vals, nx=nx, ny=ny,more=True)
            images.view(h2d['hist'], type='cont',
                        xdr=[h2d['xcenter'][0], h2d['xcenter'][-1]],
                        ydr=[h2d['ycenter'][0], h2d['ycenter'][-1]],
                        xlabel='T', ylabel='g1', levels=levels)
            h2d = eu.stat.histogram2d(Tvals, g2vals, nx=nx, ny=ny,more=True)
            images.view(h2d['hist'], type='cont',
                        xdr=[h2d['xcenter'][0], h2d['xcenter'][-1]],
                        ydr=[h2d['ycenter'][0], h2d['ycenter'][-1]],
                        xlabel='T', ylabel='g2', levels=levels)







def g1g2_to_e1e2(g1, g2):
    """
    This version without exceptions

    returns e1,e2,okflag
    """
    g = math.sqrt(g1**2 + g2**2)
    if g >= 1.:
        return LOWVAL,LOWVAL,False

    if g == 0:
        return 0.,0.,True
    e = math.tanh(2*math.atanh(g))
    if e >= 1.:
        return LOWVAL,LOWVAL,False

    fac = e/g
    e1, e2 = fac*g1, fac*g2
    return e1,e2,True

def eta1eta2_to_g1g2(eta1, eta2):
    """
    This version without exceptions

    returns e1,e2,okflag
    """
    eta = math.sqrt(eta1**2 + eta2**2)

    if eta == 0:
        return 0.,0.,True
    gtot = math.tanh(eta/2.)
    if gtot >= 1.:
        return LOWVAL,LOWVAL,False

    fac = gtot/eta
    g1, g2 = fac*eta1, fac*eta2
    return g1,g2,True

def eta1eta2_to_e1e2(eta1, eta2):
    """
    This version without exceptions

    returns e1,e2,okflag
    """
    eta = math.sqrt(eta1**2 + eta2**2)

    if eta == 0:
        return 0.,0.,True
    etot = math.tanh(eta)
    if etot >= 1.:
        return LOWVAL,LOWVAL,False

    fac = etot/eta
    e1, e2 = fac*eta1, fac*eta2
    return e1,e2,True

def check():
    import time
    import lensing.shear
    eta1=2.3
    eta2=-1.78

    n=10000
    t1=time.time()
    for i in xrange(n):
        s=lensing.shear.Shear(eta1=eta1,eta2=eta2)
        g1,g2 = s.g1,s.g2
    print time.time()-t1

    print g1,g2
    print s.e1,s.e2
    
    t1=time.time()
    for i in xrange(n):
        g1,g2,ok=eta1eta2_to_g1g2(eta1,eta2)
    print time.time()-t1

    print g1,g2
    e1,e2,ok=g1g2_to_e1e2(g1,g2)
    print e1,e2



def test(n_ggrid=19, n=1000, s2n=40, gmin=-.9, gmax=.9, show=False, clobber=False):
    """
    """
    import fitsio
    outfile=os.path.expanduser('~/tmp/test-n-ggrid%d-%06d-s2n%d.fits' % (n_ggrid,n,s2n))
    if os.path.exists(outfile) and clobber:
        os.remove(outfile)

    extra={'n_ggrid':n_ggrid,
           's2nvals':[s2n],
           'gmin':gmin,
           'gmax':gmax}
    s=BayesFitSim('bayesfit-gg06r01', extra=extra)
    print outfile

    if not os.path.exists(outfile):
        data=[]
        for i in xrange(n):
            print '-'*70
            print '%d/%d' % (i+1,n)
            out=s.process_trials_by_s2n(1, 0)
            data.append(out)

        # combine into one big struct
        alldata=eu.numpy_util.combine_arrlist(data,keep=True)

        # means for each g value
        means = zeros(n)
        for i,d in enumerate(data):
            means[i] = d['g'][:,0].sum()/d['gsens'][:,0].sum()

        print 'writing:',outfile
        with fitsio.FITS(outfile,mode='rw',clobber=True) as fits:
            fits.write(means)
            fits.write(alldata)
    else:
        print 'reading:',outfile
        with fitsio.FITS(outfile,mode='rw') as fits:
            means = fits[0].read()
            alldata=None
            if len(fits) > 1:
                alldata=fits[1].read()

    if show:
        std=means.std()
        binsize=std/4.
        eu.plotting.bhist(means, binsize=binsize, 
                          title=os.path.basename(outfile))
        
    mean = alldata['g'][:,0].sum()/alldata['gsens'][:,0].sum()
    print 'mean:',mean
    print 'mean from means:',means.mean(),"+/-",means.std()/sqrt(means.size)
    print 'error per ring:',means.std()

    return means, alldata


def randomize_e1e2(e1start,e2start, width=0.1):
    if e1start == 0 or e1start is None or e2start==0 or e2start is None:
        e1rand = 0.05*(randu()-0.5)
        e2rand = 0.05*(randu()-0.5)
    else:
        e1rand = e1start*(1 + 2*width*(randu()-0.5))
        e2rand = e2start*(1 + 2*width*(randu()-0.5))
        etot = sqrt(e1rand**2 + e2rand**2)
        if etot > 0.95:
            e1rand,e2rand=randomize_e1e2(None,None)

        """
        nmax=100
        ii=0
        while True:
            e1rand = e1start*(1 + 2*width*(randu()-0.5))
            e2rand = e2start*(1 + 2*width*(randu()-0.5))
            if etot < 0.95:
                break
            ii += 1
            if ii==nmax:
                wlog("---- hit max try on randomize e1e2, setting zero and restart")
                return randomize_e1e2(None,None)
        """
    return e1rand, e2rand


class MillerRdist:
    """
    Miller et al. 2012

    r*exp(-(r/a)^alpha)
    """
    def __init__(self, rd):
        """
        note for exponentials, Miller et al. 2012 found

            ln(rd/arcsec) ~ -1.145 - 0.269*(i_{814} - 23 )

        where rd is the median of the major axis scale length
        """
        from math import sqrt,log
        self.rd=rd
        self.alpha=4./3.
        self.a=rd/.833

        # these are only for alpha=4/3
        # note .66467=(3/4)gammainc(3/2,x**(4./3.))*gamma(3/2)
        # for x->infinity. Reaches .997 of max at x=4.3, so
        # we should generate randoms for r=a*x > 4.3*a
        self.mode=self.a*( 3.**(3./4.)/2./sqrt(2.) )
        self.norm=1./(self.a*.66467)
        self.logofnorm=log(self.norm)
        self.maxval = self.prob(self.mode)

    def lnprob(self, r):
        return self.logofnorm + log(r) - (r/self.a)**self.alpha

    def prob(self, r):
        return exp(self.lnprob(r))

    def sample(self, nrand):
        """
        generating randoms to > a*4.3 gives "3 sigma" .997
        a*5 gives .9993
        """
        minr=0
        maxr=self.a*5.
        generator=eu.random.CutGenerator(self.prob,[minr,maxr],self.maxval)
        return generator.genrand(nrand)

def test_miller_rdist(rd):
    """
    note ln(rd/arcsec) ~ -1.145 - 0.269*(i_{814} - 23 )
    """

    mrd=MillerRdist(rd)
    minr=0.01
    maxr=5.0*mrd.a
    cr=eu.random.CutGenerator(mrd.prob, [minr,maxr], mrd.maxval)
    cr.test()

    x=numpy.linspace(0.01,25.0,100000)
    p=mrd.prob(x)
    print 'predicted maxval: %.16g found on grid: %.16g' % (mrd.maxval,p.max())

