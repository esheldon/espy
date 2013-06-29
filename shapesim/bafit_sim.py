"""
MCMC with bayesian approach, using the Bernstein & Armstrong approach
"""

import os
import pprint
import numpy
from numpy import sqrt, cos, sin, exp, log, log10, pi, zeros, ones, empty, \
        random, where, array, linspace, diag, median
from numpy.random import randn
from numpy.random import random as randu
from . import shapesim
import lensing
from lensing.shear import Shear
import fimage
from fimage.convolved import NoisyConvolvedImage

import gmix_image
from gmix_image import print_pars, GMix, gmix2pars
from gmix_image.gmix_mcmc import MixMCSimple, MixMCCoellip, MixMCBD
from gmix_image.priors import GPriorBA, CenPrior

import images
import esutil as eu
from esutil.random import srandu, LogNormal, Normal
from esutil.misc import wlog

import math
import time

from sys import stderr

LOWVAL=-9999.9e9

class TryAgainError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

class BAFitSim(shapesim.BaseSim):
    def __init__(self, run, **keys):
        """
        use config files 

        can over-ride any values with extra= for quick
        tests. 

        """

        super(BAFitSim,self).__init__(run)
        self.update(keys)

        cw=self.simc.get('cen_width',None)
        if cw is not None:
            cen_dist=self.simc['cen_dist']
            if cen_dist.lower()=="normal":
                center_dist=gmix_image.priors.CenPrior([0.0]*2,[cw]*2)
            else:
                raise ValueError("implement non-normal cen dist")
            self.shapesim['center_dist']=center_dist

        if 'verbose' not in self:
            self['verbose'] = False

        self.gprior = GPriorBA(self.simc['g_width'])


    def _set_s2n_info(self, is2n):
        self._s2n_psf = self['s2n_psf']
        self._s2n = shapesim.get_s2n(self, is2n)
        self._s2n_method = self['s2n_method']
        self._s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

    def _get_npars(self):
        fitmodels=self.get_fitmodels()
        if 'coellip' in fitmodels[0]:
            ngauss=self.get_coellip_ngauss(fitmodels[0])
            npars=2*ngauss+4
        elif 'bd' in fitmodels[0]:
            npars=8
        else:
            npars=6
        return npars

   

    def process_trial_by_s2n(self, iT, is2n, isplit,
                             dowrite=False, 
                             dolog=False):

        t0=time.time()

        nellip=self.get_nellip(is2n)
        self._set_s2n_info(is2n)

        gvals = self.get_gvals(nellip)
        npars=self._get_npars()

        out = zeros(nellip*2, dtype=self.out_dtype(npars))


        i=0
        for ipair,g in enumerate(gvals):

            Tobj=self.shapesim._get_Tobj(iT)
            counts=self.shapesim._get_counts()

            ellip=lensing.util.g2e(g)

            if self['verbose'] or ( ( (ipair+1) % 10) == 0 or ipair== 0):
                stderr.write("  %s/%s pairs done g: %g\n" % ((ipair+1),nellip,g))

            while True:
                theta1 = random.random()*360.0
                theta2 = theta1 + 90.0
                ci1=self._get_one_trial(Tobj, counts, ellip, theta1)
                ci2=self._get_one_trial(Tobj, counts, ellip, theta2)

                try:
                    res1,res2=self._process_pair(ci1,ci2)
                    break
                except TryAgainError:
                    pass

            self._copy_to_output(out, i, ci1, res1)
            i += 1
            self._copy_to_output(out, i, ci2, res2)
            i += 1


        
        tm=time.time()-t0
        print 'total time:',tm
        print 'time per:',0.5*tm/nellip

        if dowrite:
            shapesim.write_output(self['run'], iT, is2n, out, itrial=isplit,
                         fs=self.fs)

        return out

    def _get_one_trial(self, Tobj, counts, ellip, theta):

        ci_nonoise = self.shapesim.get_trial(Tobj, ellip, theta, counts=counts)

        if self['retrim']:
            if 'retrim_fluxfrac' not in self:
                raise ValueError("you must set fluxfrac for a retrim")
            retrim_fluxfrac = self['retrim_fluxfrac']
            if self['verbose']:
                print >>stderr,"trimming:",retrim_fluxfrac,"before:",ci_nonoise.image.shape,
            ci_nonoise= fimage.convolved.TrimmedConvolvedImage(ci_nonoise,
                                                               fluxfrac=retrim_fluxfrac)
            if self['verbose']:
                print >>stderr,"after:",ci_nonoise.image.shape

        ci = NoisyConvolvedImage(ci_nonoise, self._s2n, self._s2n_psf,
                                 s2n_method=self._s2n_method,
                                 fluxfrac=self._s2ncalc_fluxfrac)
        return ci

    def _process_pair(self, ci1, ci2):
        reslist=[]
        for ci in [ci1,ci2]:
            res=self._run_models(ci)
            if res['flags'] != 0:
                raise TryAgainError("failed")
            reslist.append(res)
        return reslist

    def _run_models(self, ci):
        # fit models, keep the one that most looks like random error
        probrand=-9999e9
        fitmodels=self.get_fitmodels()
        for fitmodel in fitmodels:
            self._run_fitter(ci, fitmodel)
            res0 = self.fitter.get_result()

            fit_prob=res0.get('fit_prob',-9999)
            if len(fitmodels) > 1:
                print '  model:',fitmodel,'probrand:',fit_prob

            if fit_prob > probrand:
                res=res0
                probrand=fit_prob

        if len(fitmodels) > 1:
            print '    best model:',res['model']
        return res

    def _get_cen_prior(self, ci):
        cen_width = self.simc.get('cen_width',0.1)
        cen_dist = self.simc['cen_dist']
        if cen_dist.lower()=="normal":
            cen_prior=CenPrior(ci['cen'], [cen_width]*2)
        else:
            raise ValueError("implement non-normal cen dist")
        return cen_prior

    def _get_counts_prior(self, ci):
        counts_prior=None
        counts_dist = self.simc.get('counts_dist',None)
        if counts_dist is not None:
            counts_width_frac=self.simc['counts_width_frac']

            counts_mean=self.shapesim._counts_mean
            counts_width = counts_mean*counts_width_frac

            counts_prior = eu.random.get_dist(counts_dist,
                                              [counts_mean,
                                              counts_width])
        
        return counts_prior

    def _get_T_prior(self, ci):
        T_prior=None

        T_dist = self.simc.get('T_dist',None)
        if T_dist is not None:
            T_width_frac=self.simc['T_width_frac']
            T_mean = ci['Ttrue']
            T_width = T_mean*T_width_frac

            T_prior = eu.random.get_dist(T_dist, [T_mean, T_width])
 
        return T_prior


    def _run_fitter(self, ci, fitmodel):
        from gmix_image.gmix_em import GMixEMBoot

        psf_ivar=1./ci['skysig_psf']**2
        gmpsf=GMixEMBoot(ci.psf, self['ngauss_psf'], ci['cen_psf'],
                         ivar=psf_ivar,
                         maxiter=self['em_maxiter'],
                         tol=self['em_tol'])

        psf_gmix=gmpsf.get_gmix()

        Tguess = ci['Ttrue']*(1. + 0.1*srandu())
        ivar=1./ci['skysig']**2

        cen_prior=self._get_cen_prior(ci)
        counts_prior=self._get_counts_prior(ci)
        T_prior=self._get_T_prior(ci)

        if 'coellip' in fitmodel:
            ngauss=self.get_coellip_ngauss(fitmodel)
            self.fitter=MixMCCoellip(ci.image, ivar, 
                                     psf_gmix, self.gprior, ngauss,
                                     cen=ci['cen'],
                                     do_pqr=True,
                                     nwalkers=self['nwalkers'],
                                     nstep=self['nstep'], 
                                     burnin=self['burnin'],
                                     mca_a=self['mca_a'],
                                     iter=self.get('iter',False),
                                     draw_gprior=self['draw_gprior'])

        elif 'bd' in fitmodel:
            raise ValueError("fix bd")
            self.fitter=MixMCBD(ci.image, ivar, 
                                 psf_gmix, self.gprior, 
                                 cen=ci['cen'],
                                 do_pqr=True,
                                 nwalkers=self['nwalkers'],
                                 nstep=self['nstep'], 
                                 burnin=self['burnin'],
                                 mca_a=self['mca_a'],
                                 iter=self.get('iter',False),
                                 draw_gprior=self['draw_gprior'])

        else:
            sampler_type=self.get('sampler','mcmc')

            T_guess=ci['Ttrue']*(1.+0.3*srandu())
            counts_guess=ci['counts_true']*(1.0 + 0.3*srandu())
            cen_guess=ci['cen']



            if sampler_type=='mcmc':
                keys={}
                keys.update(self)
                keys['do_pqr']=True
                keys['make_plots']=self.get('make_plots',False)

                keys['cen_prior']=cen_prior
                keys['T_prior']=T_prior
                keys['counts_prior']=counts_prior
                keys['when_prior']=self['when_prior']
                keys['nwalkers']=self['nwalkers']
                keys['nstep']=self['nstep']
                keys['burnin']=self['burnin']
                keys['mca_a']=self['mca_a']
                keys['iter']=self.get('iter',False),
                keys['draw_gprior']=self['draw_gprior']

                self.fitter=MixMCSimple(ci.image,
                                        ivar, 
                                        psf_gmix,
                                        self.gprior,
                                        T_guess,
                                        counts_guess,
                                        cen_guess,
                                        fitmodel,
                                        **keys)

            elif sampler_type=='cmcmc':
                config={}
                config.update(self)

                config['model'] = fitmodel

                config['cen1_mean']=ci['cen'][0]
                config['cen1_width']=self.simc['cen_width']

                config['cen2_mean']=ci['cen'][1]
                config['cen2_width']=self.simc['cen_width']

                config['g_width']=self.simc['g_width']

                config['T_mean']=T_prior.get_mean()
                config['T_width']=T_prior.get_sigma()

                config['counts_mean']=counts_prior.get_mean()
                config['counts_width']=counts_prior.get_sigma()

                guess=self._get_cmcmc_simple_guess(ci,config)

                self.fitter=gmix_image.gmix_mcmc.MixMCC(ci.image,
                                                        ivar,
                                                        psf_gmix,
                                                        guess,
                                                        config,
                                                        gprior=self.gprior)

            elif sampler=='isample':
                keys['nsample']=self['nsample']
                g1_guess,g2_guess=self.gprior.sample2d(1)
                guess=numpy.zeros(6)
                guess[0:2] = cen_guess
                guess[2] = g1_guess[0]
                guess[3] = g2_guess[0]
                guess[4] = T_guess
                guess[5] = counts_guess

                #prior_samples=self._presample_prior_simple(cen_prior,
                #                                           self.gprior,
                #                                           T_prior,
                #                                           counts_prior)
                prior_samples=self._presample_gprior()
                self.fitter=gmix_image.gmix_isamp.GMixIsampSimple(ci.image,
                                                                  ivar,
                                                                  psf_gmix,

                                                                  cen_prior,
                                                                  self.gprior,
                                                                  T_prior,
                                                                  counts_prior,

                                                                  prior_samples,

                                                                  guess,
                                                                  fitmodel,
                                                                  **keys)

    def _get_cmcmc_simple_guess(self, ci, config):
        nwalkers=config['nwalkers']
        guess=numpy.zeros( (nwalkers, 6) )

        # cen uniform within 0.1 pixels of truth
        guess[:,0] = ci['cen'][0] + 0.1*srandu(nwalkers)
        guess[:,1] = ci['cen'][0] + 0.1*srandu(nwalkers)
        
        g_draw=config['g_draw']
        if g_draw=="prior":
            g1rand,g2rand=self.gprior.sample2d(nwalkers)
            guess[:,2]=g1rand
            guess[:,3]=g2rand
        elif g_draw=="truth":
            sh=lensing.Shear(e1=ci['e1true'],e2=ci['e2true'])
            g1=sh.g1
            g2=sh.g2

            g1rand=numpy.zeros( nwalkers )
            g2rand=numpy.zeros( nwalkers )

            nleft=nwalkers
            ngood=0
            while nleft > 0:
                g1rand_t = g1 + 0.01*srandu(nleft)
                g2rand_t = g2 + 0.01*srandu(nleft)
                g2tot=g1rand_t**2 + g2rand_t**2

                w,=numpy.where(g2tot < 0.999)
                if w.size > 0:
                    g1rand[ngood:ngood+w.size] = g1rand_t[w]
                    g2rand[ngood:ngood+w.size] = g2rand_t[w]
                    ngood += w.size
                    nleft -= w.size
            
            guess[:,2]=g1rand
            guess[:,3]=g2rand
        elif g_draw=="maxlike":
            raise ValueError("implement getting maxlike as guess")

        guess[:,4] = ci['Ttrue']*(1.+0.1*srandu(nwalkers))
        guess[:,5] = ci['counts_true']*(1.0 + 0.1*srandu(nwalkers))
        
        return guess


    def _presample_gprior(self):
        if not hasattr(self,'_gprior_samples'):
            npre=self.get('n_pre_sample',10000)
            print >>stderr,'    pre-sampling gprior:',npre
            g1,g2= self.gprior.sample2d(npre)
            self._gprior_samples= {'g1':g1, 'g2':g2}

        return self._gprior_samples

    def _presample_prior_simple(self,cen_prior,g_prior,T_prior,counts_prior):
        from gmix_image import priors
        if not hasattr(self,'_prior_samples'):

            npre=self.get('n_pre_sample',100000)
            print >>stderr,'pre-sampling the prior',npre
            npars=6
            nwalkers=20
            # per walker, 2000
            burnin=100

            start=numpy.zeros( (nwalkers,npars) )

            start[:,0:2] = cen_prior.sample(nwalkers)

            start[:,2] = 0.1*srandu(nwalkers)
            start[:,3] = 0.1*srandu(nwalkers)

            start[:,4] = T_prior.sample(nwalkers)
            start[:,5] = counts_prior.sample(nwalkers)


            comb=priors.CombinedPriorSimple(cen_prior,
                                            g_prior,
                                            T_prior,
                                            counts_prior)
            sampler = comb.sample(start,
                                  npre,
                                  burnin=burnin,
                                  nwalkers=nwalkers,
                                  get_sampler=True)
            prand = sampler.flatchain
            lnp = sampler.lnprobability
            lnp = lnp.reshape(lnp.shape[0]*lnp.shape[1])
            probs = numpy.exp(lnp)

            self._prior_samples={'samples':prand,
                                 'prob':probs}
        
        return self._prior_samples

    def get_coellip_ngauss(self, model):
        if model=='coellip1':
            return 1
        elif model=='coellip2':
            return 2
        elif model=='coellip3':
            return 3
        else:
            raise ValueError("implement '%s'" % model)
    def get_fitmodels(self):
        fitmodels=self['fitmodel']
        if not isinstance(fitmodels,list):
            fitmodels=[fitmodels]
        return fitmodels


    def _copy_to_output(self, out, i, ci, res):
        sampler_type=self.get('sampler','mcmc')

        out['flags'][i] = res['flags']
        if res['flags'] != 0:
            return

        out['s2'][i] = self.simc['Tpsf']/ci['Ttrue']
        out['sratio'][i] = sqrt(1./out['s2'][i])

        e1true=ci['e1true']
        e2true=ci['e2true']
        g1true,g2true=lensing.util.e1e2_to_g1g2(e1true,e2true)

        out['gtrue'][i,0] = g1true
        out['gtrue'][i,1] = g2true
        out['shear_true'][i,0] = ci['shear1']
        out['shear_true'][i,1] = ci['shear2']

        out['Ttrue'][i] = ci['Ttrue']

        out['s2n_admom'][i] = ci['s2n_admom']
        out['s2n_matched'][i] = ci['s2n_matched']
        out['s2n_uw'][i] = ci['s2n_uw']
        out['s2n_admom_psf'][i] = ci['s2n_admom_psf']
        out['s2n_matched_psf'][i] = ci['s2n_matched_psf']
        out['s2n_uw_psf'][i] = ci['s2n_uw_psf']

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
        for tn in ['flux','flux_err','flux_s2n']:
            if tn in res:
                out[tn][i] = res[tn]

        out['s2n_meas_w'][i] = res['s2n_w']
        out['loglike'][i] = res['loglike']
        out['chi2per'][i] = res['chi2per']
        out['dof'][i] = res['dof']
        out['fit_prob'][i] = res['fit_prob']
        if 'arate' in res:
            out['arate'][i] = res['arate']


        if 'P' in res:
            out['P'][i] = res['P']
            out['Q'][i] = res['Q']
            out['R'][i] = res['R']

        if sampler_type=='isample':
            out['lm_flux'][i] = res['lm_result']['pars'][5]
            out['lm_flux_err'][i] = res['lm_result']['perr'][5]


    def get_nellip(self, is2n):
        s2n = shapesim.get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']
        nellip = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        if nellip < self['min_gcount']:
            nellip=self['min_gcount']
        return nellip

    def get_gvals(self, nellip):
        numpy.random.seed(None)
        gvals = self.gprior.sample1d(nellip)
        return gvals


    def out_dtype(self, npars):

        dt=[('model','S20'),
            ('s2n_admom','f8'),
            ('s2n_matched','f8'),
            ('s2n_uw','f8'),
            ('s2n_admom_psf','f8'),
            ('s2n_matched_psf','f8'),
            ('s2n_uw_psf','f8'),
            ('sratio','f8'),
            ('s2','f8'),
            ('shear_true','f8',2),
            ('gtrue','f8',2),
            ('Ttrue','f8'),
            ('g','f8',2),
            ('gsens','f8',2),
            ('gcov','f8',(2,2)),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('Tmean','f8'),
            ('Terr','f8'),
            ('Ts2n','f8'),
            ('flux','f8'),
            ('flux_err','f8'),
            ('flux_s2n','f8'),
            ('s2n_meas_w','f8'),  # weighted s/n based on most likely point
            ('loglike','f8'),     # loglike of fit
            ('chi2per','f8'),     # chi^2/degree of freedom
            ('dof','i4'),         # degrees of freedom
            ('fit_prob','f8'),    # probability of the fit happening randomly
            ('arate','f8'),
            ('P','f8'),           # parameters from BA13
            ('Q','f8',2),
            ('R','f8',(2,2)),
            ('flags','i4')
           ]

        sampler_type=self.get('sampler','mcmc')
        if sampler_type=='isample':
            dt += [('lm_flux','f8'),
                   ('lm_flux_err','f8')]
        return dt




