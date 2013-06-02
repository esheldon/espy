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
from gmix_image.gmix_mcmc import MixMC, MixMCStandAlone, MixMCCoellip
from gmix_image.priors import GPriorBA, CenPrior

import images
import esutil as eu
from esutil.random import srandu
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
    def __init__(self, run, extra=None):
        """
        use config files 

        can over-ride any values with extra= for quick
        tests. 

        """
        super(BAFitSim,self).__init__(run)
        if 'verbose' not in self:
            self['verbose'] = False

        if extra:
            for k,v in extra.iteritems():
                print k,v
                self[k] = v

        self.gprior = GPriorBA(self.simc['gsigma'])


    def process_trial_by_s2n(self, is2, is2n, isplit,
                             dowrite=False, 
                             dolog=False):
        """
        its not really isplit any more, its just a split count
        """

        nellip=self.get_nellip(is2n)
        s2n_psf = self['s2n_psf']
        s2n = shapesim.get_s2n(self, is2n)
        s2n_method = self['s2n_method']
        s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

        #gvals = self.get_gvals(is2, is2n, nellip)
        gvals = self.get_gvals(nellip)

        fitmodels=self.get_fitmodels()
        if 'coellip' in fitmodels[0]:
            ngauss=self.get_coellip_ngauss(fitmodels[0])
            npars=2*ngauss+4
        else:
            npars=6

        out = zeros(nellip*2, dtype=self.out_dtype(npars))

        s2 = linspace(self.simc['mins2'],
                      self.simc['maxs2'], 
                      self.simc['nums2'])[is2]
        out['s2'] = s2
        self['s2']=s2

        i=0
        for ipair,g in enumerate(gvals):

            ellip=lensing.util.g2e(g)

            if ( (ipair+1) % 10) == 0 or ipair== 0:
                stderr.write("  %s/%s pairs done\n" % ((ipair+1),nellip))

            while True:
                theta1 = random.random()*360.0
                theta2 = theta1 + 90.0
                ci_nonoise1 = self.shapesim.get_trial(s2, ellip, theta1)
                ci_nonoise2 = self.shapesim.get_trial(s2, ellip, theta2)
                ci1 = NoisyConvolvedImage(ci_nonoise1, s2n, s2n_psf,
                                          s2n_method=s2n_method,
                                          fluxfrac=s2ncalc_fluxfrac)
                ci2 = NoisyConvolvedImage(ci_nonoise2, s2n, s2n_psf,
                                          s2n_method=s2n_method,
                                          fluxfrac=s2ncalc_fluxfrac)

                try:
                    res1,res2=self._process_pair(ci1,ci2)
                    break
                except TryAgainError:
                    pass

            self._copy_to_output(out, i, ci1, res1)
            i += 1
            self._copy_to_output(out, i, ci2, res2)
            i += 1



        if dowrite:
            shapesim.write_output(self['run'], is2, is2n, out, itrial=isplit,
                         fs=self.fs)
        return out

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

        cenprior=CenPrior(ci['cen'], [0.1]*2)

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

        else:
            self.fitter=MixMCStandAlone(ci.image, ivar, 
                                        psf_gmix, self.gprior, fitmodel,
                                        cen=ci['cen'],
                                        do_pqr=True,
                                        when_prior=self['when_prior'],
                                        nwalkers=self['nwalkers'],
                                        nstep=self['nstep'], 
                                        burnin=self['burnin'],
                                        mca_a=self['mca_a'],
                                        iter=self.get('iter',False),
                                        draw_gprior=self['draw_gprior'])

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

        e1true=ci['e1true']
        e2true=ci['e2true']
        g1true,g2true=lensing.util.e1e2_to_g1g2(e1true,e2true)

        out['gtrue'][i,0] = g1true
        out['gtrue'][i,1] = g2true
        out['shear_true'][i,0] = ci['shear1']
        out['shear_true'][i,1] = ci['shear2']

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
        if 'Fs2n' in res:
            for tn in ['Flux','Ferr','Fs2n']:
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



    def get_gvals_old(self, is2, is2n, nellip):
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
            ('Tmean','f8'),
            ('Terr','f8'),
            ('Ts2n','f8'),
            ('Flux','f8'),
            ('Ferr','f8'),
            ('Fs2n','f8'),
            ('s2n_meas_w','f8'),  # weighted s/n based on most likely point
            ('loglike','f8'),     # loglike of fit
            ('chi2per','f8'),     # chi^2/degree of freedom
            ('dof','i4'),         # degrees of freedom
            ('fit_prob','f8'),    # probability of the fit happening randomly
            ('arate','f8'),
            ('P','f8'),           # parameters from BA13
            ('Q','f8',2),
            ('R','f8',(2,2))
           ]

        return dt




