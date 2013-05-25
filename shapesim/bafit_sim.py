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
from gmix_image.gmix_mcmc import MixMC, MixMCStandAlone
from gmix_image.priors import GPriorBA, CenPrior

import images
import esutil as eu
from esutil.random import srandu
from esutil.misc import wlog

import math
import time

from sys import stderr

LOWVAL=-9999.9e9

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

            out['gtrue'][i,0] = g1true
            out['gtrue'][i,1] = g2true
            out['shear_true'][i,0] = ci_nonoise['shear1']
            out['shear_true'][i,1] = ci_nonoise['shear2']
            
            if self['verbose']:
                stderr.write('-'*70 + '\n')

            # we always write this, although slower when not verbose
            if (nellip > 1) and (( (i+1) % 10) == 0 or i== 0):
                stderr.write("  %s/%s ellip done\n" % ((i+1),nellip))

            ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                     s2n_method=s2n_method,
                                     fluxfrac=s2ncalc_fluxfrac)
            if self['verbose']:
                wlog("s2n_admom:",ci['s2n_admom'],"s2n_uw:",ci['s2n_uw'],
                     "s2n_matched:",ci['s2n_matched'])

            res=self._run_models(ci, fitmodels)

            self._copy_to_output(out, i, ci, res)

        if dowrite:
            shapesim.write_output(self['run'], is2, is2n, out, itrial=itheta,
                         fs=self.fs)
        return out

    def _copy_to_output(self, out, i, ci, res):

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

        self.fitter=MixMCStandAlone(ci.image, ivar, 
                                    psf_gmix, self.gprior, fitmodel,
                                    cen=ci['cen'],
                                    do_pqr=True,
                                    nwalkers=self['nwalkers'],
                                    nstep=self['nstep'], 
                                    burnin=self['burnin'],
                                    mca_a=self['mca_a'],
                                    iter=self.get('iter',False),
                                    draw_gprior=self['draw_gprior'])



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
            ('fit_prob','f8'),    # probability of the fit happening randomly
            ('arate','f8'),
            ('P','f8'),           # parameters from BA13
            ('Q','f8',2),
            ('R','f8',(2,2))
           ]

        return dt




