"""
Bernstein & Armstrong using the ngmix code
"""

import os
from sys import stderr
import ngmix
import numpy
import time

from .shapesim import read_config

from esutil.random import srandu

NSIGMA_RENDER=5.0

class TryAgainError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

class NGMixSim(dict):
    def __init__(self, run, s2n, npairs, **keys):
        """
        Simulate and fit the requested number of pairs at
        the specified s/n
        """

        self.run=run
        self.s2n=s2n
        self.npairs=npairs

        self.update(keys)
        self.conf = read_config(run)
        self.update(self.conf)
        self.simc = read_config(self.conf['sim'])
        self.shear=self.simc['shear']
        self.nsub=self.simc['nsub']

        self.obj_model=self.simc['obj_model']

        self.make_struct()
        self.set_priors()
        self.make_psf()
        self.set_noise()

    def run_sim(self):
        """
        Run the simulation, fitting psf and all pairs
        """
        self.fit_psf()

        i=0
        for ipair in xrange(self.npairs):
            print >>stderr,'%s/%s' % (ipair+1,self.npairs)
            while True:
                try:
                    reslist=self.process_pair()
                    break
                except TryAgainError as err:
                    print >>stderr,str(err)

            self.copy_to_output(reslist[0], i)
            i += 1
            self.copy_to_output(reslist[1], i)
            i += 1

    def process_pair(self):
        """
        Create a simulated image pair and perform the fit
        """

        imdicts = self.get_noisy_image_pair()
        reslist=[]
        for key in imdicts:
            res=self.fit_galaxy(imdicts[key])
            if res['flags'] != 0:
                raise TryAgainError("failed at %s" % key)

            reslist.append(res)
            self.print_res(res)

        return reslist

    def fit_galaxy(self, imdict):
        """
        Fit the model to the galaxy
        """

        gm=imdict['gm_pre']
        T_guess      = gm.get_T()
        counts_guess = gm.get_psum()

        fitter=ngmix.fitting.MCMCSimple(imdict['image'],
                                        imdict['wt'],
                                        imdict['jacobian'],
                                        self.obj_model,

                                        cen_prior=self.cen_prior,
                                        g_prior=self.g_prior,
                                        T_prior=self.T_prior,
                                        counts_prior=self.counts_prior,

                                        T_guess=T_guess,
                                        counts_guess=counts_guess,

                                        psf=self.psf_gmix_fit,
                                        nwalkers=self['nwalkers'],
                                        nstep=self['nstep'],
                                        burnin=self['burnin'],
                                        mca_a=self['mca_a'],
                                        do_pqr=True,
                                        do_lensfit=True)
        fitter.go()
        #fitter.make_plots(show=True)
        return fitter.get_result()

    def print_res(self,res):
        """
        print some stats
        """
        print >>stderr,'    arate:',res['arate']
        ngmix.fitting.print_pars(res['pars'],front='    pars: ',stream=stderr)
        ngmix.fitting.print_pars(res['perr'],front='    perr: ',stream=stderr)

    def fit_psf(self):
        """
        Fit the pixelized psf to a model
        """
        from ngmix.gexceptions import GMixRangeError, GMixMaxIterEM

        print >>stderr,'fitting psf'
        imsky,sky=ngmix.em.prep_image(self.psf_image)

        em=ngmix.em.GMixEM(imsky)
        guess=self.psf_gmix_true.copy()
        print 'psf guess:'
        print guess
        em.go(guess, sky, tol=1.e-5)

        self.psf_gmix_fit=em.get_gmix()
        print 'psf fit:'
        print self.psf_gmix_fit

    def set_priors(self):
        """
        Set all the priors
        """

        print >>stderr,"setting priors"
        T=self.simc['obj_T_mean']
        T_sigma = self.simc['obj_T_sigma_frac']*T
        counts=self.simc['obj_counts_mean']
        counts_sigma = self.simc['obj_counts_sigma_frac']*counts

        self.g_prior=ngmix.priors.GPriorBA(0.3)
        self.cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 0.1, 0.1)
        self.T_prior=ngmix.priors.LogNormal(T, T_sigma)
        self.counts_prior=ngmix.priors.LogNormal(counts, counts_sigma)

    def make_psf(self):
        """
        make the psf gaussian mixture model
        """

        print >>stderr,"making psf"

        self.psf_dims, self.psf_cen=self.get_dims_cen(self.simc['psf_T'])

        pars=[self.psf_cen[0],
              self.psf_cen[1],
              self.simc['psf_shape'][0],
              self.simc['psf_shape'][1],
              self.simc['psf_T'],
              1.0]
        self.psf_gmix_true=ngmix.gmix.GMixModel(pars, self.simc['psf_model'])
        
        self.psf_image=self.psf_gmix_true.make_image(self.psf_dims,
                                                     nsub=self.nsub)
    
    def set_noise(self):
        """
        Find gaussian noise that when added to the image 
        produces the requested s/n.  Use a matched filter.

         sum(pix^2)
        ------------ = S/N^2
          skysig^2

        thus
            
        sum(pix^2)
        ---------- = skysig^2
          (S/N)^2
        """
        
        from numpy.random import randn

        print >>stderr,"setting noise"

        imdict=self.get_image_pair(random=False)
        im=imdict['im1']['image']
        skysig2 = (im**2).sum()/self.s2n**2
        skysig = numpy.sqrt(skysig2)

        noise_image = skysig*randn(im.size).reshape(im.shape)
        new_im = im + noise_image

        s2n_check = numpy.sqrt( (im**2).sum()/skysig**2 )
        print >>stderr,"S/N goal:",self.s2n,"found:",s2n_check

        self.skysig=skysig
        self.ivar=1.0/skysig**2


    def get_noisy_image_pair(self, random=True):
        """
        Get an image pair, with noise added
        """
        imdict=self.get_image_pair(random=random)
        self.add_noise(imdict['im1']['image'])
        self.add_noise(imdict['im2']['image'])

        wt=numpy.zeros(imdict['im1']['image'].shape) + self.ivar
        imdict['im1']['wt']=wt
        imdict['im2']['wt']=wt
        return imdict

    def add_noise(self, im):
        """
        Add gaussian random noise
        """

        from numpy.random import randn
        im[:,:] += self.skysig*randn(im.size).reshape(im.shape)

    def get_image_pair(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors
        """

        cen_offset, shape1, shape2, T, counts=self.get_pair_pars(random=random)

        # center is just placeholder for now
        pars1=[0.0, 0.0, shape1.g1, shape1.g2, T, counts]
        pars2=[0.0, 0.0, shape2.g1, shape2.g2, T, counts]

        gm1_pre=ngmix.gmix.GMixModel(pars1, self.obj_model)
        gm2_pre=ngmix.gmix.GMixModel(pars2, self.obj_model)

        gm1  = gm1_pre.convolve(self.psf_gmix_true)
        gm2  = gm2_pre.convolve(self.psf_gmix_true)

        T = gm1.get_T()
        dims, cen = self.get_dims_cen(T)

        # jacobian is at center before offset
        j=ngmix.jacobian.UnitJacobian(cen[0], cen[1])

        cen[0] += cen_offset[0]
        cen[1] += cen_offset[1]

        gm1.set_cen(cen[0], cen[1])
        gm2.set_cen(cen[0], cen[1])

        nsub = self.nsub
        im1=gm1.make_image(dims, nsub=nsub)
        im2=gm2.make_image(dims, nsub=nsub)

        out={'im1':{'gm_pre':gm1_pre,'gm':gm1,'image':im1,'jacobian':j},
             'im2':{'gm_pre':gm2_pre,'gm':gm2,'image':im2,'jacobian':j}}
        return out

    def get_pair_pars(self, random=False):
        """
        Get pair parameters
        """
        from numpy.random import random as randu

        if random:
            cen_offset=self.cen_prior.sample()
            g = self.g_prior.sample1d(1)
            g=g[0]
            rangle1 = randu()*2*numpy.pi
            rangle2 = rangle1 + numpy.pi/2.0
            g1_1 = g*numpy.cos(rangle1)
            g2_1 = g*numpy.sin(rangle1)
            g1_2 = g*numpy.cos(rangle2)
            g2_2 = g*numpy.sin(rangle2)

            T=self.T_prior.sample()
            counts=self.counts_prior.sample()
        else:
            cen_offset=[0.0, 0.0]
            g1_1=0.0
            g2_1=0.0
            g1_2=0.0
            g2_2=0.0
            T=self.T_prior.mean
            counts=self.counts_prior.mean

        shape1=ngmix.shape.Shape(g1_1, g2_1)
        shape2=ngmix.shape.Shape(g1_2, g2_2)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        return cen_offset, shape1, shape2, T, counts

    def get_dims_cen(self, T):
        """
        Based on T, get the required dimensions and a center
        """
        sigma=numpy.sqrt(T/2.)
        dims = [2.*sigma*NSIGMA_RENDER]*2
        cen = [(dims[0]-1.)/2.]*2

        return dims, cen

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data

    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        d=self.data
        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']
        d['P'][i] = res['P']
        d['Q'][i,:] = res['Q']
        d['R'][i,:,:] = res['R']
        d['g'][i,:] = res['g']
        d['gsens'][i,:] = res['g_sens']

    def make_struct(self):
        """
        Make the output array
        """
        dt=[('pars','f8',6),
            ('pcov','f8',(6,6)),
            ('P','f8'),
            ('Q','f8',2),
            ('R','f8',(2,2)),
            ('g','f8',2),
            ('gsens','f8',2)]
        self.data=numpy.zeros(self.npairs*2, dtype=dt)

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

            center_offset=None
            while True:
                theta1 = random.random()*360.0
                theta2 = theta1 + 90.0
                if self.center_dist is not None:
                    center_offset=self.center_dist.sample()

                ci1=self._get_one_trial(Tobj, counts, ellip, theta1, center_offset=center_offset)
                ci2=self._get_one_trial(Tobj, counts, ellip, theta2, center_offset=center_offset)

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
            shapesim.write_output(self.run, iT, is2n, out, itrial=isplit,
                         fs=self.fs)

        return out

    def _get_one_trial(self, Tobj, counts, ellip, theta, center_offset=None):

        ci_nonoise = self.shapesim.get_trial(Tobj, ellip, theta, counts=counts,
                                             center_offset=center_offset)

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




