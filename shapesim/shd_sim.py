import numpy
import gmix_image
import lensing
from lensing.shear import Shear
from lensing.util import ShapeRangeError

LOWVAL=-9999.9e9

class SHDirectBase(dict):
    def __init__(self, 
                 image_list, 
                 ivar_list, 
                 psf_list, 
                 model_list, 
                 ares_list, 
                 nwalkers, 
                 burnin, 
                 nstep, 
                 **keys):
        """

        need to determine the best model for each with a pre-run of some other
        code!

        """

        self._image_list=image_list
        self._ivar_list=ivar_list
        self._psf_list=psf_list
        self._model_list=model_list

        self._ares_list=ares_list

        self._nwalkers=nwalkers
        self._burnin=burnin
        self._nstep=nstep
        self._mca_a = 3

        self._nimages=len(self._image_list)

        self.update(keys)

    def get_result(self):
        return self._result
    def get_trials(self):
        return self._trials
    def get_lnprobs(self):
        return self._lnprobs

    def _go(self):
        self.sampler=self._do_trials()

        self.trials  = self.sampler.flatchain

        lnprobs = self.sampler.lnprobability.reshape(self._nwalkers*self._nstep)
        self._lnprobs = lnprobs - lnprobs.max()

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if self.make_plots:
            self._doplots()

    def _do_trials(self):

        guess=self._get_guess()

        sampler = self._get_sampler()

        pos, prob, state = sampler.run_mcmc(guess, self._burnin)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self._nstep)

        return sampler


    def _calc_lnprob(self, pars):
        lnprob=0.0

        try:
            shear=Shear(g1=pars[0], g2=pars[1])
        except ShapeRangeError:
            return LOWVAL

        for i in xrange(self._nimages):
            image=self._image_list[i]
            ivar=self._ivar_list[i]
            psf=self._psf_list[i]
            model=self._model_list[i]
            ipars=2 + 6*i
            pars1=pars[ipars:ipars+6]

            lnprob1 = self._calc_lnprob1(image, ivar, psf, model, pars1, shear)

            lnprob += lnprob1

        return lnprob

    def _calc_lnprob1(self, image, ivar, psf, model, pars, shear):

        try:
            shape0=lensing.Shear(g1=pars[2], g2=pars[3])
            shape=shape0 + shear
        except ShapeRangeError:
            return LOWVAL

        # note these are in e space
        pars_sheared=pars.copy()
        pars_sheared[2] = shape.e1
        pars_sheared[3] = shape.e2

        gmix0 = gmix_image.GMix(pars_sheared, model=model)
        gmix = gmix.convolve(psf)

        logprob = self._get_loglike_c1(image, ivar, gmix)
        return logprob

    def _get_loglike_c(self, image, ivar, gmix):
        from gmix_image import render
        loglike,s2n,flags = render._render.loglike(image, gmix, ivar)

        if flags != 0:
            return LOWVAL
        return loglike



    def _calc_result(self):
        import mcmc

        shear=zeros(2)
        shear_cov=zeros((2,2))

        # prior is already in the distribution of
        # points.  This is simpler for most things but
        # for sensitivity we need a factor of (1/P)dP/de

        pars,pcov = mcmc.extract_stats(self.trials)

        shear[:] = pars[0:0+2]
        shear_cov[:,:] = pcov[0:0+2, 0:0+2]
 
        arates = self.sampler.acceptance_fraction
        arate = arates.mean()

        self._result={'shear':shear,
                      'shear_cov':shear_cov,
                      'pars':pars,
                      'pcov':pcov,
                      'arate':arate}


    def _get_sampler(self):
        import emcee
        sampler = emcee.EnsembleSampler(self._nwalkers, 
                                        self._npars, 
                                        self._calc_lnprob,
                                        a=self._mca_a)
        return sampler


    def _get_guess(self):
        
        nimage=len(self._image_list)
        nwalkers=self._nwalkers

        npars = 2 + 6*nimage
        self._npars=npars
        
        guess=numpy.zeros( (nwalkers,npars) )

        for i in xrange(nimage):
            ii = 2 + 6*i

            ai=self._ares_list[i]
            if ai['whyflag'] != 0:
                raise ValueError("admom with flags != 0")

            aT = ai['Irr'] + ai['icc']
            counts=self._image_list[i].sum()

            guess[:,ii+0] = ai['wrow']*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+1] = ai['wcol']*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+4] = aT*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+5] = counts*(1. + 0.1*srandu(nwalkers))

            shape=lensing.shear.Shear(e1=ai['e1'], e2=ai['e2'])
            for iw in xrange(nwalkers):
                shadd=lensing.shear.Shear(g1=0.05*srandu(), g2=0.05*srandu())

                shw = shape + shadd

                guess[iw,ii+2] = shw.g1
                guess[iw,ii+3] = shw.g2
        
        # this is the shear
        guess[:,0] = 0.05*srandu(nwalkers)
        guess[:,1] = 0.05*srandu(nwalkers)

        self._guess=guess
        return guess


 class SHDirect(SHDirectBase):
    def __init__(self, image_list, ivar_list, psf_list, model_list, ares_list, nwalkers, burnin, nstep, **keys):
        super(SHDirect,self).__init__(image_list, ivar_list, psf_list, model_list, ares_list, nwalkers, burnin, nstep, **keys)
        
        self._go()



