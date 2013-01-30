"""
One giant fit over all galaxies with two extra parameters for shear

Doesn't seem to work, completely degenerate or I have a bug
"""
import numpy
import gmix_image
import lensing
from lensing.shear import Shear
from lensing.util import ShapeRangeError
from .shapesim import read_config, ShapeSim
from . import shapesim

from esutil.random import srandu

LOWVAL=-9999.9e9

class SHDirectSim(dict):
    def __init__(self, run):
        conf=read_config(run)
        self.update(conf)

        self.simc = read_config(self['sim'])

        # for now just set this
        self.npair=30 # we do two times this, pairs at 90 degrees

        self._set_gprior()

        simpars=self.get('simpars',{})
        self.shapesim = ShapeSim(self['sim'], **simpars)

    def process_trials(self, is2, is2n, itrial=None):
        if itrial is not None:
            raise ValueError("implement itrial")

        self._make_images(is2, is2n)

        shd = SHDirect(self.image_list,
                       self.ivar_list,
                       self.psf_list,
                       self.model_list,
                       self.ares_list,
                       self['nwalkers'],
                       self['burnin'],
                       self['nstep'],
                       make_plots=self['make_plots'],
                       mca_a=self['mca_a'],
                       gprior=self.gprior)

        trials=shd.get_trials()
        shapesim.write_output(self['run'], is2, is2n, trials)

    def _make_images(self, is2, is2n):
        s2n_psf = self['s2n_psf']
        s2n = shapesim.get_s2n(self, is2n)
        s2 = numpy.linspace(self.simc['mins2'],
                            self.simc['maxs2'], 
                            self.simc['nums2'])[is2]

        self.image_list=[]
        self.ivar_list=[]
        self.psf_list=[]
        self.model_list=[]
        self.ares_list=[]

        gvals=self.gprior.sample1d(self.npair)
        for i,g in enumerate(gvals):
            print i,g
            #g=0.001 + 0.01*numpy.random.random()
            theta = 360.0*numpy.random.random()
            theta2 = theta + 90.0

            ellip=lensing.util.g2e(g)
            ci,ares,psf = self._get_ci_ares_psf(s2, ellip, theta, s2n, s2n_psf)
            ci2,ares2,psf2 = self._get_ci_ares_psf(s2, ellip, theta2, s2n, s2n_psf)

            self.image_list += [ci.image, ci2.image]
            self.ivar_list += [1./ci['skysig']**2, 1./ci2['skysig']**2]
            self.ares_list += [ares,ares2]
            self.psf_list += [psf,psf2]
            self.model_list += ['gauss','gauss']


    def _get_ci_ares_psf(self, s2, ellip, theta, s2n, s2n_psf):
        from fimage.convolved import NoisyConvolvedImage

        ci_nonoise = self.shapesim.get_trial(s2, ellip, theta)

        ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                 s2n_method=self['s2n_method'])
        # build up a fake ares 
        ares={'wrow':ci['cen_admom'][0],
              'wcol':ci['cen_admom'][1],
              'Irr':ci['cov_admom'][0],
              'Irc':ci['cov_admom'][1],
              'Icc':ci['cov_admom'][2],
              'e1':ci['e1_admom'],
              'e2':ci['e2_admom'],
              'whyflag':0}
        #print ares
        # for now only gauss
        psf=gmix_image.GMix([1.0, 
                             ci['cen_psf_admom'][0],
                             ci['cen_psf_admom'][1],
                             ci['cov_psf_admom'][0],
                             ci['cov_psf_admom'][1],
                             ci['cov_psf_admom'][2]])

        return ci, ares, psf

    def _set_gprior(self):
        import cluster_step
        exp_prior_pars=cluster_step.files.read_gprior(type='gexp')

        self.gprior=cluster_step.prior.GPriorExp(exp_prior_pars['pars'][3])

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
        self._mca_a = keys.get('mca_a', 2.0)

        self._gprior=keys.get('gprior',None)

        self._make_plots=keys.get('make_plots',False)

        self._nimages=len(self._image_list)

        self._npars = 2 + 6*self._nimages
        min_walkers=4*self._npars
        if self._nwalkers < min_walkers:
            self._nwalkers=min_walkers

        print 'nwalkers:',self._nwalkers

        self.update(keys)

    def get_result(self):
        return self._result
    def get_trials(self):
        return self._trials
    def get_lnprobs(self):
        return self._lnprobs

    def _go(self):
        self._sampler=self._do_trials()

        self._trials  = self._sampler.flatchain

        lnprobs = self._sampler.lnprobability.reshape(self._nwalkers*self._nstep)
        self._lnprobs = lnprobs - lnprobs.max()

        # get the expectation values, sensitivity and errors
        print 'calculating stats'
        self._calc_result()

        if self._make_plots:
            self._doplots()

    def _do_trials(self):

        guess=self._get_guess()

        sampler = self._get_sampler()

        print 'burning in'
        self._setup_progress(self._burnin*self._nwalkers)
        pos, prob, state = sampler.run_mcmc(guess, self._burnin)
        sampler.reset()

        print 'main run'
        self._setup_progress(self._nstep*self._nwalkers)
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

        self._update_progress()
        return lnprob

    def _calc_lnprob1(self, image, ivar, psf, model, pars, shear):

        try:
            shape0=lensing.Shear(g1=pars[2], g2=pars[3])
        except ShapeRangeError:
            return LOWVAL

        if self._gprior is not None:
            gp = self._get_lngprior(shape0.g1, shape0.g2)
        else:
            gp=0.0

        try:
            shape=shape0 + shear
        except ShapeRangeError:
            return LOWVAL

 

        # note these are in e space
        pars_sheared=pars.copy()
        pars_sheared[2] = shape.e1
        pars_sheared[3] = shape.e2

        gmix0=self._get_gmix0(pars_sheared, model)
        gmix = gmix0.convolve(psf)

        logprob = self._get_loglike_c1(image, ivar, gmix)

        logprob += gp

        return logprob

    def _get_lngprior(self, g1, g2):
        from math import sqrt, log
        g=sqrt(g1**2 + g2**2)
        gp = self._gprior.prior2d_gabs_scalar(g)
        if gp > 0:
            gp = log(gp)
        else:
            gp=LOWVAL
        return gp

    def _get_loglike_c1(self, image, ivar, gmix):
        from gmix_image import render
        loglike,s2n,flags = render._render.loglike(image, gmix, ivar)

        if flags != 0:
            return LOWVAL
        return loglike


    def _get_gmix0(self, epars, model):
        if model == 'gauss':
            type='coellip'
        else:
            raise ValueError("implement type '%s'" % model)
        gmix0 = gmix_image.GMix(epars, type=type)
        return gmix0


    def _calc_result(self):
        import mcmc

        shear=numpy.zeros(2)
        shear_cov=numpy.zeros((2,2))

        pars,pcov = mcmc.extract_stats(self._trials)

        shear[:] = pars[0:0+2]
        shear_cov[:,:] = pcov[0:0+2, 0:0+2]
 
        arates = self._sampler.acceptance_fraction
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
        
        nwalkers=self._nwalkers
        npars=self._npars

        guess=numpy.zeros( (nwalkers,npars) )

        for i in xrange(self._nimages):
            ii = 2 + 6*i

            ai=self._ares_list[i]
            if ai['whyflag'] != 0:
                raise ValueError("admom with flags != 0")

            aT = ai['Irr'] + ai['Icc']
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
        #guess[:,0] = 0.05*srandu(nwalkers)
        #guess[:,1] = 0.05*srandu(nwalkers)
        guess[:,0] = 0.038*(1.0 + 0.1*srandu(nwalkers))
        guess[:,1] = 0.001*srandu(nwalkers)

        self._guess=guess
        return guess

    def _doplots(self):
        import biggles
        import esutil as eu

        tab=biggles.Table(3,2)

        ind=numpy.arange(self._lnprobs.size)

        sh1vals = self._trials[:,0]
        sh2vals = self._trials[:,1]

        sh1p = biggles.FramedPlot()
        sh1p.add( biggles.Curve(ind, sh1vals) )
        sh1p.ylabel=r'$\gamma_1$'

        sh2p = biggles.FramedPlot()
        sh2p.add( biggles.Curve(ind, sh2vals) )
        sh2p.ylabel=r'$\gamma_2$'

        sh1std = sh1vals.std()
        binsize1=0.2*sh1std
        sh1h = eu.plotting.bhist(sh1vals, binsize=binsize1,show=False)
        sh1h.xlabel=r'$\gamma_1$'

        sh2std = sh2vals.std()
        binsize2=0.2*sh2std
        sh2h = eu.plotting.bhist(sh2vals, binsize=binsize2,show=False)
        sh2h.xlabel=r'$\gamma_2$'


        likep = biggles.FramedPlot()
        likep.add( biggles.Curve(ind, self._lnprobs) )
        likep.ylabel='ln( prob )'

        
        tab[0,0] = sh1p
        tab[0,1] = sh1h
        tab[1,0] = sh2p
        tab[1,1] = sh2h
        tab[2,0] = likep

        tab.show()

    def _setup_progress(self,ntot):
        from progressbar import ProgressBar
        self._ntot=ntot
        self._count=0
        self._progress_bar=ProgressBar(width=70)

    def _update_progress(self):
        frac=float(self._count)/self._ntot
        self._progress_bar.update(frac=frac)
        self._count += 1

class SHDirect(SHDirectBase):
    def __init__(self, image_list, ivar_list, psf_list, model_list, ares_list, nwalkers, burnin, nstep, **keys):
        super(SHDirect,self).__init__(image_list, ivar_list, psf_list, model_list, ares_list, nwalkers, burnin, nstep, **keys)
        
        self._go()

        print 'shear:',self._result['shear']
        err=numpy.sqrt( numpy.diag(self._result['shear_cov']) )
        print 'shear err:',err
        print 'arate:',self._result['arate']


