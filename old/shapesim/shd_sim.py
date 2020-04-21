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
        #self.npair=100 # we do two times this, pairs at 90 degrees
        #self.npair=10 # we do two times this, pairs at 90 degrees

        self._set_gprior()

        simpars=self.get('simpars',{})
        self.shapesim = ShapeSim(self['sim'], **simpars)

    def process_trials(self, is2, is2n, itrial=None):
        if itrial is not None:
            raise ValueError("implement itrial")

        print 'making images'
        self._make_images(is2, is2n)

        print 'running mcmc'
        if self['method']=='emcee':
            shd = SHDirect(self.image_list,
                           self.ivar_list,
                           self.psf_list,
                           self.model_list,
                           self.ares_list,
                           nwalkers=self['nwalkers'],
                           burnin=self['burnin'],
                           nstep=self['nstep'],
                           make_plots=self['make_plots'],
                           mca_a=self['mca_a'],
                           gprior=self.gprior,
                           draw_prior=False,
                           apply_gprior=False)
        elif self['method']=='mh':
            shd = SHDirectMH(self.image_list,
                             self.ivar_list,
                             self.psf_list,
                             self.model_list,
                             self.ares_list,
                             burnin=self['burnin'],
                             nstep=self['nstep'],
                             make_plots=self['make_plots'])
                             #gprior=self.gprior)

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

        self._nimages=len(self._image_list)

        self._make_plots=keys.get('make_plots',False)

        self._gprior=keys.get('gprior',None)
        # default don't actually apply it
        self._apply_gprior=keys.get('apply_gprior',False)

        # two ellipticity components for all (the shear)
        # plus centroid(2) size and flux
        self._npars = 2 + 4*self._nimages

        self.update(keys)

    def get_result(self):
        return self._result
    def get_trials(self):
        return self._trials
    def get_lnprobs(self):
        return self._lnprobs

    def get_arate(self):
        raise RuntimeError("over-ride")

    def _go(self):
        raise RuntimeError("over-ride")

    def _do_trials(self):
        raise RuntimeError("over-ride")

    def _get_guess(self):
        raise RuntimeError("over-ride")

    def _calc_lnprob(self, pars):


        try:
            shear=lensing.Shear(g1=pars[0], g2=pars[1])
        except ShapeRangeError:
            return LOWVAL

        # in e space
        pars1=numpy.zeros(6)
        pars1[2] = shear.e1
        pars1[3] = shear.e2

        lnprob=0.0
        for i in xrange(self._nimages):
            ipars=2 + 4*i
            
            image=self._image_list[i]
            ivar=self._ivar_list[i]
            psf=self._psf_list[i]
            model=self._model_list[i]

            pars1[0] = pars[ipars+0]
            pars1[1] = pars[ipars+1]
            pars1[4] = pars[ipars+2]
            pars1[5] = pars[ipars+3]

            lnprob1 = self._calc_lnprob1(image, ivar, psf, model, pars1)

            lnprob += lnprob1

            # don't want the centroids to move around
            arow=self._ares_list[i]['wrow']
            acol=self._ares_list[i]['wcol']
            w=0.1**2
            cp = -0.5*(pars1[0]-arow)**2/w - 0.5*(pars1[1]-acol)**2/w

            lnprob += cp

        # is this right?
        if self._gprior is not None and self._apply_gprior:
            gp = self._get_lngprior(shear.g1, shear.g2)
        else:
            gp=0.0

        self._update_progress()
        return lnprob

    def _calc_lnprob1(self, image, ivar, psf, model, pars):
        """
        pars in e space
        """
        gmix0=self._get_gmix0(pars, model)
        gmix = gmix0.convolve(psf)

        logprob = self._get_loglike_c1(image, ivar, gmix)

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
            gmix0=gmix_image.GMixCoellip(epars)
        else:
            raise ValueError("implement type '%s'" % model)
        return gmix0


    def _calc_result(self):
        if self._gprior is None or self._apply_gprior==False:
            print 'calculating result with no gprior'
            self._calc_result_nogprior()
        else:
            print 'calculating result with gprior'
            self._calc_result_gprior()
    
    def _calc_result_nogprior(self):
        import mcmc

        g=numpy.zeros(2)
        gcov=numpy.zeros((2,2))

        pars,pcov = mcmc.extract_stats(self._trials)

        g[:] = pars[0:0+2]
        gcov[:,:] = pcov[0:0+2, 0:0+2]
 
        arate = self.get_arate()

        self._result={'g':g,
                      'gcov':gcov,
                      'gsens':numpy.ones(2),
                      'pars':pars,
                      'pcov':pcov,
                      'arate':arate}

    def _calc_result_gprior(self):
        import mcmc

        g=numpy.zeros(2)
        gcov=numpy.zeros((2,2))
        gsens = numpy.zeros(2)

        g1vals=self._trials[:,0]
        g2vals=self._trials[:,1]

        prior = self._gprior(g1vals,g2vals)
        dpri_by_g1 = self._gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self._gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        pars,pcov = mcmc.extract_stats(self._trials)

        g[:] = pars[0:0+2]
        gcov[:,:] = pcov[0:0+2, 0:0+2]

        g1diff = g[0]-g1vals
        g2diff = g[1]-g2vals

        w,=numpy.where(prior > 0)
        if w.size == 0:
            raise ValueError("no prior values > 0!")

        gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
        gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()
 
        arate = self.get_arate()

        self._result={'g':g,
                      'gcov':gcov,
                      'gsens':gsens,
                      'pars':pars,
                      'pcov':pcov,
                      'arate':arate}



    def _doplots(self):
        import biggles
        import esutil as eu

        tab=biggles.Table(3,2)

        ind=numpy.arange(self._lnprobs.size)

        res=self._result

        sh1vals = self._trials[:,0]
        sh2vals = self._trials[:,1]
        sh1vals /= res['gsens'][0]
        sh2vals /= res['gsens'][1]

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
    def __init__(self, image_list, ivar_list, psf_list, model_list, ares_list, **keys):
        super(SHDirect,self).__init__(image_list, ivar_list, psf_list, model_list, ares_list, **keys)
        
        self._nwalkers=keys.get('nwalkers',50)
        self._burnin=keys.get('burnin',100)
        self._nstep=keys.get('nstep',100)
        self._mca_a = keys.get('mca_a', 2.0)

        # default draw from input prior
        self._draw_gprior=keys.get('draw_gprior',True)


        # two ellipticity components for all (the shear)
        # plus centroid(2) size and flux

        min_walkers=4*self._npars
        #min_walkers=8*self._npars
        if self._nwalkers < min_walkers:
            self._nwalkers=min_walkers

        print 'nwalkers:',self._nwalkers


        self._go()

        shear=self._result['g']
        err=numpy.sqrt( numpy.diag(self._result['gcov']) )
        print 'shear: ',shear
        print 'err:   ',err

        sens=self._result['gsens']
        shcorr=shear/sens
        print 'gsens:         ',sens
        print 'sens corrected:',shcorr
        print 'err            ',err/sens
        print 'arate:',self._result['arate']

    def get_arate(self):
        return self._arate

    def _go(self):
        self._sampler=self._do_trials()

        self._trials  = self._sampler.flatchain

        lnprobs = self._sampler.lnprobability.reshape(self._nwalkers*self._nstep)
        self._lnprobs = lnprobs - lnprobs.max()

        arates = self._sampler.acceptance_fraction
        self._arate = arates.mean()

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
            ii = 2 + 4*i

            ai=self._ares_list[i]
            if ai['whyflag'] != 0:
                raise ValueError("admom with flags != 0")

            aT = ai['Irr'] + ai['Icc']
            counts=self._image_list[i].sum()

            guess[:,ii+0] = ai['wrow']*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+1] = ai['wcol']*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+2] = aT*(1. + 0.1*srandu(nwalkers))
            guess[:,ii+3] = counts*(1. + 0.1*srandu(nwalkers))

        if self._draw_gprior:
            g1rand,g2rand=self._gprior.sample2d(nwalkers)
            guess[:,0] = g1rand
            guess[:,1] = g2rand
        else:
            #guess[:,0] = 0.05*srandu(nwalkers)
            #guess[:,1] = 0.05*srandu(nwalkers)
            guess[:,0] = 0.08*(1.0 + 0.1*srandu(nwalkers))
            guess[:,1] = 0.001*srandu(nwalkers)

        self._guess=guess
        return guess



class SHDirectMH(SHDirectBase):
    def __init__(self, 
                 image_list, 
                 ivar_list, 
                 psf_list, 
                 model_list, 
                 ares_list, 
                 **keys):
        """

        need to determine the best model for each with a pre-run of some other
        code!

        """
        super(SHDirectMH,self).__init__(image_list, ivar_list, psf_list, model_list, ares_list, **keys)

        self._burnin=keys.get('burnin',200)
        self._nstep=keys.get('nstep',200)
        self._draw_gprior=keys.get('draw_gprior',False)

        self._go()

        shear=self._result['g']
        err=numpy.sqrt( numpy.diag(self._result['gcov']) )
        print 'shear: ',shear
        print 'err:   ',err

        sens=self._result['gsens']
        shcorr=shear/sens
        print 'gsens:         ',sens
        print 'sens corrected:',shcorr
        print 'err            ',err/sens
        print 'arate:',self.get_arate()

    def get_arate(self):
        return self._arate


    def loglike(self,pars):
        return self._calc_lnprob(pars)
    def step(self,pars):
        from numpy.random import randn
        return pars + self._stepsize*randn(self._npars)

    def _go(self):
        self._set_guess_stepsize()
        self._do_trials()

        if self._make_plots:
            self._doplots()

    def _do_trials(self):
        import mcmc
        mh=mcmc.MCMC(self)

        self._setup_progress(self._burnin)
        burn_res=mh.run(self._burnin, self._guess)

        self._setup_progress(self._nstep)
        res=mh.run(self._nstep, burn_res['pars'][-1,:])

        self._burn_res=burn_res
        self._res=res

        self._trials=res['pars']
        self._lnprobs=res['loglike']
        self._arate=res['accepted'].astype('f8').sum()/res['accepted'].size

        self._calc_result()

    def _set_guess_stepsize(self):
        
        npars=self._npars

        guess=numpy.zeros(npars)
        # these for S/N=700
        stepsize=numpy.zeros(npars)

        for i in xrange(self._nimages):
            ii = 2 + 4*i

            ai=self._ares_list[i]
            if ai['whyflag'] != 0:
                raise ValueError("admom with flags != 0")

            aT = ai['Irr'] + ai['Icc']
            counts=self._image_list[i].sum()

            guess[ii+0] = ai['wrow']*(1. + 0.1*srandu())
            guess[ii+1] = ai['wcol']*(1. + 0.1*srandu())
            guess[ii+2] = aT*(1. + 0.1*srandu())
            guess[ii+3] = counts*(1. + 0.1*srandu())

            stepsize[ii+0] = 0.1
            stepsize[ii+1] = 0.1
            # update for s/n
            stepsize[ii+2] = 0.1
            stepsize[ii+3] = 0.05

        if self._draw_gprior:
            g1rand,g2rand=self._gprior.sample2d(1)
            guess[0] = g1rand[0]
            guess[1] = g2rand[0]
        else:
            #guess[0] = 0.1*srandu()
            #guess[1] = 0.1*srandu()
            guess[0] = 0.038
            guess[1] = -0.005

        stepsize[0]=0.1
        stepsize[1]=0.1

        self._stepsize=stepsize
        self._guess=guess
        return guess


