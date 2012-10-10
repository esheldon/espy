"""
TODO:

    - The coarse grid produces a large variance in the mean estimated shear,
    of order the grid size.

    - how to organize g1,g2 sampling?  nring=20 and then how many e values?  I
    normally just increased the number of repeats for each ring but here we
    prefer instead to get a new place in the g1,g2 plane to sample the prior.
    But this gets expensive generating the models.

    - Need to implement prior and sensitivity.


    Seeing large bias with emcee:
        - I verified the loglike calculation is the same
        - made a fixcen version; a small test shows lower bias
        - noted that low burnin did have higher bias but
        need to explore
        - re-running mcbayes-gg10r05 3 5 with tight prior on centroid (1.e-3)
        to see how well it works: better but still crap
        - ran fixcen gg10r06: still bad.

        - what about using 0.5 acceptance?

        Tried old mcmc gg10r03. Note applied prior *after*
            looks good to S/N = 10, bad at 5
        Running gg10r07 with prior *during* (like gg10r06 emcee)
            If this looks bad, it probably means I need to just calculate
            the likelihood surface and use the prior afterward.  The
            1/P dP/de is unstable.  This also probably means I will
            go back to using the affine invariant sampler as well,
            but only using it to explore the likelihood surface.
            Would have to test that too.

            If it looks good, then the problem is somehow in particulars
            of the class.  Should run with exact same numbers to see
            if agreement is found.

            It looks OK.  So, I'm going to try just plugging emcee
            into that class to see if it works. First make sure my
            faster version agrees: looks OK
            This is gg10r08. Just running s2 001 s2n 002 as a test:
                looks ok, doing the rest of this run
            Looks bad, constant offset below.  I'm guessing this
            is somehow related to the burnin period and that I
            guessed (0,0) plus noise for the ellip.
            re-running 000 012 with longer burnin just to see
                looks much better, so what's up?
                - try again with 1000 walkers and old burnin 000 011
                - try again with 10 walkers and old burnin   000 010
                    this worked best

            Doing a new run gg10r09 with 10 walkers and still
                4000 burnin, 5000 total.

            Note still choosing (0,0) centered start; with 10 walkers
            should we go back to a good guess?

        - going back to emcee now with 10 walkers.  No other changes!
            gg10r10  fixcen still.  Looks like shit, maybe because
            of the e guess.  Using same as other, zero centered.  Still shit.
                re-running 0 1 using new EmceeFitterFixCen class made
                from the mcmcfitter; looks better, re-running all now
            I've also created a new EmceeFitter class directly from
            the EmceeFitterFixCen class, so we shall see.  will run
            gg10r11 for that (testing 0 1 now s/n=15)

            Both look ok except at s/n of 10.  I wonder if increasing the
            burning would help?  Will try burn in 1000 per "in place" for
            gg10r11, s2 3 s/n 0: no difference really
                it looks biased!

            gg08r03 - old prior, emcee, fixed cen, prior after
                - doesn't look that good
            gg08r04 - old prior, emcee, fixed cen, prior during
                - running

            Idea from Anze: run with higher "temperature" to make sure we
            explore the tails.  Waiting to see if 'during' gg08r04 looks
            better, if so will try it on that, otherwise will try on
            new prior.
                mcbayes-gg10r13: temp=2, free cen, new prior
                - had wrong formula

            - along same lines, maybe need more points in after-burnin to see
            tail?  Try using 400 per in gg10r14
                - looks good
                - gg10r15 for averaging
                - gg10r16 for averaging
        

        - trying grid search bayesfit
            - gg08r07 : original prior, faster C code.
                - looks good
            - gg08r08: T marginilization
                - shows some bias at s/n=10 and 15
                - looks reasonable for higher s/n for large galaxies, but
                like total crap at high s/n for smaller galaxies
                This could be because the grid is too crude to get any
                measure of the distribution in T without some noise to
                smear things out
            - gg10r01 : new prior
                doesn't look as good as old prior.  Maybe this prior is
                just harder to get right?  Will try emcee on old prior in
                gg08r03, prior after
            - gg10r02 : new prior, T grid
                - running

    - trying exp/turb models now
        - The convolved models are coming out much more compact than true exp,
        maybe I am just making them too small or maybe my effective size is
        too small compared to real ones.

    NEED to compare mcmc to emcee version: what is the difference?
        - differences in formula?
        - burnin?
"""
import os
import numpy
from numpy import sqrt, cos, sin, exp, log, pi, zeros, ones, empty, \
        random, where, array, linspace, diag
from numpy import tanh, arctanh
from numpy.random import randn
from numpy.random import random as randu
from . import shapesim
import lensing
import fimage
from fimage.convolved import NoisyConvolvedImage
import gmix_image
import images
import esutil as eu
from esutil.misc import wlog

import math
import time

from sys import stderr

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

        self.gprior = GPrior(A=self.simc['A'],
                             B=self.simc['B'],
                             C=self.simc['C'],
                             D=self.simc['D'])
        if 'mcbayes' not in run:
            if 'n_Tgrid' in self:
                self.fitter = BayesFitterFixCen(self.gprior, 
                                                self['n_ggrid'],
                                                self['gmin'],
                                                self['gmax'],
                                                self['n_Tgrid'],
                                                self['Tmin'],
                                                self['Tmax'])
            else:
                self.fitter = BayesFitterFixTCen(self.gprior, 
                                                 self['n_ggrid'],
                                                 self['gmin'],
                                                 self['gmax'])


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
        out = zeros(nellip, dtype=self.out_dtype())
        out['s2'] = s2
        self['s2']=s2

        if dolog:
            wlog('s2n:',s2n,'s2n_psf:',s2n_psf,'s2n_method:',s2n_method)


        for i,g in enumerate(gvals):

            ellip=lensing.util.g2e(g)

            ci_nonoise = self.shapesim.get_trial(s2, ellip, theta)

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

            # self.fitter is used in here
            self._run_fitter(ci)

            res = self.fitter.get_result()

            out['g'][i,:] = res['g']
            out['gsens'][i,:] = res['gsens']
            out['gcov'][i,:,:] = res['gcov']
            if 'arate' in res:
                out['arate'][i] = res['arate']


        if dowrite:
            shapesim.write_output(self['run'], is2, is2n, out, itrial=itheta,
                         fs=self.fs)
        return out

    def _measure_gmix_psf(self, ci):
        import admom
        psfres = admom.admom(ci.psf,
                             ci['cen_uw'][0],
                             ci['cen_uw'][1], 
                             guess=2.,
                             nsub=1)

        if self['ngauss_psf']==1:
            psf=[{'p':1,
                  'row':psfres['row'],
                  'col':psfres['col'],
                  'irr':psfres['Irr'],
                  'irc':psfres['Irc'],
                  'icc':psfres['Icc']}]
        elif self['ngauss_psf']==3:

            prior=zeros(10)
            width=zeros(10) + 100

            Tpsf=psfres['Irr']+psfres['Icc']
            Tmax = Tpsf*8.3
            Tfrac1 = 1.7/8.3
            Tfrac2 = 0.8/8.3
            prior[0]=psfres['row']
            prior[1]=psfres['col']
            prior[2]=psfres['e1']
            prior[3]=psfres['e2']
            prior[4] = Tmax
            prior[5] = Tfrac1 
            prior[6] = Tfrac2

            prior[7] = 0.08
            prior[8] = 0.38
            prior[9] = 0.53

            e1start=prior[2]
            e2start=prior[3]
            prior[2],prior[3] = randomize_e1e2(e1start,e2start)

            prior[4] += prior[4]*0.05*(randu()-0.5)
            prior[5] += prior[5]*0.05*(randu()-0.5)
            prior[6] += prior[6]*0.05*(randu()-0.5)
            prior[7] += prior[7]*0.05*(randu()-0.5)
            prior[8] += prior[8]*0.05*(randu()-0.5)
            prior[9] += prior[9]*0.05*(randu()-0.5)

            gm = gmix_image.GMixFitCoellip(ci.psf, ci['skysig'],
                                           prior,width,
                                           Tpositive=True)
            psf=gm.get_gmix()
        else:
            raise ValueError("bad ngauss_psf: %s" % self['ngauss_psf'])

        return psf
    def _run_fitter(self, ci):
        """
        cheat on psf, T and cen for now
        """
        import admom

        if 'Ttrue' in ci:
            T=ci['Ttrue']
        else:
            cov=ci['covtrue']
            T = cov[0] + cov[2]
        cen = ci['cen']

        #print 's2:',self['s2'],'Tpsf/Tobj:',ci['Ttrue_psf']/ci['Ttrue']
        #print ci.image.shape
        psf=self._measure_gmix_psf(ci)


        temp=self.get('temp',None)
        if self['fixcen']:
            raise ValueError("Need to check the fixcen code for recent changes")
            self.fitter=EmceeFitterFixCen(ci.image,
                                          1./ci['skysig'],
                                          psf,
                                          ci['cen'],
                                          T,
                                          self.gprior,
                                          nwalkers=self['nwalkers'],
                                          nstep=self['nstep'], 
                                          burnin=self['burnin'],
                                          temp=temp,
                                          when_prior=self['when_prior'])

        elif True:
            cenprior=CenPrior(ci['cen'], [0.1]*2)
            self.fitter=EmceeFitter(ci.image,
                                    1./ci['skysig'],
                                    psf,
                                    cenprior,
                                    T,
                                    self.gprior,
                                    self['fitmodel'],
                                    nwalkers=self['nwalkers'],
                                    nstep=self['nstep'], 
                                    burnin=self['burnin'],
                                    temp=temp, # need to implement
                                    when_prior=self['when_prior'])

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


    def out_dtype(self):
        dt=[('s2n_admom','f8'),
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
            ('gcov','f8',(2,2))]
        if 'bayesfit' not in self['run']:
            dt += [('arate','f8')]

        return dt



class EmceeFitter:
    def __init__(self, 
                 image, 
                 ierr, 
                 psf,
                 cenprior,
                 T,
                 gprior,
                 model,
                 nwalkers=10,
                 nstep=100, 
                 burnin=400,
                 temp=None,
                 when_prior='during'):
        """
        mcmc sampling of posterior.  Amplitude is analytically marginalized.

        parameters
        ----------
        image:
            image as a numpy array
        ierr:
            Error per pixel
        psf:
            The psf gaussian mixture
        cenprior:
            The center prior object.
        T:
            Starting value for ixx+iyy of main component
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

        # cen1,cen2,e1,e2,T
        self.npars=5

        self.image=image
        self.ierr=float(ierr)
        self.T=T
        self.cenprior=cenprior

        self.model=model

        self.nwalkers=nwalkers
        self.nstep=nstep
        self.burnin=burnin
        self.gprior=gprior
        self.when_prior=when_prior

        self._set_psf(psf)

        self.Anorm = float(1)

        self.tpars=zeros(6,dtype='f8')

        if temp is not None and when_prior=='after':
            raise ValueError("don't support temp for when_prior after yet")
        self.temp=temp

        self._go()

    def get_result(self):
        return self._result

    def _go(self):
        import emcee

        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._calc_lnprob,
                                        a=2.5)
        
        guess=self._get_guess()

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
        """
        # hard priors first
        g1,g2=pars[2],pars[3]
        T=pars[4]

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0]=pars[0]
        self.tpars[1]=pars[1]
        self.tpars[2]=e1
        self.tpars[3]=e2
        self.tpars[4]=T
        self.tpars[5]=1.

        logprob = self._get_loglike_c(self.tpars)

        if self.when_prior=='during':
            gp = self._get_lngprior(g1,g2)
            logprob += gp

        cp = self.cenprior.lnprob(pars[0:2])
        logprob += cp

        if self.temp is not None:
            logprob /= self.temp
        return logprob

    def _get_loglike_c(self, pars):
        """
        pars is *full*
        """

        if self.model=='gexp':
            gmix0=gmix_image.GMixExp(pars)
        elif self.model=='gdev':
            gmix0=gmix_image.GMixDev(pars)
        elif self.model=='gauss':
            gmix0=gmix_image.GMixCoellip(pars)
        else:
            raise ValueError("bad model: '%s'" % self.model)

        gmix=gmix0.convolve(self.psf_GMix)

        loglike,flags=\
            gmix_image.render._render.loglike(self.image, 
                                              gmix,
                                              self.Anorm,
                                              self.ierr)

        if flags != 0:
            return LOWVAL
        return loglike

    def _get_lngprior(self, g1, g2):
        g=sqrt(g1**2 + g2**2)
        gp = self.gprior.prior_gabs_scalar(g)
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

        g1vals=self.trials[:,2]
        g2vals=self.trials[:,3]

        prior = self.gprior(g1vals,g2vals)
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        if self.when_prior=='after':
            # we need to multiply each by the prior
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de

            wt=None
            if self.temp is not None:
                wt=exp(self.lnprobs*(1.-1./self.temp))

            g, gcov = mcmc.extract_stats(self.trials[:,2:2+2],weights=wt)

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            if wt is None:
                gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
                gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()
            else:
                wsum=wt[w].sum()
                sum1=(wt[w]*g1diff[w]*dpri_by_g1[w]/prior[w]).sum()
                sum2=(wt[w]*g2diff[w]*dpri_by_g2[w]/prior[w]).sum()

                gsens[0]= 1.-sum1/wsum
                gsens[1]= 1.-sum2/wsum


 
        arates = self._emcee_sampler.acceptance_fraction
        arate = arates.mean()
        #print 'acceptance rate:',w.size/float(self.trials.size)
        self._result={'g':g,'gcov':gcov,'gsens':gsens,'arate':arate}
        #wlog("arate:",self._result['arate'])



    def _get_guess(self):
        guess=zeros( (self.nwalkers,self.npars) )

        guess[:,0]=self.cenprior.cen[0] + 0.01*(randu(self.nwalkers)-0.5)
        guess[:,1]=self.cenprior.cen[1] + 0.01*(randu(self.nwalkers)-0.5)

        # guess for g1,g2 is (0,0) with some scatter
        guess[:,2]=0.1*(randu(self.nwalkers)-0.5)
        guess[:,3]=0.1*(randu(self.nwalkers)-0.5)

        # guess for T is self.T with scatter
        guess[:,4] = self.T + self.T*0.1*(randu(self.nwalkers)-0.5)

        return guess

    def _set_psf(self, psf):
        self.psf_gmix = psf

        if not isinstance(psf[0],dict):
            raise ValueError("psf must be list of dicts")
        self.psf_pars = gmix_image.gmix2pars(psf)
        self.psf_GMix=gmix_image.GMix(psf)

    def _doplots(self):

        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(5,2)

        cen1vals=self.trials[:,0]
        cen2vals=self.trials[:,1]
        g1vals=self.trials[:,2]
        g2vals=self.trials[:,3]
        Tvals=self.trials[:,4]

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

        likep = biggles.FramedPlot()
        likep.add( biggles.Curve(ind, self.lnprobs) )
        likep.ylabel='ln( prob )'


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))

        res=self.get_result()
        print 'acceptance rate:',res['arate']
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]

        cenw = cen1vals.std()
        cen_bsize=cenw*0.2
        hplt_cen0 = eu.plotting.bhist(cen1vals,binsize=cen_bsize,
                                      color='blue',
                                      show=False)
        hplt_cen = eu.plotting.bhist(cen2vals,binsize=cen_bsize,
                                     color='red',
                                     show=False, plt=hplt_cen0)
        hplt_cen.add(key)

        bsize1=errs[0]*0.2
        bsize2=errs[1]*0.2
        hplt_g1 = eu.plotting.bhist(g1vals,binsize=bsize1,
                                  show=False)
        hplt_g2 = eu.plotting.bhist(g2vals,binsize=bsize2,
                                  show=False)

        Tsdev = Tvals.std()
        Tbsize=Tsdev*0.2
        hplt_T = eu.plotting.bhist(Tvals,binsize=Tbsize,
                                  show=False)

        hplt_g1.xlabel=r'$g_1$'
        hplt_g2.xlabel=r'$g_2$'
        hplt_T.xlabel='T'
        tab[0,0] = burn_cen
        tab[1,0] = burn_g1
        tab[2,0] = burn_g2
        tab[3,0] = burn_T

        tab[0,1] = hplt_cen
        tab[1,1] = hplt_g1
        tab[2,1] = hplt_g2
        tab[3,1] = hplt_T
        tab[4,0] = likep
        tab.show()

class EmceeFitterFixCen:
    def __init__(self, 
                 image, 
                 ierr, 
                 psf,
                 cen,
                 T,
                 gprior,
                 nwalkers=10,
                 nstep=100, 
                 burnin=400,
                 temp=None,
                 when_prior='during'):
        """
        parameters
        ----------
        image:
            image as a numpy array
        ierr:
            Error per pixel
        psf:
            The psf gaussian mixture
        cen:
            The center.
        T:
            Starting value for ixx+iyy of main component
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

        # e1,e2,T
        self.npars=3

        self.image=image
        self.ierr=float(ierr)
        self.T=T
        self.cen=cen

        self.nwalkers=nwalkers
        self.nstep=nstep
        self.burnin=burnin
        self.gprior=gprior
        self.when_prior=when_prior

        self._set_psf(psf)

        self.Anorm = float(1)

        self.tpars=zeros(6,dtype='f8')

        if temp is not None and when_prior=='after':
            raise ValueError("don't support temp for when_prior after yet")
        self.temp=temp

        self._go()

    def get_result(self):
        return self._result

    def _go(self):
        import emcee

        guess=zeros( (self.nwalkers,self.npars) )

        # guess for g1,g2 is (0,0) with some scatter
        guess[:,0]=0.1*(randu(self.nwalkers)-0.5)
        guess[:,1]=0.1*(randu(self.nwalkers)-0.5)

        # guess for T is self.T with scatter
        guess[:,2] = self.T + self.T*0.1*(randu(self.nwalkers)-0.5)

        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._calc_lnprob,
                                        a=4.)
        

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
        """
        # hard priors first
        g1,g2=pars[0],pars[1]
        T=pars[2]

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0]=self.cen[0]
        self.tpars[1]=self.cen[1]
        self.tpars[2]=e1
        self.tpars[3]=e2
        self.tpars[4]=T
        self.tpars[5]=1.

        logprob = self._get_loglike_c(self.tpars)

        if self.when_prior=='during':
            gp = self._get_lngprior(g1,g2)
            logprob += gp

        if self.temp is not None:
            logprob /= self.temp
        return logprob




    def _get_loglike_c(self, pars):
        #t0=time.time()

        loglike,flags=\
            gmix_image.render._render.loglike_coellip(self.image, 
                                                      pars, 
                                                      self.psf_pars, 
                                                      self.Anorm,
                                                      self.ierr)
        if flags != 0:
            return LOWVAL
        #self.logl_tm += time.time()-t0
        return loglike

    def _get_lngprior(self, g1, g2):
        g=sqrt(g1**2 + g2**2)
        gp = self.gprior.prior_gabs_scalar(g)
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

        g1vals=self.trials[:,0]
        g2vals=self.trials[:,1]

        prior = self.gprior(g1vals,g2vals)
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        if self.when_prior=='after':
            # we need to multiply each by the prior
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de
            wt=None
            if self.temp is not None:
                wt=exp(self.lnprobs*(1.-1./self.temp))
                #print wt
            g, gcov = mcmc.extract_stats(self.trials[:,0:2], weights=wt)
            """
            g2, gcov2 = mcmc.extract_stats(self.trials[:,0:2])
            print g
            print g2
            print gcov
            print gcov2
            """

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            if wt is None:
                gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
                gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()
            else:
                wsum=wt[w].sum()
                sum1=(wt[w]*g1diff[w]*dpri_by_g1[w]/prior[w]).sum()
                sum2=(wt[w]*g2diff[w]*dpri_by_g2[w]/prior[w]).sum()

                gsens[0]= 1.-sum1/wsum
                gsens[1]= 1.-sum2/wsum

 
        arates = self._emcee_sampler.acceptance_fraction
        arate = arates.mean()
        #print 'acceptance rate:',w.size/float(self.trials.size)
        self._result={'g':g,'gcov':gcov,'gsens':gsens,'arate':arate}
        #wlog("arate:",self._result['arate'])




    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _doplots(self):
        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(4,2)

        ind=numpy.arange(self.trials[:,0].size)

        g1vals=self.trials[:,0]
        g2vals=self.trials[:,1]
        Tvals=self.trials[:,2]

        burn_g1=biggles.FramedPlot()
        burn_g1.add(biggles.Curve(ind, g1vals))
        burn_g1.ylabel=r'$g_1$'

        burn_g2=biggles.FramedPlot()
        burn_g2.add(biggles.Curve(ind, g2vals))
        burn_g2.ylabel=r'$g_2$'

        burn_T=biggles.FramedPlot()
        burn_T.add(biggles.Curve(ind, Tvals))
        burn_T.ylabel='T'

        likep = biggles.FramedPlot()
        likep.add( biggles.Curve(ind, self.lnprobs) )
        likep.ylabel='ln( prob )'


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))

        res=self.get_result()
        print 'acceptance rate:',res['arate']
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]

        bsize1=errs[0]*0.2
        bsize2=errs[1]*0.2
        hplt_g1 = eu.plotting.bhist(g1vals,binsize=bsize1,
                                    show=False)
        hplt_g2 = eu.plotting.bhist(g2vals,binsize=bsize2,
                                  show=False)

        Tsdev = Tvals.std()
        Tbsize=Tsdev*0.2
        wlog("T binsize:",Tbsize)
        hplt_T = eu.plotting.bhist(Tvals,binsize=Tbsize, show=False)

        hplt_g1.xlabel=r'$g_1$'
        hplt_g2.xlabel=r'$g_2$'
        hplt_T.xlabel='T'
        tab[0,0] = burn_g1
        tab[1,0] = burn_g2
        tab[2,0] = burn_T

        tab[0,1] = hplt_g1
        tab[1,1] = hplt_g2
        tab[2,1] = hplt_T
        tab[3,0] = likep
        tab.show()




class MCMCFitter:
    """
    likelihood is exp(0.5*A*B^2)

    A is 
        sum( (model/err)^2 ) which we fix to sum ((ydata/err)^2)
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that the
    relative height of the likelihood is valid given the other simplifications
    we have made. We choose 1 because, using the re-normalization below,
    the value actually cancels.

    We must normalize the model according to the A condition.  create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S) 

    If using "exp", 3 gaussian mixture model, you can examine run
    gmix-fit-et10r99

        Moving parts

            - Want cen with tight prior
            - Want e1,e2. prior applied *after* exploring likelihood surface
            because the prior can be cuspy and this confuses the MCMC
            - Want a T that represents the size, no prior at all but a limited
            range

        Fixed parts of model

            - Want a delta function; fixed at zero for now: can only detect
            non-zero at *very* high S/N
            - Want a fixed F that is a multiple of T for the second gauss.  For
            s2 ~ 1, this appears to be 

                ~3.8 for s2 ~ 1

            - Want fixed p values. these  are only meaningful relative since we
            marginalize over amplitude. For s2=1

                0.379 0.563, 0.0593
    """

    def __init__(self, 
                 image, 
                 ierr, 
                 psf,
                 guess,
                 s2n, 
                 gprior,
                 cenwidth,
                 when_prior='during',
                 nstep=2000, 
                 burnin=200):
        """
        parameters
        ----------
        image: 2D numpy array
            image to be fit
        ierr: scalar or array

            1/(Error per pixel) or an array of inverse sqrt(variance) same
            shape as image.  sqrt(var) is used instead of var because of the
            way the analytic normalization marginalization occurs.
        psf: gaussian mixture
            A list of dictionaries
        guess:
            Starting guess for pars

            [e1,e2,cen1,cen2,T]
        s2n:
            S/N value; will be used to create step sizes.
        gprior:
            The prior on the g1,g2 surface.  A GPrior object
            or similar.
        cenwidth: scalar
            The prior width on the center.
        nstep: optional
            Number of steps in MCMC chain, default 2000
        burnin:
            Number of burn in steps, default 200.
        when_prior: string
            'after' to search the pure likelihood surface and apply
            the prior afterward, otherwise during the seach.
        """

        # e1,e2,cen1,cen2,T
        self.npars=5

        self.image=image
        self.ierr=ierr
        self._set_psf(psf)
        self.guess=array(guess,dtype='f8')
        self.s2n=s2n
        self.gprior=gprior
        self.cenwidth=cenwidth

        self.nstep=nstep
        self.burnin=burnin
        self.when_prior=when_prior

        gstepsize=self._get_gstepsize(s2n)
        Tstepsize=self._get_Tstepsize(s2n)
        cen_stepsize=1.6*cenwidth

        steps=[cen_stepsize,cen_stepsize,gstepsize,gstepsize,Tstepsize]
        #print 'guess:',guess
        #print 'steps:',steps
        #stop
        self.stepsize=array(steps,dtype='f8')


        if isinstance(ierr,numpy.ndarray):
            if ierr.shape != image.shape:
                raise ValueError("ierr image shape %s does not match image "
                                 "%s " % (ierr.shape,image.shape))
        if self.guess.size != self.npars:
            raise ValueError("expected guess with length %d" % self.npars)

        self.A = 1


        # for the p
        self.tpars=zeros(self.npars+1)
        self.model=zeros(self.image.shape)

        self._go()

    def get_result(self):
        return self._result

    def _go(self):
        import mcmc
        mcmc=mcmc.MCMC(self)
        
        #print self.guess
        self.trials = mcmc.run(self.nstep, self.guess)

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if False:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop

    def loglike(self, pars):
        """
        Packaged up for the MCMC code
        """
        # hard priors first
        try:
            # make sure we pass shape tests
            self.shape=lensing.Shear(g1=pars[2], g2=pars[3])
        except:
            return LOWVAL

        if pars[4] < 0:
            return LOWVAL

        if self.when_prior=='after':
            return self._get_loglike(pars)
        else:
            return self._get_logprob(pars)

    def step(self, pars):
        """
        public for the MCMC code
        """
        nextpars=zeros(self.npars)
        """
        f=random.random()
        if f > 0.5:
            nextpars = pars + self.small_stepsize*randn(self.npars)
        else:
            nextpars = pars + self.stepsize*randn(self.npars)
        """
        nextpars = pars + self.stepsize*randn(self.npars)
        #print 'nextpars:',nextpars
        return nextpars

    def _get_logprob(self, pars):
        """
        get log(like) + log(prior)
        """
        loglike = self._get_loglike(pars)
        gp = log( self.gprior(pars[2], pars[3]) )
        gp=gp[0]
        logprob = gp + loglike
        return logprob

    def _get_loglike(self, pars):
        # like is exp(0.5*A*B^2) where
        # A is sum((model/err)^2) and is fixed
        # and
        #   B = sum(model*image/err^2)/A
        #     = sum(model/err * image/err)/A

        # build up B
        # this is model/err
        Btmp=self._get_normalized_model(pars)

        # Now multiply by image/err
        Btmp *= self.image
        Btmp *= self.ierr
        # now do the sum and divide by A
        B = Btmp.sum()/self.A

        # A actually fully cancels here when we perform the
        # renormalization to make sum( (model/err)^2 ) == A
        loglike = self.A * B**2/2
        return loglike


    def _calc_result(self):
        """
        if when_prior='after', We apply the shape prior here

        We marginalize over all parameters but g1,g2, which
        are index 0 and 1 in the pars array
        """
        import mcmc
        burn=self.burnin

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g1vals=self.trials['pars'][burn:,2]
        g2vals=self.trials['pars'][burn:,3]

        prior = self.gprior(self.trials['pars'][burn:,2],
                            self.trials['pars'][burn:,3])
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        if self.when_prior=='after':
            # we need to multiply each by the prior
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de
            means, cov = mcmc.extract_stats(self.trials['pars'][burn:,:])
            #mcmc.print_stats(means,cov,names=['cen1','cen2','g1','g2','T'])

            g=means[2:4]
            gcov=cov[2:4,2:4]

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
            gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()

 
        w,=where(self.trials['accepted']==1)
        arate = w.size/float(self.trials.size)
        print 'acceptance rate:',w.size/float(self.trials.size)
        self._result={'g':g,'gcov':gcov,'gsens':gsens,'arate':arate}


    def _get_normalized_model(self, pars):
        """
        Model/err.  Normalized such that sum(model/err)^2 ) = A

        ymod is a temporary
        """
        model = self._get_model(pars)

        ymod = model*self.ierr

        ymodsum = ymod.sum()
        ymod2sum = (ymod**2).sum()

        # this ensures sum( (model/err)^2 ) == A
        norm = sqrt(ymodsum**2*self.A/ymod2sum)

        # also divide again by err since we actually need model/err^2
        ymod *= (norm/ymodsum)

        return ymod



    def _get_gstepsize(self, s2n):
        # these are pretty good for gauss-gauss and just e1,e2
        # going for about 0.3-0.35 acceptance rate for 2 pars
        # for higher dims should tend closer to .25
        s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,   60,   70,   80,   90,  100]
        steps = [.40, .22, .13, .10, .08, .07, .05, .04, .034, .030, .026, .023, .021]

        s2ns=array(s2ns,dtype='f8')
        steps=array(steps,dtype='f8')

        return numpy.interp(s2n, s2ns, steps)

    def _get_Tstepsize(self, s2n):
        s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,  60,  70,  80,  90, 100]
        steps = [2.0, 1.8, 1.5, 1.1, .80, .60, .45, .35, .28, .20, .20, .20, .20]

        s2ns=array(s2ns,dtype='f8')
        steps=array(steps,dtype='f8')

        return numpy.interp(s2n, s2ns, steps)


    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _doplots(self):
        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(3,2)
        # g1
        burn_g1=mcmc.plot_burnin(self.trials['pars'][:,2], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_1$',
                               show=False)
        # g2
        burn_g2=mcmc.plot_burnin(self.trials['pars'][:,3], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_2$',
                               show=False)
        # T
        burn_T=mcmc.plot_burnin(self.trials['pars'][:,4], 
                               self.trials['loglike'],
                               self.burnin, ylabel='T',
                               show=False)


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))
        w,=where(self.trials['accepted']==1)
        print 'acceptance rate:',w.size/float(self.trials.size)
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]
        bsize1=errs[0]*0.2
        bsize2=errs[1]*0.2
        hplt_g1 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 2],binsize=bsize1,
                                  show=False)
        hplt_g2 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 3],binsize=bsize2,
                                  show=False)

        Tsdev = self.trials['pars'][self.burnin:,4].std()
        Tbsize=Tsdev*0.2
        hplt_T = eu.plotting.bhist(self.trials['pars'][self.burnin:, 4],binsize=Tbsize,
                                  show=False)

        hplt_g1.xlabel=r'$g_1$'
        hplt_g2.xlabel=r'$g_2$'
        hplt_T.xlabel='T'
        tab[0,0] = burn_g1
        tab[1,0] = burn_g2
        tab[2,0] = burn_T
        tab[0,1] = hplt_g1
        tab[1,1] = hplt_g2
        tab[2,1] = hplt_T
        tab.show()




class MCMCFitterFixCen:
    """
    Fixed center but T floating

    likelihood is exp(0.5*A*B^2)

    A is 
        sum( (model/err)^2 ) which we fix to sum ((ydata/err)^2)
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that the
    relative height of the likelihood is valid given the other simplifications
    we have made. We choose 1 because, using the re-normalization below,
    the value actually cancels.

    We must normalize the model according to the A condition.  create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S) 

    If using "exp", 3 gaussian mixture model, you can examine run
    gmix-fit-et10r99

        Moving parts

            - Want cen with tight prior
            - Want e1,e2. prior applied *after* exploring likelihood surface
            because the prior can be cuspy and this confuses the MCMC
            - Want a T that represents the size, no prior at all but a limited
            range

        Fixed parts of model

            - Want a delta function; fixed at zero for now: can only detect
            non-zero at *very* high S/N
            - Want a fixed F that is a multiple of T for the second gauss.  For
            s2 ~ 1, this appears to be 

                ~3.8 for s2 ~ 1

            - Want fixed p values. these  are only meaningful relative since we
            marginalize over amplitude. For s2=1

                0.379 0.563, 0.0593
    """

    def __init__(self, 
                 image=None, 
                 pixerr=None, 
                 guess=[0.,0.],
                 gprior=None,
                 psf=None,
                 cen=None,
                 T=None,
                 s2n=None, 
                 nstep=6400, 
                 burnin=1000,
                 when_prior='during',
                 sampler='mcmc'):
        """
        parameters
        ----------
        image:
            image as a numpy array
        pixerr:
            Error per pixel
        guess:
            Starting guess for pars
        psf:
            The psf gaussian mixture
        cen:
            The center.
        T:
            Starting value for ixx+iyy of main component
        nstep:
            Number of steps in MCMC chain.
        burnin:
            Number of burn in steps.
        gprior:
            The prior on the g1,g2 surface.
        s2n:
            S/N value; will be used to create step sizes.
        """

        self.make_plots=False
        self._set_data(image=image, 
                       pixerr=pixerr, 
                       guess=guess,
                       gprior=gprior,
                       psf=psf,
                       cen=cen,
                       T=T,
                       s2n=s2n, 
                       nstep=nstep, 
                       burnin=burnin,
                       sampler=sampler,
                       when_prior=when_prior)

        self._go()

    def get_result(self):
        return self._result

    def _go(self):
        if self.sampler=='mcmc':
            self._go_mcmc()
        else:
            self._go_emcee()

    def _go_mcmc(self):
        import mcmc
        mcmc=mcmc.MCMC(self)

        guess=zeros(self.npars)
        # guess for g1,g2 is (0,0 with some scatter
        guess[0:2]=0.1*(randu(2)-0.5)
        # guess for T is self.T with scatter
        guess[2] = self.T + self.T*0.1*(randu()-0.5)
        
        self.trials = mcmc.run(self.nstep, guess)

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if self.make_plots:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop

    def _go_emcee(self):
        import emcee
        self.nwalkers=10
        guess=zeros( (self.nwalkers,self.npars) )


        # guess for g1,g2 is (0,0) with some scatter
        guess[:,0]=0.1*(randu(self.nwalkers)-0.5)
        guess[:,1]=0.1*(randu(self.nwalkers)-0.5)

        # guess for T is self.T with scatter
        guess[:,2] = self.T + self.T*0.1*(randu(self.nwalkers)-0.5)

        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self.loglike,
                                        a=4.)
        
        burn_per=self.burnin/self.nwalkers
        nstep_per=(self.nstep-self.burnin)/self.nwalkers

        pos, prob, state = sampler.run_mcmc(guess, burn_per)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, nstep_per)
        
        # dummy
        self.burnin=10

        dtype = [('pars',('f8',self.npars)), ('loglike','f8')]
        trials = zeros(self.nwalkers*nstep_per, dtype=dtype)

        lnprobs = sampler.lnprobability.reshape(self.nwalkers*nstep_per)
        lnprobs = lnprobs - lnprobs.max()

        trials['pars'] = sampler.flatchain
        trials['loglike'] = lnprobs

        self.trials=trials

        self._emcee_sampler=sampler

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if self.make_plots:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop


    def loglike(self, pars):
        """
        pars are [g1,g2,T]
        Packaged up for the MCMC code
        """
        # hard priors first
        g1,g2=pars[0],pars[1]
        T=pars[2]

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0]=self.cen[0]
        self.tpars[1]=self.cen[1]
        self.tpars[2]=e1
        self.tpars[3]=e2
        self.tpars[4]=T
        self.tpars[5]=1.

        logprob = self._get_loglike_c(self.tpars)
        if self.when_prior=='during':
            gp = self._get_lngprior(g1,g2)
            logprob += gp

        return logprob


    def loglike_old(self, pars):
        """
        Packaged up for the MCMC code
        """
        # hard priors first
        try:
            shape=lensing.Shear(g1=pars[0], g2=pars[1])
        except:
            return LOWVAL

        if pars[2] < 0:
            return LOWVAL

        if self.when_prior=='after':
            return self._get_loglike(shape, pars[2], self.cen)
        else:
            return self._get_logprob(shape, pars[2], self.cen)


    def step(self, pars):
        """
        public for the MCMC code
        """
        #nextpars=zeros(self.npars)
        """
        f=random.random()
        if f > 0.5:
            nextpars = pars + self.small_stepsize*randn(self.npars)
        else:
            nextpars = pars + self.stepsize*randn(self.npars)
        """
        nextpars = pars + self.stepsize*randn(self.npars)
        return nextpars

    def _get_logprob(self, shape, T, cen):
        """
        get log(like) + log(prior)
        """
        loglike = self._get_loglike(shape, T, cen)
        gp = self._get_lngprior(shape.g1, shape.g2)
        logprob = gp + loglike
        return logprob


    def _get_loglike_c(self, pars):
        #t0=time.time()

        loglike,flags=\
            gmix_image.render._render.loglike_coellip(self.image, 
                                                      pars, 
                                                      self.psf_pars, 
                                                      self.Anorm,
                                                      self.ierr)
        if flags != 0:
            return LOWVAL
        #self.logl_tm += time.time()-t0
        return loglike

    def _get_lngprior(self, g1, g2):
        g=sqrt(g1**2 + g2**2)
        gp = self.gprior.prior_gabs_scalar(g)
        if gp > 0:
            gp = log(gp)
        else:
            gp=LOWVAL
        return gp

    def _get_loglike(self, shape, T, cen):
        # like is exp(0.5*A*B^2) where
        # A is sum((model/err)^2) and is fixed
        # and
        #   B = sum(model*image/err^2)/A
        #     = sum(model/err * image/err)/A

        # build up B
        # this is model/err
        Btmp=self._get_normalized_model(shape, T, cen)

        # Now multiply by image/err
        Btmp *= self.image
        Btmp *= 1./self.pixerr
        # now do the sum and divide by A
        B = Btmp.sum()/self.Anorm

        # A actually fully cancels here when we perform the
        # renormalization to make sum( (model/err)^2 ) == A
        loglike = self.Anorm * B**2/2

        return loglike



    def _get_normalized_model(self, shape, T, cen):
        """
        Model/err.  Normalized such that sum(model/err)^2 ) = A
        """
        model = self._get_model(shape, T, cen)

        ymod = model/self.pixerr

        ymodsum = ymod.sum()
        ymod2sum = (ymod**2).sum()

        # this ensures sum( (model/err)^2 ) == A
        norm = sqrt(ymodsum**2*self.Anorm/ymod2sum)

        # also divide again by err since we actually need model/err^2
        ymod *= (norm/ymodsum)

        return ymod


    def _calc_result(self):
        """
        if when_prior='after', We apply the shape prior here

        We marginalize over all parameters but g1,g2, which
        are index 0 and 1 in the pars array
        """
        import mcmc
        burn=self.burnin

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g1vals=self.trials['pars'][burn:,0]
        g2vals=self.trials['pars'][burn:,1]

        prior = self.gprior(self.trials['pars'][burn:,0],
                            self.trials['pars'][burn:,1])
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        if self.when_prior=='after':
            # we need to multiply each by the prior
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de
            g, gcov = mcmc.extract_stats(self.trials['pars'][burn:,0:2])

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
            gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()

 
        if self.sampler=='emcee':
            arates = self._emcee_sampler.acceptance_fraction
            arate = arates.mean()
        else:
            w,=where(self.trials['accepted']==1)
            arate = w.size/float(self.trials.size)
        #print 'acceptance rate:',w.size/float(self.trials.size)
        self._result={'g':g,'gcov':gcov,'gsens':gsens,'arate':arate}
        #wlog("arate:",self._result['arate'])



    def _get_gstepsize(self, s2n):
        # these are pretty good for gauss-gauss and just e1,e2
        # going for about 0.3-0.35 acceptance rate for 2 pars
        # for higher dims should tend closer to .25
        if self.when_prior=='during':
            s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,   60,   70,   80,   90,  100]
            steps = [.20, .14, .10, .08, .06, .05, .04, .034, .030, .026, .021, .017, .015]
        else:
            s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,   60,   70,   80,   90,  100]
            steps = [.40, .22, .13, .10, .08, .07, .05, .04, .034, .030, .026, .023, .021]

        s2ns=array(s2ns,dtype='f8')
        steps=array(steps,dtype='f8')

        return numpy.interp(s2n, s2ns, steps)

    def _get_Tstepsize(self, s2n):
        if self.when_prior=='during':
            s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,  60,  70,  80,  90, 100]
            steps = [1.0, 1.2, 1.0, .80, .50, .45, .35, .28, .20, .18, .16, .15, .13]
        else:
            s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,  60,  70,  80,  90, 100]
            steps = [2.0, 1.8, 1.5, 1.1, .80, .60, .45, .35, .28, .20, .20, .20, .20]

        s2ns=array(s2ns,dtype='f8')
        steps=array(steps,dtype='f8')

        return numpy.interp(s2n, s2ns, steps)


    def _set_data(self, 
                  image=None, 
                  pixerr=None, 
                  guess=None,
                  gprior=None,
                  psf=None,
                  cen=None,
                  T=None,
                  s2n=None, 
                  nstep=6400,
                  burnin=1000,
                  sampler='mcmc',
                  when_prior='during'):

        # e1,e2,T
        self.npars=3

        self.image=image
        self.pixerr=pixerr
        self.ierr=float(1./pixerr)
        self.T=T
        self.cen=cen
        self.nstep=nstep
        self.burnin=burnin
        self.s2n=s2n
        self.gprior=gprior
        self.when_prior=when_prior
        self.sampler=sampler

        gstepsize=self._get_gstepsize(s2n)
        Tstepsize=self._get_Tstepsize(s2n)
        self.stepsize=array([gstepsize,gstepsize,Tstepsize], dtype='f8')

        self._set_psf(psf)

        if gprior is None or s2n is None:
            raise ValueError("BayesFitter(gprior=, nstep=6400, burnin=1000, s2n=)")
        if image is None or pixerr is None:
            raise ValueError("send at least image and pixerr")
        if T is None or cen is None:
            raise ValueError("for now send T and cen to fix them")
        #if self.guess.size != self.npars:
        #    raise ValueError("expected guess with length %d" % self.npars)

        self.Anorm = float(1)

        # This will make sure make sure the guess is valid
        #self.guess_shape = lensing.Shear(guess[0],guess[1])

        self.tpars=zeros(6,dtype='f8')

        self._model = zeros( image.shape )

    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _doplots(self):
        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(3,2)
        burn1=mcmc.plot_burnin(self.trials['pars'][:,0], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_1$',
                               show=False)
        burn2=mcmc.plot_burnin(self.trials['pars'][:,1], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_2$',
                               show=False)
        burn3=mcmc.plot_burnin(self.trials['pars'][:,2], 
                               self.trials['loglike'],
                               self.burnin, ylabel='T',
                               show=False)


        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))
        print 'acceptance rate:',self._result['arate']
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]
        bsize1=errs[0]*0.2
        bsize2=errs[1]*0.2
        hplt1 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 0],binsize=bsize1,
                                  show=False)
        hplt2 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 1],binsize=bsize2,
                                  show=False)

        Tsdev = self.trials['pars'][self.burnin:,2].std()
        Tbsize=Tsdev*0.2
        hplt3 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 2],binsize=Tbsize,
                                  show=False)

        hplt1.xlabel=r'$g_1$'
        hplt2.xlabel=r'$g_2$'
        hplt3.xlabel='T'
        tab[0,0] = burn1
        tab[1,0] = burn2
        tab[2,0] = burn3
        tab[0,1] = hplt1
        tab[1,1] = hplt2
        tab[2,1] = hplt3
        tab.show()



class MCMCFitterFixTCen:
    """
    Fixing T and cen

    likelihood is exp(0.5*A*B^2)

    A is 
        sum( (model/err)^2 ) which we fix to sum ((ydata/err)^2)
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that the
    relative height of the likelihood is valid given the other simplifications
    we have made. We choose 1 because, using the re-normalization below,
    the value actually cancels.

    We must normalize the model according to the A condition.  create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S) 

    If using "exp", 3 gaussian mixture model, you can examine run
    gmix-fit-et10r99

        Moving parts

            - Want cen with tight prior
            - Want e1,e2. prior applied *after* exploring likelihood surface
            because the prior can be cuspy and this confuses the MCMC
            - Want a T that represents the size, no prior at all but a limited
            range

        Fixed parts of model

            - Want a delta function; fixed at zero for now: can only detect
            non-zero at *very* high S/N
            - Want a fixed F that is a multiple of T for the second gauss.  For
            s2 ~ 1, this appears to be 

                ~3.8 for s2 ~ 1

            - Want fixed p values. these  are only meaningful relative since we
            marginalize over amplitude. For s2=1

                0.379 0.563, 0.0593
    """

    def __init__(self, 
                 image=None, 
                 pixerr=None, 
                 guess=[0.,0.],
                 gprior=None,
                 psf=None,
                 cen=None, # send to fix cen
                 T=None,  # send to fix T
                 s2n=None, 
                 nstep=6400, 
                 when_prior='after',
                 burnin=1000):
        """
        parameters
        ----------
        image:
            image as a numpy array
        pixerr:
            Error per pixel
        guess:
            Starting guess for pars
        psf:
            The psf gaussian mixture
        cen:
            The center, to be held fixed
        T:
            ixx+iyy to be held fixed.
        nstep:
            Number of steps in MCMC chain.
        burnin:
            Number of burn in steps.
        gprior:
            The prior on the g1,g2 surface.
        s2n:
            S/N value; will be used to create step sizes.
            Currently step size is fixed actually
        """

        self._set_data(image=image, 
                       pixerr=pixerr, 
                       guess=guess,
                       gprior=gprior,
                       psf=psf,
                       cen=cen, # send to fix cen
                       T=T,  # send to fix T
                       s2n=s2n, 
                       nstep=nstep, 
                       burnin=burnin)

        self._go()

    def get_result(self):
        return self._result

    def _go(self):
        import mcmc
        mcmc=mcmc.MCMC(self)

        # need to improve this when we add other parameters
        guess=0.1*random.random(self.npars)
        
        self.trials = mcmc.run(self.nstep, guess)

        # get the expectation values, sensitivity and errors
        self._calc_result()

        if False:
            self._doplots()
            key=raw_input('hit a key (q to quit): ')
            if key=='q':
                stop


    def loglike(self, pars):
        """
        Packaged up for the MCMC code
        """
        try:
            shape=lensing.Shear(g1=pars[0], g2=pars[1])
        except:
            return LOWVAL

        if self.when_prior=='after':
            return self._get_loglike(shape, self.T, self.cen)
        else:
            return self._get_logprob(shape, self.T, self.cen)

    def step(self, pars):
        """
        public for the MCMC code
        """
        nextpars=zeros(self.npars)
        """
        f=random.random()
        if f > 0.5:
            nextpars = pars + self.small_stepsize*randn(self.npars)
        else:
            nextpars = pars + self.stepsize*randn(self.npars)
        """
        nextpars = pars + self.stepsize*randn(self.npars)
        return nextpars

    def _get_logprob(self, shape, T, cen):
        """
        get log(like) + log(prior)
        """
        loglike = self._get_loglike(shape, T, cen)
        gp = log( self.gprior(pars[2], pars[3]) )
        gp=gp[0]
        logprob = gp + loglike
        return logprob

    def _get_loglike(self, shape, T, cen):
        # like is exp(0.5*A*B^2) where
        # A is sum((model/err)^2) and is fixed
        # and
        #   B = sum(model*image/err^2)/A
        #     = sum(model/err * image/err)/A

        # build up B
        # this is model/err
        Btmp=self._get_normalized_model(shape, T, cen)

        # Now multiply by image/err
        Btmp *= self.image
        Btmp *= 1./self.pixerr
        # now do the sum and divide by A
        B = Btmp.sum()/self.A

        # A actually fully cancels here when we perform the
        # renormalization to make sum( (model/err)^2 ) == A
        loglike = self.A * B**2/2
        return loglike


    def _calc_result(self):
        """
        if when_prior='after', We apply the shape prior here
        """
        burn=self.burnin

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g1vals=self.trials['pars'][burn:,0]
        g2vals=self.trials['pars'][burn:,1]

        prior = self.gprior(self.trials['pars'][burn:,0],
                            self.trials['pars'][burn:,1])
        dpri_by_g1 = self.gprior.dbyg1(g1vals,g2vals)
        dpri_by_g2 = self.gprior.dbyg2(g1vals,g2vals)

        psum = prior.sum()

        if self.when_prior=='after':
            # we need to multiply each by the prior
            g[0] = (g1vals*prior).sum()/psum
            g[1] = (g2vals*prior).sum()/psum

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            gcov[0,0] = (g1diff**2*prior).sum()/psum
            gcov[0,1] = (g1diff*g2diff*prior).sum()/psum
            gcov[1,0] = gcov[0,1]
            gcov[1,1] = (g2diff**2*prior).sum()/psum

            # now the sensitivity is 
            #  sum( (<g>-g) L*dP/dg )
            #  ----------------------
            #        sum(L*P)
            #
            # the likelihood is already in the points

            gsens[0] = 1. - (g1diff*dpri_by_g1).sum()/psum
            gsens[1] = 1. - (g2diff*dpri_by_g2).sum()/psum
        else:
            # prior is already in the distribution of
            # points.  This is simpler for most things but
            # for sensitivity we need a factor of (1/P)dP/de
            g, gcov = mcmc.extract_stats(self.trials['pars'][burn:,:])

            g1diff = g[0]-g1vals
            g2diff = g[1]-g2vals

            w,=where(prior > 0)
            if w.size == 0:
                raise ValueError("no prior values > 0!")

            gsens[0]= 1.-(g1diff[w]*dpri_by_g1[w]/prior[w]).mean()
            gsens[1]= 1.-(g2diff[w]*dpri_by_g2[w]/prior[w]).mean()

 
        w,=where(self.trials['accepted']==1)
        arate = w.size/float(self.trials.size)
        #print 'acceptance rate:',w.size/float(self.trials.size)
        self._result={'g':g,'gcov':gcov,'gsens':gsens,'arate':arate}


    def _get_normalized_model(self, shape, T, cen):
        """
        Model/err.  Normalized such that sum(model/err)^2 ) = A
        """
        model = self._get_model(shape, T, cen)

        ymod = model/self.pixerr

        ymodsum = ymod.sum()
        ymod2sum = (ymod**2).sum()

        # this ensures sum( (model/err)^2 ) == A
        norm = sqrt(ymodsum**2*self.A/ymod2sum)

        # also divide again by err since we actually need model/err^2
        ymod *= (norm/ymodsum)


        return ymod


    def _get_stepsize(self, s2n):
        # these are pretty good for gauss-gauss and just e1,e2
        # going for about 0.3-0.35 acceptance rate for 2 pars
        # for higher dims should tend closer to .25
        s2ns  = [  5,  10,  15,  20,  25,  30,  40,  50,   60,   70,   80,   90,  100]
        steps = [.40, .22, .13, .10, .08, .07, .05, .04, .034, .030, .026, .023, .021]

        s2ns=array(s2ns,dtype='f8')
        steps=array(steps,dtype='f8')

        return numpy.interp(s2n, s2ns, steps)

    def _set_data(self, 
                  image=None, 
                  pixerr=None, 
                  guess=None,
                  gprior=None,
                  psf=None,
                  cen=None,
                  T=None,
                  s2n=None, 
                  nstep=6400,
                  when_prior='after',
                  burnin=1000):

        self.npars=2

        self.image=image
        self.pixerr=pixerr
        #self.guess=0.1*random.random(self.npars)
        self.T=T
        self.cen=cen
        self.nstep=nstep
        self.burnin=burnin
        self.s2n=s2n
        self.gprior=gprior

        stepsize=self._get_stepsize(s2n)
        #print 'stepsize:',stepsize

        self.stepsize=array([stepsize,stepsize], dtype='f8')
        #self.stepsize=array([0.1,0.1], dtype='f8')
        #self.small_stepsize=array([0.01,0.01], dtype='f8')
        #self.stepsize=array([0.8,0.8], dtype='f8')
        #self.small_stepsize=array([0.01,0.01], dtype='f8')

        self._set_psf(psf)

        if gprior is None or s2n is None:
            raise ValueError("BayesFitter(gprior=, nstep=6400, burnin=1000, s2n=)")
        if image is None or pixerr is None:
            raise ValueError("send at least image and pixerr")
        if T is None or cen is None:
            raise ValueError("for now send T and cen to fix them")
        #if self.guess.size != self.npars:
        #    raise ValueError("expected guess with length %d" % self.npars)

        self.A = 1

        # This will make sure make sure the guess is valid
        #self.guess_shape = lensing.Shear(guess[0],guess[1])


        self.when_prior=when_prior

    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _doplots(self):
        import mcmc
        import biggles
        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(2,2)
        burn1=mcmc.plot_burnin(self.trials['pars'][:,0], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_1$',
                               show=False)
        burn2=mcmc.plot_burnin(self.trials['pars'][:,1], 
                               self.trials['loglike'],
                               self.burnin, ylabel=r'$g_2$',
                               show=False)

        g = self._result['g']
        gcov = self._result['gcov']
        errs = sqrt(diag(gcov))
        w,=where(self.trials['accepted']==1)
        print 'acceptance rate:',w.size/float(self.trials.size)
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'g1sens:',self._result['gsens'][0]
        print 'g2sens:',self._result['gsens'][1]
        bsize1=errs[0]*0.2
        bsize2=errs[1]*0.2
        hplt1 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 0],binsize=bsize1,
                                  show=False)
        hplt2 = eu.plotting.bhist(self.trials['pars'][self.burnin:, 1],binsize=bsize2,
                                  show=False)

        tab[0,0] = burn1
        tab[1,0] = burn2
        tab[0,1] = hplt1
        tab[1,1] = hplt2
        tab.show()


class BayesFitterFixCen:
    """
    likelihood is exp(0.5*A*B^2)

    This analytically marginalizes over the amplitide.  We
    actually check a grid over the parameters.  There is also
    the MCMC Fitter

    A is 
        sum( (model/err)^2 ) which we fix to sum ((ydata/err)^2)
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that the
    relative height of the likelihood is valid given the other simplifications
    we have made. We arbitrarily choose A=1 because, using the re-normalization
    below, the value actually cancels.

    To normalize the model according to the A condition, create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S)
    """
    def __init__(self, prior, n_ggrid, gmin, gmax, n_Tgrid, Tmin, Tmax):
        self.prior=prior

        # range for search
        self.n_ggrid=n_ggrid # in both g1 and g2
        self.gmin=gmin
        self.gmax=gmax

        self.n_Tgrid=n_Tgrid
        self.Tmin=Tmin
        self.Tmax=Tmax

        self.g1vals = linspace(self.gmin, self.gmax, self.n_ggrid)
        self.g2vals = self.g1vals.copy()

        self.Tvals = linspace(self.Tmin, self.Tmax, self.n_Tgrid)

        self.g1matrix = self.g1vals.reshape(self.n_ggrid,1)*ones(self.n_ggrid)
        self.g2matrix = self.g2vals*ones(self.n_ggrid).reshape(self.n_ggrid,1)

        self._set_prior_matrices()

        self.tpars = zeros(6, dtype='f8')

        self.Anorm = float(1)

        self.models={}

    def get_like(self):
        return self._like
    def get_result(self):
        return self._result

    def process_image(self, 
                      image, 
                      pixerr, 
                      cen,
                      psf):

        self.image=image
        self.pixerr=float(pixerr)
        self.ierr = float(1./pixerr)
        self.cen=cen

        self._set_psf(psf)


        self._go()


    def _go(self):
        cen   = self.cen
        dims = self.image.shape

        loglike=LOWVAL + zeros((self.n_Tgrid,self.n_ggrid, self.n_ggrid))

        for iT,T in enumerate(self.Tvals):
            for i1,g1 in enumerate(self.g1vals):
                for i2,g2 in enumerate(self.g2vals):
                    loglike[T,i1,i2] = self._calc_loglike(cen, g1, g2, T)

        loglike -= loglike.max()
        liketot = exp(loglike)

        # marginalize over T
        like = liketot.sum(0)

        self._like=like
        #images.multiview(like)
        #key=raw_input('q to quit: ')
        #if key == 'q':
        #    stop

        # get the expectation values, sensitiviey and errors
        self._calc_result()

    def _calc_loglike(self, cen, g1, g2, T):
        """
        pars are [g1,g2,T]
        """

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0]=cen[0]
        self.tpars[1]=cen[1]
        self.tpars[2]=e1
        self.tpars[3]=e2
        self.tpars[4]=T
        self.tpars[5]=1.

        loglike,flags=\
            gmix_image.render._render.loglike_coellip(self.image, 
                                                      self.tpars, 
                                                      self.psf_pars, 
                                                      self.Anorm,
                                                      self.ierr)
        if flags != 0:
            return LOWVAL

        return loglike 



    def _calc_result(self):
        lp = self._like*self.prior_matrix
        lpsum = lp.sum()

        g1 = (lp*self.g1matrix).sum()/lpsum
        g2 = (lp*self.g2matrix).sum()/lpsum

        # order matters for sensitivity below
        g1diff = g1-self.g1matrix
        g2diff = g2-self.g2matrix
        g11var = (lp*g1diff**2).sum()/lpsum
        g22var = (lp*g2diff**2).sum()/lpsum
        g12var = (lp*g1diff*g2diff).sum()/lpsum

        
        ldp1 = self._like*self.prior_d1_matrix
        ldp2 = self._like*self.prior_d2_matrix
        
        g1sens = 1.-(g1diff*ldp1).sum()/lpsum
        g2sens = 1.-(g2diff*ldp2).sum()/lpsum

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g[0] = g1
        g[1] = g2
        gsens[0] = g1sens
        gsens[1] = g2sens
        gcov[0,0] = g11var
        gcov[0,1] = g12var
        gcov[1,0] = g12var
        gcov[1,1] = g22var


        self._result={'g':g,'gcov':gcov,'gsens':gsens}


    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _set_prior_matrices(self):
        self.prior_matrix = zeros( (self.n_ggrid,self.n_ggrid) )
        self.prior_d1_matrix = zeros( (self.n_ggrid,self.n_ggrid) )
        self.prior_d2_matrix = zeros( (self.n_ggrid,self.n_ggrid) )

        for i1,g1 in enumerate(self.g1vals):
            for i2,g2 in enumerate(self.g2vals):
                self.prior_matrix[i1,i2] = self.prior(g1,g2)
                self.prior_d1_matrix[i1,i2] = self.prior.dbyg1(g1,g2)
                self.prior_d2_matrix[i1,i2] = self.prior.dbyg2(g1,g2)



class BayesFitterFixTCen:
    """
    likelihood is exp(0.5*A*B^2)

    This analytically marginalizes over the amplitide.  We
    actually check a grid over the parameters.  There is also
    the MCMC Fitter

    A is 
        sum( (model/err)^2 ) which we fix to sum ((ydata/err)^2)
    B is
        sum( model*data/err^2 )/A

    We must fix A as a constant for every model we generate, so that the
    relative height of the likelihood is valid given the other simplifications
    we have made. We arbitrarily choose A=1 because, using the re-normalization
    below, the value actually cancels.

    To normalize the model according to the A condition, create
        ymod = model/err
    which will have a normalization S. Then the renormalization is given by
        N = sqrt( S*A/sum(ymod^2) )
    (i.e. set ymod = ymod*N/S)
    """
    def __init__(self, prior, n_ggrid, gmin, gmax):
        self.prior=prior

        # range for search
        self.n_ggrid=n_ggrid # in both g1 and g2
        self.gmin=gmin
        self.gmax=gmax

        self.g1vals = linspace(self.gmin, self.gmax, self.n_ggrid)
        self.g2vals = self.g1vals.copy()

        self.g1matrix = self.g1vals.reshape(self.n_ggrid,1)*ones(self.n_ggrid)
        self.g2matrix = self.g2vals*ones(self.n_ggrid).reshape(self.n_ggrid,1)

        self._set_prior_matrices()

        self.tpars = zeros(6, dtype='f8')
        self.Anorm = float(1)

        self.models={}

    def get_like(self):
        return self._like
    def get_result(self):
        return self._result

    def process_image(self, 
                      image, 
                      pixerr, 
                      cen,
                      T,
                      psf):

        self.image=image
        self.pixerr=float(pixerr)
        self.ierr = float(1./pixerr)
        self.T=T
        self.cen=cen

        self._set_psf(psf)


        self._go()


    def _go(self):
        T     = self.T
        cen   = self.cen
        dims = self.image.shape

        loglike=LOWVAL + zeros((self.n_ggrid, self.n_ggrid))

        for i1,g1 in enumerate(self.g1vals):
            for i2,g2 in enumerate(self.g2vals):
                loglike[i1,i2] = self._calc_loglike(cen, g1, g2, T)

        loglike -= loglike.max()
        like = exp(loglike)

        self._like=like
        #images.multiview(like)
        #key=raw_input('q to quit: ')
        #if key == 'q':
        #    stop

        # get the expectation values, sensitiviey and errors
        self._calc_result()

    def _calc_loglike(self, cen, g1, g2, T):
        """
        pars are [g1,g2,T]
        """

        e1,e2,ok = g1g2_to_e1e2(g1,g2)
        if not ok:
            return LOWVAL

        if T < 0:
            return LOWVAL

        self.tpars[0]=cen[0]
        self.tpars[1]=cen[1]
        self.tpars[2]=e1
        self.tpars[3]=e2
        self.tpars[4]=T
        self.tpars[5]=1.

        loglike,flags=\
            gmix_image.render._render.loglike_coellip(self.image, 
                                                      self.tpars, 
                                                      self.psf_pars, 
                                                      self.Anorm,
                                                      self.ierr)
        if flags != 0:
            return LOWVAL

        return loglike 



    def _calc_result(self):
        lp = self._like*self.prior_matrix
        lpsum = lp.sum()

        g1 = (lp*self.g1matrix).sum()/lpsum
        g2 = (lp*self.g2matrix).sum()/lpsum

        # order matters for sensitivity below
        g1diff = g1-self.g1matrix
        g2diff = g2-self.g2matrix
        g11var = (lp*g1diff**2).sum()/lpsum
        g22var = (lp*g2diff**2).sum()/lpsum
        g12var = (lp*g1diff*g2diff).sum()/lpsum

        
        ldp1 = self._like*self.prior_d1_matrix
        ldp2 = self._like*self.prior_d2_matrix
        
        g1sens = 1.-(g1diff*ldp1).sum()/lpsum
        g2sens = 1.-(g2diff*ldp2).sum()/lpsum

        g=zeros(2)
        gcov=zeros((2,2))
        gsens = zeros(2)

        g[0] = g1
        g[1] = g2
        gsens[0] = g1sens
        gsens[1] = g2sens
        gcov[0,0] = g11var
        gcov[0,1] = g12var
        gcov[1,0] = g12var
        gcov[1,1] = g22var


        self._result={'g':g,'gcov':gcov,'gsens':gsens}


    def _set_psf(self, psf):
        self.psf_gmix = psf
        self.psf_pars = None

        if psf is not None:
            if not isinstance(psf[0],dict):
                raise ValueError("psf must be list of dicts")
            self.psf_pars = gmix_image.gmix2pars(psf)

    def _set_prior_matrices(self):
        self.prior_matrix = zeros( (self.n_ggrid,self.n_ggrid) )
        self.prior_d1_matrix = zeros( (self.n_ggrid,self.n_ggrid) )
        self.prior_d2_matrix = zeros( (self.n_ggrid,self.n_ggrid) )

        for i1,g1 in enumerate(self.g1vals):
            for i2,g2 in enumerate(self.g2vals):
                self.prior_matrix[i1,i2] = self.prior(g1,g2)
                self.prior_d1_matrix[i1,i2] = self.prior.dbyg1(g1,g2)
                self.prior_d2_matrix[i1,i2] = self.prior.dbyg2(g1,g2)


    def _go_old(self):
        T     = self.T
        cen   = self.cen
        dims = self.image.shape

        loglike=-9999.0e9 + zeros((self.n_ggrid, self.n_ggrid))

        for i1,g1 in enumerate(self.g1vals):
            for i2,g2 in enumerate(self.g2vals):

                try:
                    # will raise exception if g > 1 or e > 1
                    sh=self.get_shape(i1,i2)
                except lensing.ShapeRangeError:
                    continue

                # like is exp(0.5*A*B^2) where
                # A is sum((model/err)^2) and is fixed
                # and
                #   B = sum(model*image/err^2)/A
                #     = sum(model/err * image/err)/A

                # build up B
                # this is model/err
                Btmp=self._get_normalized_model(i1, i2, T, dims, cen)

                # Now multiply by image/err
                Btmp *= self.image
                Btmp *= 1./self.pixerr
                # now do the sum and divide by A
                B = Btmp.sum()/self.A

                # A actually fully cancels here when we perform the
                # renormalization to make sum( (model/err)^2 ) == A
                arg = self.A * B**2/2

                loglike1 = arg
                loglike2 = self._calc_loglike(cen, g1, g2, T)
                wlog(loglike1-loglike2)
                loglike[i1,i2] = arg

                #B = (mod_over_err2*self.image).sum()/self.A
                
                # A actually fully cancels here when we perform the
                # renormalization to make sum( (model/err)^2 ) == A
                #arg = self.A*B**2/2
                #loglike[i1,i2] = arg

        loglike -= loglike.max()
        like = exp(loglike)

        self._like=like
        #images.multiview(like)
        #key=raw_input('q to quit')
        #if key == 'q':
        #    stop

        # get the expectation values, sensitiviey and errors
        self._calc_result()



class CenPrior:
    def __init__(self, cen, sigma):
        self.cen=cen
        self.sigma=sigma
        self.sigma2=[s**2 for s in sigma]

    def lnprob(self, pos):
        lnprob0 = -(self.cen[0]-pos[0])**2/self.sigma2[0]
        lnprob1 = -(self.cen[1]-pos[1])**2/self.sigma2[1]
        return lnprob0 + lnprob1

class GPrior:
    """
    This is in g1,g2 space

    2D
    Prob = A cos(|g| pi/2) exp( - [ 2 |g| / B / (1 + |g|^D) ]^C )
    d/dE(  A cos(sqrt(E^2+q^2) pi/2) exp( - ( 2 sqrt(E^2 + q^2) / B / (1 + sqrt(E^2 + q^2)^D) )^C ) )

    For 1D prob, you need to multiply by 2*pi*|g|
    """
    def __init__(self, A=12.25, B=0.03, C=0.45, D=13.):
        # A actually depends on norm when doing full thing
        self.A = A
        self.B = B
        self.C = C
        self.D = D

        self.maxval = self(0., 0.)

    def __call__(self, g1, g2):
        """
        Get the prior for the input g1,g2 value(s)

        This is the 2 dimensional prior.  If you want to just
        generate the |g| values, use prior1d
        """
        g = sqrt(g1**2 + g2**2)
        return self.prior_gabs(g)

    def dbyg1(self, g1, g2, h=1.e-6):
        """
        Derivative with respect to g1 at the input g1,g2 location

        Uses central difference and a small enough step size
        to use just two points
        """
        ff = self(g1+h/2, g2)
        fb = self(g1-h/2, g2)

        return (ff - fb)/h

    def dbyg2(self, g1, g2, h=1.e-6):
        """
        Derivative with respect to g2 at the input g1,g2 location

        Uses central difference and a small enough step size
        to use just two points
        """
        ff = self(g1, g2+h/2)
        fb = self(g1, g2-h/2)
        return (ff - fb)/h


    def prior_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        g = array(g, ndmin=1, copy=False)
        prior = zeros(g.size)

        w,=where(g < 1)
        if w.size > 0:
            prior[w] = self.A * cos(g[w]*pi/2)*exp( - ( 2*g[w] / self.B / (1 + g[w]**self.D) )**self.C )
        return prior

    def prior_gabs_scalar(self, g):
        return self.A * math.cos(g*pi/2)*math.exp( - ( 2*g / self.B / (1 + g**self.D) )**self.C )

    def sample(self, nrand, as_shear=False):
        """
        Get random g1,g2 values

        parameters
        ----------
        nrand: int
            Number to generate
        as_shear: bool, optional
            If True, get a list of Shear objects
        """
        g1 = zeros(nrand)
        g2 = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand2 = random.random(nleft)
            grand = sqrt(grand2)
            # now uniform angles
            rangle = random.random(nleft)*2*pi

            # now get cartesion locations in g1,g2 plane
            g1rand = grand*cos(rangle)
            g2rand = grand*sin(rangle)

            # now finally the height from [0,maxval)
            h = self.maxval*random.random(nleft)

            pvals = self(g1rand, g2rand)

            w,=where(h < pvals)
            if w.size > 0:
                g1[ngood:ngood+w.size] = g1rand[w]
                g2[ngood:ngood+w.size] = g2rand[w]
                ngood += w.size
                nleft -= w.size

        if as_shear:
            from lensing.shear import Shear
            shlist=[]
            for g1i,g2i in zip(g1,g2):
                shlist.append(Shear(g1=g1i,g2=g2i))
            return shlist
        else:
            return g1, g2


    def prior1d(self, g):
        """
        Get the 1d prior for an input |g| value(s).

        To generate 2-d g1,g2, use prior()
        """
        return 2*pi*g*self.prior_gabs(g)

    def sample1d(self, nrand):
        """
        Get random |g| from the 1d distribution

        parameters
        ----------
        nrand: int
            Number to generate
        """

        if not hasattr(self,'maxval1d'):
            self.set_maxval1d()

        g = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand = random.random(nleft)

            # now finally the height from [0,maxval)
            h = self.maxval1d*random.random(nleft)

            pvals = self.prior1d(grand)

            w,=where(h < pvals)
            if w.size > 0:
                g[ngood:ngood+w.size] = grand[w]
                ngood += w.size
                nleft -= w.size
   
        return g

    def set_maxval1d(self):
        import scipy.optimize
        
        (minvalx, fval, iterations, fcalls, warnflag) \
                = scipy.optimize.fmin(self.prior1dneg, 0.1, full_output=True, 
                                      disp=False)
        if warnflag != 0:
            raise ValueError("failed to find min: warnflag %d" % warnflag)
        self.maxval1d = -fval

    def prior1dneg(self, g, *args):
        """
        So we can use the minimizer
        """
        return -self.prior1d(g)

def g1g2_to_e1e2(g1, g2):
    """
    This version without exceptions

    returns e1,e2,okflag
    """
    g = sqrt(g1**2 + g2**2)
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

def randomize_e1e2(e1start,e2start):
    if e1start == 0 and e2start==0:
        e1rand = 0.05*(randu()-0.5)
        e2rand = 0.05*(randu()-0.5)
    else:
        nmax=100
        ii=0
        while True:
            e1rand = e1start*(1 + 0.2*(randu()-0.5))
            e2rand = e2start*(1 + 0.2*(randu()-0.5))
            etot = sqrt(e1rand**2 + e2rand**2)
            if etot < 0.95:
                break
            ii += 1
            if ii==nmax:
                wlog("---- hit max try on randomize e1e2, setting zero and restart")
                return randomize_e1e2(0.0,0.0)

    return e1rand, e2rand


