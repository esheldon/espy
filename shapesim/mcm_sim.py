import numpy
from numpy import sqrt, diag, log10, exp
import gmix_image
from gmix_image.util import get_estyle_pars, calculate_some_stats, print_pars
import lensing
from lensing.shear import Shear
from lensing.util import ShapeRangeError
from .shapesim import read_config, ShapeSim
from . import shapesim

import copy

from esutil.random import srandu

LOWVAL=-9999.9e9

class MCMSim(dict):
    def __init__(self, run):
        conf=read_config(run)
        self.update(conf)

        self.simc = read_config(self['sim'])

        # for now just set this
        #self.npair=6000 # we do two times this, pairs at 90 degrees

        self._set_gprior()

        simpars=self.get('simpars',{})
        self.shapesim = ShapeSim(self['sim'], **simpars)

    def get_npair(self, is2n):
        s2n = shapesim.get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']
        nellip = shapesim.get_s2n_nrepeat(s2n, fac=s2n_fac)

        if nellip < self['min_gcount']:
            nellip=self['min_gcount']
        return nellip


    def process_trial(self, is2, is2n, itrial):

        gsum=numpy.zeros(2)
        gsum_means=numpy.zeros(2)
        nsum_means=numpy.zeros(2)

        wsum=numpy.zeros(2)
        nsum=0

        npair=self.get_npair(is2n)
        gvals=self.gprior.sample1d(npair)

        gmm_list=[]
        for i,g in enumerate(gvals):
            imd1,imd2=self._make_ring_pair(g, is2, is2n)
            print '%d/%d' % (i+1,npair)

            mcm1=self._run_mcm(imd1)
            mcm2=self._run_mcm(imd2)

            gmm1=mcm1.get_like_mixture(ngauss=3)
            gmm2=mcm2.get_like_mixture(ngauss=3)

            gmm_list += [gmm1,gmm2]

            res1=mcm1.get_result()
            res2=mcm2.get_result()

            gsum_means += res1['g']
            nsum_means += 1

            wsum[0] += (1./res1['gcov'][0,0])
            wsum[1] += (1./res1['gcov'][1,1])

            gsum_means += res2['g']
            nsum_means += 1

            wsum[0] += (1./res2['gcov'][0,0])
            wsum[1] += (1./res2['gcov'][1,1])

            err=sqrt(1./wsum)
            mn=gsum_means/nsum_means
            print '  g1m: %s +/- %s    g2m: %s +/- %s' % (mn[0],err[0],mn[1],err[1])



        self._ng1=200
        self._ng2=200
        g1range=[mn[0]-8.*err[0], mn[0]+8.*err[0]]
        g2range=[mn[1]-8.*err[1], mn[1]+8.*err[1]]
        lnp_surf = self._get_cum_loglike_surf(gmm_list, 
                                              g1range=g1range, 
                                              g2range=g2range)

        self.lnp_surf = lnp_surf - lnp_surf.max()

        g1max,g2max=self._get_maxlike_loc(lnp_surf)
        print 'maxlike:     ',g1max,g2max

        self._test_gmm_mean(gmm_list)
        if False:
            import images
            plt=images.multiview(exp(self.lnp_surf),
                                 xdr=g2range, 
                                 ydr=g1range,
                                 xlabel='g2', ylabel='g1',
                                 show=False)
            plt_log=images.multiview(self.lnp_surf,
                                 xdr=g2range, 
                                 ydr=g1range,
                                 xlabel='g2', ylabel='g1',
                                 show=False)

            fname='/tmp/test.png'
            fname_log='/tmp/test-log.png'
            print fname
            print fname_log
            plt.write_img(1024,1024,fname)
            plt_log.write_img(1024,1024,fname_log)

        #trials=shd.get_trials()
        #shapesim.write_output(self['run'], is2, is2n, trials)

        #shapesim.write_output(self['run'], is2, is2n, out, itrial=itrial)

    def _test_gmm_mean(self, gmm_list):
        wsum=0.0
        g1censum = 0.0
        g2censum = 0.0
        for i,gmm in enumerate(gmm_list):
            wsum += gmm.weights_.sum()
            g1censum += (gmm.means_[:,0]*gmm.weights_).sum()
            g2censum += (gmm.means_[:,1]*gmm.weights_).sum()


        g1mean=g1censum/wsum
        g2mean=g2censum/wsum

        print 'From gauss mix:'
        print '  g1m: %s  g2m: %s ' % (g1mean,g2mean)

    def _get_cum_loglike_surf_slow(self, gmm_list, g1range=[-1,1], g2range=[-1,1]):
        g1range=numpy.array(g1range).clip(-1.,1.)
        g2range=numpy.array(g2range).clip(-1.,1.)

        lnp_surf = numpy.zeros( (self._ng1,self._ng2) )

        g1vals = numpy.linspace(g1range[0], g1range[1], self._ng1)
        g2vals = numpy.linspace(g2range[0], g2range[1], self._ng2)

        g1g2=numpy.zeros( (1,2) )

        print 'accumulating loglike grid'

        max_lnprob=-9999.e40
        for i1 in xrange(self._ng1):
            g1=g1vals[i1]
            print i1
            for i2 in xrange(self._ng2):
                g2=g2vals[i2]

                g1g2[0,0] = g1
                g1g2[0,1] = g2
                lnp=0.0
                for gmm in gmm_list:
                    lnp += gmm.score(g1g2)

                lnp_surf[i1,i2] = lnp

                if lnp > max_lnprob:
                    max_lnprob = lnp
                    g1max = g1
                    g2max = g2
        return lnp_surf, g1max, g2max

    def _get_cum_loglike_surf(self, gmm_list, g1range=[-1,1], g2range=[-1,1]):
        import esutil as eu

        g1range=numpy.array(g1range).clip(-1.,1.)
        g2range=numpy.array(g2range).clip(-1.,1.)

        self._g1range=g1range
        self._g2range=g2range
        mg=numpy.mgrid[0:self._ng1, 0:self._ng2]
        self._rows = mg[0].copy()
        self._cols = mg[1].copy()
        self._g1grid = eu.numpy_util.arrscl(self._rows, g1range[0], g1range[1])
        self._g2grid = eu.numpy_util.arrscl(self._cols, g2range[0], g2range[1])
        
        self._g1g2_points = numpy.zeros( (self._ng1*self._ng2, 2) )
        self._g1g2_points[:,0] = self._g1grid.ravel()
        self._g1g2_points[:,1] = self._g2grid.ravel()

        print 'accumulating loglike'
        for i,gmm in enumerate(gmm_list):
            lnp= self._get_gmm_loglike_surface(gmm)
            if i==0:
                lnp_surf=lnp
            else:
                lnp_surf += lnp

        print
        return lnp_surf

    def _get_gmm_loglike_surface(self, gmm):
        surf=gmm.score(self._g1g2_points)
        surf=surf.reshape(self._ng1, self._ng2)
        return surf

    def _get_maxlike_loc(self, lnprob):

        w=lnprob.argmax()
        #row=(self._rows.ravel())[w]
        #col=(self._cols.ravel())[w]

        g1=( self._g1grid.ravel() )[w]
        g2=( self._g2grid.ravel() )[w]

        return g1,g2

    def _run_mcm(self, imd):
        mcm=MCM(imd['image'],imd['ivar'],imd['psf'],self.gprior,
                self['fitmodel'], 
                self['nwalkers'],self['burnin'],self['nstep'],
                ares=imd['ares'],
                make_plots=self['make_plots'],
                mca_a=self['mca_a'])
        return mcm

    def _make_ring_pair(self, g, is2, is2n):
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
        
        theta = 360.0*numpy.random.random()
        theta2 = theta + 90.0

        ellip=lensing.util.g2e(g)

        ci,ares,psf = self._get_ci_ares_psf(s2, ellip, theta, s2n, s2n_psf)
        ci2,ares2,psf2 = self._get_ci_ares_psf(s2, ellip, theta2, s2n, s2n_psf)

        imd1={'image':ci.image,
              'ivar':1./ci['skysig']**2,
              'ares':ares,
              'psf':psf,
              'model':'gauss'}
        imd2={'image':ci2.image,
              'ivar':1./ci2['skysig']**2,
              'ares':ares2,
              'psf':psf2,
              'model':'gauss'}

        return imd1,imd2

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


class MCM:
    def __init__(self, image, ivar, psf, gprior, model, nwalkers, burnin, nstep, **keys):
        """
        mcmc sampling of posterior.

        Two modes of operation - send a center guess and admom will
        be run internally, or send ares=, with wrow,wcol,Irr,Irc,Icc

        parameters
        ----------
        image:
            sky subtracted image as a numpy array
        ivar:
            1/(Error per pixel)**2
        psf:
            The psf gaussian mixture as a GMix object
        cen:
            The center guess.  Ignored if ares= is sent.
        gprior:
            The prior on the g1,g2 surface, just for drawing
        model:
            Type of model, gexp, gdev, gauss

        nwalkers: 
            Number of walkers, default 20
        nstep: 
            Number of steps in MCMC chain, default 200
        burnin: 
            Number of burn in steps
        mca_a: optional
            For affine invariant chain, default 2
        ares: optional
            The output from a run of admom.  The whyflag
            field must be zero.
        cen_width: bool
            Use this as a width on the prior,
            with the center set the adaptive moments solution.
            Default is broad, 1.0
        """
        
        self.make_plots=keys.get('make_plots',False)

        # cen1,cen2,e1,e2,T,p
        self.npars=6

        self.image=image
        self.ivar=float(ivar)
        self.model=model

        self.psf_gmix=psf

        self.gprior=gprior

        self.Tprior=keys.get('Tprior',None)

        self.nwalkers=nwalkers
        self.burnin=burnin
        self.nstep=nstep
        self.mca_a=keys.get('mca_a',2.0)
        
        self.cen_guess=keys.get('cen',None)
        self.ares=keys.get('ares',None)

        self.cen_width=keys.get('cen_width',1.0)

        if self.cen_guess is None and self.ares is None:
            raise ValueError("send cen= or ares=")
        if self.ares is not None and self.ares['whyflag']!=0:
            raise ValueError("If you enter ares it must have "
                             "whyflag==0")

        self.counts=self.image.sum()

        self._go()

    def get_result(self):
        return self._result

    def get_gmix(self):
        epars=get_estyle_pars(self._result['pars'])
        return self._get_convolved_gmix(epars)

    def get_trials(self):
        self._trials

    def _get_convolved_gmix(self, epars):
        """
        epars must be in e1,e2 space
        """
        if self.model=='gauss':
            type='coellip'
        else:
            type=self.model
        gmix0=gmix_image.GMix(epars, type=type)
        gmix=gmix0.convolve(self.psf_gmix)
        return gmix


    def _go(self):
        self.sampler=self._do_trials()

        self._trials  = self.sampler.flatchain

        lnprobs = self.sampler.lnprobability.reshape(self.nwalkers*self.nstep)
        self.lnprobs = lnprobs - lnprobs.max()

        self._calc_result()

        if self.make_plots:
            self._doplots()

    def _get_sampler(self):
        import emcee
        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._calc_lnprob,
                                        a=self.mca_a)
        return sampler

    def _do_trials(self):

        tau_max=0.10
        while True:
            sampler = self._get_sampler()
            guess=self._get_guess()
            pos, prob, state = sampler.run_mcmc(guess, self.burnin)
            sampler.reset()
            pos, prob, state = sampler.run_mcmc(pos, self.nstep)

            try:
                acor=sampler.acor
                tau = (sampler.acor/self.nstep).max()
                if tau > tau_max:
                    print "tau",tau,"greater than",tau_max
                    self.nwalkers = self.nwalkers*2
                    sampler = self._get_sampler()
                else:
                    break
            except:
                # something went wrong with acor, run some more
                pass
        self._tau=tau
        return sampler

    def _calc_lnprob(self, pars):
        epars=get_estyle_pars(pars)
        if epars is None:
            return LOWVAL

        logprob = self._get_loglike_c(epars)

        cp = self.cenprior.lnprob(pars[0:2])
        logprob += cp

        if self.Tprior is not None:
            Tp = self.Tprior.lnprob(pars[4])
            logprob += Tp

        return logprob

 
    def _get_loglike_c(self, epars):
        """
        These epars are in e space
        """
        from gmix_image import render

        gmix=self._get_convolved_gmix(epars)

        loglike,s2n,flags=\
            render._render.loglike(self.image, gmix, self.ivar)

        if flags != 0:
            return LOWVAL
        return loglike


    def _calc_result(self):
        """
        We marginalize over all parameters but g1,g2, which
        are index 0 and 1 in the pars array
        """
        import mcmc

        g=numpy.zeros(2)
        gsum=numpy.zeros(2)
        gcov=numpy.zeros((2,2))

        gsum[0] = self._trials[:,2].sum()
        gsum[1] = self._trials[:,3].sum()

        pars,pcov = mcmc.extract_stats(self._trials)

        g[:] = pars[2:4]
        gcov[:,:] = pcov[2:4, 2:4]
 
        arates = self.sampler.acceptance_fraction
        arate = arates.mean()

        max_epars=self._get_maxprob_epars()
        gmix=self._get_convolved_gmix(max_epars)

        stats=calculate_some_stats(self.image, 
                                   self.ivar, 
                                   gmix,
                                   self.npars)

        Tmean=pars[4]
        Terr=sqrt(pcov[4,4])
        Ts2n=pars[4]/sqrt(pcov[4,4])

        self._result={'model':self.model,
                      'g':g,
                      'gcov':gcov,
                      'gsum':gsum,
                      'nsum':self._trials.shape[0],
                      'pars':pars,
                      'perr':sqrt(diag(pcov)),
                      'pcov':pcov,
                      'Tmean':Tmean,
                      'Terr':Terr,
                      'Ts2n':Ts2n,
                      'arate':arate}

        self._result.update(stats)

    def _get_maxprob_epars(self):
        wmax=self.lnprobs.argmax()
        max_pars = self._trials[wmax,:].copy()
        max_epars=get_estyle_pars(max_pars)
        return max_epars

    def get_maxprob_model(self):
        from gmix_image import render
        max_epars=self._get_maxprob_epars()
        gmix=self._get_convolved_gmix(max_epars)
        model=render.gmix2image(gmix,self.image.shape)
        return model

    def _run_admom(self, image, ivar, cen, Tguess):
        import admom

        ntry=10
        for i in xrange(ntry):
            ares = admom.admom(image,
                                    cen[0],
                                    cen[1],
                                    sigsky=sqrt(1/ivar),
                                    guess=Tguess/2,
                                    nsub=1)
            if ares['whyflag']==0:
                break
        if i==(ntry-1):
            raise ValueError("admom failed %s times" % ntry)

        return ares

    def get_like_mixture(self, ngauss=3, niter=100, min_covar=1.0e-06):
        #from scikits.learn import mixture
        from sklearn import mixture

        g1vals=self._trials[:,2]
        g2vals=self._trials[:,3]
        data=numpy.zeros((g1vals.size,2))
        data[:,0] = g1vals
        data[:,1] = g2vals

        #gmm=mixture.GMM(n_states=ngauss,cvtype='full')
        #gmm.fit(data, n_iter=niter, min_covar=min_covar)
        gmm=mixture.GMM(n_components=ngauss,covariance_type='full',min_covar=min_covar)
        #gmm=mixture.DPGMM(n_components=ngauss,covariance_type='full',min_covar=min_covar)
        gmm.fit(data)

        if not hasattr('gmm','covars_'):
            gmm.covars_ = gmm._get_covars()
            """
            for i,prec in enumerate(gmm.precs_):
                cov = numpy.linalg.inv(prec)
                gmm.covars_.append(cov)
            """


        return gmm

    def get_like_hist(self, 
                      g1range=[-1.0,1.0], 
                      ng1=2000, 
                      g2range=[-1.0,1.0],
                      ng2=2000):
        import esutil as eu

        g1vals=self._trials[:,2]
        g2vals=self._trials[:,3]
        h2d = eu.stat.histogram2d(g1vals, g2vals, 
                                  nx=ng1, xmin=g1range[0], xmax=g1range[1], 
                                  ny=ng2, ymin=g2range[0], ymax=g2range[1], 
                                  more=True)
        return h2d

    def _get_guess(self):
        from gmix_image.priors import CenPrior
        if self.ares is None:
            self.ares=self._run_admom(self.image, self.ivar, 
                                      self.cen_guess, 8.0)

        
        cen=[self.ares['wrow'],self.ares['wcol']]
        self.cenprior=CenPrior(cen, [self.cen_width]*2)

        Tadmom=self.ares['Irr'] + self.ares['Icc']

        guess=numpy.zeros( (self.nwalkers,self.npars) )

        guess[:,0]=self.cenprior.cen[0] + 0.01*srandu(self.nwalkers)
        guess[:,1]=self.cenprior.cen[1] + 0.01*srandu(self.nwalkers)

        g1rand,g2rand=self.gprior.sample2d(self.nwalkers)
        guess[:,2] = g1rand
        guess[:,3] = g2rand

        guess[:,4] = Tadmom*(1 + 0.1*srandu(self.nwalkers))
        guess[:,5] = self.counts*(1 + 0.1*srandu(self.nwalkers))

        self._guess=guess
        return guess


    def _doplots(self):
        import mcmc
        import biggles
        import esutil as eu

        biggles.configure("default","fontsize_min",1.2)
        tab=biggles.Table(6,2)

        cen1vals=self._trials[:,0]
        cen2vals=self._trials[:,1]
        Tvals=self._trials[:,4]
        g1vals=self._trials[:,2]
        g2vals=self._trials[:,3]
        g1lab=r'$g_1$'
        g2lab=r'$g_2$'

        ampvals=self._trials[:,5]

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
        print 's2n weighted:',res['s2n_w']
        print 'acceptance rate:',res['arate'],'mca_a',self.mca_a
        print 'T:  %.16g +/- %.16g' % (Tvals.mean(), Tvals.std())

        print_pars(self._result['pars'])
        print_pars(sqrt(diag(self._result['pcov'])))
        print 'g1: %.16g +/- %.16g' % (g[0],errs[0])
        print 'g2: %.16g +/- %.16g' % (g[1],errs[1])
        print 'chi^2/dof: %.3f/%i = %f' % (res['chi2per']*res['dof'],res['dof'],res['chi2per'])
        print 'probrand:',res['fit_prob']

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

        if True: 
            import images
            nx = ny = 20
            levels=8
            h2d = eu.stat.histogram2d(g1vals, g2vals, nx=nx, ny=ny,more=True)
            images.view(h2d['hist'],
                             xdr=[h2d['xcenter'][0], h2d['xcenter'][-1]],
                             ydr=[h2d['ycenter'][0], h2d['ycenter'][-1]],
                             xlabel='g1', ylabel='g2', levels=levels)
            if True:
                gmm=self.get_like_mixture(ngauss=3)
                print gmm.weights_
                print gmm.means_

                covmean=0*gmm.covars_[0].copy()
                for i,cov in enumerate(gmm.covars_):
                    covmean += gmm.weights_[i]*cov
                """
                covmean=0*gmm.precs_[0].copy()
                for i,prec in enumerate(gmm.precs_):
                    cov = numpy.linalg.inv(prec)
                    covmean += gmm.weights_[i]*cov
                """

                wsum=gmm.weights_.sum()
                covmean /= wsum

                print 'weighted means'
                print (gmm.means_[:,0]*gmm.weights_).sum()/wsum,
                print (gmm.means_[:,1]*gmm.weights_).sum()/wsum
                
                print 'weighted error'
                print sqrt(diag(covmean))

                rnd=gmm.sample(g1vals.size*10)
                print rnd.shape
                rh2d = eu.stat.histogram2d(rnd[:,0], rnd[:,1], nx=nx, ny=ny,
                                           xmin=h2d['xlow'][0], xmax=h2d['xhigh'][-1],
                                           ymin=h2d['ylow'][0], ymax=h2d['yhigh'][-1],
                                           more=True)
                rh2d['hist'] = rh2d['hist']*float(g1vals.size)/float(rh2d['hist'].sum())
                images.compare_images(h2d['hist'], rh2d['hist'])

        if False: 
            import images
            model=self.get_maxprob_model()
            images.compare_images(self.image,model,
                                  label1='image [%d,%d]' % self.image.shape,
                                  label2='model')

        key=raw_input('hit a key (q to quit): ')
        if key=='q':
            stop
        print





