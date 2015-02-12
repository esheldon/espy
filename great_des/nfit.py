from __future__ import print_function
import os
from sys import stderr,stdout
import time
import numpy
from numpy import array, sqrt, zeros, log, exp

import ngmix
from ngmix.fitting import print_pars
from ngmix.gexceptions import GMixMaxIterEM, GMixRangeError
from ngmix.observation import Observation

from gmix_meds.util import FromPSFGuesser, FixedParsGuesser, FromParsGuesser

import meds

# starting new values for these
DEFVAL=-9999
PDEFVAL=9999
BIG_DEFVAL=-9.999e9
BIG_PDEFVAL=9.999e9


NO_CUTOUTS=2**0
PSF_FIT_FAILURE=2**1
PSF_LARGE_OFFSETS=2**2
EXP_FIT_FAILURE=2**3
DEV_FIT_FAILURE=2**4

BOX_SIZE_TOO_BIG=2**5
S2N_TOO_LOW=2**6

NO_ATTEMPT=2**30

#PSF_S2N=1.e6
PSF_OFFSET_MAX=0.25
PSF_TOL=1.0e-5
EM_MAX_TRY=3
EM_MAX_ITER=100

ISAMP_BAD_COV=2**7

_CHECKPOINTS_DEFAULT_MINUTES=[10,30,60,90]

class MedsFit(dict):
    def __init__(self, meds_file, truth_file, psf_file, **keys):
        """
        This differs from nfit in gmix_meds in that there is no coadd

        parameters
        ----------
        meds_file: string
            meds path
        truth_file: string
            truth info to get the PSF file
        psf_file: string
            file holding the psf images

        The following are sent through keywords

        obj_range: optional
            a 2-element sequence or None.  If not None, only the objects in the
            specified range are processed and get_data() will only return those
            objects. The range is inclusive unlike slices.
        psf_model: string, int
            e.g. "em2"
        psf_offset_max: optional
            max offset between multi-component gaussians in psf models

        checkpoints: number, optional
            Times after which to checkpoint, seconds
        checkpoint_file: string, optional
            File which will hold a checkpoint.
        checkpoint_data: dict, optional
            The data representing a previous checkpoint, object and
            psf fits
        """

        self['meds_file']=meds_file
        self['truth_file']=truth_file
        self['psf_file']=psf_file

        self._load_meds()
        self._load_truth()
        self._load_psf_fobj()

        self.update(keys)

        self['psf_ngauss']=get_em_ngauss(self['psf_model'])
        self['psf_sigma_guess']=self['psf_fwhm_guess']/2.3548200450309493

        # make sure some optional ones are set
        self['shear_expand']=self.get('shear_expand',None)
        self['obj_range']=self.get('obj_range',None)
        self['make_plots']=self.get('make_plots',False)

        self._setup_checkpoints()
        self._set_index_list()

        if self._checkpoint_data is None:
            self._make_struct()

    def get_data(self):
        """
        Get the data structure.  If a subset was requested, only those rows are
        returned.
        """
        return self.data

    def get_meds_meta(self):
        """
        get copies of the meta data
        """
        return self.meds_meta.copy()

    def get_magzp(self):
        """
        Get the magnitude zero point.
        """
        return self.meds_meta['magzp_ref'][0]

    def do_fits(self):
        """
        Fit all objects in our list
        """

        t0=time.time()

        last=self.index_list[-1]
        num=len(self.index_list)

        for dindex in xrange(num):

            if self.data['processed'][dindex]==1:
                # was checkpointed
                continue

            self.dindex = dindex
            self.mindex = self.index_list[dindex]

            self._load_psf_image()
            self._load_image_and_weight()

            print( 'index: %d:%d' % (self.mindex,last), )

            self._fit_all()

            tm=time.time()-t0

            self._try_checkpoint(tm) # only at certain intervals

        tm=time.time()-t0
        print("time:",tm)
        print("time per:",tm/num)


    def _fit_all(self):
        """
        Process the indicated object through the requested fits
        """

        t0=time.time()

        self.res={'flags':0}

        # for checkpointing
        self.data['processed'][self.dindex]=1
        self.data['nimage_use'][self.dindex] = 1
        self.data['number'][self.dindex] = self.meds['number'][self.mindex]

        self._fit_psf()

        if self.res['flags'] != 0:
            print('not fitting object due to psf failure')
            return

        self._fit_galaxy()

        self._copy_to_output()

        self.data['time'][self.dindex] = time.time()-t0

    def _fit_psf(self):
        """
        Fit a psf to a single gauss and more complex model
        """
        from ngmix import GMixMaxIterEM

        self.fitting_galaxy=False

        # first 1 gaussian
        jacobian=self._get_jacobian(self.psf_cen_guess_pix)

        psf_obs = Observation(self.psf_image,
                              jacobian=jacobian)

        fitter1=self._fit_em_1gauss(psf_obs,
                                    self['psf_sigma_guess'])

        if fitter1 is None:
            self.res['flags'] = PSF_FIT_FAILURE
            print("psf em1 failure at object",self.mindex)
            return


        psf_gmix1=fitter1.get_gmix()
        print("psf fwhm1:",2.35*sqrt( psf_gmix1.get_T()/2. ))

        self.psf_gmix1=psf_gmix1
        self.res['psf_gmix_em1']=psf_gmix1

        sigma_guess_new = sqrt( psf_gmix1.get_T()/2. )
        fitter=self._fit_em_ngauss(psf_obs,
                                   sigma_guess_new,
                                   self['psf_ngauss'])

        if fitter is None:
            self.res['flags'] = PSF_FIT_FAILURE
            print("psf failure at object",index)
        else:
            psf_gmix = fitter.get_gmix()
            print("psf fit:")
            print(psf_gmix)
            print("psf fwhm:",2.35*sqrt( psf_gmix.get_T()/2. ))

            self.psf_gmix=psf_gmix
            self.res['psf_gmix']=psf_gmix

            psf_obs.set_gmix(psf_gmix)
            self.res['psf_obs'] = psf_obs

            if self['make_plots']:
                self._compare_psf(fitter)

    def _fit_galaxy(self):
        """
        Fit psf flux and other models
        """

        self.fitting_galaxy=True

        self.res['jacobian'] = self._get_jacobian(self.gal_cen_guess_pix)

        if self['trim_image']:
            self._make_trimmed_image()

        obs = self._get_observation()
        obs.set_psf( self.res['psf_obs'] )

        self._fit_galaxy_psf_flux(obs)
        self._fit_galaxy_model(obs)

    def _fit_galaxy_psf_flux(self, obs):
        """
        Get flux fitting the psf model

        Call this after fitting galaxy with em 1
        """

        print("    fitting galaxy with psf")

        fitter=ngmix.fitting.TemplateFluxFitter(obs, do_psf=True)
        fitter.go()

        res=fitter.get_result()


        self.res['psf_flux_flags'] = res['flags']
        self.res['psf_flux'] = res['flux']
        self.res['psf_flux_err'] = res['flux_err']
        self.res['psf_flux_s2n'] = res['flux']/res['flux_err']

        print("        psf flux: %g +/- %g" % (res['flux'],res['flux_err']))


    def _fit_galaxy_model(self, obs):
        """
        Run through and fit all the models
        """
        import gmix_meds.nfit

        res=self.res
        model=self['fit_model']

        ps2n=self.res['psf_flux_s2n']
        if ps2n < self['min_s2n']:
            print("    psf s/n too low:",ps2n)
            self.res['flags'] |= S2N_TOO_LOW
            return

        print('    fitting',model,'using maxlike')
        max_guesser=self._get_guesser_from_psf()
        max_fitter=self._fit_simple_max(obs,model,max_guesser) 

        self._print_galaxy_res(max_fitter)

        # faking errors, as they are not needed
        print('    fitting',model,'using mcmc')
        max_res=max_fitter.get_result()
        model_guesser = FromParsGuesser(max_res['pars'],max_res['pars']*0.1)
        fitter=self._fit_simple_mcmc(obs,
                                     model,
                                     model_guesser)

        #res['gauss_fitter'] = gauss_fitter
        res['max_fitter'] = max_fitter
        res['galaxy_fitter'] = fitter
        res['galaxy_res'] = fitter.get_result()
        res['max_res'] = max_fitter.get_result()
        res['flags'] = res['galaxy_res']['flags']

        self._add_shear_info(res['galaxy_res'], fitter)

        self._print_galaxy_res(fitter)

        if self['make_plots']:
            self._do_gal_plots(res['galaxy_fitter'])


    def _fit_simple_mcmc(self, obs, model, guesser):
        """
        Fit the simple model
        """

        from ngmix.fitting import MCMCSimple

        # note flat on g!
        prior=self['search_prior']

        epars=self['emcee_pars']

        guess=guesser(n=epars['nwalkers'])

        fitter=MCMCSimple(obs,
                          model,
                          prior=prior,
                          use_logpars=self['use_logpars'],
                          nwalkers=epars['nwalkers'],
                          mca_a=epars['a'])

        pos=fitter.run_mcmc(guess,epars['burnin'])
        pos=fitter.run_mcmc(pos,epars['nstep'],thin=epars['thin'])

        g_prior = self['prior'].g_prior
        trials = fitter.get_trials()
        weights = g_prior.get_prob_array2d(trials[:,2], trials[:,3])

        fitter.calc_result(weights=weights)

        self.weights=weights

        return fitter


    def _get_guess_simple(self):
        """
        width is relative for T and counts
        """

        width = 0.01


        res=self.res
        gmix = res['em_gmix']
        g1,g2,T = gmix.get_g1g2T()
        F = res['em_gauss_flux']


        nwalkers = self['nwalkers']
        guess=numpy.zeros( (nwalkers, 6) )

        guess[:,0] = width*srandu(nwalkers)
        guess[:,1] = width*srandu(nwalkers)

        guess_shape=get_shape_guess(g1,g2,nwalkers,width=width)
        guess[:,2]=guess_shape[:,0]
        guess[:,3]=guess_shape[:,1]

        #guess[:,4] = log10( get_positive_guess(T,nwalkers,width=width) )
        guess[:,4] = ( get_positive_guess(T,nwalkers,width=width) )

        # got anything better?
        guess[:,5] = log10( get_positive_guess(F,nwalkers,width=width) )

        return guess


    def _fit_simple_max(self, obs, model, guesser):
        from ngmix.fitting import MaxSimple

        max_pars=self['max_pars']
        if max_pars['method']=='lm':
            return self._fit_simple_lm(obs, model, guesser)

        prior=self['search_prior']
        print("search prior:",prior)
        fitter=MaxSimple(obs,
                         model,
                         prior=prior,
                         use_logpars=self['use_logpars'],
                         **max_pars)

        for i in xrange(max_pars['ntry']):
            guess=guesser(prior=prior)
            print_pars(guess,front='    guess: ')

            fitter.run_max(guess, **max_pars)
            res=fitter.get_result()
            if res['flags']==0:
                break
        return fitter

    def _fit_simple_lm(self, obs, model, guesser):
        from ngmix.fitting import LMSimple

        max_pars=self['max_pars']

        prior=self['search_prior']
        fitter=LMSimple(obs,
                        model,
                        prior=prior,
                        use_logpars=self['use_logpars'],
                        lm_pars=max_pars['lm_pars'])

        for i in xrange(max_pars['ntry']):
            guess=guesser(prior=prior)
            print_pars(guess,front='    guess: ')

            fitter.run_max(guess)
            res=fitter.get_result()
            if res['flags']==0:
                break
        return fitter


    def _get_guesser_from_psf(self):
        """
        take flux guesses from psf take canonical center (0,0)
        and near zero ellipticity.  Size is taken from around the
        expected psf size, which is about 0.9''

        The size will often be too big

        """
        print('        getting guess from psf')

        dindex=self.dindex
        data=self.data

        psf_flux=self.res['psf_flux']
        psf_flux=psf_flux.clip(min=0.1, max=1.0e9)

        # arbitrary
        T = 2*(0.9/2.35)**2

        if self['use_logpars']:
            scaling='log'
        else:
            scaling='linear'

        guesser=FromPSFGuesser(T, psf_flux, scaling=scaling)
        return guesser


    def _load_image_and_weight(self):
        """
        Load the image, weight map
        """

        dindex=self.dindex
        mindex=self.mindex
        
        self.image = self.meds.get_cutout(mindex,0).astype('f8')
        if self['noisefree']:
            self.wt = 0*self.image + (1.0/self['skynoise']**2)
        else:
            self.wt = self.meds.get_cutout(mindex,0,type='weight')

        self.gal_cen_guess_pix=(array(self.image.shape)-1)/2.
        self.data['nimage_tot'][dindex] = 1

    def _load_psf_image(self):
        """
        Get psf images for the SE images
        associated with the cutouts
        """

        psf_id = self.truth['id_psf'][self.mindex]
        self.psf_image = self.psf_fobj[psf_id].read().astype('f8')
        
        self.psf_cen_guess_pix=(array(self.psf_image.shape)-1)/2.

    def _get_jacobian(self, cen):
        jdict = self.meds.get_jacobian(self.mindex,0)

        jacobian=ngmix.Jacobian(cen[0],
                                cen[1],
                                jdict['dudrow'],
                                jdict['dudcol'],
                                jdict['dvdrow'],
                                jdict['dvdcol'])

        return jacobian


    def _fit_em_1gauss(self, obs, sigma_guess):
        """
        Just run the fitter

        cen is the global cen not within jacobian
        """
        return self._fit_with_em(obs, sigma_guess, 1)

    def _fit_em_ngauss(self, obs, sigma_guess, ngauss):
        """
        Start with fit from using 1 gauss, which must be entered

        cen is the global cen not within jacobian
        """

        fitter=self._fit_with_em(obs, sigma_guess, ngauss)

        return fitter



    def _fit_with_em(self, obs_in, sigma_guess, ngauss):
        """
        Fit the image using EM
        """

        im_with_sky, sky = ngmix.em.prep_image(obs_in.image)

        new_obs=Observation(im_with_sky, jacobian=obs_in.jacobian)

        ntry,maxiter,tol = self._get_em_pars()
        for i in xrange(ntry):
            guess = self._get_em_guess(sigma_guess, ngauss)
            try:
                fitter=self._do_fit_em_with_full_guess(new_obs,
                                                       sky,
                                                       guess)
                tres=fitter.get_result()
                #print("em numiter:",tres['numiter'])
                break
            except GMixMaxIterEM:
                fitter=None
            except GMixRangeError:
                fitter=None

        return fitter

    def _do_fit_em_with_full_guess(self,
                                   obs,
                                   sky,
                                   guess):

        ntry,maxiter,tol = self._get_em_pars()

        fitter=ngmix.em.GMixEM(obs)
        fitter.go(guess, sky, maxiter=maxiter, tol=tol)

        return fitter


    def _get_em_guess(self, sigma, ngauss):
        """
        Guess for the EM algorithm
        """

        if ngauss==1:
            return self._get_em_guess_1gauss(sigma)
        elif ngauss==2:
            return self._get_em_guess_2gauss(sigma)
        elif ngauss==3:
            return self._get_em_guess_3gauss(sigma)
        else:
            return self._get_em_guess_ngauss(sigma,ngauss)

    def _get_em_guess_1gauss(self, sigma):

        sigma2 = sigma**2
        pars=array( [1.0 + 0.1*srandu(),
                     0.2*srandu(),
                     0.2*srandu(), 
                     sigma2*(1.0 + 0.5*srandu()),
                     0.2*sigma2*srandu(),
                     sigma2*(1.0 + 0.5*srandu())] )

        return ngmix.gmix.GMix(pars=pars)

    def _get_em_guess_2gauss(self, sigma):

        sigma2 = sigma**2

        pars=array( [_em2_pguess[0],
                     0.1*srandu(),
                     0.1*srandu(),
                     _em2_fguess[0]*sigma2*(1.0 + 0.1*srandu()),
                     0.0,
                     _em2_fguess[0]*sigma2*(1.0 + 0.1*srandu()),

                     _em2_pguess[1],
                     0.1*srandu(),
                     0.1*srandu(),
                     _em2_fguess[1]*sigma2*(1.0 + 0.1*srandu()),
                     0.0,
                     _em2_fguess[1]*sigma2*(1.0 + 0.1*srandu())] )


        return ngmix.gmix.GMix(pars=pars)

    def _get_em_guess_3gauss(self, sigma):

        sigma2 = sigma**2

        glist=[]

        pars=array( [_em3_pguess[0]*(1.0+0.1*srandu()),
                     0.1*srandu(),
                     0.1*srandu(),
                     _em3_fguess[0]*sigma2*(1.0 + 0.1*srandu()),
                     0.01*srandu(),
                     _em3_fguess[0]*sigma2*(1.0 + 0.1*srandu()),

                     _em3_pguess[1]*(1.0+0.1*srandu()),
                     0.1*srandu(),
                     0.1*srandu(),
                     _em3_fguess[1]*sigma2*(1.0 + 0.1*srandu()),
                     0.01*srandu(),
                     _em3_fguess[1]*sigma2*(1.0 + 0.1*srandu()),

                     _em3_pguess[2]*(1.0+0.1*srandu()),
                     0.1*srandu(),
                     0.1*srandu(),
                     _em3_fguess[2]*sigma2*(1.0 + 0.1*srandu()),
                     0.01*srandu(),
                     _em3_fguess[2]*sigma2*(1.0 + 0.1*srandu())]

                  )


        return ngmix.gmix.GMix(pars=pars)


    def _get_em_guess_ngauss(self, sigma, ngauss):

        sigma2 = sigma**2

        glist=[]

        for i in xrange(ngauss):
            if i > 2:
                ig=2
            else:
                ig=i

            glist += [_em3_pguess[ig]*(1.0 + 0.1*srandu()),
                      0.1*srandu(),
                      0.1*srandu(),
                      _em3_fguess[ig]*sigma2*(1.0 + 0.1*srandu()),
                      0.1*srandu(),
                      _em3_fguess[ig]*sigma2*(1.0 + 0.1*srandu())]

        pars=array(glist)

        return ngmix.gmix.GMix(pars=pars)


    def _get_em_pars(self):
        if self.fitting_galaxy:
            return self['gal_em_ntry'], self['gal_em_maxiter'], self['gal_em_tol']
        else:
            return self['psf_em_ntry'], self['psf_em_maxiter'], self['psf_em_tol']

    def _make_trimmed_image(self):
        res=self.res

        cen_jacob=array( res['jacobian'].get_cen() )

        em_gmix=res['em_gmix']

        T=em_gmix.get_T()
        cen=array(em_gmix.get_cen())/self['pixel_scale']

        cen = cen + cen_jacob

        sigma=sqrt(T/2.)/self['pixel_scale']
        #print("cen:",cen)
        #print("from gauss:",res['em_gmix'])
        #print("got T:",T,"sigma (pixels):",sigma)

        radius = sigma*self['trim_nsigma']

        minrow=int(cen[0]-radius)
        maxrow=int(cen[0]+radius+1)
        mincol=int(cen[1]-radius)
        maxcol=int(cen[1]+radius+1)

        sh=self.image.shape

        if minrow < 0:
            minrow=0
        if maxrow > sh[0]:
            maxrow=sh[0]

        if mincol < 0:
            mincol=0
        if maxcol > sh[1]:
            maxcol=sh[1]

        self.trimmed_image = self.image[minrow:maxrow, mincol:maxcol]
        self.trimmed_wt    = self.wt[minrow:maxrow, mincol:maxcol]
        self.trimmed_cen   = array(cen)-array([minrow,mincol])
        res['trimmed_jacobian'] = self._get_jacobian(self.trimmed_cen)

        print("    trimmed image:",self.trimmed_image.shape,
              "cen:  ",self.trimmed_cen)

    def _get_observation(self):
        """
        Get the appropriate observation
        """
        im, wt, jacob=self._get_image_data()
        return Observation(im, weight=wt, jacobian=jacob)

    def _get_image_data(self):
        """
        Get the appropriate image data for the current object.
        """

        if self['trim_image']:
            image=self.trimmed_image
            wt=self.trimmed_wt
            jacobian=self.res['trimmed_jacobian']
        else:
            image=self.image
            wt=self.wt
            jacobian=self.res['jacobian']
        return image, wt, jacobian

    def _compare_psf(self, fitter):
        """
        compare psf image to best fit model
        """
        import images

        model=self['psf_model']

        model_image = fitter.make_image(counts=self.psf_image.sum())

        plt=images.compare_images(self.psf_image,
                                  model_image,
                                  label1='psf',
                                  label2=model,
                                  show=False)

        pname='psf-resid-%s-%06d.png' % (model, self.mindex)
        print("          ",pname)
        plt.write_img(1400,800,pname)

    def _do_gal_plots(self, fitter):
        """
        Make residual plot and trials plot
        """
        self._compare_gal(fitter)
        self._make_trials_plot(fitter)
        self._plot_autocorr(fitter)

    def _compare_gal(self, fitter):
        """
        compare psf image to best fit model
        """
        import images

        model=self['fit_model']
        title = '%d %s' % (self.mindex, model)

        gmix = fitter.get_gmix()

        obs = self._get_observation()

        res=self.res
        psf_gmix = res['psf_gmix']
        gmix_conv = gmix.convolve(psf_gmix)

        image=obs.image
        model_image = gmix_conv.make_image(image.shape,
                                           jacobian=obs.jacobian)

        plt=images.compare_images(image,
                                  model_image,
                                  label1='galaxy',
                                  label2=model,
                                  show=False)
        plt.title=title
        pname='gal-resid-%06d-%s.png' % (self.mindex,model)

        resid_std = (image-model_image).std()
        print("    residual std:",resid_std)
        print("          ",pname)
        plt.write_img(1400,800,pname)


    def _make_trials_plot(self, fitter):
        """
        Plot the trials
        """

        model=self['fit_model']

        title = '%d %s' % (self.mindex, model)

        width,height=800,800
        weights=self.weights
        pdict=fitter.make_plots(title=title, weights=weights, do_triangle=True)
        pdict['trials'].title=title

        if 'wtrials' in pdict:
            wtrials_pname='wtrials-%06d-%s.png' % (self.mindex,model)
            print("          ",wtrials_pname)

            pdict['wtrials'].title=title
            pdict['wtrials'].write_img(width,height,wtrials_pname)

        trials_pname='trials-%06d-%s.png' % (self.mindex,model)
        print("          ",trials_pname)
        pdict['trials'].write_img(width,height,trials_pname)

        if 'triangle' in pdict:
            pname='triangle-%06d-%s.png' % (self.mindex,model)
            print("          ",pname)
            pdict['triangle'].savefig(pname)



    def _plot_autocorr(self, fitter):
        """
        Plot the trials
        """

        model=self['fit_model']

        trials=fitter.get_trials()

        plot_arr = plot_autocorr(trials)
        plot_arr.title = '%d %s' % (self.mindex, model)

        width,height=800,800

        pname='autocorr-%06d-%s.png' % (self.mindex,model)

        print("          ",pname)
        plot_arr.write_img(width, height, pname)



    def _print_galaxy_res(self, fitter):
        res=fitter.get_result()

        print_pars(res['pars'], front="    pars: ")

        if 'pars_err' in res:
            print_pars(res['pars_err'], front="    err:  ")
            if 'arate' in res:
                mess="            s/n: %.1f  arate: %.2f  tau: %.1f"
                tup = (res['s2n_w'],res['arate'],res['tau'])
                mess=mess % tup
                print(mess)
            elif 'efficiency' in res:
                mess="            s/n: %.1f  neff: %.1f  efficiency: %.2f"
                tup = (res['s2n_w'],res['neff'],res['efficiency'])
                mess=mess % tup
                print(mess)
            elif 's2n_w' in res:
                mess="            s/n: %.1f" % res['s2n_w']
                if 'nfev' in res:
                    mess="%s nfev: %d" % (mess,res['nfev'])

                print(mess)
        else:
            print("    NO COV PRESENT")


    def _copy_to_output(self):
        """
        Copy the galaxy fits
        """

        dindex=self.dindex
        res=self.res
        data=self.data

        # overall flags, model flags copied below
        data['flags'][dindex] = res['flags']

        if 'psf_gmix' in res:
            self._copy_galaxy_pars()

    def _copy_galaxy_pars(self):
        """
        Copy from the result dict to the output array
        """

        dindex=self.dindex
        allres = self.res

        data=self.data

        # allres flags is or'ed with the galaxy fitting flags
        if allres['flags'] != 0:
            if self['fitter']=='mcmc':
                print("        Not copying pars due to failure")
                return

            # for max fitter, we might "fail" but still have something
            # worth copying.  But if there are no errors, there is no
            # point in continuing
            if 'pars_err' not in allres:
                print("        Not copying pars due to failure")
                return

        res=allres['galaxy_res']
        res_max=allres['max_res']

        pars=res['pars']
        pars_cov=res['pars_cov']

        T=pars[4]
        T_err=sqrt(pars_cov[4, 4])

        flux=pars[5]
        flux_err=sqrt(pars_cov[5, 5])

        data['pars'][dindex,:] = pars
        data['pars_cov'][dindex,:,:] = pars_cov

        data['max_flags'][dindex] = res_max['flags']
        data['pars_max'][dindex,:] = res_max['pars']
        if 'pars_cov' in res_max:
            data['pars_max_cov'][dindex,:,:] = res_max['pars_cov']

        data['flux'][dindex] = flux
        data['flux_err'][dindex] = flux_err
        data['T'][dindex] = T
        data['T_err'][dindex] = T_err

        if T_err != 0:
            T_s2n = T/T_err
        else:
            T_s2n=-9999.0
        data['T_s2n'][dindex]=T_s2n


        data['g'][dindex,:] = res['g']
        data['g_cov'][dindex,:,:] = res['g_cov']

        if self['fitter']=='mcmc':
            data['arate'][dindex] = res['arate']
            data['tau'][dindex] = res['tau']
        elif self['fitter']=='isample':
            data['efficiency'][dindex] = res['efficiency']
            data['neff'][dindex] = res['neff']
            #data['niter'][dindex] = res['niter']

        for sn in _stat_names:
            if sn in res:
                data[sn][dindex] = res[sn]

        if 'P' in res:
            data['P'][dindex] = res['P']
            data['Q'][dindex,:] = res['Q']
            data['R'][dindex,:,:] = res['R']
        if 'g_sens' in res:
            data['g_sens'][dindex,:] = res['g_sens']

    def _load_meds(self):
        """
        Load all listed meds files
        """

        print(self['meds_file'])
        self.meds = meds.MEDS(self['meds_file'])
        self.meds_meta=self.meds.get_meta()
        self.nobj_tot = self.meds.size

        self.mindex=0
        jtmp=self._get_jacobian([0.0, 0.0])

        self['pixel_scale'] = jtmp.get_scale()
        print("pixel scale:",self['pixel_scale'])

    def _load_truth(self):
        """
        load the truth file for getting the psf index
        """
        import fitsio
        print(self['truth_file'])
        self.truth = fitsio.read(self['truth_file'],lower=True)
        
    def _load_psf_fobj(self):
        """
        Load the psf file as a FITS object
        """
        import fitsio
        print(self['psf_file'])
        self.psf_fobj = fitsio.FITS(self['psf_file'])

    def _set_index_list(self):
        """
        set the list of indices to be processed
        """
        if self['obj_range'] is None:
            start=0
            end=self.nobj_tot-1
        else:
            start=self['obj_range'][0]
            end=self['obj_range'][1]

        self.index_list = numpy.arange(start,end+1)


    def _setup_checkpoints(self):
        """
        Set up the checkpoint times in minutes and data
        """
        self['checkpoints'] = self.get('checkpoints',_CHECKPOINTS_DEFAULT_MINUTES)
        self['n_checkpoint']    = len(self['checkpoints'])
        self['checkpointed']    = [0]*self['n_checkpoint']
        self['checkpoint_file'] = self.get('checkpoint_file',None)

        self._set_checkpoint_data()

        if self['checkpoint_file'] is not None:
            self['do_checkpoint']=True
        else:
            self['do_checkpoint']=False

    def _set_checkpoint_data(self):
        """
        See if checkpoint data was sent
        """
        self._checkpoint_data=self.get('checkpoint_data',None)
        if self._checkpoint_data is not None:
            self.data=self._checkpoint_data['data']

    def _try_checkpoint(self, tm):
        """
        Checkpoint at certain intervals.  
        Potentially modified self['checkpointed']
        """

        should_checkpoint, icheck = self._should_checkpoint(tm)

        if should_checkpoint:
            self._write_checkpoint(tm)
            self['checkpointed'][icheck]=1

    def _should_checkpoint(self, tm):
        """
        Should we write a checkpoint file?
        """

        should_checkpoint=False
        icheck=-1

        if self['do_checkpoint']:
            tm_minutes=tm/60

            for i in xrange(self['n_checkpoint']):

                checkpoint=self['checkpoints'][i]
                checkpointed=self['checkpointed'][i]

                if tm_minutes > checkpoint and not checkpointed:
                    should_checkpoint=True
                    icheck=i

        return should_checkpoint, icheck

    def _write_checkpoint(self, tm):
        """
        Write out the current data structure to a temporary
        checkpoint file.
        """
        import fitsio

        print('checkpointing at',tm/60,'minutes')
        print(self['checkpoint_file'])

        with fitsio.FITS(self['checkpoint_file'],'rw',clobber=True) as fobj:
            fobj.write(self.data, extname="model_fits")

    def _add_shear_info(self, res, fitter):
        """
        Add pqr or lensfit info
        """

        trials=fitter.get_trials()
        g=trials[:,2:2+2]

        g_prior=self['prior'].g_prior

        ls=ngmix.lensfit.LensfitSensitivity(g, g_prior)
        res['g_sens'] = ls.get_g_sens()
        res['nuse'] = ls.get_nuse()

        pqrobj=ngmix.pqr.PQR(g, g_prior)
        P,Q,R = pqrobj.get_pqr()
        res['P']=P
        res['Q']=Q
        res['R']=R



    def _make_struct(self):
        """
        make the output structure
        """

        np=ngmix.gmix.get_model_npars(self['fit_model'])

        dt=[('number','i4'),
            ('processed','i1'),
            ('flags','i4'),
            ('nimage_tot','i4'),
            ('nimage_use','i4'),
            ('time','f8'),

            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),

            ('max_flags','i4'),
            ('pars_max','f8',np),
            ('pars_max_cov','f8',(np,np)),

            ('flux','f8'),
            ('flux_err','f8'),
            ('T','f8'),
            ('T_err','f8'),
            ('T_s2n','f8'),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),

            ('s2n_w','f8'),
            ('chi2per','f8'),
            ('dof','f8'),
           ]

        if self['fitter'] == 'mcmc':
            dt += [('arate','f4'),('tau','f4')]
        elif self['fitter']=='isample':
            dt += [('efficiency','f4'),('neff','f4')]#,('niter','i2')]

        if self['do_shear']:
            dt += [('P', 'f8'),
                   ('Q', 'f8', 2),
                   ('R', 'f8', (2,2)),
                   ('g_sens','f8',2)]


        num=self.index_list.size
        data=zeros(num, dtype=dt)

        #data['psf_flux'] = DEFVAL
        #data['psf_flux_err'] = PDEFVAL

        #data['em_gauss_flux'] = DEFVAL
        #data['em_gauss_flux_err'] = PDEFVAL
        #data['em_gauss_cen'] = DEFVAL


        data['pars'] = DEFVAL
        data['pars_cov'] = PDEFVAL
        data['pars_max'] = DEFVAL
        data['pars_max_cov'] = PDEFVAL
        data['flux'] = DEFVAL
        data['flux_err'] = PDEFVAL
        data['T'] = DEFVAL
        data['T_err'] = PDEFVAL
        data['g'] = DEFVAL
        data['g_cov'] = PDEFVAL

        data['s2n_w'] = DEFVAL
        data['chi2per'] = PDEFVAL

        if self['do_shear']:
            data['P'] = DEFVAL
            data['Q'] = DEFVAL
            data['R'] = DEFVAL
            data['g_sens'] = DEFVAL

        if self['fitter']=='isample':
            data['efficiency'] = DEFVAL
            data['neff']       = DEFVAL
     
        self.data=data

class MedsFitMax(MedsFit):
    def _fit_galaxy_model(self, obs):
        """
        Run through and fit all the models
        """
        import gmix_meds.nfit

        res=self.res
        model=self['fit_model']

        ps2n=self.res['psf_flux_s2n']
        if ps2n < self['min_s2n']:
            print("    psf s/n too low:",ps2n)
            self.res['flags'] |= S2N_TOO_LOW
            return

        print('    fitting',model,'using maxlike')
        max_guesser=self._get_guesser_from_psf()
        max_fitter=self._fit_simple_max(obs,model,max_guesser) 
        fitres=max_fitter.get_result()

        # fool _copy_galaxy_pars by setting the fitters to the same
        # fitter

        res['max_fitter'] = max_fitter
        res['galaxy_fitter'] = max_fitter
        res['galaxy_res'] = fitres
        res['max_res'] = fitres
        res['flags'] = res['galaxy_res']['flags']

        self._print_galaxy_res(max_fitter)


class MedsFitISample(MedsFit):
    def _fit_galaxy_model(self, obs):
        """
        Run through and fit all the models
        """
        import gmix_meds.nfit
        
        assert self['g_prior_during']==True,"g prior must be during"

        ipars=self['isample_pars']

        res=self.res
        model=self['fit_model']

        ps2n=self.res['psf_flux_s2n']
        if ps2n < self['min_s2n']:
            print("    psf s/n too low:",ps2n)
            self.res['flags'] |= S2N_TOO_LOW
            return

        print('    fitting',model,'using maxlike')
        max_guesser=self._get_guesser_from_psf()
        max_fitter=self._fit_simple_max(obs,model,max_guesser) 
        max_res=max_fitter.get_result()

        res['max_fitter'] = max_fitter
        res['max_res'] = max_res

        if max_res['flags'] != 0:
            sampler=None
        else:            
            self._try_replace_cov(max_fitter)
            self._print_galaxy_res(max_fitter)

            use_fitter = max_fitter
            niter=len(ipars['nsample'])
            for i,nsample in enumerate(ipars['nsample']):
                sampler=self._make_sampler(use_fitter)
                if sampler is None:
                    break

                sampler.make_samples(nsample)

                sampler.set_iweights(max_fitter.calc_lnprob)
                sampler.calc_result()

                tres=sampler.get_result()

                print("    eff iter %d: %.2f" % (i,tres['efficiency']))
                use_fitter = sampler

        if sampler is None:
            gres={'flags':ISAMP_BAD_COV}
            res['flags']=ISAMP_BAD_COV
        else:
            self._add_shear_info(sampler)
            gres=sampler.get_result()
            res['flags'] = gres['flags']
            self._print_galaxy_res(sampler)

        res['galaxy_fitter'] = sampler
        res['galaxy_res'] = gres


        if self['make_plots']:
            self._do_gal_plots(res['galaxy_fitter'])


    def _add_shear_info(self, sampler):
        """
        lensfit and pqr

        call calc_result *before* calling this method
        """

        maxres = self.res['max_res']

        # this is the full prior
        prior=self['search_prior']
        g_prior=prior.g_prior

        iweights = sampler.get_iweights()
        samples = sampler.get_samples()
        g_vals=samples[:,2:2+2]

        remove_prior=True

        res=sampler.get_result()

        # keep for later if we want to make plots
        self.weights=iweights

        # we are going to mutate the result dict owned by the sampler
        res['s2n_w'] = maxres['s2n_w']

        ls=ngmix.lensfit.LensfitSensitivity(g_vals,
                                            g_prior,
                                            weights=iweights,
                                            remove_prior=remove_prior)
        g_sens = ls.get_g_sens()
        g_mean = ls.get_g_mean()

        res['g_sens'] = g_sens
        res['nuse'] = ls.get_nuse()

        # not able to use extra weights yet
        '''
        pqrobj=ngmix.pqr.PQR(g, g_prior,
                             shear_expand=self.shear_expand,
                             remove_prior=remove_prior)


        P,Q,R = pqrobj.get_pqr()
        res['P']=P
        res['Q']=Q
        res['R']=R
        '''

    def _make_sampler(self, fitter):
        from ngmix.fitting import GCovSampler, GCovSamplerT
        from numpy.linalg import LinAlgError

        ipars=self['isample_pars']

        res=fitter.get_result()
        icov = res['pars_cov']*ipars['ifactor']**2

        try:
            if ipars['sampler']=='T':
                sampler=GCovSamplerT(res['pars'],
                                     icov,
                                     ipars['df'],
                                     min_err=ipars['min_err'],
                                     max_err=ipars['max_err'])
            else:
                sampler=GCovSampler(res['pars'],
                                    icov,
                                    min_err=ipars['min_err'],
                                    max_err=ipars['max_err'])
        except LinAlgError:
            print("        bad cov")
            sampler=None

        return sampler


    def _try_replace_cov(self, fitter):
        """
        the lm cov sucks, try to replace it
        """

        # reference to res
        res=fitter.get_result()

        print("        replacing cov")
        max_pars=self['max_pars']
        fitter.calc_cov(max_pars['cov_h'], max_pars['cov_m'])

        if res['flags'] != 0:
            print("        replacement failed")
            res['flags']=0


_em2_fguess=array([0.5793612389470884,1.621860687127999])
_em2_pguess=array([0.596510042804182,0.4034898268889178])

_em3_pguess = array([0.596510042804182,0.4034898268889178,1.303069003078001e-07])
_em3_fguess = array([0.5793612389470884,1.621860687127999,7.019347162356363],dtype='f8')


_stat_names=['s2n_w',
             'chi2per',
             'dof']


def get_model_names(model):
    names=['flags',
           'pars',
           'pars_cov',
           'flux',
           'flux_err',
           'flux_cov',
           'g',
           'g_cov',
           'g_sens',
           'e',
           'e_cov',
           'e_sens',
           'P',
           'Q',
           'R',
           'iter',
           'tries',
           'arate',
           'tau']
    names += _stat_names

    ndict={}
    for n in names:
        ndict[n] = '%s_%s' % (model,n)

    return ndict

def get_em_ngauss(name):
    ngauss=int( name[2:] )
    return ngauss

def srandu(num=None):
    """
    Generate random numbers in the symmetric distribution [-1,1]
    """
    return 2*(numpy.random.random(num)-0.5)

def get_positive_guess(val, n, width=0.01):
    """
    Get guess, making sure positive
    """

    if val <= 0.0:
        print("val <= 0: %s" % val)
        print("using arbitrary value of 0.1")
        val=0.1
        #raise GMixRangeError("val <= 0: %s" % val)

    vals=zeros(n)-9999.0
    while True:
        w,=numpy.where(vals <= 0)
        if w.size == 0:
            break
        else:
            vals[w] = val*(1.0 + width*srandu(w.size))

    return vals

def get_shape_guess(g1, g2, n, width=0.01):
    """
    Get guess, making sure in range
    """

    gtot = sqrt(g1**2 + g2**2)
    if gtot > 0.98:
        g1 = g1*0.9
        g2 = g2*0.9

    guess=zeros( (n, 2) )
    shape=ngmix.Shape(g1, g2)

    for i in xrange(n):

        while True:
            try:
                g1_offset = width*srandu()
                g2_offset = width*srandu()
                shape_new=shape.copy()
                shape_new.shear(g1_offset, g2_offset)
                break
            except GMixRangeError:
                pass

        guess[i,0] = shape_new.g1
        guess[i,1] = shape_new.g2

    return guess

def plot_autocorr(trials, window=100, show=False, **kw):
    import biggles
    import emcee

    arr=biggles.FramedArray(trials.shape[1], 1)
    arr.uniform_limits=True

    func=emcee.autocorr.function(trials)
    tau2 = emcee.autocorr.integrated_time(trials, window=window)

    xvals=numpy.arange(func.shape[0])
    zc=biggles.Curve( [0,func.shape[0]-1],[0,0] )

    for i in xrange(trials.shape[1]):
        pts=biggles.Curve(xvals,func[:,i],color='blue')
        
        lab=biggles.PlotLabel(0.9,0.9,
                              r'$%s tau\times 2: %s$' % (i,tau2[i]),
                              halign='right')
        arr[i,0].add(pts,zc,lab)

    if show:
        arr.show(**kw)

    return arr


