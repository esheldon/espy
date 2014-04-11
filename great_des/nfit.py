from __future__ import print_function
import os
from sys import stderr,stdout
import time
import numpy
from numpy import array, sqrt, zeros

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

NO_ATTEMPT=2**30

#PSF_S2N=1.e6
PSF_OFFSET_MAX=0.25
PSF_TOL=1.0e-5
EM_MAX_TRY=3
EM_MAX_ITER=100

_CHECKPOINTS_DEFAULT_MINUTES=[10,30,60,90]

class MedsFit(object):
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

        self.meds_file=meds_file
        self.truth_file=truth_file
        self.psf_file=psf_file

        self._load_meds()
        self._load_truth()
        self._load_psf_fobj()

        self.conf={}
        self.conf.update(keys)

        self.fitter = keys['fitter']
        self.fit_model = keys['fit_model']

        self.trim_image=keys.get('trim_image',False)
        self.nwalkers=keys.get('nwalkers',80)
        self.burnin=keys.get('burnin',400)
        self.nstep=keys.get('nstep',800)
        self.mca_a=keys.get('mca_a',2.0)

        self.do_pqr=keys.get('do_pqr',True)

        self.shear_expand=keys.get('shear_expand',None)

        self._unpack_priors()

        self._setup_checkpoints()

        self.obj_range=keys.get('obj_range',None)
        self._set_index_list()

        self.psf_model=keys['psf_model']
        self.psf_ngauss=get_em_ngauss(self.psf_model)

        self.psf_fwhm_guess=keys.get('psf_fwhm_guess')
        self.psf_sigma_guess=self.psf_fwhm_guess/2.3548200450309493

        self.region=keys.get('region','seg_and_sky')
        self.max_box_size=keys.get('max_box_size',2048)

        self.make_plots=keys.get('make_plots',False)
        self.prompt=keys.get('prompt',True)

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
        import ngmix
        from ngmix import GMixMaxIterEM

        self.fitting_galaxy=False

        # first 1 gaussian
        jacobian=self._get_jacobian(self.psf_cen_guess_pix)

        fitter1=self._fit_em_1gauss(self.psf_image,
                                    jacobian,
                                    self.psf_sigma_guess)

        if fitter1 is None:
            self.res['flags'] = PSF_FIT_FAILURE
            print("psf em1 failure at object",self.mindex)
            return


        psf_gmix1=fitter1.get_gmix()
        self.psf_gmix1=psf_gmix1
        self.res['psf_gmix_em1']=psf_gmix1

        fitter=self._fit_em_ngauss(self.psf_image,
                                   jacobian,
                                   psf_gmix1,
                                   self.psf_ngauss)

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

            if self.make_plots:
                self._compare_psf(fitter)

    def _fit_galaxy(self):
        """
        Fit psf flux and other models
        """

        self.fitting_galaxy=True

        print('    fitting gal em 1gauss')

        self._fit_galaxy_em()
        self._fit_galaxy_psf_flux()

        if self.trim_image:
            self._make_trimmed_image()
        self._fit_galaxy_model()

    def _fit_galaxy_em(self):
        """

        Fit a single gaussian with em to find a decent center.  We don't get a
        flux out of that, but we can get flux using the _fit_flux routine

        """

        print("    fitting galaxy with em")
        # first the structural fit
        sigma_guess = sqrt( self.res['psf_gmix'].get_T()/2.0 )

        jacobian = self._get_jacobian(self.gal_cen_guess_pix)

        print('      sigma guess:',sigma_guess)
        print('      cen guess:  ',self.gal_cen_guess_pix)
        fitter=self._fit_em_1gauss(self.image,
                                   jacobian, 
                                   sigma_guess)

        em_gmix = fitter.get_gmix()
        print("      em gmix:",em_gmix)

        # now get a flux
        print('    fitting robust gauss flux')
        flux, flux_err = self._fit_flux(self.image,
                                        self.wt,
                                        jacobian,
                                        em_gmix)

        self.res['em_gauss_flux'] = flux
        self.res['em_gauss_flux_err'] = flux_err
        self.res['em_gauss_cen'] = em_gmix.get_cen()
        self.res['em_gmix'] = em_gmix
        self.res['jacobian'] = jacobian


    def _fit_galaxy_psf_flux(self):
        """
        Get flux fitting the psf model

        Call this after fitting galaxy with em 1
        """

        print("    fitting galaxy with psf")
        res=self.res

        gmix=res['psf_gmix'].copy()

        # jacobian is not on the best center from
        # the robust single gaussian fit
        flux, flux_err = self._fit_flux(self.image,
                                        self.wt,
                                        res['jacobian'],
                                        gmix)
        res['psf_flux'] = flux
        res['psf_flux_err'] = flux_err

    def _fit_flux(self, image, weight_image, jacob, gmix):
        """
        Fit the flux from a fixed gmix model.  This is linear and always
        succeeds.
        """
        import ngmix

        fitter=ngmix.fitting.PSFFluxFitter(image,
                                           weight_image,
                                           jacob,
                                           gmix)
        fitter.go()
        res=fitter.get_result()

        flux=res['flux']
        flux_err=res['flux_err']
        mess='         %s +/- %s' % (flux,flux_err)
        print(mess)

        return flux, flux_err

    def _fit_galaxy_model(self):
        """
        Run through and fit all the models
        """

        model=self.fit_model
        fitter_type=self.fitter
        print('    fitting',model,"using",fitter_type)

        if model=='sersic':
            if fitter_type=='mcmc':
                if self.conf['iter']:
                    fitter=self._fit_sersic_mcmc_iter()
                else:
                    fitter=self._fit_sersic_mcmc()
            elif fitter_type=='lm':
                raise RuntimeError("don't use lm, it is shit")
                fitter=self._fit_sersic_lm()
            else:
                raise ValueError("bad fitter type: '%s'" % fitter_type)
        else:
            raise ValueError("bad model: '%s'" % model)

        self.res['galaxy_fitter'] = fitter
        self.res['galaxy_res'] = fitter.get_result()
        self.res['flags'] |= self.res['galaxy_res']['flags']

        self._print_galaxy_res()

        if fitter_type=='mcmc' and self.make_plots:
            self._do_gal_plots(self.res['galaxy_fitter'])


    def _fit_sersic_mcmc(self):
        """
        Fit a sersic model, taking guesses from our previous em fits
        """
        import ngmix
        from ngmix.fitting import print_pars

        res=self.res

        image,wt,jacobian=self._get_image_data()

        full_guess=self._get_guess_sersic_maxlike()

        fitter=ngmix.fitting.MCMCSersic(image,
                                        wt,
                                        jacobian,
                                        psf=res['psf_gmix'],

                                        nwalkers=self.nwalkers,
                                        burnin=self.burnin,
                                        nstep=self.nstep,
                                        mca_a=self.mca_a,

                                        full_guess=full_guess,

                                        shear_expand=self.shear_expand,

                                        cen_prior=self.cen_prior,
                                        T_prior=self.T_prior,
                                        counts_prior=self.counts_prior,
                                        g_prior=self.g_prior,
                                        n_prior=self.n_prior,
                                        do_pqr=self.do_pqr)
        fitter.go()

        return fitter

    def _fit_sersic_mcmc_iter(self):
        """
        Fit a sersic model, taking guesses from our previous em fits
        """
        import ngmix
        from ngmix.fitting import print_pars

        res=self.res

        image,wt,jacobian=self._get_image_data()

        ntry=self.conf['ntry']
        for i in xrange(ntry):
            print("    try: %s/%s" % (i+1,ntry))
            if i==0:
                full_guess=self._get_guess_sersic(self.nwalkers)
            else:
                best_pars=fitter.best_pars
                print_pars(best_pars,front="    first pars: ")
                full_guess=self._get_guess_sersic_frompars(best_pars)

            fitter=ngmix.fitting.MCMCSersic(image,
                                            wt,
                                            jacobian,
                                            psf=res['psf_gmix'],

                                            nwalkers=self.nwalkers,
                                            burnin=self.burnin,
                                            nstep=self.nstep,
                                            mca_a=self.mca_a,

                                            full_guess=full_guess,

                                            shear_expand=self.shear_expand,

                                            cen_prior=self.cen_prior,
                                            T_prior=self.T_prior,
                                            counts_prior=self.counts_prior,
                                            g_prior=self.g_prior,
                                            n_prior=self.n_prior,
                                            do_pqr=self.do_pqr)
            fitter.go()

        return fitter


    def _get_guess_sersic_maxlike(self):
        """
        fit with lm and use as guess.  If it fails, fall back
        to the guess based on em
        """
        from ngmix.fitting import print_pars

        fitter=self._fit_sersic_lm()
        lmres=fitter.get_result()

        if lmres['flags'] != 0:
            print("    LM failed, using standard guess instead")
            full_guess=self._get_guess_sersic(self.nwalkers)
        else:
            pars=lmres['pars']
            #print_pars(pars,front="    LM pars: ")
            full_guess=self._get_guess_sersic_frompars(pars)

        return full_guess

    def _fit_sersic_lm(self):
        """
        Fit a sersic model, taking guesses from our previous em fits
        """
        import ngmix
        from ngmix.fitting import print_pars

        image,wt,jacobian=self._get_image_data()

        res=self.res

        chi2per_min=1.e9

        random_n=False
        start_pars=None
        ntry=self.conf['ntry_lm']
        for i in xrange(ntry):

            guess=self._get_guess_sersic(1, random_n=random_n, pars=start_pars)
            guess=guess[0,:]

            if i==0:
                print_pars(guess,front="    LM start: ")

            tfitter=ngmix.fitting.LMSersic(image,
                                           wt,
                                           jacobian,
                                           guess,

                                           psf=res['psf_gmix'],

                                           lm_pars=self.conf['lm_pars'],

                                           cen_prior=self.cen_prior,
                                           T_prior=self.T_prior,
                                           counts_prior=self.counts_prior,
                                           g_prior=self.g_prior,
                                           n_prior=self.n_prior)
            tfitter.go()

            tres=tfitter.get_result()
            if tres['flags']==0:
                chi2per=tres['chi2per']
                if chi2per < chi2per_min:
                    print_pars(tres['pars'],front="    best LM pars %02s: " % (i+1))
                    chi2per_min=chi2per
                    fitter=tfitter

                start_pars=tres['pars']
            else:
                print("    lm fail:",tres['flags'])
                start_pars=None
                if i==0:
                    fitter=tfitter

            random_n=True

        return fitter


    def _get_guess_sersic(self,
                          num,
                          pars=None,
                          random_n=True,
                          widths=[0.05, 0.05, 0.05, 0.05]):
        """
        Guess based on the em fit, with n drawn from prior
        """

        res=self.res

        if pars is None:
            gmix = res['em_gmix']
            cen1,cen2=gmix.get_cen()
            g1,g2,T = gmix.get_g1g2T()
            flux = res['em_gauss_flux']
        else:
            cen1=pars[0]
            cen2=pars[1]
            g1=pars[2]
            g2=pars[3]
            T=pars[4]
            flux=pars[5]

        shapes=get_shape_guess(g1, g2, num, width=widths[1])

        n_prior=self.n_prior

        if random_n:
            nvals = self.n_prior.sample(num)
        else:
            nvals=zeros(num)
            nvals[:] = 0.5*(n_prior.minval+n_prior.maxval)

        guess=zeros( (num, 7) )

        guess[:,0] = cen1 + widths[0]*srandu(num)
        guess[:,1] = cen2 + widths[0]*srandu(num)

        guess[:,2]=shapes[:,0]
        guess[:,3]=shapes[:,1]

        guess[:,4] = get_positive_guess(T,num,width=widths[2])
        guess[:,5] = get_positive_guess(flux,num,width=widths[3])

        guess[:,6] = nvals

        return guess


    def _get_guess_sersic_frompars(self, pars):
        """
        Guess centered on the input pars
        """
        nwalkers = self.nwalkers

        nmin=self.n_prior.minval
        nmax=self.n_prior.maxval

        guess=zeros( (nwalkers, 7) )
        guess[:,0] = pars[0] + 0.01*srandu(nwalkers)
        guess[:,1] = pars[1] + 0.01*srandu(nwalkers)

        shapes=get_shape_guess(pars[2],pars[3], nwalkers, width=0.01)
        guess[:,2] = shapes[:,0]
        guess[:,3] = shapes[:,1]

        guess[:,4] = get_positive_guess(pars[4],nwalkers,width=0.01)
        guess[:,5] = get_positive_guess(pars[5],nwalkers,width=0.01)

        nleft=nwalkers
        ngood=0
        while nleft > 0:
            vals = pars[6]*(1.0 + 0.01*srandu(nleft))
            w,=numpy.where( (vals > nmin) & (vals < nmax) )
            nkeep=w.size
            if nkeep > 0:
                guess[ngood:ngood+nkeep,6] = vals[w]
                nleft -= w.size
                ngood += w.size

        return guess




    def _load_image_and_weight(self):
        """
        Load the image, weight map
        """

        dindex=self.dindex
        mindex=self.mindex
        
        self.image = self.meds.get_cutout(mindex,0).astype('f8')
        if self.conf['noisefree']:
            self.wt = 0*self.image + (1.0/self.conf['skynoise']**2)
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
        import ngmix
        jdict = self.meds.get_jacobian(self.mindex,0)

        jacobian=ngmix.Jacobian(cen[0],
                                cen[1],
                                jdict['dudrow'],
                                jdict['dudcol'],
                                jdict['dvdrow'],
                                jdict['dvdcol'])

        return jacobian


    def _fit_em_1gauss(self, im, jacobian, sigma_guess):
        """
        Just run the fitter

        cen is the global cen not within jacobian
        """
        return self._fit_with_em(im, jacobian, sigma_guess, 1)

    def _fit_em_ngauss(self, im, jacobian, gmix1, ngauss):
        """
        Start with fit from using 1 gauss, which must be entered

        cen is the global cen not within jacobian
        """

        sigma_guess = sqrt( gmix1.get_T()/2. )
        fitter=self._fit_with_em(im, jacobian, sigma_guess, ngauss)

        return fitter



    def _fit_with_em(self, im, jacobian, sigma_guess, ngauss):
        """
        Fit the image using EM
        """
        import ngmix
        from ngmix.gexceptions import GMixMaxIterEM, GMixRangeError
        conf=self.conf

        im_with_sky, sky = ngmix.em.prep_image(im)

        ntry,maxiter,tol = self._get_em_pars()
        for i in xrange(ntry):
            guess = self._get_em_guess(sigma_guess, ngauss)
            #print("em guess:")
            #print(guess)
            try:
                fitter=self._do_fit_em_with_full_guess(im_with_sky,
                                                       sky,
                                                       guess,
                                                       jacobian)
                tres=fitter.get_result()
                #print("em numiter:",tres['numiter'])
                break
            except GMixMaxIterEM:
                fitter=None
            except GMixRangeError:
                fitter=None

        return fitter

    def _do_fit_em_with_full_guess(self,
                                   image,
                                   sky,
                                   guess,
                                   jacobian):
        import ngmix

        ntry,maxiter,tol = self._get_em_pars()

        fitter=ngmix.em.GMixEM(image, jacobian=jacobian)
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
        import ngmix

        sigma2 = sigma**2
        pars=array( [1.0 + 0.1*srandu(),
                     0.2*srandu(),
                     0.2*srandu(), 
                     sigma2*(1.0 + 0.5*srandu()),
                     0.2*sigma2*srandu(),
                     sigma2*(1.0 + 0.5*srandu())] )

        return ngmix.gmix.GMix(pars=pars)

    def _get_em_guess_2gauss(self, sigma):
        import ngmix

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
        import ngmix

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
        import ngmix

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
        conf=self.conf
        if self.fitting_galaxy:
            return conf['gal_em_ntry'], conf['gal_em_maxiter'], conf['gal_em_tol']
        else:
            return conf['psf_em_ntry'], conf['psf_em_maxiter'], conf['psf_em_tol']

    def _make_trimmed_image(self):
        res=self.res

        cen_jacob=array( res['jacobian'].get_cen() )

        em_gmix=res['em_gmix']

        T=em_gmix.get_T()
        cen=array(em_gmix.get_cen())/self.pixel_scale

        cen = cen + cen_jacob

        sigma=sqrt(T/2.)/self.pixel_scale
        #print("cen:",cen)
        #print("from gauss:",res['em_gmix'])
        #print("got T:",T,"sigma (pixels):",sigma)

        radius = 5.0*sigma

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

    def _get_image_data(self):
        """
        Get the appropriate image data for the current object.
        """

        if self.trim_image:
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

        model=self.conf['psf_model']

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

    def _compare_gal(self, fitter):
        """
        compare psf image to best fit model
        """
        import images

        model=self.fit_model

        gmix = fitter.get_gmix()

        image, wt, jacobian=self._get_image_data()

        res=self.res
        psf_gmix = res['psf_gmix']
        gmix_conv = gmix.convolve(psf_gmix)

        model_image = gmix_conv.make_image(image.shape,
                                           jacobian=jacobian)

        plt=images.compare_images(image,
                                  model_image,
                                  label1='galaxy',
                                  label2=model,
                                  show=False)
        pname='gal-resid-%06d-%s.png' % (self.mindex,model)

        resid_std = (image-model_image).std()
        print("    residual std:",resid_std)
        print("          ",pname)
        plt.write_img(1400,800,pname)


    def _make_trials_plot(self, fitter):
        """
        Plot the trials
        """

        model=self.fit_model
        width,height=800,800

        tup=fitter.make_plots(title=model)

        if isinstance(tup, tuple):
            p,wp = tup
            wtrials_pname='wtrials-%06d-%s.png' % (self.mindex,model)
            print("          ",wtrials_pname)
            wp.write_img(width,height,wtrials_pname)
        else:
            p = tup

        trials_pname='trials-%06d-%s.png' % (self.mindex,model)
        print("          ",trials_pname)
        p.write_img(width,height,trials_pname)

    def _print_galaxy_res(self):
        from ngmix.fitting import print_pars
        res=self.res['galaxy_res']

        print_pars(res['pars'], front="    pars: ")
        print_pars(res['pars_err'], front="    err:  ")

        #self._print_galaxy_cen(res)
        #self._print_galaxy_shape(res)
        #self._print_galaxy_T(res)
        #self._print_galaxy_flux(res)

        #if self.fit_model=='sersic':
        #    self._print_galaxy_n(res)
        
        if 'arate' in res:
            print('        arate:',res['arate'])

    def _print_galaxy_n(self, res):
        """
        print the center
        """
        n=res['pars'][6]
        n_err=res['pars_err'][6]
        mess='        n: %.4g +/- %.4g'
        mess = mess % (n, n_err)
        print(mess)


    def _print_galaxy_cen(self, res):
        """
        print the center
        """
        pars=res['pars']
        perr=res['pars_err']
        mess='        cen: %.4g +/- %.4g %.4g +/- %.4g'
        mess = mess % (pars[0],perr[0],pars[1],perr[1])
        print(mess)

    def _print_galaxy_shape(self, res):
        """
        print shape info
        """
        g1=res['pars'][2]
        g1err=sqrt(res['pars_cov'][2,2])
        g2=res['pars'][3]
        g2err=sqrt(res['pars_cov'][3,3])

        mess='        g1: %.4g +/- %.4g g2: %.4g +/- %.4g'
        mess = mess % (g1,g1err,g2,g2err)
        print(mess)

    def _print_galaxy_flux(self, res):
        """
        print in a nice format
        """

        flux = res['pars'][5:].sum()
        flux_err = sqrt( res['pars_cov'][5:, 5:].sum() )
        s2n=flux/flux_err

        print('        flux: %s +/- %s Fs2n: %s' % (flux,flux_err,s2n))

    def _print_galaxy_T(self, res):
        """
        print T, Terr, Ts2n and sigma
        """

        T = res['pars'][4]
        Terr = sqrt( res['pars_cov'][4,4] )

        if Terr > 0:
            Ts2n=T/Terr
        else:
            Ts2n=-9999.0
        if T > 0:
            sigma=sqrt(T/2.)
        else:
            sigma=-9999.0

        tup=(T,Terr,Ts2n,sigma)
        print('        T: %s +/- %s Ts2n: %s sigma: %s' % tup)


    def _copy_to_output(self):
        """
        Copy the galaxy fits
        """

        dindex=self.dindex
        res=self.res
        data=self.data
        conf=self.conf

        # overall flags, model flags copied below
        data['flags'][dindex] = res['flags']

        if 'psf_gmix' in res:
            self._copy_galaxy_pars()

    def _copy_galaxy_pars(self):
        """
        Copy from the result dict to the output array
        """

        dindex=self.dindex
        conf=self.conf
        allres = self.res

        data=self.data

        # em fit to convolved galaxy
        if 'em_gauss_flux' in allres:
            data['em_gauss_flux'][dindex] = allres['em_gauss_flux']
            data['em_gauss_flux_err'][dindex] = allres['em_gauss_flux_err']
            data['em_gauss_cen'][dindex] = allres['em_gauss_cen']
        else:
            print("        em gauss not found, not copying pars")
            return

        if 'psf_flux' in allres:
            data['psf_flux'][dindex] = allres['psf_flux']
            data['psf_flux_err'][dindex] = allres['psf_flux_err']
        else:
            print("        em gauss not found, not copying pars")
            return

        # allres flags is or'ed with the galaxy fitting flags
        if allres['flags'] != 0:
            print("        Not copying pars due to failure")
            return

        res=allres['galaxy_res']

        pars=res['pars']
        pars_cov=res['pars_cov']

        flux=pars[5]
        flux_err=sqrt(pars_cov[5, 5])

        data['pars'][dindex,:] = pars
        data['pars_cov'][dindex,:,:] = pars_cov

        data['flux'][dindex] = flux
        data['flux_err'][dindex] = flux_err

        data['g'][dindex,:] = res['g']
        data['g_cov'][dindex,:,:] = res['g_cov']

        if 'arate' in res:
            data['arate'][dindex] = res['arate']
            data['tau'][dindex] = res['tau']

        for sn in _stat_names:
            if sn in res:
                data[sn][dindex] = res[sn]

        if conf['do_pqr']:
            data['P'][dindex] = res['P']
            data['Q'][dindex,:] = res['Q']
            data['R'][dindex,:,:] = res['R']

    def _load_meds(self):
        """
        Load all listed meds files
        """

        print(self.meds_file)
        self.meds = meds.MEDS(self.meds_file)
        self.meds_meta=self.meds.get_meta()
        self.nobj_tot = self.meds.size

        self.mindex=0
        jtmp=self._get_jacobian([0.0, 0.0])

        self.pixel_scale = jtmp.get_scale()
        print("pixel scale:",self.pixel_scale)

    def _load_truth(self):
        """
        load the truth file for getting the psf index
        """
        import fitsio
        print(self.truth_file)
        self.truth = fitsio.read(self.truth_file,lower=True)
        
    def _load_psf_fobj(self):
        """
        Load the psf file as a FITS object
        """
        import fitsio
        print(self.psf_file)
        self.psf_fobj = fitsio.FITS(self.psf_file)

    def _set_index_list(self):
        """
        set the list of indices to be processed
        """
        if self.obj_range is None:
            start=0
            end=self.nobj_tot-1
        else:
            start=self.obj_range[0]
            end=self.obj_range[1]

        self.index_list = numpy.arange(start,end+1)


    def _setup_checkpoints(self):
        """
        Set up the checkpoint times in minutes and data
        """
        self.checkpoints = self.conf.get('checkpoints',_CHECKPOINTS_DEFAULT_MINUTES)
        self.n_checkpoint    = len(self.checkpoints)
        self.checkpointed    = [0]*self.n_checkpoint
        self.checkpoint_file = self.conf.get('checkpoint_file',None)

        self._set_checkpoint_data()

        if self.checkpoint_file is not None:
            self.do_checkpoint=True
        else:
            self.do_checkpoint=False

    def _set_checkpoint_data(self):
        """
        See if checkpoint data was sent
        """
        self._checkpoint_data=self.conf.get('checkpoint_data',None)
        if self._checkpoint_data is not None:
            self.data=self._checkpoint_data['data']

    def _try_checkpoint(self, tm):
        """
        Checkpoint at certain intervals.  
        Potentially modified self.checkpointed
        """

        should_checkpoint, icheck = self._should_checkpoint(tm)

        if should_checkpoint:
            self._write_checkpoint(tm)
            self.checkpointed[icheck]=1

    def _should_checkpoint(self, tm):
        """
        Should we write a checkpoint file?
        """

        should_checkpoint=False
        icheck=-1

        if self.do_checkpoint:
            tm_minutes=tm/60

            for i in xrange(self.n_checkpoint):

                checkpoint=self.checkpoints[i]
                checkpointed=self.checkpointed[i]

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
        print(self.checkpoint_file)

        with fitsio.FITS(self.checkpoint_file,'rw',clobber=True) as fobj:
            fobj.write(self.data, extname="model_fits")

    def _unpack_priors(self):
        conf=self.conf

        self.T_prior=conf['T_prior']
        self.counts_prior=conf['counts_prior']
        self.g_prior=conf['g_prior']

        self.n_prior = conf.get('n_prior',None)

        # in arcsec (or units of jacobian)
        self.cen_prior=conf.get("cen_prior",None)


    def _make_struct(self):
        """
        make the output structure
        """
        import ngmix

        np=ngmix.gmix.get_model_npars(self.fit_model)

        dt=[('number','i4'),
            ('processed','i1'),
            ('flags','i4'),
            ('nimage_tot','i4'),
            ('nimage_use','i4'),
            ('time','f8'),

            ('em_gauss_flux','f8'),
            ('em_gauss_flux_err','f8'),
            ('em_gauss_cen','f8',2),

            ('psf_flux',    'f8'),
            ('psf_flux_err','f8'),

            ('fit_flags','i4'),
            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),

            ('flux','f8'),
            ('flux_err','f8'),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),

            ('s2n_w','f8'),
            ('chi2per','f8'),
            ('dof','f8'),
           ]

        if self.fitter == 'mcmc':
            dt += [('arate','f8'), ('tau','f8')]

        if self.do_pqr:
            dt += [('P', 'f8'),
                   ('Q', 'f8', 2),
                   ('R', 'f8', (2,2))]


        num=self.index_list.size
        data=zeros(num, dtype=dt)

        data['psf_flux'] = DEFVAL
        data['psf_flux_err'] = PDEFVAL

        data['em_gauss_flux'] = DEFVAL
        data['em_gauss_flux_err'] = PDEFVAL
        data['em_gauss_cen'] = DEFVAL

        data['fit_flags'] = NO_ATTEMPT

        data['pars'] = DEFVAL
        data['pars_cov'] = PDEFVAL
        data['flux'] = DEFVAL
        data['flux_err'] = PDEFVAL
        data['g'] = DEFVAL
        data['g_cov'] = PDEFVAL

        data['s2n_w'] = DEFVAL
        data['chi2per'] = PDEFVAL

        if self.fitter=='mcmc':
            data['tau'] = PDEFVAL

        if self.do_pqr:
            data['P'] = DEFVAL
            data['Q'] = DEFVAL
            data['R'] = DEFVAL

     
        self.data=data

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
    from ngmix.gexceptions import GMixRangeError

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
    import ngmix
    from ngmix.gexceptions import GMixRangeError

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


