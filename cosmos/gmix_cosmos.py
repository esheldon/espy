from sys import stderr
import time
import numpy
import galsim

from . import files

NO_ATTEMPT=2**30
DEFVAL=-9999
PDEFVAL=9999
BIG_DEFVAL=-9.999e9
BIG_PDEFVAL=9.999e9


def fit_cosmos(models, output_file, **keys):
    """
    Do some fits and write an output file.
    """
    import fitsio
    fitter=CosmosFitter(models, **keys)
    fitter.do_fits()

    data=fitter.get_data()

    print >>stderr,'writing:',output_file
    with fitsio.FITS(output_file, 'rw', clobber=True) as fobj:
        fobj.write(data, extension="model_fits")

class CosmosFitter(object):
    """
    Fit models to the cosmos galaxies released as part of great3.
    
    Convolve with a gaussian PSF, and add a little noise.
    """
    def __init__(self,
                 models,
                 obj_range=None,
                 psf_fwhm=0.3,
                 pixel_scale=0.1,
                 sky_sigma0=0.005, # before whitening, ~0.01 after
                 rotate=True,
                 make_plots=False,
                 plot_base=None, # not yet used
                 nwalkers=80,
                 burnin=1000,
                 nstep=200,
                 mca_a=3.0,
                 min_arate=0.25,
                 ntry=1):
        self._models      = models
        self._psf_fwhm    = float(psf_fwhm)    # arcsec
        self._sky_sigma0  = float(sky_sigma0)
        self._rotate      = rotate
        self._pixel_scale = float(pixel_scale)
        self._pixel_area  = self._pixel_scale**2
        self._make_plots  = make_plots
        self._plot_dir    = plot_base

        self._nwalkers=nwalkers
        self._burnin=burnin
        self._nstep=nstep
        self._mca_a=mca_a
        self._min_arate=min_arate
        self._ntry=ntry

        self.read_data()
        self.set_index_list(obj_range)

        self.setup_rand()

        self.make_pixel()
        self.make_psf()

        self.make_struct()

    def get_data(self):
        """
        Get the fit data
        """
        return self._data

    def set_index_list(self, obj_range):
        """
        Set the range of objects to process
        """
        index_list=numpy.arange(self._nobj)
        if obj_range is not None:
            self._index_list=index_list[obj_range[0]:obj_range[1]+1]
        else:
            self._index_list=index_list

    def do_fits(self):
        """
        Perform all the fits
        """

        self.fit_psf()
        self.fit_galaxies()

    def make_pixel(self):
        """
        Make a pixel object for the convolutions
        """
        self._pix_obj = galsim.Pixel(self._pixel_scale)

    def make_psf(self, show=False):
        """
        Make the psf image
        """
        psf_obj_prepix = galsim.Gaussian(fwhm=self._psf_fwhm,
                                         flux=1.0)
        self._psf_obj = galsim.Convolve([psf_obj_prepix, self._pix_obj])

        self._psf_image_obj = self._psf_obj.draw(dx=self._pixel_scale)
        self._psf_image=self._psf_image_obj.array.astype('f8')

        self._psf_cen = [ 0.5*self._psf_image.shape[0] ]*2

        if show:
            import images
            images.multiview(self.psf_image)

    def fit_psf(self):
        """
        Fit the PSF image
        """
        import ngmix

        sigma = self._psf_fwhm/2.35
        Tguess = 2*sigma**2

        j=ngmix.Jacobian(self._psf_cen[0],
                         self._psf_cen[1],
                         self._pixel_scale,
                         0.,
                         0.,
                         self._pixel_scale)

        imsky,sky=ngmix.em.prep_image(self._psf_image)
        em=ngmix.em.GMixEM(imsky, jacobian=j)

        gmix_guess=ngmix.gmix.GMixModel([0.0, 0.0, 0.0, 0.0, Tguess, 1.0],'gauss')
        print >>stderr,'psf guess:'
        print >>stderr,gmix_guess

        em.go(gmix_guess, sky)

        self._psf_gmix = em.get_gmix()

        print >>stderr,'psf fit:'
        print >>stderr,self._psf_gmix
        print >>stderr,""

    def fit_galaxies(self):
        """
        Convolve the real galaxy images with a psf to make it a
        simple gaussian
        """

        tm0=time.time()

        self.set_priors()
        last=self._index_list[-1]
        ndo=self._index_list.size

        for dindex in xrange(ndo):

            rindex=self._index_list[dindex]
            self.rindex=rindex

            print >>stderr,'%s:%s' % (rindex, last)
            gal_im_obj,wt = self.make_galaxy_image(dindex)
            self.fit_galaxy_models(gal_im_obj,wt,dindex)

        tm=time.time()-tm0
        print >>stderr,'time per galaxy:',tm/ndo

    def fit_galaxy_models(self, gal_im_obj, wt, dindex):
        """
        Fit all the models
        """

        im=gal_im_obj.array.astype('f8')
        print >>stderr," ",im.shape
        counts_guess,T_guess0,row,col=self.get_guesses(gal_im_obj,im,wt)
        j=self.get_jacobian(row,col)

        print >>stderr,'  T_guess0:     ',T_guess0
        print >>stderr,'  counts_guess: ',counts_guess

        for model in self._models:
            print >>stderr,'    model:',model
            if model=='bdc':
                res=self.fit_bdc(im,wt,j,model,counts_guess,T_guess0)
            elif model=='bdf':
                res=self.fit_bdf(im,wt,j,model,counts_guess,T_guess0)
            else:
                res=self.fit_simple_model(im,wt,j,model,counts_guess,T_guess0)

            self.copy_result(dindex, res)

            self.print_some_stats(res)

    def get_guesses(self, gal_im_obj, im, wt):
        """
        Try adaptive moments, otherwise fall back on simple
        guesses
        """
        import ngmix
        try:
            moms=gal_im_obj.FindAdaptiveMom()
            counts_guess = moms.moments_amp * self._pixel_area
            # extra factor of 2 because this is the weight sigma
            T_guess0     = 2 * 2*moms.moments_sigma**2 * self._pixel_area
            row_guess    = moms.moments_centroid.y -1
            col_guess    = moms.moments_centroid.x -1
        except:
            print >>stderr,'Warning: adaptive moments failed, trying EM'
            counts_guess = im.sum() * self._pixel_area
            try:

                row_guess=0.5*im.shape[0] 
                col_guess=0.5*im.shape[1] 
                gm_guess=ngmix.gmix.GMixModel([row_guess,col_guess,0.0, 0.0, 8.0, 1.0],
                                              'gauss')
                imsky,sky=ngmix.em.prep_image(im)
                fitter=ngmix.em.GMixEM(imsky)
                fitter.go(gm_guess, sky, maxiter=1000, tol=1.e-3)

                gm=fitter.get_gmix()
                T_guess0 = gm.get_T() * self._pixel_area
                row_guess,col_guess = gm.get_cen()
            except:
                # will use row_guess from above
                print >>stderr,'Warning: EM also failed!'
                T_guess0=0.25

        return counts_guess, T_guess0, row_guess, col_guess

    def fit_simple_model(self, im, wt, j, model, counts_guess, T_guess0):
        """
        fit a single model
        """

        import ngmix

        T_guess=_guess_Tfac[model]*T_guess0
        print >>stderr,'    T_guess:      ',T_guess

        fitter=ngmix.fitting.MCMCSimple(im, wt, j, model,
                                        psf=self._psf_gmix,
                                        T_guess=T_guess,
                                        counts_guess=counts_guess,
                                        T_prior=self._T_prior,
                                        counts_prior=self._counts_prior,
                                        cen_prior=self._cen_prior,
                                        nwalkers=self._nwalkers,
                                        burnin=self._burnin,
                                        nstep=self._nstep,
                                        mca_a=self._mca_a,
                                        ntry=self._ntry,
                                        min_arate=self._min_arate)
        fitter.go()

        if self._make_plots:
            plt,plt_resid_list=fitter.make_plots(show=True,
                                                 do_residual=True,
                                                 title=model)

        res=fitter.get_result()
        res['T_guess0']=T_guess0
        res['T_guess']=T_guess
        res['counts_guess']=counts_guess
        return res

    def fit_bdf(self, im, wt, j, model, counts_guess, T_guess0):
        """
        fit bdc

        We need to generate a full guess for each walker
        """
        import ngmix
        npars=7

        T_guess=_guess_Tfac[model]*T_guess0
        print >>stderr,'    T_guess:      ',T_guess

        full_guess=self.get_full_bdf_guess(T_guess,counts_guess)

        nwalkers = self._nwalkers
        for i in xrange(self._ntry):
            fitter=ngmix.fitting.MCMCBDF(im, wt, j, model,
                                         psf=self._psf_gmix,
                                         full_guess=full_guess,
                                         T_prior=self._T_prior, # flat
                                         counts_prior=self._counts_prior,
                                         bfrac_prior=self._bfrac_prior,
                                         cen_prior=self._cen_prior,
                                         nwalkers=nwalkers,
                                         burnin=self._burnin,
                                         nstep=self._nstep,
                                         mca_a=self._mca_a,
                                         #ntry=self._ntry,
                                         min_arate=self._min_arate
                                        )
            fitter.go()

            res=fitter.get_result()

            if res['arate'] >= self._min_arate:
                break
            else:
                #print >>stderr,"        arate",res['arate'],"<",self._min_arate
                nwalkers = nwalkers*2
                print >>stderr,'        trying nwalkers:',nwalkers
                trials=fitter.trials
                full_guess=numpy.zeros( (nwalkers,npars) )
                full_guess[:,:] = trials[-nwalkers:, :]

        if self._make_plots:
            fitter.make_plots(show=True, do_residual=True, title=model)

        res['T_guess0']=T_guess0
        res['T_guess']=T_guess
        res['counts_guess']=counts_guess
        return res


    def fit_bdc(self, im, wt, j, model, counts_guess, T_guess0):
        """
        fit bdc

        We need to generate a full guess for each walker
        """
        import ngmix
        npars=8

        T_guess=_guess_Tfac[model]*T_guess0
        print >>stderr,'    T_guess:      ',T_guess

        full_guess=self.get_full_bdc_guess(T_guess,counts_guess)

        nwalkers = self._nwalkers
        for i in xrange(self._ntry):
            fitter=ngmix.fitting.MCMCBDC(im, wt, j, model,
                                         psf=self._psf_gmix,
                                         full_guess=full_guess,
                                         T_b_prior=self._T_prior, # flat
                                         T_d_prior=self._T_prior,
                                         #counts_b_prior=self._counts_prior,
                                         #counts_d_prior=self._counts_prior,
                                         cen_prior=self._cen_prior,
                                         nwalkers=nwalkers,
                                         burnin=self._burnin,
                                         nstep=self._nstep,
                                         mca_a=self._mca_a,
                                         #ntry=self._ntry,
                                         min_arate=self._min_arate
                                        )
            fitter.go()

            res=fitter.get_result()

            if res['arate'] >= self._min_arate:
                break
            else:
                #print >>stderr,"        arate",res['arate'],"<",self._min_arate
                nwalkers = nwalkers*2
                print >>stderr,'        trying nwalkers:',nwalkers
                trials=fitter.trials
                full_guess=numpy.zeros( (nwalkers,npars) )
                full_guess[:,:] = trials[-nwalkers:, :]

        if self._make_plots:
            fitter.make_plots(show=True, do_residual=True, title=model)

        res['T_guess0']=T_guess0
        res['T_guess']=T_guess
        res['counts_guess']=counts_guess
        return res

    def get_full_bdf_guess(self, T_guess, counts_guess):
        """
        Get a full guess with an entry for each walker
        """
        n=self._nwalkers
        npars=7

        guess=numpy.zeros( (n, npars) )
        cen1,cen2=self._cen_prior.sample(n)

        guess[:,0]=cen1
        guess[:,1]=cen2
        guess[:,2]=0.1*srandu(n)
        guess[:,3]=0.1*srandu(n)

        # need to tune these guesses
        guess[:,4] = 2*T_guess*(1.0 + 0.1*srandu(n))

        # bulge flux
        guess[:,5] = 2*counts_guess*(1.0 + 0.1*srandu(n))
        # disk flux
        guess[:,6] = 2*counts_guess*(1.0 + 0.1*srandu(n))

        return guess

    def get_full_bdc_guess(self, T_guess, counts_guess):
        """
        Get a full guess with an entry for each walker
        """
        npars=8
        n=self._nwalkers

        guess=numpy.zeros( (n, npars) )
        cen1,cen2=self._cen_prior.sample(n)

        guess[:,0]=cen1
        guess[:,1]=cen2
        guess[:,2]=0.1*srandu(n)
        guess[:,3]=0.1*srandu(n)

        # need to tune these guesses
        # bulge T
        guess[:,4] = T_guess*2.0*(1.0 + 0.1*srandu(n))
        # disk T
        guess[:,5] = T_guess*0.7*(1.0 + 0.1*srandu(n))

        # bulge flux
        guess[:,6] = 2*counts_guess*(1.0 + 0.1*srandu(n))
        # disk flux
        guess[:,7] = 2*counts_guess*(1.0 + 0.1*srandu(n))

        return guess

    def get_jacobian(self, rowcen, colcen):
        """
        Get a jacobian for the input center
        """
        import ngmix
        j=ngmix.Jacobian(rowcen,
                         colcen,
                         self._pixel_scale,
                         0.,
                         0.,
                         self._pixel_scale)
        return j

    def print_some_stats(self, res):
        """
        Print some stats from the fitter result
        """
        from ngmix.fitting import print_pars
        print_pars(res['pars'],front='      pars: ')
        print_pars(res['pars_err'],front='      perr: ')
        if res['model'] in ['exp','dev']:
            flux_s2n=res['pars'][5]/res['pars_err'][5]
            T_s2n=res['pars'][4]/res['pars_err'][4]
            print >>stderr,'      T_s2n: %g flux_s2n: %g' % (T_s2n, flux_s2n)

        mess='      arate: %(arate)g chi2per: %(chi2per)g aic: %(aic)g bic: %(bic)s'
        mess = mess % res
        print >>stderr,mess


    def make_galaxy_image(self, dindex, show=False, get_obj=False):
        """
        Make the galaxy image object
        """

        rindex=self._index_list[dindex]
        gal_obj=self.make_galaxy_obj(rindex)

        # apply center offset in pixel coords
        dx = self._randu() - 0.5
        dy = self._randu() - 0.5

        gal_im_obj=gal_obj.draw(dx=self._pixel_scale, offset=(dx,dy))

        gal_im_obj.addNoise(self._rand_gauss_obj)
        #sky_var=gal_obj.noise.applyWhiteningTo(gal_im_obj)
        #sky_sigma = numpy.sqrt(sky_var)

        sh=gal_im_obj.array.shape
        #wt=numpy.zeros(sh) + 1./sky_sigma**2
        wt=numpy.zeros(sh) + 1./self._sky_sigma0**2
        #print 'sky_sigma0:',self._sky_sigma0,'sky_sigma:',sky_sigma

        if show:
            import images
            print 'max image:',gal_im_obj.array.max()
            images.multiview(gal_im_obj.array)

        if get_obj:
            return gal_obj, gal_im_obj, wt
        else:
            return gal_im_obj, wt

    def make_galaxy_obj(self, rindex):
        """
        Make the galaxy object
        """

        gal_pre = galsim.RealGalaxy(self._real_cat, index=rindex)

        # Rotate by a random angle
        if self._rotate:
            theta = 2.*numpy.pi * self._randu() * galsim.radians
            gal_pre.applyRotation(theta)

        # Make the pixel and psf convolved image
        gal_obj = galsim.Convolve([self._psf_obj, self._pix_obj, gal_pre])
        return gal_obj


    def read_data(self):
        """
        Read the catalogs and set pixel scale and nobj
        """
        self._cat=files.read_cat()
        self._fit_cat=files.read_fits_cat()
        self._real_cat = files.get_galsim_catalog() # RealGalaxyCatalog
        self._nobj=self._real_cat.getNObjects()

    def set_priors(self):
        """
        Set some priors
        """
        import ngmix

        # all pretty broad
        self._T_prior=ngmix.priors.FlatPrior(0.0001, 10000.0)
        self._counts_prior=ngmix.priors.FlatPrior(1.0e-6, 1.0e6)
        self._bfrac_prior=ngmix.priors.TruncatedGaussian(0.0,0.1,0.0,1.0)

        # width arcsec
        self._cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 0.1, 0.1)

    def setup_rand(self):
        """
        Set up the uniform random number generator and
        the gaussian generator
        """
        self._randu  = galsim.UniformDeviate()
        self._rand_gauss_obj = galsim.GaussianNoise(self._randu,
                                                    sigma=self._sky_sigma0)

    def copy_result(self, dindex, res):
        """
        Copy a result dict from the fitter to the output
        """
        model=res['model']
        n=get_model_names(model)

        data=self._data

        data['T_guess0'][dindex] = res['T_guess0']
        data['T_guess'][dindex] = res['T_guess']
        data['counts_guess'][dindex] = res['counts_guess']

        data[n['flags']][dindex] = res['flags']
        if res['flags']==0:
            pars=res['pars']
            pars_cov=res['pars_cov']

            data[n['pars']][dindex,:] = pars
            data[n['pars_cov']][dindex,:,:] = pars_cov

            for sn in _stat_names:
                data[n[sn]][dindex] = res[sn]

            data[n['arate']][dindex] = res['arate']
            if res['tau'] is not None:
                data[n['tau']][dindex] = res['tau']


    def make_struct(self):
        """
        make the output structure
        """
        import ngmix

        dt=[('id','i4'),
            ('T_guess0','f8'),
            ('T_guess','f8'),
            ('counts_guess','f8')]

        for model in self._models:
            n=get_model_names(model)

            np=ngmix.gmix.get_model_npars(model)

            dt+=[(n['flags'],'i4'),
                 (n['pars'],'f8',np),
                 (n['pars_cov'],'f8',(np,np)),
                
                 (n['s2n_w'],'f8'),
                 (n['chi2per'],'f8'),
                 (n['dof'],'f8'),
                 (n['aic'],'f8'),
                 (n['bic'],'f8'),
                 (n['arate'],'f8'),
                 (n['tau'],'f8'),
                ]


        num=self._index_list.size
        data=numpy.zeros(num, dtype=dt)

        for model in self._models:
            n=get_model_names(model)

            data[n['flags']] = NO_ATTEMPT

            data[n['pars']] = DEFVAL
            data[n['pars_cov']] = PDEFVAL

            data[n['s2n_w']] = DEFVAL
            data[n['chi2per']] = PDEFVAL
            data[n['aic']] = BIG_PDEFVAL
            data[n['bic']] = BIG_PDEFVAL

            data[n['tau']] = PDEFVAL
     
        data['id'] = self._index_list
        self._data=data

def get_model_names(model):
    names=['rfc_flags',
           'rfc_tries',
           'rfc_iter',
           'rfc_pars',
           'rfc_pars_cov',
           'flags',
           'pars',
           'pars_cov',
           'flux',
           'flux_err',
           'flux_cov',
           'g',
           'g_cov',
           'g_sens',
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


_guess_Tfac={'exp':1.0,
             'dev':25.0, # can be larger for some galaxies!
             'bdc':10.0,
             'bdf':1.0} # need to calibrate this
_stat_names=['s2n_w',
             'chi2per',
             'dof',
             'aic',
             'bic']


def srandu(n=None):
    return 2*(numpy.random.random(n)-0.5)
