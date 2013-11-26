from sys import stderr
import numpy
import galsim

from . import files

_guess_Tfac={'exp':2.0,
             'dev':50.0} # can be larger for some galaxies!
class CosmosFitter(object):
    def __init__(self,
                 models,
                 obj_range=None,
                 psf_fwhm=0.3,
                 pixel_scale=0.1,
                 sky_sigma0=0.005, # before whitening, ~0.01 after
                 rotate=True,
                 make_plots=False):
        self._models      = models
        self._psf_fwhm    = float(psf_fwhm)    # arcsec
        self._sky_sigma0  = float(sky_sigma0)
        self._rotate      = rotate
        self._pixel_scale = float(pixel_scale)
        self._make_plots  = make_plots


        self._nwalkers=40
        self._burnin=400
        self._nstep=100


        self.read_data()
        self.set_index_list(obj_range)

        self.setup_rand()

        self.make_pixel()
        self.make_psf()

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

        self.set_priors()
        last=self._index_list[-1]

        for index in self._index_list:
            print >>stderr,'%s:%s' % (index, last)
            gal_im_obj,wt = self.make_galaxy_image(index)
            self.fit_galaxy_models(gal_im_obj,wt)


    def fit_galaxy_models(self, gal_im_obj, wt):
        """
        Fit all the models
        """
        import ngmix
        from ngmix.fitting import print_pars

        moms=gal_im_obj.FindAdaptiveMom()

        counts_guess = moms.moments_amp * self._pixel_scale**2
        T_guess0    = 2*moms.moments_sigma**2 * self._pixel_scale**2
        print >>stderr,'  T_guess0:     ',T_guess0
        print >>stderr,'  counts_guess: ',counts_guess

        j=ngmix.Jacobian(moms.moments_centroid.y-1.,
                         moms.moments_centroid.x-1.,
                         self._pixel_scale,
                         0.,
                         0.,
                         self._pixel_scale)

        im=gal_im_obj.array.astype('f8')
        for model in self._models:
            print >>stderr,'    model:',model
            T_guess=_guess_Tfac[model]*T_guess0
            print >>stderr,'    T_guess:      ',T_guess
            fitter=ngmix.fitting.MCMCSimple(im, wt, j, model,
                                            psf=self._psf_gmix,
                                            T_guess=T_guess,
                                            counts_guess=counts_guess,
                                            T_prior=self._T_prior,
                                            cen_prior=self._cen_prior,
                                            nwalkers=self._nwalkers,
                                            burnin=self._burnin,
                                            nstep=self._nstep)
            fitter.go()
            res=fitter.get_result()
            print_pars(res['pars'],front='      pars: ')
            print_pars(res['perr'],front='      perr: ')
            flux_s2n=res['pars'][5]/res['perr'][5]
            T_s2n=res['pars'][4]/res['perr'][4]
            print >>stderr,'      T_s2n: %g flux_s2n: %g' % (T_s2n, flux_s2n)
            print >>stderr,'      arate: %(arate)g chi2per: %(chi2per)g' % res

            if self._make_plots:
                fitter.make_plots(show=True, do_residual=True, title=model)


    def make_galaxy_image(self, index, show=False, get_obj=False):
        """
        Make the galaxy image object
        """
        gal_obj=self.make_galaxy_obj(index)

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

    def make_galaxy_obj(self, index):
        """
        Make the galaxy object
        """

        gal_pre = galsim.RealGalaxy(self._real_cat, index = index)

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
        # width arcsec
        self._cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 0.1, 0.1)
        #self._cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 1.0, 1.0)

    def setup_rand(self):
        """
        Set up the uniform random number generator and
        the gaussian generator
        """
        self._randu  = galsim.UniformDeviate()
        self._rand_gauss_obj = galsim.GaussianNoise(self._randu,
                                                    sigma=self._sky_sigma0)

    def make_struct(self):
        """
        make the output structure
        """

        dt=[('id','i4')]

        simple_npars=6
        cov_shape=(simple_npars,simple_npars)

        for model in self._models:
            n=get_model_names(model)

            np=simple_npars

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
     
        self.data=data

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


