from sys import stderr
import numpy
import galsim

from . import files

class CosmosFitter(object):
    def __init__(self, models, psf_fwhm=0.2, pixel_scale=0.15, sky_sigma=0.0001):
        self._models      = models
        self._psf_fwhm    = psf_fwhm    # arcsec
        self._sky_sigma   = sky_sigma
        self._pixel_scale = pixel_scale

        self.read_data()
        self.setup_rand()

    def do_fits(self):
        """
        Perform all the fits
        """
        self.make_pixel()
        self.make_psf()

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
                                         flux=100.0)
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

        self.psf_gmix = em.get_gmix()

        print 'psf fit:'
        print self.psf_gmix

    def fit_galaxies(self):
        """
        Convolve the real galaxy images with a psf to make it a
        simple gaussian
        """
        for i in xrange(self._nobj):
            gal_im_obj,wt = self.make_galaxy_image(i)

            self.fit_galaxy_models(gal_im_obj,wt)

    def fit_galaxy_models(self, gal_im_obj):
        """
        Fit all the models
        """
        import ngmix

        moms=gal_im_obj.FindAdaptiveMoments()
        flux_guess = moms.moments_amp
        T_guess    = 4 * 2*moms.moments_sigma**2
        for model in self._models:
            fitter=ngmix.fitting.MCMCSimple

    def make_galaxy_image(self, index, add_noise=True, show=False):
        """
        Make the galaxy image object
        """
        gal_obj=self.make_galaxy_obj(index)

        # apply center offset in pixel coords
        dx = self._randu() - 0.5
        dy = self._randu() - 0.5

        #print >>stderr,"drawing"
        gal_im_obj=gal_obj.draw(dx=self._pixel_scale, offset=(dx,dy))

        wt=0*gal_im_obj.array + 1./self._sky_sigma**2
        if add_noise:
            gal_im_obj.addNoise(self._rand_gauss_obj)

        if show:
            import images
            print 'max image:',gal_im_obj.array.max()
            images.multiview(gal_im_obj.array)

        return gal_im_obj, wt

    def make_galaxy_obj(self, index):
        """
        Make the galaxy object
        """

        #print >>stderr,"getting real gal_pre"
        gal_pre = galsim.RealGalaxy(self._real_cat, index = index)
        #gal.setFlux(gal_flux)

        # Rotate by a random angle
        theta = 2.*numpy.pi * self._randu() * galsim.radians
        gal_pre.applyRotation(theta)

        # Make the pixel and psf convolved image
        #print >>stderr,"convolving"
        gal_obj = galsim.Convolve([self._psf_obj, self._pix_obj, gal_pre])
        return gal_obj
        #return gal_pre


    def read_data(self):
        """
        Read the catalogs and set pixel scale and nobj
        """
        self._cat=files.read_cat()
        self._fit_cat=files.read_fits_cat()
        self._real_cat = files.get_galsim_catalog() # RealGalaxyCatalog
        self._nobj=self._real_cat.getNObjects()

    def setup_rand(self):
        """
        Set up the uniform random number generator and
        the gaussian generator
        """
        self._randu  = galsim.UniformDeviate()
        self._rand_gauss_obj = galsim.GaussianNoise(self._randu,
                                                    sigma=self._sky_sigma)
