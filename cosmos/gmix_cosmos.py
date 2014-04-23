from __future__ import print_function
from sys import stderr
import time
import numpy
from numpy import sqrt, array, zeros, where
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

    print('writing:',output_file)
    with fitsio.FITS(output_file, 'rw', clobber=True) as fobj:
        fobj.write(data, extension="model_fits")

class CosmosFitter(object):
    """
    Fit models to the cosmos galaxies released as part of great3.
    
    Convolve with a gaussian PSF, and add a little noise.
    """
    def __init__(self,
                 model,
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
                 ntry=1,
                 **keys):

        self._conf=keys
        self._model       = model
        self._psf_fwhm    = float(psf_fwhm)    # arcsec
        self._sky_sigma0  = float(sky_sigma0)
        self._rotate      = rotate
        self._pixel_scale = float(pixel_scale)
        self._pixel_area  = self._pixel_scale**2
        self._make_plots  = make_plots
        self._plot_base   = plot_base

        self._max_image_size = keys.get('max_size',96)

        self._shear = keys.get('shear',None)
        self._shear_expand=self._shear
        self._do_pqr = (self._shear is not None)

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

        self.set_priors()
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
        print('psf guess:')
        print(gmix_guess)

        em.go(gmix_guess, sky, maxiter=50000, tol=1.0e-6)

        self._psf_gmix = em.get_gmix()

        print('psf fit:')
        print(self._psf_gmix)
        print("")

    def fit_galaxies(self):
        """
        Convolve the real galaxy images with a psf to make it a
        simple gaussian
        """

        tm0=time.time()

        last=self._index_list[-1]
        ndo=self._index_list.size

        for dindex in xrange(ndo):

            self.res={'flags':0}

            rindex=self._index_list[dindex]
            self.rindex=rindex

            print('%s:%s' % (rindex, last))
            gal_im_obj,wt = self.make_galaxy_image(dindex)

            self.gal_image=gal_im_obj.array.astype('f8')
            shape=self.gal_image.shape
            print(shape)

            if (shape[0] > self._max_image_size
                    or shape[1] > self._max_image_size):
                print("    Skipping large image")
                self._data['flags'][dindex] = NO_ATTEMPT
            else:
                self.weight_image = wt
                self.gal_cen_guess=(array(shape)-1.0)/2.0

                self.fit_galaxy_em()
                self.fit_simple_model()

                self.copy_result(dindex)
                self.print_some_stats()

                if self._make_plots:
                    self._do_gal_plots()

        tm=time.time()-tm0
        print('time per galaxy:',tm/ndo)


    def fit_galaxy_em(self):
        """

        Fit a single gaussian with em to find a decent center.  We don't get a
        flux out of that, but we can get flux using the _fit_flux routine

        """

        # first the structural fit
        sigma_guess = sqrt( self._psf_gmix.get_T()/2.0 )
        print('      sigma guess:',sigma_guess)
        fitter=self._fit_em_1gauss(self.gal_image,
                                   self.gal_cen_guess,
                                   sigma_guess)

        em_gmix = fitter.get_gmix()
        print("      em gmix:",em_gmix)

        row_rel, col_rel = em_gmix.get_cen()
        em_cen = self.gal_cen_guess + array([row_rel,col_rel])

        print("      em gauss cen:",em_cen)

        jacob=self.get_jacobian(em_cen)

        # now get a flux
        print('    fitting robust gauss flux')
        flux, flux_err = self._fit_flux(self.gal_image,
                                        self.weight_image,
                                        jacob,
                                        em_gmix)

        self.res['em_gauss_flux'] = flux
        self.res['em_gauss_flux_err'] = flux_err
        self.res['em_gauss_cen'] = jacob.get_cen()
        self.res['em_gmix'] = em_gmix
        self.res['jacob'] = jacob


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




    def _fit_em_1gauss(self, im, cen_guess, sigma_guess):
        """
        Fit the image using EM

        cen is the global cen not within jacobian
        """
        import ngmix
        from ngmix.gexceptions import GMixMaxIterEM, GMixRangeError

        im_with_sky, sky = ngmix.em.prep_image(im)

        jacob = self.get_jacobian(cen_guess)

        ntry=200
        maxiter=5000
        tol=5.0e-4
        for i in xrange(ntry):
            guess = self._get_em_guess_1gauss(sigma_guess)
            print("em guess:")
            print(guess)
            try:
                fitter=ngmix.em.GMixEM(im_with_sky, jacobian=jacob)
                fitter.go(guess, sky, maxiter=maxiter, tol=tol)


                tres=fitter.get_result()
                print("em numiter:",tres['numiter'])
                break
            except GMixMaxIterEM:
                fitter=None
            except GMixRangeError:
                fitter=None

        return fitter


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



    def fit_simple_model(self):
        """
        fit a single model
        """

        import ngmix
        
        # possible value errors
        nretry=10
        for i in xrange(nretry):
            try:
                full_guess=self._get_guess_simple()
                fitter=ngmix.fitting.MCMCSimple(self.gal_image,
                                                self.weight_image,
                                                self.res['jacob'],
                                                self._model,
                                                psf=self._psf_gmix,

                                                T_prior=self._T_prior,
                                                counts_prior=self._counts_prior,
                                                cen_prior=self._cen_prior,
                                                g_prior=self._g_prior,
                                                g_prior_during=self._g_prior_during,

                                                shear_expand=self._shear_expand,

                                                nwalkers=self._nwalkers,
                                                burnin=self._burnin,
                                                nstep=self._nstep,
                                                mca_a=self._mca_a,

                                                ntry=self._ntry,
                                                min_arate=self._min_arate,
                                                do_pqr=self._do_pqr)
                fitter.go()
                break
            except ValueError as err:
                print("got value error: %s" % str(err))


        self.res['fitter'] = fitter
        self.res.update( fitter.get_result() )



    def _get_guess_simple(self,
                          widths=[0.01, 0.01, 0.01, 0.01, 0.01, 0.01]):
        """
        Get a guess centered on the truth

        width is relative for T and counts
        """

        

        nwalkers = self._nwalkers

        res=self.res
        gmix = res['em_gmix']
        g1,g2,T = gmix.get_g1g2T()
        flux = res['em_gauss_flux']

        guess=numpy.zeros( (nwalkers, 6) )

        # centers relative to jacobian center
        guess[:,0] = widths[0]*srandu(nwalkers)
        guess[:,1] = widths[1]*srandu(nwalkers)

        guess_shape=get_shape_guess(g1, g2, nwalkers, width=widths[2])
        guess[:,2]=guess_shape[:,0]
        guess[:,3]=guess_shape[:,1]

        guess[:,4] = get_positive_guess(T,nwalkers,width=widths[4])
        guess[:,5] = get_positive_guess(flux,nwalkers,width=widths[5])

        return guess




    def do_make_plots(self, fitter, model):
        trials_plot=self._plot_base+'-%s-%06d-trials.png' % (model,self.rindex)
        diff_plot=self._plot_base+'-%s-%06d-diff.png' % (model,self.rindex)

        tp,dplist=fitter.make_plots(show=False, do_residual=True, title=model)
        dp=dplist[0]

        print(trials_plot)
        tp.write_img(1100,1100,trials_plot)
        print(diff_plot)
        dp.write_img(1100,1100,diff_plot)


    def _do_gal_plots(self):
        """
        Make residual plot and trials plot
        """
        self._compare_gal()
        self._make_trials_plot()

    def _make_trials_plot(self):
        """
        Plot the trials
        """
        width,height=800,800
        tup=self.res['fitter'].make_plots(title=self._model)

        if isinstance(tup, tuple):
            p,wp = tup
            wtrials_pname='wtrials-%06d-%s.png' % (self.rindex,self._model)
            print("          ",wtrials_pname)
            wp.write_img(width,height,wtrials_pname)
        else:
            p = tup

        trials_pname='trials-%06d-%s.png' % (self.rindex,self._model)
        print("          ",trials_pname)
        p.write_img(width,height,trials_pname)

    def _compare_gal(self):
        """
        compare psf image to best fit model
        """
        import images

        res=self.res
        gmix = res['fitter'].get_gmix()

        psf_gmix = self._psf_gmix
        gmix_conv = gmix.convolve(psf_gmix)

        model_image = gmix_conv.make_image(self.gal_image.shape,
                                           jacobian=res['jacob'])

        plt=images.compare_images(self.gal_image,
                                  model_image,
                                  label1='galaxy',
                                  label2=self._model,
                                  show=False)
        pname='gal-resid-%06d-%s.png' % (self.rindex,self._model)
        print("          ",pname)
        plt.write_img(1400,800,pname)



    def get_jacobian(self, cen):
        """
        Get a jacobian for the input center
        """
        import ngmix
        j=ngmix.Jacobian(cen[0],
                         cen[1],
                         self._pixel_scale,
                         0.,
                         0.,
                         self._pixel_scale)
        return j

    def print_some_stats(self):
        """
        Print some stats from the fitter result
        """
        from ngmix.fitting import print_pars

        res=self.res

        print_pars(res['pars'],front='      pars: ')
        print_pars(res['pars_err'],front='      perr: ')

        T_s2n=res['pars'][4]/res['pars_err'][4]
        flux_s2n=res['pars'][5]/res['pars_err'][5]
        print('      T_s2n: %g flux_s2n: %g' % (T_s2n, flux_s2n))

        mess='      arate: %(arate)g chi2per: %(chi2per)g'
        mess = mess % res
        print(mess)


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
            print('max image:',gal_im_obj.array.max())
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

        if self._shear is not None:
            gal_pre.applyShear(g1=self._shear[0], g2=self._shear[1])

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

        # width arcsec
        self._cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 0.1, 0.1)

        self._g_prior=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
        self._g_prior_during=self._conf.get('g_prior_during',False)

    def setup_rand(self):
        """
        Set up the uniform random number generator and
        the gaussian generator
        """
        self._randu  = galsim.UniformDeviate()
        self._rand_gauss_obj = galsim.GaussianNoise(self._randu,
                                                    sigma=self._sky_sigma0)

    def copy_result(self, dindex):
        """
        Copy a result dict from the fitter to the output
        """

        res=self.res

        data=self._data

        data['flags'][dindex] = res['flags']
        if res['flags']==0:
            pars=res['pars']
            pars_cov=res['pars_cov']

            data['pars'][dindex,:] = pars
            data['pars_cov'][dindex,:,:] = pars_cov

            for sn in _stat_names:
                data[sn][dindex] = res[sn]

            data['arate'][dindex] = res['arate']
            if res['tau'] is not None:
                data['tau'][dindex] = res['tau']

            if self._do_pqr:
                data['P'][dindex] = res['P']
                data['Q'][dindex] = res['Q']
                data['R'][dindex] = res['R']

    def make_struct(self):
        """
        make the output structure
        """
        import ngmix

        np=ngmix.gmix.get_model_npars(self._model)
        dt=[('id','i4'),
            ('flags','i4'),
            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),
            ('s2n_w','f8'),
            ('chi2per','f8'),
            ('dof','f8'),
            ('aic','f8'),
            ('bic','f8'),
            ('arate','f8'),
            ('tau','f8') ]

        if self._do_pqr:
            dt += [('P','f8'),
                   ('Q','f8',2),
                   ('R','f8',(2,2))]

        num=self._index_list.size
        data=numpy.zeros(num, dtype=dt)


        data['flags'] = NO_ATTEMPT

        data['pars'] = DEFVAL
        data['pars_cov'] = PDEFVAL

        data['s2n_w'] = DEFVAL
        data['chi2per'] = PDEFVAL
        data['aic'] = BIG_PDEFVAL
        data['bic'] = BIG_PDEFVAL

        data['tau'] = PDEFVAL
     
        data['id'] = self._index_list

        if self._do_pqr:
            data['P'] = DEFVAL
            data['Q'] = DEFVAL
            data['R'] = DEFVAL

        self._data=data


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

def get_shape_guess(g1, g2, n, width=0.01):
    """
    Get guess, making sure in range
    """
    import ngmix
    from ngmix.gexceptions import GMixRangeError

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
        w,=where(vals <= 0)
        if w.size == 0:
            break
        else:
            vals[w] = val*(1.0 + width*srandu(w.size))

    return vals


