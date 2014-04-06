"""
Create some sersic profiles using galsim and fit gaussian
mixtures to them
"""
from __future__ import print_function

import os
from numpy import zeros, sqrt, array
from pprint import pprint

def fit_sersic(n=1.0, hlr=4.0, ngauss=6, **keys):
    fitter=SersicFitter(n, hlr, ngauss)
    fitter.dofit(**keys)

    fitter.make_plots()
    fitter.write_fit()

    return fitter

def get_dir():
    dir=os.environ['LENSDIR']
    dir=os.path.join(dir,'sersic-gmix')
    return dir

def get_root(n, ngauss):
    dir=get_dir()
    root='sersic-%.2f-ngauss-%d' % (n, ngauss)
    root=os.path.join(dir, root)
    return root

def get_fname(n, ngauss, type, ext='fits'):
    root=get_root(n, ngauss)
    fname='%s-%s.%s' % (root, type, ext)
    return fname

def get_spline_fname(ngauss):
    d=get_dir()
    pfile='sersic-ngauss-%d-splines.pickle' % ngauss
    pfile=os.path.join(d, pfile)
    return pfile

def read_spline_data(ngauss):
    import pickle
    fname=get_spline_fname(ngauss)
    print("reading:",fname)
    with open(fname) as fobj:
        spline_data=pickle.load(fobj)

    return spline_data

def get_combined_fname(ngauss):
    d=get_dir()
    pfile='sersic-ngauss-%d-combined.fits' % ngauss
    pfile=os.path.join(d, pfile)
    return pfile


class SersicFitter(dict):
    def __init__(self, n=1.0, hlr=4.0, ngauss=6):
        """
        hlr in pixels
        """
        import ngmix

        self['n']=n
        self['hlr']=hlr
        self['ngauss']=ngauss

        self['pixel_scale']=1.0
        self['flux']=100.0
        self['sky_noise']=1.0e-6

        # need some psf
        self['psf_sigma']=1.414

        self['npars']=ngmix.gmix.get_coellip_npars(ngauss)

        self._make_names()

        self._make_image()
        self._make_jacobian()

    def dofit(self, **keys):
        self['nwalkers'] = keys.get('nwalkers',320)
        self['burnin'] = keys.get('burnin',5000)
        self['nstep'] = keys.get('nstep',1000)
        pprint(self)

        self._do_fit()

    def _do_fit(self):
        import ngmix

        self._fit_psf()
        self._fit_galaxy()

    def _fit_galaxy(self):
        import ngmix

        first_T_mean=3.0e-3*(self['hlr']/4.0)
        first_T_sigma=first_T_mean*1.0e-4
        first_T_prior=ngmix.priors.LogNormal(first_T_mean, first_T_sigma)

        guess=self._get_guess(first_T_prior)

        fitter=ngmix.fitting.MCMCCoellip(self.image,
                                         self.wt,
                                         self.jacobian,
                                         psf=self.psf_gmix,
                                         nwalkers=self['nwalkers'],
                                         burnin=self['burnin'],
                                         nstep=self['nstep'],
                                         first_T_prior=first_T_prior,
                                         full_guess=guess)
        fitter.go()

        self.fitter=fitter
        self.res=fitter.get_result()

        self._add_pars_norm()

    def _add_pars_norm(self):
        import ngmix

        ngauss=self['ngauss']

        pars_norm=self.res['pars'].copy()
        Fvals=pars_norm[4+ngauss:].copy()
        ftot=Fvals.sum()

        Tvals=pars_norm[4:4+ngauss].copy()
        Ttot=(Tvals*Fvals).sum()/ftot
        Tnorm=Tvals/Ttot

        Fnorm=Fvals/ftot

        pars_norm[4+ngauss:] = Fnorm
        pars_norm[4:4+ngauss] = Tnorm

        self.res['pars_norm'] = pars_norm

        ngmix.fitting.print_pars(pars_norm,front='pars norm: ')


    def _get_guess(self, first_T_prior):
        import ngmix
        from ngmix.priors import srandu

        nwalkers=self['nwalkers']

        ngauss=self['ngauss']


        if os.path.exists(self.fits_name):
            import fitsio
            print("reading guess from:",self.fits_name)
            data=fitsio.read(self.fits_name)
            pars0 = data['pars'][0,4:]

            T=1.0
            F=1.0
        else:
            if ngauss==4:
                pars0=array([0.01183116, 0.06115546,  0.3829298 ,  2.89446939,
                             0.19880675,  0.18535747, 0.31701891,  0.29881687])
            elif ngauss==5:
                pars0=array([1.e-4,    0.01183116, 0.06115546,  0.3829298 ,  2.89446939,
                             0.2,  0.2,        0.2,         0.2,         0.2])
            elif ngauss==6:
                #pars0=array([1.e-4,    0.01183116, 0.1,  0.3829298 ,  2.89446939,    8.0,
                #             0.16,     0.16,       0.16,        0.16,         0.16,         0.16])
                pars0=array([0,            0.00259716,  0.0595377,    0.45871,    1.0, 4.0, 
                             3.85681e-05,  0.0167603,  0.0513565,   0.515547,    0.1, 0.3])

            # this is a terrible guess
            T = 2*self['hlr']**2
            F = self['flux']

        full_guess=zeros( (nwalkers, self['npars']) )
        full_guess[:,0] = 0.1*srandu(nwalkers)
        full_guess[:,1] = 0.1*srandu(nwalkers)
        full_guess[:,2] = 0.01*srandu(nwalkers)
        full_guess[:,3] = 0.01*srandu(nwalkers)

        for i in xrange(ngauss):
            if i==0 and first_T_prior is not None:
                print("sampling first_T_prior")
                full_guess[:,4+i] = first_T_prior.sample(nwalkers)
                full_guess[:,4+ngauss+i] = F*pars0[ngauss+i]*(1.0 + 0.01*srandu(nwalkers))
            else:
                full_guess[:,4+i] = T*pars0[i]*(1.0 + 0.01*srandu(nwalkers))
                full_guess[:,4+ngauss+i] = F*pars0[ngauss+i]*(1.0 + 0.01*srandu(nwalkers))

        return full_guess

    def _fit_psf(self):
        import ngmix

        # psf using EM
        im_psf_sky,sky=ngmix.em.prep_image(self.psf_image)

        fitter=ngmix.em.GMixEM(im_psf_sky)

        cen=(self.psf_image.shape[0]-1)/2.0

        parguess=[1.0, cen, cen, self['psf_sigma']**2, 0.0, self['psf_sigma']**2]
        guess=ngmix.gmix.GMix(pars=parguess)

        fitter.go(guess, sky, maxiter=5000)

        res=fitter.get_result()
        self.psf_gmix = fitter.get_gmix()

        print('psf numiter:',res['numiter'],'fdiff:',res['fdiff'])
        print("psf_gmix:",self.psf_gmix)


    def _make_image(self):
        import galsim 

        gal = galsim.Sersic(self['n'],
                            half_light_radius=self['hlr'],
                            flux=self['flux'])

        psf = galsim.Gaussian(sigma=self['psf_sigma'], flux=1.0)
        pixel=galsim.Pixel(self['pixel_scale'])

        gal_final = galsim.Convolve([gal, psf, pixel])
        psf_final = galsim.Convolve([psf, pixel])

        # deal with massive unannounced api changes
        try:
            image_obj = gal_final.draw(scale=self['pixel_scale'])
            psf_obj   = psf_final.draw(scale=self['pixel_scale'])
        except:
            image_obj = gal_final.draw(dx=self['pixel_scale'])
            psf_obj   = psf_final.draw(dx=self['pixel_scale'])

        image_obj.addNoise(galsim.GaussianNoise(sigma=self['sky_noise']))

        image = image_obj.array.astype('f8')

        print("image sum:",image.sum())

        # never use more than 200x200
        #if image.shape[0] > 200:
        #    cen=int(0.5*image.shape[0])
        #    xmin=cen-100
        #    xmax=cen+100
        #    if xmin < 0:
        #        xmin=0
        #    if xmax > image.shape[1]:
        #        xmax=image.shape[1]
        #
        #    image=image[xmin:xmax, xmin:xmax]

        psf_image = psf_obj.array.astype('f8')

        wt = image*0 + ( 1.0/self['sky_noise']**2 )

        self.image=image
        self.psf_image = psf_image
        self.wt=wt
        print("image dims:",image.shape)
        print("psf image dims:",psf_image.shape)

    def _make_jacobian(self):
        import ngmix
        cen0=(self.image.shape[0]-1)/2.
        self.jacobian=ngmix.jacobian.UnitJacobian(cen0, cen0)

    def show_image(self):
        import images
        images.multiview(self.image,title='object')

    def show_psf_image(self):
        import images
        images.multiview(self.psf_image,title='psf')

    def make_plots(self):
        import images

        burnp, histp=self.fitter.make_plots(separate=True,
                                            show=False,
                                            fontsize_min=0.2)

        gmix0=self.fitter.get_gmix()
        gmix=gmix0.convolve(self.psf_gmix)
        model_image=gmix.make_image(self.image.shape,
                                    jacobian=self.jacobian)

        residp = images.compare_images(self.image, model_image,
                                       label1='image',
                                       label2='model',
                                       show=False)

        print(self.burn_png)
        burnp.write_img(1800,1000,self.burn_png)

        print(self.hist_png)
        histp.write_img(1800,1000,self.hist_png)

        print(self.resid_png)
        residp.write_img(1800,1000,self.resid_png)

    def _make_names(self):
        self.fits_name=get_fname(self['n'],self['ngauss'],'fit',ext='fits')
        self.burn_png=get_fname(self['n'],self['ngauss'],'steps',ext='png')
        self.hist_png=get_fname(self['n'],self['ngauss'],'hist',ext='png')
        self.resid_png=get_fname(self['n'],self['ngauss'],'resid',ext='png')

    def write_fit(self):
        import fitsio


        npars=self['npars']
        output=zeros(1, dtype=[('n','f8'),
                               ('hlr','f8'),
                               ('arate','f8'),
                               ('chi2per','f8'),
                               ('dof','f8'),
                               ('pars','f8',npars),
                               ('pars_norm','f8',npars),
                               ('pars_err','f8',npars),
                               ('pars_cov','f8',(npars,npars))])


        output['n'] = self['n']
        output['hlr'] = self['hlr']
        output['arate'] = self.res['arate']
        output['chi2per'] = self.res['chi2per']
        output['dof'] = self.res['dof']
        output['pars'][0,:] = self.res['pars']
        output['pars_norm'][0,:] = self.res['pars_norm']
        output['pars_err'][0,:] = self.res['pars_err']
        output['pars_cov'][0,:,:] = self.res['pars_cov']

        print(self.fits_name)
        fitsio.write(self.fits_name, output, clobber=True)


