from __future__ import print_function
import os
from sys import stderr,stdout
import time
import numpy
from numpy import array, sqrt, zeros, log, exp, arange

import ngmix
from ngmix.fitting import print_pars
from ngmix.gexceptions import GMixMaxIterEM, GMixRangeError
from ngmix.observation import Observation
from ngmix.jacobian import Jacobian

from .nfit import srandu

import fitsio
import meds

# starting new values for these
DEFVAL=-9999
PDEFVAL=9999

NO_CUTOUTS=2**0
PSF_FIT_FAILURE=2**1
PSF_LARGE_OFFSETS=2**2
GAL_FIT_FAILURE=2**3

NO_ATTEMPT=2**30

ISAMP_BAD_COV=2**7

_CHECKPOINTS_DEFAULT_MINUTES=[10,30,60,90]

class MedsFitBase(dict):
    def __init__(self, meds_file, truth_file, psf_file, **keys):
        self['meds_file']=meds_file
        self['truth_file']=truth_file
        self['psf_file']=psf_file

        self.update(keys)

        self.set_defaults()
        self.load_data()
        self.set_indices()

        self.make_struct()

    def set_defaults(self):
        """
        deal with default parameters and conversions
        """
        # in arcsec
        sigma_guess=self['psf_pars']['fwhm_guess']/2.35
        self['psf_Tguess'] = 2*sigma_guess**2

        self['make_plots']=self.get('make_plots',False)

        self['randomize_psf'] = self.get('randomize_psf',False)

    def get_data(self):
        """
        get a reference to the data structure
        """
        return self.data

    def do_fits(self):
        """
        loop and fit all objects
        """

        t0=time.time()

        last=self.indices.size-1
        for dindex,mindex in enumerate(self.indices):
            self.dindex=dindex
            self.mindex=mindex

            print("%d:%d  %d:%d" % (dindex, last, mindex, self['end']))

            self.data['number'][dindex] = self.meds['number'][mindex]

            self.make_psf_observation()
            self.make_galaxy_observation()

            self.run_fitters()

            if self['make_plots']:
                self.compare_psf()
                self.do_gal_plots()

        tm=time.time()-t0
        num=len(self.indices)
        print("time:",tm)
        print("time per:",tm/num)

    def run_fitters(self):
        from great3.generic import PSFFailure,GalFailure

        dindex=self.dindex
        boot=self.get_bootstrapper()

        # find the center and reset the jacobian
        boot.find_cen()

        self.boot=boot

        try:

            self.fit_psf()
            self.fit_psf_flux()

            try:

                self.fit_galaxy()
                self.copy_galaxy_result()
                self.print_galaxy_result()

            except GalFailure:
                print("    galaxy fitting failed")
                self.data['flags'][dindex] = GAL_FIT_FAILURE

        except PSFFailure:
            print("    psf fitting failed")
            self.data['flags'][dindex] = PSF_FIT_FAILURE


    def fit_psf(self):

        dindex=self.dindex
        boot=self.boot

        psf_pars=self['psf_pars']

        boot.fit_psf(psf_pars['model'],
                     Tguess=self['psf_Tguess'],
                     ntry=psf_pars['ntry'])

        self.psf_fitter=self.boot.get_psf_fitter()

        self.copy_psf_result()

    def fit_psf_flux(self):
        """
        fit psf model to galaxy with one free parameter for flux
        """
        boot=self.boot
        dindex=self.dindex

        boot.fit_gal_psf_flux()

        data=self.data
        data['psf_flux'][dindex] = boot.psf_flux
        data['psf_flux_err'][dindex] = boot.psf_flux_err

    def fit_galaxy(self):
        """
        over-ride for different fitters
        """
        raise RuntimeError("over-ride me")

    def fit_max(self):
        """
        do a maximum likelihood fit
        """
        boot=self.boot

        model=self['model_pars']['model']
        max_pars=self['max_pars']

        boot.fit_max(model,
                     max_pars,
                     prior=self['prior'],
                     ntry=max_pars['ntry'])

    def get_bootstrapper(self):
        """
        get the bootstrapper for fitting psf through galaxy
        """
        from great3.sfit import get_bootstrapper
        boot = get_bootstrapper(self.psf_obs, self.obs)
        return boot

    def make_psf_observation(self):
        """
        read the image and weight data
        """

        ext=self.get_psf_ext()

        self['psf_id'] = ext
        self.data['psf_id'][self.dindex] = ext

        image0 = self.psf_obj[ext][:,:]

        noise=self['psf_pars']['addnoise']
        image = image0 + numpy.random.normal(scale=noise,
                                             size=image0.shape)

        weight = image.copy()
        weight *= 0
        weight += 1.0/noise**2

        jacob=self.get_jacobian()
        self.psf_obs = Observation(image, weight=weight, jacobian=jacob)

    def get_psf_ext(self):
        if self['randomize_psf']:
            ext=numpy.random.randint(0,self.npsf)
            print("getting random psf:",ext)
        else:
            ext=self.truth['id_psf'][self.mindex]

        return ext

    def make_galaxy_observation(self):
        """
        read the image and weight data
        """

        image = self.meds.get_cutout(self.mindex, 0)

        calc_weight=self.get('calc_weight',False)
        if calc_weight:
            import esutil as eu
            border=5
            nrow,ncol=image.shape
            row,col=numpy.mgrid[0:nrow, 0:ncol]
            w=numpy.where((row < 5) | (row > (nrow-5-1)) 
                          | (col < 5) | (col > (ncol-5-1))  )

            mn, skysig_meas = eu.stat.sigma_clip(image[w].ravel())
            print("    skysig_meas: %g" % skysig_meas)

            weight = 0*image + 1.0/skysig_meas**2

        else:
            if self['noisefree']:
                weight = image*0 + 1.0/self['skynoise']**2
            else:
                weight = self.meds.get_cutout(self.mindex, 0, type='weight')

        jacob=self.get_jacobian()
        self.obs = Observation(image, weight=weight, jacobian=jacob)


    def get_psf_guesser(self):
        """
        get guesser based of size of psf and psf flux
        """
        from gmix_meds.util import FromPSFGuesser
        data=self.data
        T=data['psf_T'][self.dindex]
        flux=data['psf_flux'][self.dindex]

        if self['use_logpars']:
            scaling='log'
        else:
            scaling='linear'

        guesser=FromPSFGuesser(T, flux, scaling=scaling)

        return guesser

    def get_flux_and_prior_guesser(self):
        """
        from the psf flux and the prior
        """

        from gmix_meds.util import FluxAndPriorGuesser

        psf_flux=self.data['psf_flux'][self.dindex]
        psf_flux=psf_flux.clip(min=0.1, max=1.0e9)

        if self['use_logpars']:
            scaling='log'
        else:
            scaling='linear'

        guesser=FluxAndPriorGuesser(psf_flux, self['prior'],scaling=scaling)
        return guesser

    def get_T_flux_and_prior_guesser(self):
        """
        from the psf flux and the prior
        """
        from gmix_meds.util import TFluxAndPriorGuesser

        psf_gmix = self.psf_obs.get_gmix()
        psf_T = psf_gmix.get_T()

        psf_flux=self.data['psf_flux'][self.dindex]
        psf_flux=psf_flux.clip(min=0.1, max=1.0e9)

        if self['use_logpars']:
            scaling='log'
        else:
            scaling='linear'

        guesser=TFluxAndPriorGuesser(psf_T,
                                     psf_flux,
                                     self['prior'],
                                     scaling=scaling)
        return guesser



    def try_replace_cov(self, fitter):
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

    def compare_psf(self):
        """
        compare psf image to best fit model
        """
        import images

        fitter=self.psf_fitter

        model=self['psf_pars']['model']

        obs=self.psf_obs
        if 'em' in model:
            model_image = fitter.make_image(counts=obs.image.sum())
        else:
            gm=fitter.get_gmix()
            j=obs.get_jacobian()
            model_image = gm.make_image(obs.image.shape,
                                        jacobian=j)

        plt=images.compare_images(obs.image,
                                  model_image,
                                  label1='psf',
                                  label2=model,
                                  show=False)

        pname='psf-resid-%s-%06d.png' % (model, self.mindex)
        print("          ",pname)
        plt.write_img(1400,800,pname)

    def do_gal_plots(self):
        """
        Make residual plot and trials plot
        """
        res=self.gal_fitter.get_result()
        if res['flags'] != 0:
            return

        self.compare_gal()
        #self.make_trials_plot()
        #self.plot_autocorr()

    def compare_gal(self):
        """
        compare psf image to best fit model
        """
        import images

        fitter=self.gal_fitter

        model=self['model_pars']['model']
        title = '%d %s' % (self.mindex, model)

        gmix = fitter.get_gmix()

        obs = self.obs
        res=fitter.get_result()

        psf_gmix = self.psf_obs.get_gmix()
        gmix_conv = gmix.convolve(psf_gmix)

        image=obs.image
        model_image = gmix_conv.make_image(image.shape,
                                           jacobian=obs.get_jacobian())

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


    def get_jacobian(self):
        """
        get the jacobian and return a Jacobian object
        """
        jdict = self.meds.get_jacobian(self.mindex,0)
        jacob = Jacobian(jdict['row0']-1,
                         jdict['col0']-1,
                         jdict['dudrow'],
                         jdict['dudcol'],
                         jdict['dvdrow'],
                         jdict['dvdcol'])

        return jacob

    def load_data(self):
        """
        read or load all data
        """
        print("loading:",self['meds_file'])
        self.meds=meds.MEDS(self['meds_file'])

        print("reading:",self['truth_file'])
        self.truth=fitsio.read(self['truth_file'])

        print("loading:",self['psf_file'])
        self.psf_obj = fitsio.FITS(self['psf_file'])
        self.npsf = len(self.psf_obj)

    def set_indices(self):
        """
        this version we don't support work dir
        """

        obj_range = self.get('obj_range',None)

        if obj_range is not None:
            self['start'] = obj_range[0]
            self['end'] = obj_range[1]
            self.indices = arange(obj_range[0], obj_range[1]+1)
        else:
            self['start']=0
            self['end']=self.meds.size-1
            self.indices = arange(self.meds.size)

    def copy_psf_result(self):
        """
        copy some subset of the psf parameters
        """

        ppars=self['psf_pars']

        data=self.data
        fitter=self.psf_fitter

        res=fitter.get_result()

        data['psf_flags'][self.dindex] = res['flags']

        if 'nfev' in res:
            data['psf_nfev'][self.dindex] = res['nfev']
        elif 'numiter' in res:
            data['psf_nfev'][self.dindex] = res['numiter']

        if res['flags'] != 0:
            return

        psf_gmix=fitter.get_gmix()
        g1,g2,T=psf_gmix.get_g1g2T()

        print("    psf_id: %d psf_fwhm: %.3f g: %.3g %.3g" % (self['psf_id'],sqrt(T/2)*2.35,g1,g2) )

        if 'em' in ppars['model']:
            print("    niter: %d fdiff: %g" % (res['numiter'],res['fdiff']))
        else:
            print_pars(res['pars'],    front='    psf_pars: ')
            print_pars(res['pars_err'],front='    psf_perr: ')

        data['psf_g'][self.dindex, 0] = g1
        data['psf_g'][self.dindex, 1] = g2
        data['psf_T'][self.dindex] = T

    def get_namer(self):
        from great3.generic import Namer

        if self['use_logpars']:
            n=Namer('log')
        else:
            n=Namer()

        return n

    def copy_galaxy_result(self):
        """
        copy some subset of the psf parameters
        """
        from pprint import pprint
        n=self.get_namer()

        data=self.data
        dindex=self.dindex

        fitter=self.gal_fitter

        res=fitter.get_result()
        #pprint(res)

        if res['flags'] != 0:
            print("    galaxy fit failure")
            data['flags'][dindex] = GAL_FIT_FAILURE
            return
        else:
            data['flags'][dindex]=0

        jacob=self.boot.gal_obs.jacobian
        jrow, jcol = jacob.get_cen()
        scale = jacob.get_scale()
        row = jrow + res['pars'][0]/scale
        col = jcol + res['pars'][1]/scale

        data['pars'][dindex] = res['pars']
        data['pars_cov'][dindex] = res['pars_cov']

        data['cen_pix'][dindex] = array([row,col])

        data['g'][dindex] = res['g']
        data['g_cov'][dindex] = res['g_cov']

        data[n('flux')][dindex] = res['pars'][5]
        data[n('flux_err')][dindex] = sqrt(res['pars_cov'][5,5])
        data[n('T')][dindex] = res['pars'][4]
        data[n('T_err')][dindex] = sqrt(res['pars_cov'][4,4])

        if self['use_logpars']:
            Ts2n = 1.0/data[n('T_err')][dindex]
        else:
            Ts2n = data[n('T')][dindex]/data[n('T_err')][dindex]
        
        data['T_s2n'][dindex] = Ts2n

        data['s2n_w'][dindex] = res['s2n_w']
        data['chi2per'][dindex] = res['chi2per']
        data['dof'][dindex] = res['dof']

    def print_galaxy_result(self):
        res=self.gal_fitter.get_result()

        if 'pars' in res:
            print_pars(res['pars'],    front='    gal_pars: ')
        if 'pars_err' in res:
            print_pars(res['pars_err'],front='    gal_perr: ')

    def make_dtype(self):
        """
        make the output data type
        """

        n=self.get_namer()

        np=ngmix.gmix.get_model_npars(self['model_pars']['model'])

        dt=[
            ('number','i4'),
            ('flags','i4'),

            ('psf_id','i4'),
            ('psf_flags','i4'),
            ('psf_nfev','i4'),
            ('psf_g','f8',2),
            ('psf_T','f8'),

            ('psf_flux','f8'),
            ('psf_flux_err','f8'),

            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),

            ('cen_pix','f8',2),
            (n('flux'),'f8'),
            (n('flux_err'),'f8'),
            (n('T'),'f8'),
            (n('T_err'),'f8'),
            ('T_s2n','f8'),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),

            ('s2n_w','f8'),
            ('chi2per','f8'),
            ('dof','f8'),
           ]

        self.dtype=dt

    def make_struct(self):
        self.make_dtype()

        n=self.get_namer()

        num=self.indices.size
        data=zeros(num, dtype=self.dtype)

        data['flags'] = NO_ATTEMPT
        data['psf_g'] = DEFVAL
        data['psf_T'] = DEFVAL

        data['psf_flux'] = DEFVAL
        data['psf_flux_err'] = PDEFVAL

        data['pars'] = DEFVAL
        data['pars_cov'] = PDEFVAL

        data[n('flux')] = DEFVAL
        data[n('flux_err')] = PDEFVAL
        data[n('T')] = DEFVAL
        data[n('T_err')] = PDEFVAL
        data['g'] = DEFVAL
        data['g_cov'] = PDEFVAL

        data['s2n_w'] = DEFVAL
        data['chi2per'] = PDEFVAL

    
        self.data=data

    def set_fracdev_stuff(self):
        self['fracdev_grid']=self.get('fracdev_grid',None)
        self['fracdev_prior'] = self.get('fracdev_prior',None)


class MedsFitMax(MedsFitBase):
    def fit_galaxy(self):
        """
        fit with max like, using a MaxRunner object
        """
        self.fit_max()
        self.gal_fitter=self.boot.get_max_fitter()


    def copy_galaxy_result(self):
        """
        extra copies beyond the default
        """
        super(MedsFitMax,self).copy_galaxy_result()
        res=self.gal_fitter.get_result()
        if 'nfev' in res:
            self.data['nfev'][self.dindex] = res['nfev']

    def print_galaxy_result(self):
        super(MedsFitMax,self).print_galaxy_result()
        res=self.gal_fitter.get_result()

        if 's2n_w' in res:
            tup=(res['s2n_w'],res['nfev'],res['chi2per'])
            print("    s2n: %.1f nfev: %d chi2per: %.3f" % tup)


    def make_dtype(self):
        super(MedsFitMax,self).make_dtype()

        self.dtype += [
            ('nfev','i4'),
        ]

    def make_struct(self):
        super(MedsFitMax,self).make_struct()
        self.data['nfev'] = PDEFVAL



class CompositeMedsFitMax(MedsFitMax):
    def __init__(self, *args, **keys):
        super(CompositeMedsFitMax,self).__init__(*args, **keys)
        self.set_fracdev_stuff()

    def copy_galaxy_result(self):
        """
        extra copies beyond the default
        """
        super(CompositeMedsFitMax,self).copy_galaxy_result()
        res=self.gal_fitter.get_result()
        if 'fracdev' in res:
            self.data['fracdev'][self.dindex] = res['fracdev']
            self.data['fracdev_err'][self.dindex] = res['fracdev_err']

    def make_dtype(self):
        super(CompositeMedsFitMax,self).make_dtype()

        self.dtype += [
            ('fracdev','f8'),
            ('fracdev_err','f8'),
        ]

    def make_struct(self):
        super(CompositeMedsFitMax,self).make_struct()
        self.data['fracdev'] = PDEFVAL
        self.data['fracdev_err'] = PDEFVAL

    def get_bootstrapper(self):
        """
        get the bootstrapper for fitting psf through galaxy
        """
        from great3.sfit import get_bootstrapper
        
        # keywords pass on fracdev_prior and fracdev_grid
        boot = get_bootstrapper(self.psf_obs,
                                self.obs,
                                type='composite',
                                **self)
        return boot



class MedsFitShearBase(MedsFitBase):
    def fit_galaxy(self):
        pass

    def make_dtype(self):
        super(MedsFitShearBase,self).make_dtype()

        np=ngmix.gmix.get_model_npars(self['model_pars']['model'])
        self.dtype += [
            ('max_flags','i4'),
            ('max_pars','f8',np),
            ('max_cov','f8',(np,np)),
            ('P', 'f8'),
            ('Q', 'f8', 2),
            ('R', 'f8', (2,2)),
            ('g_sens','f8',2)
        ]

    def make_struct(self):
        super(MedsFitShearBase,self).make_struct()

        data=self.data
        data['max_flags'] = NO_ATTEMPT
        data['max_pars']  = DEFVAL
        data['max_cov']   = PDEFVAL

        data['P'] = DEFVAL
        data['Q'] = DEFVAL
        data['R'] = DEFVAL
        data['g_sens'] = DEFVAL

 
class MedsFitISample(MedsFitShearBase):
    def fit_galaxy(self):
        """
        call super to fit max like
        """
        super(MedsFitISample,self).fit_galaxy()
        self.fit_max()
        self.do_isample()

        self.add_shear_info()

        self.gal_fitter=self.boot.get_isampler()

    def do_isample(self):
        """
        run isample on the bootstrapper
        """
        ipars=self['isample_pars']
        self.boot.isample(ipars, prior=self['prior'])

    def add_shear_info(self):
        """
        add shear information based on the gal_fitter
        """

        boot=self.boot
        max_fitter=boot.get_max_fitter()
        sampler=boot.get_isampler()

        # this is the full prior
        prior=self['prior']
        g_prior=prior.g_prior

        iweights = sampler.get_iweights()
        samples = sampler.get_samples()
        g_vals=samples[:,2:2+2]

        res=sampler.get_result()

        # keep for later if we want to make plots
        self.weights=iweights

        # we are going to mutate the result dict owned by the sampler
        stats = max_fitter.get_fit_stats(res['pars'])
        res.update(stats)

        ls=ngmix.lensfit.LensfitSensitivity(g_vals,
                                            g_prior,
                                            weights=iweights,
                                            remove_prior=True)
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

    def print_galaxy_result(self):
        super(MedsFitISample,self).print_galaxy_result()
        res=self.gal_fitter.get_result()

        if 's2n_w' in res:
            tup=(res['s2n_w'],res['chi2per'])
            print("    s2n: %.1f chi2per: %.3f" % tup)

    def copy_galaxy_result(self):
        """
        extra copies beyond the default
        """
        super(MedsFitISample,self).copy_galaxy_result()
        res=self.gal_fitter.get_result()
        if 'g_sens' in res:
            data=self.data
            dindex=self.dindex
            data['g_sens'][dindex,:]   = res['g_sens']
            data['efficiency'][dindex] = res['efficiency']
            data['neff'][dindex]       = res['neff']


    def make_dtype(self):
        super(MedsFitISample,self).make_dtype()

        self.dtype += [
            ('efficiency','f4'),
            ('neff','f4'),
        ]

class CompositeMedsFitISample(MedsFitISample):

    def copy_galaxy_result(self):
        """
        extra copies beyond the default
        """
        super(CompositeMedsFitISample,self).copy_galaxy_result()
        res=self.gal_fitter.get_result()
        if 'fracdev' in res:
            data=self.data
            dindex=self.dindex
            data['fracdev'][dindex]        = res['fracdev']
            data['fracdev_noclip'][dindex] = res['fracdev_noclip']
            data['fracdev_err'][dindex]    = res['fracdev_err']
            data['TdByTe'][dindex]         = res['TdByTe']

    def make_dtype(self):
        super(CompositeMedsFitISample,self).make_dtype()

        self.dtype += [
            ('fracdev','f8'),
            ('fracdev_noclip','f8'),
            ('fracdev_err','f8'),
            ('TdByTe','f8'),
        ]

    def make_struct(self):
        super(CompositeMedsFitISample,self).make_struct()
        self.data['fracdev'] = PDEFVAL
        self.data['fracdev_noclip'] = PDEFVAL
        self.data['fracdev_err'] = PDEFVAL
        self.data['TdByTe'] = PDEFVAL

    def get_bootstrapper(self):
        """
        get the bootstrapper for fitting psf through galaxy
        """
        from great3.sfit import get_bootstrapper
        
        # keywords pass on fracdev_prior and fracdev_grid
        boot = get_bootstrapper(self.psf_obs,
                                self.obs,
                                type='composite',
                                **self)
        return boot



class PSFRunner(object):
    """
    wrapper to generate guesses and run the psf fitter a few times
    """
    def __init__(self, obs, model, Tguess, lm_pars):
        self.obs=obs

        mess="psf model should be turb or gauss,got '%s'" % model
        assert model in ['turb','gauss'],mess

        self.model=model
        self.lm_pars=lm_pars
        self.set_guess0(Tguess)

    def go(self, ntry=1):
        from ngmix.fitting import LMSimple

        for i in xrange(ntry):
            guess=self.get_guess()
            fitter=LMSimple(self.obs,self.model,lm_pars=self.lm_pars)
            fitter.go(guess)

            res=fitter.get_result()
            if res['flags']==0:
                break

        self.fitter=fitter

    def get_guess(self):
        guess=self.guess0.copy()

        guess[0:0+2] + 0.01*srandu(2)
        guess[2:2+2] + 0.1*srandu(2)
        guess[4] = guess[4]*(1.0 + 0.1*srandu())
        guess[5] = guess[5]*(1.0 + 0.1*srandu())

        return guess

    def set_guess0(self, Tguess):
        Fguess = self.obs.image.sum()
        self.guess0=array( [0.0, 0.0, 0.0, 0.0, Tguess, Fguess] )

class CompositeMaxRunner(object):
    """
    wrapper to generate guesses and run the psf fitter a few times
    """
    def __init__(self, obs, pars, guesser, fracdev, TdByTe,
                 use_logpars=True, prior=None):
        pass


class MaxRunner(object):
    """
    wrapper to generate guesses and run the psf fitter a few times
    """
    def __init__(self, obs, model, pars, guesser, use_logpars=True, prior=None):
        self.obs=obs

        self.pars=pars
        self.method=pars['method']
        if self.method == 'lm':
            self.send_pars=pars['lm_pars']
        else:
            self.send_pars=pars

        mess="model should be exp or dev,got '%s'" % model
        assert model in ['exp','dev'],mess

        self.model=model
        self.use_logpars=use_logpars
        self.prior=prior

        self.bestof = pars.get('bestof',1)

        self.guesser=guesser

    def go(self, ntry=1):
        if self.method=='lm':
            method=self._go_lm
        elif self.method=='nm':
            method=self._go_nm
        else:
            raise ValueError("bad method '%s'" % self.method)

        lnprob_max=-numpy.inf
        fitter_best=None
        for i in xrange(self.bestof):
            method(ntry=ntry)

            res=self.fitter.get_result()
            if res['flags']==0:
                if res['lnprob'] > lnprob_max:
                    lnprob_max = res['lnprob']
                    fitter_best=self.fitter
        
        if fitter_best is not None:
            self.fitter=fitter_best

    def _go_lm(self, ntry=1):
        from ngmix.fitting import LMSimple

        for i in xrange(ntry):
            guess=self.guesser()
            fitter=LMSimple(self.obs,
                            self.model,
                            lm_pars=self.send_pars,
                            use_logpars=self.use_logpars,
                            prior=self.prior)

            fitter.go(guess)

            res=fitter.get_result()
            if res['flags']==0:
                break

        self.fitter=fitter

    def _go_nm(self, ntry=1):
        from ngmix.fitting import MaxSimple

        for i in xrange(ntry):
            guess=self.guesser()
            fitter=MaxSimple(self.obs, self.model,
                             method='Nelder-Mead',
                             use_logpars=self.use_logpars,
                             prior=self.prior)

            fitter.run_max(guess, **self.send_pars)

            res=fitter.get_result()
            if res['flags']==0:
                break

        self.fitter=fitter


