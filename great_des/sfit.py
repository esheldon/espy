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

from gmix_meds.util import FromPSFGuesser, FixedParsGuesser, FromParsGuesser, FromFullParsGuesser

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
        sigma_guess=self['psf_fwhm_guess']/2.35
        self['psf_Tguess'] = 2*sigma_guess**2

    def get_data(self):
        """
        get a reference to the data structure
        """
        return self.data

    def do_fits(self):
        """
        loop and fit all objects
        """

        last=self.indices.size-1
        for dindex,mindex in enumerate(self.indices):
            self.dindex=dindex
            self.mindex=mindex

            print("%d:%d  %d:%d" % (dindex, last, self['start'],self['end']))

            self.data['number'][dindex] = mindex

            self.make_psf_observation()
            self.make_galaxy_observation()

            psf_flags=self.fit_psf()

            if psf_flags != 0:
                continue

            self.fit_psf_flux()

            self.fit_galaxy()

            self.copy_galaxy_result()
            self.print_galaxy_result()

    def fit_psf(self):
        """
        load the observation and do the fit
        """

        psf_runner=PSFRunner(self.psf_obs,
                             self['psf_model'],
                             self['psf_Tguess'],
                             self['psf_max_pars'])

        psf_runner.go(ntry=self['psf_ntry'])

        self.psf_fitter = psf_runner.fitter

        res=self.psf_fitter.get_result()

        if res['flags']==0:
            #print("    psf success, setting psf in obs")
            gmix=self.psf_fitter.get_gmix()
            self.psf_obs.set_gmix(gmix)
            self.obs.set_psf(self.psf_obs)

            self.copy_psf_result()
        else:
            print("    psf fitting failed")
            self.data['flags'][dindex] = PSF_FIT_FAILURE

        return res['flags']

    def fit_psf_flux(self):
        """
        use psf as a template, measure flux (linear)
        """
        fitter=ngmix.fitting.TemplateFluxFitter(self.obs, do_psf=True)
        fitter.go()

        res=fitter.get_result()

        data=self.data
        data['psf_flux'][self.dindex] = res['flux']
        data['psf_flux_err'][self.dindex] = res['flux_err']

        print("    psf flux: %.3g +/- %.3g" % (res['flux'],res['flux_err']))

        self.psf_flux_fitter=fitter

    def make_psf_observation(self):
        """
        read the image and weight data
        """

        ext=self.truth['id_psf'][self.mindex]

        self['psf_id'] = ext

        image0 = self.psf_obj[ext][:,:]

        image = image0 + numpy.random.normal(scale=self['psf_addnoise'],
                                             size=image0.shape)

        weight = image.copy()
        weight *= 0
        weight += 1.0/self['psf_addnoise']**2

        jacob=self.get_jacobian()
        self.psf_obs = Observation(image, weight=weight, jacobian=jacob)

    def make_galaxy_observation(self):
        """
        read the image and weight data
        """

        image = self.meds.get_cutout(self.mindex, 0)
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
        data=self.data
        T=data['psf_T'][self.dindex]
        flux=data['psf_flux'][self.dindex]

        if self['use_logpars']:
            scaling='log'
        else:
            scaling='linear'

        guesser=FromPSFGuesser(T, flux, scaling=scaling)

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

    def set_indices(self):
        """
        this version we don't support work dir
        """

        obj_range = self.get('obj_range',None)
        self['start'] = obj_range[0]
        self['end'] = obj_range[1]

        if obj_range is not None:
            self.indices = arange(obj_range[0], obj_range[1]+1)
        else:
            self.indices = arange(self.meds.size)

    def copy_psf_result(self):
        """
        copy some subset of the psf parameters
        """

        data=self.data
        fitter=self.psf_fitter

        res=fitter.get_result()

        data['psf_flags'][self.dindex] = res['flags']

        if 'nfev' in res:
            data['psf_nfev'][self.dindex] = res['nfev']

        if res['flags'] != 0:
            return

        psf_gmix=fitter.get_gmix()
        g1,g2,T=psf_gmix.get_g1g2T()

        print("    psf_id: %d psf_fwhm: %.3f g: %.3g,%.3g" % (self['psf_id'],sqrt(T/2)*2.35,g1,g2) )
        print_pars(res['pars'],    front='    psf_pars: ')
        print_pars(res['pars_err'],front='    psf_perr: ')

        data['psf_g'][self.dindex, 0] = g1
        data['psf_g'][self.dindex, 1] = g2
        data['psf_T'][self.dindex] = T

    def get_namer(self):
        from gmix_meds.util import Namer

        if self['use_logpars']:
            n=Namer('log')
        else:
            n=Namer()

        return n

    def copy_galaxy_result(self):
        """
        copy some subset of the psf parameters
        """

        n=self.get_namer()

        data=self.data
        dindex=self.dindex

        fitter=self.gal_fitter

        res=fitter.get_result()

        if res['flags'] != 0:
            print("    galaxy fit failure")
            data['flags'][dindex] = GAL_FIT_FAILURE
            return

        data['pars'][dindex] = res['pars']
        data['pars_cov'][dindex] = res['pars_cov']

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
            print_pars(res['pars_err'],front='    gal_perr: ')

    def make_dtype(self):
        """
        make the output data type
        """

        n=self.get_namer()

        np=ngmix.gmix.get_model_npars(self['fit_model'])

        dt=[
            ('number','i4'),
            ('flags','i4'),

            ('psf_flags','i4'),
            ('psf_nfev','i4'),
            ('psf_g','f8',2),
            ('psf_T','f8'),

            ('psf_flux','f8'),
            ('psf_flux_err','f8'),

            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),

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


class MedsFitMax(MedsFitBase):
    def fit_galaxy(self):
        """
        fit with max like, using a MaxRunner object
        """
        data=self.data

        max_pars=self['max_pars']

        guesser=self.get_psf_guesser()
        runner=MaxRunner(self.obs,
                         self['fit_model'],
                         max_pars['lm_pars'],
                         guesser,
                         use_logpars=self['use_logpars'],
                         prior=self['prior'])

        runner.go(ntry=max_pars['ntry'])

        self.gal_fitter=runner.fitter

        res=self.gal_fitter.get_result()
        if res['flags']==0:
            self.try_replace_cov(self.gal_fitter)

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

class MedsFitShearBase(MedsFitBase):
    def fit_galaxy(self):
        pass

    def make_dtype(self):
        super(MedsFitShearBase,self).make_dtype()

        np=ngmix.gmix.get_model_npars(self['fit_model'])
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
    def make_dtype(self):
        super(MedsFitIsample,self).make_dtype()

        self.dtype += [
            ('efficiency','f4'),
            ('neff','f4'),
        ]

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


class MaxRunner(object):
    """
    wrapper to generate guesses and run the psf fitter a few times
    """
    def __init__(self, obs, model, lm_pars, guesser, use_logpars=True, prior=None):
        self.obs=obs

        mess="model should be exp or dev,got '%s'" % model
        assert model in ['exp','dev'],mess

        self.model=model
        self.lm_pars=lm_pars
        self.use_logpars=use_logpars
        self.prior=prior

        self.guesser=guesser

    def go(self, ntry=1):
        from ngmix.fitting import LMSimple

        for i in xrange(ntry):
            guess=self.guesser()
            fitter=LMSimple(self.obs,
                            self.model,
                            lm_pars=self.lm_pars,
                            use_logpars=self.use_logpars,
                            prior=self.prior)

            fitter.go(guess)

            res=fitter.get_result()
            if res['flags']==0:
                break

        self.fitter=fitter


