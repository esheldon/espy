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

class StackSim(dict):
    def __init__(self, run, is2, is2n, itrial):
        conf=read_config(run)
        self.update(conf)

        self.simc = read_config(self['sim'])

        self.is2=is2
        self.is2n=is2n
        self.itrial=itrial

        self.s2n = shapesim.get_s2n(self, is2n)
        self.s2 = numpy.linspace(self.simc['mins2'],
                                 self.simc['maxs2'], 
                                 self.simc['nums2'])[is2]
        self.sratio=numpy.sqrt(1./self.s2)
        self.npair=self.get_npair()

        self._set_gprior()

        simpars=self.get('simpars',{})
        self.shapesim = ShapeSim(self['sim'], **simpars)


    def get_npair(self):
        s2n_fac = self['s2n_fac']
        nellip = shapesim.get_s2n_nrepeat(self.s2n, fac=s2n_fac)

        if nellip < self['min_gcount']:
            nellip=self['min_gcount']
        return nellip


    def go(self):

        gsum=numpy.zeros(2)
        gsum_means=numpy.zeros(2)
        nsum_means=numpy.zeros(2)

        wsum=numpy.zeros(2)
        nsum=0

        gvals=self.gprior.sample1d(self.npair)

        for i,g in enumerate(gvals):
            imd1,imd2=self._make_ring_pair(g)
            if (((i+1) % 100) == 0):
                print '%d/%d' % (i+1,self.npair)

            if i==0:
                image_stack=0*imd1['image'].copy()
                psf_stack=0*imd1['psf_image'].copy()
                image_skyvar=0.0
                psf_skyvar=0.0

            image_stack += imd1['image']
            image_stack += imd2['image']

            image_skyvar += imd1['image_skyvar']
            image_skyvar += imd2['image_skyvar']

            psf_stack += imd1['psf_image']
            psf_stack += imd2['psf_image']

            psf_skyvar += imd1['psf_skyvar']
            psf_skyvar += imd2['psf_skyvar']
            

        self._image_stack=image_stack
        self._image_skyvar=image_skyvar
        self._psf_stack=psf_stack
        self._psf_skyvar=psf_skyvar

        if False:
            import images
            images.multiview(image_stack)

        self._fit_stacks()
        self._set_result()

    def write(self):
        shapesim.write_output(self['run'], self.is2, self.is2n, self._data, 
                              itrial=self.itrial)

    def _set_result(self):

        pars     = self._fitter.get_pars()
        pcov     = self._fitter.get_pcov()
        psf_pars = self._psf_fitter.get_pars()
        psf_pcov = self._psf_fitter.get_pcov()

        res=self._fitter.get_result()

        s2n=self._fitter.get_weighted_s2n()
        psf_s2n=self._psf_fitter.get_weighted_s2n()

        data=self._get_struct(pars.size, psf_pars.size)

        data['model'] = self['fitmodel']
        data['npair'] = self.npair
        data['nimage'] = self.npair*2
        data['s2n'] = self.s2n
        data['s2'] = self.s2
        data['sratio'] = self.sratio

        data['shear_true'] = self.simc['shear']


        data['e'] = pars[2:2+2]
        data['ecov'] = pcov[2:2+2, 2:2+2]
        data['pars'] = pars
        data['pcov'] = pcov
        data['Tmean'] = pars[4]
        data['Terr'] = numpy.sqrt(pcov[4,4])
        data['Ts2n'] = data['Tmean']/data['Terr']

        data['psf_pars'] = psf_pars
        data['psf_pcov'] = psf_pcov

        data['loglike'] = res['loglike']
        data['chi2per'] = res['chi2per']
        data['dof'] = res['dof']
        data['fit_prob'] = res['fit_prob']

        data['s2n_stack'] = s2n
        data['image_stack'] = self._image_stack
        data['image_skyvar'] = self._image_skyvar

        data['psf_s2n_stack'] = psf_s2n
        data['psf_stack'] = self._psf_stack
        data['psf_skyvar'] = self._psf_skyvar

        self._data=data


    def _get_struct(self, npars, psf_npars):
        shape=self._image_stack.shape
        psf_shape=self._psf_stack.shape

        dt=[('model','S5'),
            ('npair','i4'),
            ('nimage','i4'),
            ('s2n','f8'),  # per object
            ('s2','f8'),
            ('sratio','f8'),
            ('shear_true','f8',2),

            ('e','f8',2),
            ('ecov','f8',(2,2)),
            ('Tmean','f8'),
            ('Terr','f8'),
            ('Ts2n','f8'),

            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('psf_pars','f8',psf_npars),
            ('psf_pcov','f8',(psf_npars,psf_npars)),

            ('loglike','f8'),     # loglike of fit
            ('chi2per','f8'),     # chi^2/degree of freedom
            ('dof','i4'),         # degrees of freedom
            ('fit_prob','f8'),    # probability of the fit happening randomly

            ('s2n_stack','f8'),
            ('image_stack','f8',shape),
            ('image_skyvar','f8'),

            ('psf_s2n_stack','f8'),
            ('psf_stack','f8',psf_shape),
            ('psf_skyvar','f8')
           ]

        data=numpy.zeros(1, dtype=dt)
        return data

    def _fit_stacks(self):

        ares = self._run_admom(self._image_stack, self._image_skyvar)
        psf_ares = self._run_admom(self._psf_stack, self._psf_skyvar)

        psf_fitter = self._fit_stack(self._psf_stack, self._psf_skyvar, psf_ares, name='PSF')

        psf_gmix=psf_fitter.get_gmix()

        fitter=self._fit_stack(self._image_stack, self._image_skyvar, ares, psf=psf_gmix, name='Images')

        self._fitter=fitter
        self._psf_fitter=psf_fitter

    def _fit_stack(self, image_stack, skyvar, ares, psf=None, name=''):

        row_guess=ares['wrow'] 
        col_guess=ares['wcol'] 
        e1guess=ares['e1']
        e2guess=ares['e2']
        Tguess=ares['Irr']+ares['Icc']
        counts=image_stack.sum()

        while True:
            if False:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*8.0*(1.+0.1*srandu()),
                                    Tguess*2.42*(1.+0.1*srandu()),
                                    Tguess*0.20*(1.+0.1*srandu()),
                                    counts*0.28*(1.+0.1*srandu()), 
                                    counts*0.34*(1.+0.1*srandu()), 
                                    counts*0.38*(1.+0.1*srandu())])
            elif True:
                # appropriate for stacking gaussians
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*3.0*(1.+0.1*srandu()),
                                    Tguess*0.9*(1.+0.1*srandu()),
                                    counts*0.6*(1.+0.1*srandu()), 
                                    counts*0.4*(1.+0.1*srandu())])

            elif False:
                prior=numpy.array( [row_guess+0.01*srandu(),
                                    col_guess+0.01*srandu(),
                                    e1guess+0.05*srandu(),
                                    e2guess+0.05*srandu(),
                                    Tguess*(1.+0.1*srandu()),
                                    counts*(1.+0.1*srandu())] )

            width=numpy.abs(prior)*1.e6

            fitter=gmix_image.gmix_fit.GMixFitCoellip(image_stack, sqrt(skyvar), prior, width, 
                                                      psf=psf, verbose=False)
            flags=fitter.get_flags()

            if flags==0:
                break
            else:
                print flags

        res=fitter.get_result()
        print_pars(prior,front='prior: ')
        print_pars(res['pars'], front='pars: ')
        print_pars(res['perr'], front='perr: ')

        model=fitter.get_model()
        #images.multiview(im/im.max(), nonlinear=1)
        if False:
            import images
            images.compare_images(image_stack/image_stack.max(), 
                                  model/model.max(),
                                  nonlinear=1,
                                  title=name)

        return fitter


    def _run_admom(self, image_stack, skyvar):
        import admom

        cen_guess=numpy.array(image_stack.shape,dtype='f8')/2.
        while True:
            ares = admom.admom(image_stack, 
                               cen_guess[0], 
                               cen_guess[1],
                               sigsky=sqrt(skyvar),
                               guess=4.0,
                               nsub=1)

            if ares['whyflag'] == 0:
                break
        return ares 

    def _make_ring_pair(self, g):
        s2n_psf = self['s2n_psf']

        self.image_list=[]
        self.ivar_list=[]
        self.psf_list=[]
        self.model_list=[]
        self.ares_list=[]
        
        theta = 360.0*numpy.random.random()
        theta2 = theta + 90.0

        ellip=lensing.util.g2e(g)

        ci,ares,psf = self._get_ci_ares_psf(self.s2, ellip, theta, self.s2n, s2n_psf)
        ci2,ares2,psf2 = self._get_ci_ares_psf(self.s2, ellip, theta2, self.s2n, s2n_psf)

        imd1={'image':ci.image,
              'psf_image':ci.psf,
              'image_skyvar':ci['skysig']**2,
              'ivar':1./ci['skysig']**2,
              'psf_skyvar':ci['skysig_psf']**2,
              'psf_ivar':1./ci['skysig_psf']**2,
              'ares':ares,
              'psf':psf,
              'model':'gauss'}
        imd2={'image':ci2.image,
              'psf_image':ci2.psf,
              'image_skyvar':ci2['skysig']**2,
              'ivar':1./ci2['skysig']**2,
              'psf_skyvar':ci2['skysig_psf']**2,
              'psf_ivar':1./ci2['skysig_psf']**2,
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

