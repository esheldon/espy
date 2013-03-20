"""
Generate image simulations and process them with the
deswl pipeline
"""

import os
from sys import stderr
import numpy
from numpy import zeros, sqrt, tanh, arctanh, random

import fimage
from fimage import mom2sigma
import deswl
from esutil.numpy_util import ahelp
from esutil.misc import wlog

from lensing.util import e2gamma, gamma2e, e1e2_to_g1g2, g1g2_to_e1e2
from . import shapesim

import images

SHAPELETS_EDGE = 0x40

class DESWLSim(shapesim.BaseSim):
    def __init__(self, run):
        super(DESWLSim,self).__init__(run)

        cf = os.environ['WL_DIR']
        cf = os.path.join(cf,'etc','wl.config')
        if not os.path.isfile(cf):
            raise ValueError("config file not found: "+cf)

        self.config_fname = cf


    def run(self, ci):
        """
        Run the psf and object through deswl
        """
        dt = [('flags','i8'),
              ('sigma_psf','f8'),
              ('gal_prepsf_sigma','f8'),
              ('sigma0','f8'),
              ('s2','f8'),
              ('nu','f8'),
              ('gamma1','f8'),
              ('gamma2','f8'),
              ('gamma','f8'),
              ('gcov11','f8'),
              ('gcov12','f8'),
              ('gcov22','f8')]
        out=zeros(1, dtype=dt)
        out['sigma_psf'] = -9999
        out['gal_prepsf_sigma'] = -9999
        out['sigma0'] = -9999
        out['s2'] = 9999
        out['nu'] = -9999
        out['gamma1'] = -9999
        out['gamma2'] = -9999
        out['gcov11'] = 9999
        out['gcov12'] = 9999
        out['gcov22'] = 9999

        sky=0.0
        skysig=ci['skysig']
        skysig_psf=ci['skysig_psf']

        psf_sigma_guess=\
            fimage.mom2sigma(ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2])
        # randomizing might help with the gamma==0 bug?
        psf_sigma_guess += 0.1*random.random()
        sigma0_guess=\
            fimage.mom2sigma(ci['cov_admom'][0]+ci['cov_admom'][2])


        sigma_obj = fimage.mom2sigma(ci['cov_uw'][0]+ci['cov_uw'][2])

        #psf_max_aperture = self['maxaper_nsig']*psf_sigma_guess
        #shear_max_aperture = self['maxaper_nsig']*sigma_obj
        psf_max_aperture = ci.psf.shape[0]
        shear_max_aperture = ci.image.shape[0]

        wlpsfobj = deswl.cwl.WLObject(ci.psf,
                                      float(ci['cen'][0]), 
                                      float(ci['cen'][1]),
                                      float(sky),
                                      float(skysig_psf**2),
                                      float(psf_max_aperture), 
                                      psf_sigma_guess)
        out['flags'] = wlpsfobj.get_flags()
        if out['flags'] != 0:
            self.wlog('psf wlobj flags:',out['flags'])
            return out
        out['sigma_psf'] = wlpsfobj.get_sigma0()

        wlobj = deswl.cwl.WLObject(ci.image,
                                   float(ci['cen'][0]), 
                                   float(ci['cen'][1]),
                                   float(sky),
                                   float(skysig**2),
                                   float(shear_max_aperture), 
                                   sigma0_guess)
        out['flags'] = wlobj.get_flags()
        if out['flags'] != 0:
            self.wlog('wlobj flags:',out['flags'])
            return out
        out['sigma0'] = wlobj.get_sigma0()

        wlshear = deswl.cwl.WLShear(self.config_fname,
                                    wlobj,
                                    wlpsfobj,
                                    self['psf_order'], self['gal_order'],
                                    self['maxaper_nsig'])

        out['flags'] = wlshear.get_flags()
        if out['flags'] != 0:
            self.wlog('wlshear flags:',out['flags'])
            return out

        out['s2'] = out['sigma_psf']**2/(out['sigma0']**2-out['sigma_psf']**2)
        out['gamma1'] = wlshear.get_shear1()
        if out['gamma1'] == 0:
            self.wlog("found gamma==0 bug")
            out['flags'] = -2
            return out

        out['gamma2'] = wlshear.get_shear2()
        out['gamma'] = sqrt(out['gamma1']**2 + out['gamma2']**2)

        # usually not useful
        out['gal_prepsf_sigma'] = wlshear.get_prepsf_sigma() 

        out['nu'] = wlshear.get_nu()
        out['gcov11'] = wlshear.get_cov11()
        out['gcov12'] = wlshear.get_cov12()
        out['gcov22'] = wlshear.get_cov22()
        return out


    def copy_output(self, s2, ellip, s2n, ci, res):
        st = numpy.zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n'] = s2n
        st['ellip'] = ellip
        st['e1true'] = ci['e1true']
        st['e2true'] = ci['e2true']
        st['etrue']  = ci['etrue']
        st['e1_uw'] = ci['e1_image0_uw']
        st['e2_uw'] = ci['e2_image0_uw']
        st['e_uw']  = ci['e_image0_uw']
        st['gamma'] = e2gamma(st['etrue'])
        st['gamma1'],st['gamma2'] = e1e2_to_g1g2(st['e1true'],st['e2true'])

        size2psf = ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2]
        size2obj = ci['cov_image0_uw'][0]+ci['cov_image0_uw'][2]
        st['s2noweight'] = size2psf/size2obj

        s2psf_am = ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2]
        s2obj_am = ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2]
        st['s2admom'] = s2psf_am/s2obj_am
        st['sigma_psf_admom'] = \
            mom2sigma(ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2])
        st['sigma_admom'] = \
            mom2sigma(ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2])
        st['sigma0_admom'] = \
            mom2sigma(ci['cov_admom'][0]+ci['cov_admom'][2])


        st['sigma_psf_meas'] = res['sigma_psf']
        #st['sigma_meas'] = res['sigma']
        st['gal_prepsf_sigma_meas'] = res['gal_prepsf_sigma']
        st['sigma0_meas'] = res['sigma0']
        #st['s2_meas'] = res['sigma_psf']**2/(res['sigma0']**2-res['sigma_psf']**2)
        st['s2_meas'] = res['s2']
        st['gamma1_meas'] = res['gamma1']
        st['gamma2_meas'] = res['gamma2']
        st['gamma_meas'] = res['gamma']
        e = gamma2e(res['gamma'])
        e1,e2 = g1g2_to_e1e2(res['gamma1'],res['gamma2'])
        st['e1_meas'] = e1
        st['e2_meas'] = e2
        st['e_meas'] = e

        st['s2n_meas'] = res['nu']
        st['gcov11'] = res['gcov11']
        st['gcov12'] = res['gcov12']
        st['gcov22'] = res['gcov22']
        st['flags'] = res['flags']
        return st

    def out_dtype(self):
        dt = [('s2n','f8'),
              ('ellip','f8'),
              
              ('s2','f8'),         # requested (spsf/sobj)**2
              ('s2noweight','f8'), # unweighted s2 of object before noise
              ('sigma_psf_admom','f8'),
              ('sigma_admom','f8'),
              ('sigma0_admom','f8'),
              ('s2admom','f8'),    # s2 from admom, generally different
              ('etrue','f8'),
              ('e1true','f8'),
              ('e2true','f8'),
              ('e_uw','f8'),
              ('e1_uw','f8'),
              ('e2_uw','f8'),
              ('gamma','f8'),
              ('gamma1','f8'),
              ('gamma2','f8'),

              ('s2_meas','f8'),
              ('s2n_meas','f8'),    # same as nu
              ('sigma_psf_meas','f8'),
              ('sigma0_meas','f8'),
              ('gal_prepsf_sigma_meas','f8'),
              ('gamma_meas','f8'),
              ('gamma1_meas','f8'),
              ('gamma2_meas','f8'),
              ('e_meas','f8'),
              ('e1_meas','f8'),
              ('e2_meas','f8'),
              ('gcov11','f8'),
              ('gcov12','f8'),
              ('gcov22','f8'),
              ('flags','i8')]
        return dt

