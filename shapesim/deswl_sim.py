"""
Generate image simulations and process them with the
deswl pipeline
"""

import numpy
from numpy import zeros, sqrt, tanh, arctanh

import fimage
from fimage import mom2sigma
import deswl
from esutil.numpy_util import ahelp
from esutil.misc import wlog
from . import shapesim

def e2gamma(e):
    return tanh(0.5*arctanh(e))
def e1e2_to_g1g2(e1, e2):
    e = sqrt(e1**2 + e2**2)
    g = e2gamma(e)
    fac = g/e
    g1, g2 = fac*e1, fac*e2
    return g1,g2


class DESWLSim(dict):
    def __init__(self, run, verbose=False):
        conf=shapesim.read_config('deswl',run)
        for k,v in conf.iteritems():
            self[k] = v
        self.verbose = verbose
 
    def process_trials(self, is2, ie):
        """
        Generate random realizations of a particular element in the s2 and
        ellip sequences.

            is2 is a number between 0 and self['nums2']-1
            ie  is a number between 0 and self['nume']-1
        """
        
        out = numpy.zeros(self['ntrial'], dtype=self.out_dtype())
        ss = shapesim.ShapeSim(self['sim'])

        s2n = self['s2n']
        s2,ellip = self.get_s2_e(is2, ie)

        for i in xrange(self['ntrial']):
            iter=0
            while iter < self['itmax']:
                ci=ss.get_trial(s2,ellip,s2n)
                res = self.run_wl(ci)
                if res['flags'][0] == 0:
                    st = self.copy_output(s2, ellip, s2n, ci, res)

                    #ahelp(res)
                    out[i] = st
                    #ahelp(out[i])
                    break
                else:
                    iter += 1
            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        return out

    def run_wl(self, ci):
        """
        Run the psf and object through deswl
        """
        dt = [('flags','i8'),
              ('sigma_psf','f8'),
              ('sigma','f8'),
              ('s2','f8'),
              ('nu','f8'),
              ('gamma1','f8'),
              ('gamma2','f8'),
              ('gamma','f8'),
              ('gcov11','f8'),
              ('gcov12','f8'),
              ('gcov22','f8')]
        out=zeros(1, dtype=dt)

        sky=0.0
        skysig=ci['skysig']
        if skysig == 0:
            skysig=1
        psf_sigma_guess=\
            fimage.mom2sigma(ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2])
        psf_aperture = 4*psf_sigma_guess
        sigma_obj = fimage.mom2sigma(ci['cov_uw'][0]+ci['cov_uw'][2])
        shear_aperture = 4.*sigma_obj

        wlq = deswl.cwl.WLQuick(self['psf_order'],
                                self['gal_order'])

        wlq.set_psf(ci.psf, 
                    float(ci['cen'][0]),
                    float(ci['cen'][1]), 
                    float(sky), 
                    float(psf_aperture))
        wlq.set_image(ci.image,
                      float(ci['cen'][0]), 
                      float(ci['cen'][1]),
                      float(sky),
                      float(skysig**2),
                      float(shear_aperture))

        out['flags'] += wlq.calculate_psf_sigma(psf_sigma_guess)
        if out['flags'] != 0:
            wlog('psf shapelets flags:',out['flags'])
            return out
        out['sigma_psf'] = wlq.get_psf_sigma()

        out['flags'] += wlq.calculate_psf_shapelets()
        if out['flags'] != 0:
            wlog('psf shapelets flags:',out['flags'])
            return out

        out['flags'] += wlq.calculate_shear()
        if out['flags'] != 0:
            wlog('shear flags:',out['flags'])
            return out

        out['sigma'] = wlq.get_sigma()
        out['s2'] = out['sigma_psf']/out['sigma']
        out['nu'] = wlq.get_nu()
        out['gamma1'] = wlq.get_shear1()
        out['gamma2'] = wlq.get_shear2()
        out['gamma'] = sqrt(out['gamma1']**2 + out['gamma2']**2)
        out['gcov11'] = wlq.get_cov11()
        out['gcov12'] = wlq.get_cov12()
        out['gcov22'] = wlq.get_cov22()
        return out

    def get_s2_e(self, is2, ie):
        self.check_is2_ie(is2, ie)
        s2 = numpy.linspace(self['mins2'],self['maxs2'], self['nums2'])[is2]
        ellip = numpy.linspace(self['mine'],self['maxe'], self['nume'])[ie]

        return s2, ellip

    def check_is2_ie(self, is2, ie):
        max_is2 = self['nums2']-1
        max_ie  = self['nume']-1
        if (is2 < 0) or (is2 > max_is2):
            raise ValueError("is2 must be within [0,%d], "
                             "got %d" % (max_is2,is2))
        if (ie < 0) or (ie > max_ie):
            raise ValueError("ie must be within [0,%d], "
                             "got %d" % (max_ie,ie))


    def copy_output(self, s2, ellip, s2n, ci, res):
        st = numpy.zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n'] = s2n
        st['ellip'] = ellip
        st['e1true'] = ci['e1true']
        st['e2true'] = ci['e2true']
        st['etrue']  = ci['etrue']
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
        st['sigma_tot_admom'] = \
            mom2sigma(ci['cov_admom'][0]+ci['cov_admom'][2])


        st['sigma_psf_meas'] = res['sigma_psf']
        st['sigma_meas'] = res['sigma']
        # this gives almost exactly 0.5 of what it should be
        st['s2_meas'] = res['sigma_psf']**2/(res['sigma']**2-res['sigma_psf']**2)
        #st['s2_meas'] = res['sigma_psf']**2/res['sigma']**2
        st['gamma1_meas'] = res['gamma1']
        st['gamma2_meas'] = res['gamma2']
        st['gamma_meas'] = res['gamma']
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
              ('sigma_tot_admom','f8'),
              ('s2admom','f8'),    # s2 from admom, generally different
              ('etrue','f8'),
              ('e1true','f8'),
              ('e2true','f8'),
              ('gamma','f8'),
              ('gamma1','f8'),
              ('gamma2','f8'),

              ('s2_meas','f8'),
              ('s2n_meas','f8'),    # same as nu
              ('sigma_psf_meas','f8'),
              ('sigma_meas','f8'),
              ('gamma_meas','f8'),
              ('gamma1_meas','f8'),
              ('gamma2_meas','f8'),
              ('gcov11','f8'),
              ('gcov12','f8'),
              ('gcov22','f8'),
              ('flags','i8')]
        return dt

