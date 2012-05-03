"""
Generate image simulations and process them with the
deswl pipeline
"""

from sys import stderr
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
        conf=shapesim.read_config(run)
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
        
        maxwrite_ci=1
        nwrite_ci=0

        out = numpy.zeros(self['ntrial'], dtype=self.out_dtype())
        ss = shapesim.ShapeSim(self['sim'])

        s2n = self['s2n']
        s2,ellip = self.get_s2_e(is2, ie)

        for i in xrange(self['ntrial']):
            stderr.write(".")
            iter=0
            while iter < self['itmax']:
                ci=ss.get_trial(s2,ellip,s2n)
                res = self.run_wl(ci)
                if res['flags'][0] == 0:
                    st = self.copy_output(s2, ellip, s2n, ci, res)
                    out[i] = st
                    #self.write_ci(ci, is2, ie, res=st)
                    break
                else:
                    if res['flags'][0] == -2 and nwrite_ci < maxwrite_ci:
                        self.write_ci(ci, is2, ie,error=True)
                        nwrite_ci += 1
                    iter += 1
            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        stderr.write("\n")
        shapesim.write_output(self['run'], is2, ie, out)
        return out

    def write_ci(self, ci, is2, ie, error=False, res=None):
        """
        Write the ci to a file in the outputs directory
        """
        import fitsio
        import tempfile
        if error:
            suffix='-err'
        else:
            suffix='-good'
        rand=tempfile.mktemp(dir='',suffix=suffix)
        url=shapesim.get_output_url(self['run'], is2, ie)
        url = url.replace('.rec','-'+rand+'.fits')
        h = {}
        for k,v in self.iteritems():
            h[k] = v
        for k,v in ci.iteritems():
            h[k] = v

        with fitsio.FITS(url, mode='rw', clobber=True) as fobj:
            fobj.write(ci.image, header=h, extname='image')
            fobj.write(ci.psf, extname='psf')
            fobj.write(ci.image0, extname='image0')

            if res is not None:
                fobj.write(res, extname='results')


    def run_wl(self, ci):
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
        if self['s2n'] <= 0:
            skysig=0.001
        else:
            skysig=ci['skysig']

        psf_sigma_guess=\
            fimage.mom2sigma(ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2])
        sigma0_guess=\
            fimage.mom2sigma(ci['cov_uw'][0]+ci['cov_uw'][2])
        psf_aperture = 4*psf_sigma_guess
        sigma_obj = fimage.mom2sigma(ci['cov_uw'][0]+ci['cov_uw'][2])
        shear_aperture = 4.*sigma_obj

        wlpsfobj = deswl.cwl.WLObject(ci.psf,
                                      float(ci['cen'][0]), 
                                      float(ci['cen'][1]),
                                      float(sky),
                                      float(1),
                                      float(psf_aperture), 
                                      psf_sigma_guess)
        out['flags'] = wlpsfobj.get_flags()
        if out['flags'] != 0:
            wlog('psf wlobj flags:',out['flags'])
            return out
        out['sigma_psf'] = wlpsfobj.get_sigma0()

        wlobj = deswl.cwl.WLObject(ci.image,
                                   float(ci['cen'][0]), 
                                   float(ci['cen'][1]),
                                   float(sky),
                                   float(skysig**2),
                                   float(shear_aperture), 
                                   sigma0_guess)
        out['flags'] = wlobj.get_flags()
        if out['flags'] != 0:
            wlog('wlobj flags:',out['flags'])
            return out
        out['sigma0'] = wlobj.get_sigma0()

        wlshear = deswl.cwl.WLShear(wlobj,wlpsfobj,
                                    self['psf_order'], self['gal_order'])

        out['flags'] = wlshear.get_flags()
        if out['flags'] != 0:
            wlog('wlshear flags:',out['flags'])
            return out

        out['s2'] = out['sigma_psf']**2/(out['sigma0']**2-out['sigma_psf']**2)
        out['gamma1'] = wlshear.get_shear1()
        if out['gamma1'] == 0:
            wlog("found gamma==0 bug")
            out['flags'] = -2
            return out

        out['gamma2'] = wlshear.get_shear2()
        out['gamma'] = sqrt(out['gamma1']**2 + out['gamma2']**2)

        # usually not useful
        out['gal_prepsf_sigma'] = wlshear.get_prepsf_sigma() 

        if self['s2n'] > 0:
            out['nu'] = wlshear.get_nu()
            out['gcov11'] = wlshear.get_cov11()
            out['gcov12'] = wlshear.get_cov12()
            out['gcov22'] = wlshear.get_cov22()
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
              ('gcov11','f8'),
              ('gcov12','f8'),
              ('gcov22','f8'),
              ('flags','i8')]
        return dt

