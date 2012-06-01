"""
Generate image simulations and process them with the
gmix fitting pipeline
"""

from numpy import random, zeros, sqrt, array, ceil
import sys
from sys import stderr
from lensing.util import e2gamma, e1e2_to_g1g2
from . import shapesim
from fimage import mom2sigma, cov2sigma
from pprint import pprint 
import copy
import images

try:
    import gmix_image
    from gmix_image import pars2gmix_coellip
except ImportError:
    stderr.write("could not import gmix_image")

try:
    import admom
except ImportError:
    stderr.write("could not import admom")

class GMixFitSim(shapesim.BaseSim):
    """
    We only override

        .run()
        .out_dtype()
        .copy_output()

    """
    def __init__(self, run):
        super(GMixFitSim,self).__init__(run)
        if 'verbose' not in self:
            self['verbose'] = False

    def run(self, ci):
        """
        Process the input convolved image

        Output will be a dict with
        --------------------------
        flags:
            Flags of the last processing
        psf_res:
            Result of psf processing.
        res:
            Result of image processing, if psf processing succeeded.
        """
        show=False
        out={}

        coellip_psf=self['coellip_psf']
        coellip_obj=self['coellip_obj']
        out['psf_res'] = self.process_image(ci.psf, 
                                            self['ngauss_psf'],
                                            ci['cen_psf_admom'],
                                            ci['cov_psf_admom'],
                                            coellip=coellip_psf)
        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:
            out['res'] = self.process_image(ci.image, 
                                            self['ngauss_obj'],
                                            ci['cen_admom'],
                                            ci['cov_admom'],
                                            psf=out['psf_res']['gmix'],
                                            coellip=coellip_obj)
            out['flags'] = out['res']['flags']
            if show and out['flags'] == 0:
                pprint(out['res'])
                self.show_residual(ci, out['psf_res']['gmix'], 
                                   objmix=out['res']['gmix'])
            elif show:
                self.show_residual(ci, out['psf_res']['gmix'])
        if out['flags'] != 0 and self['verbose']:
            print 'flags:',gmix_image.flagname(out['flags'])
        return out

    def process_image(self, image, ngauss, cen, cov, psf=None,
                      coellip=True):
        if not coellip:
            raise ValueError("must use coellip for now")

        counts = image.sum()
        guess = self.get_guess_coellip(counts, ngauss, cen, cov, psf=psf)
        gm = gmix_image.GMixFitCoellip(image,guess,psf=psf)
        #images.multiview(image)

        if gm.ier > 4:
            flags = gm.ier
        else:
            flags = 0
        out={'gmix':    gm.gmix,
             'pcov':    gm.pcov,
             'flags':   flags,
             'ier':     gm.ier,
             'numiter': gm.numiter,
             'coellip': coellip}
        return out

    def get_guess_coellip(self, counts, ngauss, cen, cov, psf=None):
        npars = 2*ngauss+4
        guess=zeros(npars)
        guess[0] = cen[0]
        guess[1] = cen[1]

        guess[2] = cov[0]
        guess[3] = cov[1]
        guess[4] = cov[2]

        if psf is not None:
            psfmoms = gmix_image.total_moms(psf)
            guess[2] -= psfmoms['irr']
            guess[3] -= psfmoms['irc']
            guess[4] -= psfmoms['icc']
        
        # If psf sent, this is an object. If ngauss==3, 
        # make guesses good for an exp disk
        if psf is not None and ngauss == 3:
            guess[5:5+3] = array([0.419696,0.0725887,0.499471])
            guess[5:5+3] = counts/guess[5:5+3].sum()
            guess[8] = 0.227659
            guess[9] = 3.57138
        else:
            # generic guesses
            guess[5:5+ngauss] = counts/ngauss

            if ngauss > 1:
                if ngauss == 2:
                    guess[5+ngauss] = 3.0
                elif ngauss == 3:
                    guess[5+ngauss] = 0.5
                    guess[5+ngauss+1] = 3.0
                else:
                    # 4 or mor
                    guess[5+ngauss] = 0.5
                    guess[5+ngauss+1] = 3.0
                    guess[5+ngauss+2:] = 4.0

        return guess

    def show_residual(self, ci, psfmix, objmix=None):
        """
        Show plots of the input compared with the fit gaussian mixtures.
        """
        
        psfmodel = gmix_image.gmix2image(psfmix,ci.psf.shape,
                                         counts=ci.psf.sum()) 
        images.compare_images(ci.psf,psfmodel,
                              label1='psf',label2='gmix')

        if objmix is not None:
            skysig=None
            if ci['skysig'] > 0:
                skysig=ci['skysig']
            model0 = gmix_image.gmix2image(objmix,ci.image0.shape,
                                           counts=ci.image0.sum()) 
            model = gmix_image.gmix2image(objmix,ci.image.shape,
                                          psf=psfmix,
                                          counts=ci.image.sum()) 

            images.compare_images(ci.image0,model0,
                                  label1='object0',label2='gmix',
                                  skysig=skysig)
            images.compare_images(ci.image,model,
                                  label1='object+psf',label2='gmix',
                                  skysig=skysig)
        stop

    def copy_output(self, s2, ellip, s2n, ci, res):
        st = zeros(1, dtype=self.out_dtype())

        # first copy inputs and data from the CI
        st['s2'] = s2
        st['s2n'] = s2n
        st['ellip'] = ellip
        st['e1true'] = ci['e1true']
        st['e2true'] = ci['e2true']
        st['etrue']  = ci['etrue']
        st['gamma'] = e2gamma(st['etrue'])
        st['gamma1'],st['gamma2'] = e1e2_to_g1g2(st['e1true'],st['e2true'])

        st['irr_uw'] = ci['cov_uw'][0]
        st['irc_uw'] = ci['cov_uw'][1]
        st['icc_uw'] = ci['cov_uw'][2]

        st['irr_psf_uw'] = ci['cov_psf_uw'][0]
        st['irc_psf_uw'] = ci['cov_psf_uw'][1]
        st['icc_psf_uw'] = ci['cov_psf_uw'][2]

        size2psf = ci['cov_psf_uw'][0]+ci['cov_psf_uw'][2]
        size2obj = ci['cov_image0_uw'][0]+ci['cov_image0_uw'][2]
        st['s2_uw'] = size2psf/size2obj

        s2psf_am = ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2]
        s2obj_am = ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2]
        st['s2admom'] = s2psf_am/s2obj_am
        st['sigma_psf_admom'] = \
            mom2sigma(ci['cov_psf_admom'][0]+ci['cov_psf_admom'][2])
        st['sigma_admom'] = \
            mom2sigma(ci['cov_image0_admom'][0]+ci['cov_image0_admom'][2])
        st['sigma0_admom'] = \
            mom2sigma(ci['cov_admom'][0]+ci['cov_admom'][2])

        if 'psf_res' in res:
            for s,r in zip( st['gmix_psf'], res['psf_res']['gmix']):
                for k in ['p','row','col','irr','irc','icc']:
                    s[k] = r[k]
            psf_moms = gmix_image.total_moms(res['psf_res']['gmix'])
            st['irr_psf_meas'] = psf_moms['irr']
            st['irc_psf_meas'] = psf_moms['irc']
            st['icc_psf_meas'] = psf_moms['icc']
            st['sigma_psf_meas'] = 0.5*(psf_moms['irr']+psf_moms['icc'])

            st['numiter_psf'] = res['psf_res']['numiter']

        if 'res' in res:
            for s,r in zip( st['gmix'], res['res']['gmix']):
                for k in ['p','row','col','irr','irc','icc']:
                    s[k] = r[k]

            moms = gmix_image.total_moms(res['res']['gmix'])
            st['irr_meas'] = moms['irr']
            st['irc_meas'] = moms['irc']
            st['icc_meas'] = moms['icc']
            st['s2_meas'] = \
                (psf_moms['irr']+psf_moms['icc'])/(moms['irr']+moms['icc'])
            st['sigma_meas'] = 0.5*(moms['irr']+moms['icc'])

            st['e1_meas'] = (moms['icc']-moms['irr'])/(moms['icc']+moms['irr']) 
            st['e2_meas'] = 2*moms['irc']/(moms['icc']+moms['irr']) 
            st['e_meas'] = sqrt(st['e1_meas']**2 + st['e2_meas']**2)


            st['gamma_meas'] = e2gamma(st['e_meas'])
            st['gamma1_meas'],st['gamma2_meas'] = \
                    e1e2_to_g1g2(st['e1_meas'],st['e2_meas'])

            st['flags'] = res['res']['flags']
            st['numiter'] = res['res']['numiter']

        else:
            st['s2_meas'] = -9999



        # figure out how to measure this
        st['s2n_meas'] = st['s2n']


        return st


    def out_dtype(self):
        gmix_dt = [('p','f8'),('row','f8'),('col','f8'),
                   ('irr','f8'),('irc','f8'),('icc','f8')]
        dt=[('s2n','f8'),
            ('ellip','f8'),

            ('s2','f8'),         # requested (spsf/sobj)**2
            ('s2_uw','f8'), # unweighted s2 of object before noise
            ('sigma_psf_admom','f8'),
            ('sigma_admom','f8'),
            ('sigma0_admom','f8'),
            ('s2admom','f8'),    # s2 from admom, generally different

            ('irr_uw','f8'),
            ('irc_uw','f8'),
            ('icc_uw','f8'),
            ('irr_psf_uw','f8'),
            ('irc_psf_uw','f8'),
            ('icc_psf_uw','f8'),

            ('etrue','f8'),
            ('e1true','f8'),
            ('e2true','f8'),
            ('gamma','f8'),
            ('gamma1','f8'),
            ('gamma2','f8'),

            ('numiter','i8'),
            ('numiter_psf','i8'),

            ('flags','i8'),

            ('s2n_meas','f8'),    # use admom s2n

            ('s2_meas','f8'),
            ('irr_psf_meas','f8'),
            ('irc_psf_meas','f8'),
            ('icc_psf_meas','f8'),
            ('irr_meas','f8'),
            ('irc_meas','f8'),
            ('icc_meas','f8'),
            ('sigma_meas','f8'),
            ('sigma_psf_meas','f8'),
            ('e_meas','f8'),
            ('e1_meas','f8'),
            ('e2_meas','f8'),
            ('gamma_meas','f8'),
            ('gamma1_meas','f8'),
            ('gamma2_meas','f8'),

            ('gmix_psf',gmix_dt,self['ngauss_psf']),
            ('gmix',gmix_dt,self['ngauss_obj'])]

        return dt


