"""
Generate image simulations and process them with the
gmix EM pipeline
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
    from gmix_image import GMIXEM_ERROR_NEGATIVE_DET
except ImportError:
    stderr.write("could not import gmix_image")

try:
    import admom
except ImportError:
    stderr.write("could not import admom")

def corr_e(e, R, s2n):
    bias = 1.+4./s2n**2*(1.-3./R + 1./R**2 + e**2)
    ecorr = e/bias
    return ecorr



class GMixEMSim(shapesim.BaseSim):
    """
    We only override

        .run()
        .out_dtype()
        .copy_output()

    """
    def __init__(self, run):
        super(GMixEMSim,self).__init__(run)
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

        out['psf_res'] = self.process_image(ci.psf, 
                                            self['ngauss_psf'],
                                            ci['cen_psf_admom'],
                                            ci['cov_psf_admom'],
                                            show=False)
        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:
            coellip=self.get('coellip',False)
            cocenter=self.get('cocenter',False)
            out['res'] = self.process_image(ci.image, 
                                            self['ngauss'],
                                            ci['cen_admom'],
                                            ci['cov_admom'],
                                            psf=out['psf_res']['gmix'],
                                            coellip=coellip,
                                            cocenter=cocenter,
                                            show=False)
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
                      coellip=False, cocenter=False, show=False):
        im=image.copy()

        # need no zero pixels and sky value
        im_min = im.min()
        if im_min <= 0:
            im -= im_min
            sky=0.001*im.max()
            im += sky
        else:
            sky = im_min


        # In the iteration, we can sometimes run into negative determinants.
        # we will retry a few times with different random offsets in that case


        flags = GMIXEM_ERROR_NEGATIVE_DET
        ntry=0
        while flags == GMIXEM_ERROR_NEGATIVE_DET and ntry < self['max_retry']:
            guess = self.get_guess(ngauss, cen, cov)
            gm = gmix_image.GMixEM(im,guess,
                                   sky=sky,
                                   maxiter=self['gmix_maxiter'],
                                   tol=self['gmix_tol'],
                                   coellip=coellip,
                                   cocenter=cocenter,
                                   psf=psf,
                                   verbose=self['verbose'])
            flags = gm.flags
            #stderr.write("gmix niter: %d\n" % gm.numiter)
            ntry += 1
        #stderr.write("ntry: %d " % ntry)
        out={'gmix': gm.pars,
             'flags': gm.flags,
             'numiter':gm.numiter,
             'fdiff':gm.fdiff,
             'ntry':ntry}
        if flags == 0 and show:
            print 'psf'
            pprint(psf)
            print 'gmix'
            pprint(out['gmix'])
            imrec=gmix_image.gmix2image(out['gmix'], image.shape, 
                                        psf=psf,counts=image.sum())
            images.compare_images(image,imrec,
                                  label1='image',label2='%d gauss' % ngauss)
        return out

    def get_guess(self, ngauss, cen, cov):
        # We use the input moments as guesses
        # cannot use [{...}]*ngauss, uses same objects!
        cenoff=0.2
        covoff=0.5
        guess=[]
        for i in xrange(ngauss):
            g={'p':1./ngauss,
               'row':cen[0],
               'col':cen[1],
               'irr':cov[0],
               'irc':cov[1],
               'icc':cov[2]}
            guess.append(g)
        # perturb them
        for g in guess:
            #g['p'] += 0.2*g['p']*random.random()
            g['row'] += cenoff*(random.random()-0.5)
            g['col'] += cenoff*(random.random()-0.5)
            g['irr'] += covoff*(random.random()-0.5)
            g['irc'] += covoff*(random.random()-0.5)
            g['icc'] += covoff*(random.random()-0.5)
        
        return guess

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
            st['ntry_psf'] = res['psf_res']['ntry']
            st['fdiff_psf'] = res['psf_res']['fdiff']

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
            st['ntry'] = res['res']['ntry']
            st['fdiff'] = res['res']['fdiff']
            
            # only makes sense for ngauss==1, need to adapt
            if s2n > 0:
                Tpsf = psf_moms['irr']+psf_moms['icc']
                Tobj = moms['irr']+moms['icc']
                R = 1.-Tpsf/(Tpsf+Tobj)

                st['e1_corr'] = corr_e(st['e1_meas'], R, s2n)
                st['e2_corr'] = corr_e(st['e2_meas'], R, s2n)
                st['e_corr'] = sqrt(st['e1_corr']**2 + st['e2_corr']**2)

                st['gamma_corr'] = e2gamma(st['e_corr'])
                st['gamma1_corr'],st['gamma2_corr'] = \
                        e1e2_to_g1g2(st['e1_corr'],st['e2_corr'])

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
            ('ntry','i8'),
            ('fdiff','f8'),
            ('numiter_psf','i8'),
            ('ntry_psf','i8'),
            ('fdiff_psf','f8'),

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

            ('e_corr','f8'),
            ('e1_corr','f8'),
            ('e2_corr','f8'),
            ('gamma_corr','f8'),
            ('gamma1_corr','f8'),
            ('gamma2_corr','f8'),

            ('gmix_psf',gmix_dt,self['ngauss_psf']),
            ('gmix',gmix_dt,self['ngauss']),
            ('gmix_corr',gmix_dt,self['ngauss'])]

        return dt


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
