"""
Generate image simulations and process them with the
gmix pipeline
"""

from numpy import random, zeros
import sys
import lensing
from . import shapesim

try:
    import gmix_image
    from gmix_image import GMIX_ERROR_NEGATIVE_DET
except ImportError:
    sys.stderr("could not import gmix_image")

try:
    import admom
except ImportError:
    sys.stderr("could not import admom")

class GMixSim(shapesim.BaseSim):
    """
    We only override

        .run()
        .out_dtype()
        .copy_output()

    """
    def __init__(self, run):
        super(GMixSim,self).__init__(run)

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
        
        out={}

        out['psf_res'] = self.process_image(ci.psf, 
                                            self['psf_ngauss'],
                                            ci['cen_psf_admom'],
                                            ci['cov_psf_admom'])
        out['flags'] = out['psf_res']['flags']
        if out['flags'] == 0:

            out['res'] = self.process_image(ci.image, 
                                            self['ngauss'],
                                            ci['cen_admom'],
                                            ci['cov_admom'],
                                            psf=out['psf_res']['gmix'])

            out['flags'] = out['res']['flags']

        return out

    def process_image(self, image, ngauss, cen, cov, psf=None):
        im=image.copy()

        im_min = im.min()
        if im_min <= 0:
            im -= im_min
            sky=0.001*im.max()
            im += sky
        else:
            sky = im_min


        # In the iteration, we can sometimes run into negative determinants.
        # we will retry a few times with different random offsets in that case

        flags = GMIX_ERROR_NEGATIVE_DET
        ntry=0
        while flags == GMIX_ERROR_NEGATIVE_DET and ntry < self['max_retry']:
            guess = self.get_guess(ngauss, cen, cov)
            gm = gmix_image.GMix(im,guess,
                                 sky=sky,
                                 maxiter=self['maxiter'],
                                 tol=self['tol'],
                                 psf=psf)
            flags = gm.flags
            ntry += 1
        out={'gmix': gm.pars,
             'flags': gm.flags,
             'numiter':gm.numiter,
             'fdiff':gm.fdiff,
             'ntry':ntry}
        return out

    def get_guess(self, ngauss, cen, cov):
        # We use the input moments as guesses
        guess=ngauss*[{'p':1./ngauss,
                       'row':cen[0],
                       'col':cen[1],
                       'irr':cov[0],
                       'irc':cov[1],
                       'icc':cov[2]}]
        # perturb them
        for g in guess:
            g['row'] += 0.2*random.random()
            g['col'] += 0.2*random.random()
            g['irr'] += 0.2*random.random()
            g['irc'] += 0.2*random.random()
            g['icc'] += 0.2*random.random()
        return guess

    def copy_output(self, s2, ellip, s2n, ci, res):
        st = zeros(1, dtype=self.out_dtype())
        return st
    def out_dtype(self):
        dt=[('s2n','f8')]
        return dt
