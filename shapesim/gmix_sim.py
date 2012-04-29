"""
Generate image simulations and process them with the
gmix pipeline
"""

import sys
import lensing
from . import shapesim

try:
    import gmix_image
    from gmix_image import GMIX_ERROR_NEGATIVE_DET
except:
    sys.stderr("could not import gmix_image")

try:
    import admom
except:
    sys.stderr("could not import admom")

class GMixSim(shapesim.BaseSim):
    """
    We only override

        .run()
        .out_dtype()
        .copy_output()

    """
    def __init__(self, run):
        super(self,GMixSim).__init__(run)

    def run(self, ci):
        """
        Process the input convolved image
        """

        # need to ensure there is some sky
        
        im_min = ci.image.min()
        im = ci.image.copy()
        if im_min <= 0:
            im -= im_min
            sky=0.001*im.max()
            im += sky
        else:
            sky = im_min

        # In the iteration, we can sometimes run into negative determinants.
        # we will retry a few times in that case

        gm = gmix_image.GMix(im,gd,
                             sky=sky,
                             maxiter=self['maxiter'],
                             tol=self['tol'])
        out={'gmix': gm.pars,
             'flags': gm.flags,
             'numiter':gm.numiter,
             'fdiff':gm.fdiff}
        return None

    def process_image(self, image, ngauss):
        im=image.copy()

        im_min = im.min()
        if im_min <= 0:
            im -= im_min
            sky=0.001*im.max()
            im += sky
        else:
            sky = im_min

        # to generate guesses, we run adaptive moments
        # In the iteration, we can sometimes run into negative determinants.
        # we will retry a few times with different random offsets in that case

        flags = GMIX_ERROR_NEGATIVE_DET
        while flags == GMIX_ERROR_NEGATIVE_DET:
            gm = gmix_image.GMix(im,guess,
                                 sky=sky,
                                 maxiter=self['maxiter'],
                                 tol=self['tol'])
            flags = gm.flags
        out={'gmix': gm.pars,
             'flags': gm.flags,
             'numiter':gm.numiter,
             'fdiff':gm.fdiff}

    def copy_output(self, s2, ellip, s2n, ci, res):
        st = numpy.zeros(1, dtype=self.out_dtype())
        return st
    def out_dtype(self):
        dt=[('s2n','f8')]
        return dt
