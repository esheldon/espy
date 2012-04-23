import os
from os.path import join as path_join
import esutil as eu
from esutil.misc import wlog
import numpy
from numpy import ogrid, array, sqrt, where
from numpy.random import standard_normal
import fimage
import admom

class ShapeSim(dict):
    """
    The config file defines the PSF model and size as well as the
    galaxy model and signal to noise ratio (negative for no noise)
    
    Then use the run(s2, ellip) method to generate
    a single realization.

    This runs a single s2 and ellip

    Note there are more things in the config than we use here.
    """
    def __init__(self, run, verbose=False):
        conf=read_config('sim',run)
        for k,v in conf.iteritems():
            self[k] = v
        self.verbose = verbose
    
    def get_trial(self, s2, ellip, s2n):
        """
        Genereate a realization of the input size ratio squared and total
        ellipticity.

            - a random orientation is chosen
            - if s2n is > 0 in the config, appropriate noise is also added.

        """
        theta = 180.0*numpy.random.random()
        ci = self.new_convolved_image(s2, ellip, theta)

        if s2n > 0:
            ci.image_nonoise = ci.image
            self.add_noise(ci, s2n)
        else:
            ci['skysig'] = 0.0

        return ci

    def new_convolved_image(self, s2, obj_ellip, obj_theta):
        """
        Generate a convolved image with the input parameters and the psf and
        object models listed in the config.
        """
        psf_cov=fimage.ellip2mom(2*self['psf_sigma']**2,
                              e=self['psf_ellip'],theta=0)
        if self['psfmodel'] == 'dgauss':
            psf_cov1=psf_cov
            psf_cov2=psf_cov*self['psf_sigrat']**2
            b=self['psf_cenrat']
            psfpars = dict(model = 'dgauss',
                           cov1 = psf_cov1,
                           cov2 = psf_cov2,
                           cenrat=b)
        else:
            psfpars = dict(model = 'gauss', cov = psf_cov)

        sigma = self['psf_sigma']/sqrt(s2)
        cov=fimage.ellip2mom(2*sigma**2,e=obj_ellip,theta=obj_theta)
        objpars = dict(model = self['objmodel'], cov=cov)

        ci = fimage.convolved.ConvolvedImageFFT(objpars,psfpars)

        return ci

    def add_noise(self, ci, s2n):
        """
        Add noise to a convolved image based on requested S/N.  We only add
        background noise so the S/N is

              sum(pix)
        -------------------   = S/N
        sqrt(npix*skysig**2)

        thus
            
            sum(pix)
        ----------------   = skysig
        sqrt(npix)*(S/N)

        
        We use an aperture of 3sigma but if no pixels
        are returned we increase the aperture until some
        are returned.

        Side effects:
            ci.image is set to the noisy image
            ci['skysig'] is set to the noise level
        """

        cen = ci['cen_uw']
        T = ci['cov_uw'][0] + ci['cov_uw'][2]
        sigma = fimage.mom2sigma(T)
        imagenn = ci.image_nonoise
        shape = imagenn.shape

        row,col=ogrid[0:shape[0], 0:shape[1]]
        rm = array(row - cen[0], dtype='f8')
        cm = array(col - cen[1], dtype='f8')

        radpix = sqrt(rm**2 + cm**2)
        # sigfac is 4 in admom
        sigfac = 4.0
        step=0.1
        w = where(radpix <= sigfac*sigma)
        npix = w[0].size + w[1].size
        while npix == 0:
            sigfac += step
            w = where(radpix <= sigfac*sigma)
            npix = w[0].size + w[1].size

        pix = imagenn[w]
        wt = sqrt(pix)
        wsignal = (wt*pix).sum()
        wsum = wt.sum()

        skysig = wsignal/sqrt(wsum)/s2n

        noise_image = \
            skysig*standard_normal(imagenn.size).reshape(shape)
        ci.image = imagenn + noise_image

        out = admom.admom(ci.image, cen[0], cen[1], guess=T/2., sigsky=skysig)

        # fix up skysig based on measurement
        skysig = out['s2n']/s2n*skysig
        noise_image = \
            skysig*standard_normal(imagenn.size).reshape(shape)
        ci.image = imagenn + noise_image
        ci['skysig'] = skysig


        if False:
            out = admom.admom(ci.image, cen[0], cen[1], guess=T/2., sigsky=skysig)
            wlog("    target S/N:            ",s2n)
            wlog("    meas S/N after noise:  ",out['s2n'])


def get_config_dir():
    d=os.environ['ESPY_DIR']
    return path_join(d,'shapesim','config')
def get_config_file(type,run):
    d=get_config_dir()
    name='%s-%s.yaml' % (type,run)
    return path_join(d, name)
def read_config(type,run):
    f=get_config_file(type,run)
    return eu.io.read(f)


