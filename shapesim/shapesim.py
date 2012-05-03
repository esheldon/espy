from sys import stderr
import os
from os.path import join as path_join
import esutil as eu
from esutil.misc import wlog
import numpy
from numpy import ogrid, array, sqrt, where, linspace, median, zeros
from numpy.random import standard_normal
import fimage
import admom
import fitsio
import time

class ShapeSim(dict):
    """
    The config file defines the PSF model and size as well as the
    galaxy model but not it's size or ellipticity or noise
    properties.
    
    Then use the get_trial(s2, ellip, s2n) method to generate a single
    realization.
    """
    def __init__(self, run, **keys):
        conf=read_config(run)
        for k,v in conf.iteritems():
            self[k] = v
        for k,v in keys.iteritems():
            self[k] = v

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
        psfmodel = self['psfmodel']
        objmodel = self['objmodel']

        if psfmodel in ['gauss','dgauss']:
            psfpars, psf_sigma_tot = self._get_gauss_psf_pars()
        elif psfmodel == 'turb':
            psfpars = {'model':'turb','psf_fwhm':self['psf_fwhm']}
            psf_sigma_tot = self['psf_fwhm']/fimage.convolved.TURB_SIGMA_FAC
        else:
            raise ValueError("unknown psf model: '%s'" % psfmodel)

        sigma = psf_sigma_tot/sqrt(s2)
        cov=fimage.ellip2mom(2*sigma**2,e=obj_ellip,theta=obj_theta)
        objpars = {'model':objmodel, 'cov':cov}

        if psfmodel in ['gauss','dgauss']:
            if objmodel == 'gauss':
                ci = fimage.convolved.ConvolverAllGauss(objpars,psfpars, **self)
            else:
                ci = fimage.convolved.ConvolverGaussFFT(objpars,psfpars, **self)
            #ci = fimage.convolved.ConvolvedImageFFT(objpars,psfpars, **self)
        else:
            ci = fimage.convolved.ConvolverTurbulence(objpars,psfpars, **self)
            #ci = fimage.convolved.ConvolvedTurbulentPSF(objpars,psfpars, **self)

        return ci

    def _get_gauss_psf_pars(self):
        """
        for gauss or double gauss psf
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
            psum = 1+self['psf_cenrat']
            cov11 = (psf_cov1[0] + psf_cov2[0]*self['psf_cenrat'])/psum
            cov22 = (psf_cov1[2] + psf_cov2[2]*self['psf_cenrat'])/psum
            psf_sigma_tot = sqrt( (cov11+cov22)/2)
        else:
            psfpars = dict(model = 'gauss', cov = psf_cov)
            psf_sigma_tot = self['psf_sigma']

        return psfpars, psf_sigma_tot

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


class BaseSim(dict):
    def __init__(self, run):
        conf=read_config(run)
        for k,v in conf.iteritems():
            self[k] = v
        numpy.random.seed(self['seed'])

    def run(self, ci):
        """
        Process the input convolved image.

        over-ride this
        """
        raise RuntimeError("Override the .run() method")

    def copy_output(self, s2, ellip, s2n, ci, res):
        """
        Copy the result structure and convolved image
        to the array for output

        over-ride this
        """
        raise RuntimeError("Override the .copy_output() method")

    def out_dtype(self):
        """
        The output dtype

        over-ride this
        """
        raise RuntimeError("Override the .out_dtype() method")

    def process_trials(self, is2, ie, max_write_ci=1):
        """
        Generate random realizations of a particular element in the s2 and
        ellip sequences.

        parameters
        ----------

        is2: integer
            A number between 0 and self['nums2']-1
        ie: integer
            A number is a number between 0 and self['nume']-1
        max_write_ci: optional
            Max number of failed ci to write to disk. Default 1
        """
        import images 
        nwrite_ci=0

        out = numpy.zeros(self['ntrial'], dtype=self.out_dtype())
        ss = ShapeSim(self['sim'])

        s2n = self['s2n']
        s2,ellip = self.get_s2_e(is2, ie)

        for i in xrange(self['ntrial']):
            stderr.write(".")
            iter=0
            while iter < self['itmax']:

                ci=ss.get_trial(s2,ellip,s2n)
                res = self.run(ci)

                if res['flags'] == 0:
                    st = self.copy_output(s2, ellip, s2n, ci, res)
                    out[i] = st
                    break
                else:
                    images.multiview(ci.image,title='image')
                    images.multiview(ci.psf,title='psf')
                    stop
                    if nwrite_ci < max_write_ci:
                        self.write_ci(ci, is2, ie)
                        nwrite_ci += 1
                    iter += 1

            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        stderr.write("\n")
        write_output(self['run'], is2, ie, out)
        return out

    def write_ci(self, ci, is2, ie):
        """
        Write the ci to a file in the outputs directory
        """
        import tempfile
        rand=tempfile.mktemp(dir='')
        url=get_output_url(self['run'], is2, ie)
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

    def get_s2_e(self, is2, ie):
        """
        Extract the s2 and e corresponding to the input indices
        """
        self.check_is2_ie(is2, ie)
        s2 = linspace(self['mins2'],self['maxs2'], self['nums2'])[is2]
        ellip = linspace(self['mine'],self['maxe'], self['nume'])[ie]

        return s2, ellip

    def check_is2_ie(self, is2, ie):
        """
        Verify the is2 and ie are within range
        """
        max_is2 = self['nums2']-1
        max_ie  = self['nume']-1
        if (is2 < 0) or (is2 > max_is2):
            raise ValueError("is2 must be within [0,%d], "
                             "got %d" % (max_is2,is2))
        if (ie < 0) or (ie > max_ie):
            raise ValueError("ie must be within [0,%d], "
                             "got %d" % (max_ie,ie))




def get_config_dir():
    d=os.environ['ESPY_DIR']
    return path_join(d,'shapesim','config')
def get_config_file(run):
    d=get_config_dir()
    name='%s.yaml' % run
    return path_join(d, name)
def read_config(run):
    """
    run could be 'name' in sim
    """
    f=get_config_file(run)
    c = eu.io.read(f)
    if 'run' in c:
        n='run'
    else:
        n='name'
    if c[n] != run:
        raise ValueError("%s in config does not match "
                         "itself: '%s' instead of '%s'" % (n,c[n],run))
    return c
def get_simdir():
    dir=os.environ.get('LENSDIR')
    return path_join(dir, 'shapesim')

def get_run_dir(run):
    dir=get_simdir()
    return path_join(dir,run)


def get_wq_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'wq')
    return dir

def get_wq_url(run, is2, ie):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    """
    dir=get_wq_dir(run)
    f='%s-%03i-%03i.yaml' % (run,is2,ie)

    return path_join(dir, f)

def get_plot_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'plots')
    return dir

def get_plot_file(run, s2max=None, yrange=None):
    d=get_plot_dir(run)
    f='%s' % run

    if s2max is not None:
        f += '-s2max%0.3f' % s2max

    if yrange is not None:
        f += '-yr%0.3f-%0.3f' % tuple(yrange)
    f += '.eps'
    f = path_join(d, f)
    return f

def get_output_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'outputs')
    return dir

def get_output_url(run, is2, ie):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    """
    dir=get_output_dir(run)
    f='%s-%03i-%03i.rec' % (run,is2,ie)
    return path_join(dir, f)

def write_output(run, is2, ie, data):
    f=get_output_url(run, is2, ie)
    wlog("Writing output:",f)
    eu.io.write(f, data)

def read_output(run, is2, ie, verbose=False):
    f=get_output_url(run, is2, ie)
    if verbose:
        wlog("reading output:",f)
    return eu.io.read(f)

def read_all_outputs(run, average=False):
    """
    Data are grouped as a list by is2 and then sublists by ie
    """
    data=[]
    c=read_config(run)
    for is2 in xrange(c['nums2']):
        s2data=[]
        for ie in xrange(c['nume']):
            try:
                edata = read_output(run, is2, ie)
                s2data.append(edata)
            except:
                pass
        data.append(s2data)

    if average:
        return average_outputs(data)
    else:
        return data

def average_outputs(data):
    """
    Take the results from read_all_outputs and average the
    trials for each ellip value.  The result will be a list
    of arrays, one for each s2.  Each array will have one
    entry for each ellipticity
    """
    out=[]
    dt = data[0][0].dtype
    for s2data in data:
        # s2data is a list of arrays, one for each ellip.  the array has an
        # entry for each trial.  Average over the trials

        d=zeros(len(s2data),dtype=dt)
        for i,edata in enumerate(s2data):
            for n in d.dtype.names:
                if edata[n].dtype.names is None:
                    d[n][i] = median(edata[n])
        out.append(d)
    return out
