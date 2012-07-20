from sys import stderr
import os
from os.path import join as path_join
import esutil as eu
from esutil.misc import wlog
from esutil.stat import wmom
import numpy
from numpy import ogrid, array, sqrt, where, linspace, median, zeros
from numpy.random import standard_normal
import fimage
from fimage import add_noise_admom, add_noise_dev, add_noise_uw
from fimage.conversions import mom2sigma, cov2sigma, etheta2e1e2, ellip2mom
from fimage.convolved import NoisyConvolvedImage
import lensing

import pprint

import admom
import fitsio
import time

class ShapeSim(dict):
    """
    The config file defines the PSF model and size as well as the
    galaxy model but not it's size or ellipticity or noise
    properties.
    
    """
    def __init__(self, simname, **keys):
        conf=read_config(simname)
        for k,v in conf.iteritems():
            self[k] = v

        # over-ride things
        for k,v in keys.iteritems():
            self[k] = v

        self.fs = 'hdfs'
        self.cache_list={}

        self['verbose'] = self.get('verbose',False)

        wlog("sim self:")
        pprint.pprint(self, stream=stderr)

    def write_trial(self, is2, ie, itheta=None):
        """
        Write simulate image/psf to the cache, to be used later.

        Unless this is a ring test, q random orientation is chosen.  A pickle
        of the entire convolved image object is written as well as a fits file
        with the images.

        parameters
        ----------
        is2: integer
            A number between 0 and self['nums2']-1
        ie: integer
            A number is a number between 0 and self['nume']-1

        All keys are written
        """

        s2,ellip = get_s2_e(self, is2, ie)
        theta = get_theta(self, itheta=itheta)

        ci=self.get_trial(s2,ellip,theta)

        wlog("dims: [%s,%s]" % ci.image.shape)
        wlog("theta:",theta)

        if itheta is not None:
            fits_file = get_theta_cache_url(self['name'],is2,ie,itheta,fs=self.fs)
            clobber=True
        else:
            fits_file = get_random_cache_url(self['name'], is2, ie, fs=self.fs)
            clobber=False

        wlog("writing cache file:",fits_file)
        shear=self.get_shear()
        if shear:
            h={'delta1':shear.e1,'delta2':shear.e2,
               'shear1':shear.g1,'shear2':shear.g2}
        else:
            h=None

        with eu.hdfs.HDFSFile(fits_file,verbose=True) as hdfs_file:
            ci.write_fits(hdfs_file.localfile, extra_keys=h)
            hdfs_file.put(clobber=clobber)
        return ci


    def get_trial(self, s2, ellip, theta):
        """
        Genereate a realization of the input size ratio squared and total
        ellipticity.

            - Unless dotrim is false, The image is trimmed to where 0.999937 is
            contained that is 4-sigma for a gaussian, but will have a different
            meaning for other profiles.

        parameters
        ----------
        s2: s2 value
        ellip: ellipticity value
        s2n: S/N ratio for object after convolution
        s2n_psf: S/N ratio for psf
        """

        ci_full = self.new_convolved_image(s2, ellip, theta)
        if self['dotrim']:
            if self['verbose']:
                wlog("trimming")
            ci = fimage.convolved.TrimmedConvolvedImage(ci_full)
        else:
            ci=ci_full

        return ci

    def read_random_cache(self, is2, ie):
        """
        Read data from a random cache file
        """
        f=self.get_random_cache_file(is2,ie)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile)
        return ci

    def read_random_ring_cache(self, is2, ie):
        """
        Read data from a random cache file
        """
        f=self.get_random_cache_file(is2,ie)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci1 = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile,
                                                          hid=1)
            ci2 = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile,
                                                          hid=2)
        return ci1, ci2



    def read_theta_cache(self, is2, ie, itheta):
        """
        Read data from a random cache file
        """
        f=get_theta_cache_url(self['name'],is2,ie,itheta,fs=self.fs)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile)
        return ci


    def get_random_cache_file(self, is2, ie):
        """
        Get random file from cache, with replacement.
        """
        key = '%s-%03i-%03i' % (self['name'],is2,ie)
        flist = self.cache_list.get(key,None)

        pattern=get_cache_pattern(self['name'],is2,ie,fs=self.fs)

        if flist is None:
            if self.fs == 'hdfs':
                flist=eu.hdfs.ls(pattern)
            else:
                import glob
                flist=glob.glob(pattern)
        self.cache_list[key] = flist

        if len(flist) < 100:
            raise ValueError("less than 100 files in cache for "
                             "%s %s %s" % (self['name'],is2,ie))
        i=eu.numpy_util.randind(len(flist),1)
        return flist[i]

    def show_ci(self, ci):
        import images
        images.multiview(ci.image0,title='pre-psf image')
        images.multiview(ci.psf,title='psf')
        images.multiview(ci.image,title='convolved image')
        key=raw_input('hit a key:')
        if key == 'q':
            stop

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
        if psfmodel == 'turb' and objmodel == 'dev':
            pass
            #sigma *= 1.75
            # might want to try this, but really seems we need something 
            # size dependent.  Why?
            #sigma *= 1.5

        shear=self.get_shear()
        cov=self.get_cov(sigma, obj_ellip, obj_theta, shear=shear)
        objpars = {'model':objmodel, 'cov':cov}

        if psfmodel in ['gauss','dgauss']:
            if objmodel == 'gauss':
                ci = fimage.convolved.ConvolverAllGauss(objpars,psfpars, **self)
            else:
                ci = fimage.convolved.ConvolverGaussFFT(objpars,psfpars, **self)
        else:
            ci = fimage.convolved.ConvolverTurbulence(objpars,psfpars, **self)

        ci['obj_theta'] = obj_theta
        return ci

    def get_shear(self):
        sh = self['shear']
        if sh is not None:
            if len(sh) != 2:
                raise ValueError("shear in config should have the "
                                 "form [g1,g2]")
            shear=lensing.Shear(g1=sh[0],g2=sh[1])
        else:
            shear=None
        return shear

    def get_cov(self, sigma, e, theta, shear=None):
        if shear:
            e1,e2 = etheta2e1e2(e, theta)
            shape=lensing.Shear(e1=e1,e2=e2)
            sheared_shape = shape + shear
            cov = ellip2mom(e1=sheared_shape.e1,
                            e2=sheared_shape.e2,
                            T=2*sigma**2)
        else:
            cov=ellip2mom(2*sigma**2,e=e,theta=theta)
        return cov

    def _get_gauss_psf_pars(self):
        """
        for gauss or double gauss psf
        """
        e1 = self['psf_e1']
        e2 = self['psf_e2']
        psf_cov=fimage.ellip2mom(2*self['psf_sigma']**2,
                                 e1=e1, e2=e2)
        #wlog("psf_cov:",psf_cov)
        #wlog("psf e1:",(psf_cov[2]-psf_cov[0])/(psf_cov[2]+psf_cov[0]))
        #wlog("psf e2:",2*psf_cov[1]/(psf_cov[2]+psf_cov[0]))
        if self['psfmodel'] == 'dgauss':
            psf_cov1=psf_cov
            psf_cov2=psf_cov*self['psf_sigrat']**2

            b=self['psf_cenrat']
            psfpars = {'model':'dgauss',
                       'cov1':psf_cov1,
                       'cov2':psf_cov2,
                       'cenrat':b}

            psf_cov = (psf_cov1 + b*psf_cov2)/(1+b)
            #psum = 1+self['psf_cenrat']
            #cov11 = (psf_cov1[0] + psf_cov2[0]*self['psf_cenrat'])/psum
            #cov22 = (psf_cov1[2] + psf_cov2[2]*self['psf_cenrat'])/psum
            psf_sigma_tot = cov2sigma(psf_cov)
        else:
            psfpars = {'model':'gauss', 'cov':psf_cov}
            psf_sigma_tot = self['psf_sigma']

        return psfpars, psf_sigma_tot




class BaseSim(dict):
    def __init__(self, run):
        conf=read_config(run)
        for k,v in conf.iteritems():
            self[k] = v
        numpy.random.seed(self['seed'])

        self.simc = read_config(self['sim'])

        self.fs='hdfs'

        self['verbose'] = self.get('verbose',False)

        simpars=self.get('simpars',{})
        simpars['verbose'] = self['verbose']
        self.shapesim = ShapeSim(self['sim'], **simpars)

        wlog("run self:")
        pprint.pprint(self, stream=stderr)

        orient=self.simc['orient']
        if orient != 'ring':
            raise ValueError("no longer support anything but ring")
    
    def wlog(self, *args):
        if self['verbose']:
            wlog(*args)

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


    def process_trials(self, is2, ie_or_is2n, itrial=None):
        runtype=self['runtype']
        if runtype == 'byellip':
            if itrial is not None:
                self.process_trial_by_e(is2, ie_or_is2n, itrial,
                                        dowrite=True, dolog=True)
            else:
                self.process_trials_by_e(is2, ie_or_is2n)
        else:
            if itrial is not None:
                self.process_trial_by_s2n(is2, ie_or_is2n, itrial,
                                          dowrite=True, dolog=True)
            else:
                self.process_trials_by_s2n(is2, ie_or_is2n)

    def process_trial_by_s2n(self, is2, is2n, itheta, 
                             dowrite=False, 
                             dolog=False):
        """
        Process a singe element in the ring, with nrepeat
        possible noise realizations
        """
        # fixed ie
        ie = self['ie']
        s2,ellip = get_s2_e(self.simc, is2, ie)
        s2n_psf = self['s2n_psf']
        s2n = get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']
        s2n_method = self['s2n_method']
        s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        # do this once and get noise realizations.
        ci_nonoise = self.get_a_trial(self.shapesim, is2, ie, 
                                      itheta=itheta)

        if dolog:
            wlog("ring theta: %s/%s" % (itheta+1,self.simc['nring']))
            wlog('ellip:',ellip,'s2n:',s2n,
                 's2n_psf:',s2n_psf,'s2n_method:',s2n_method)

        out = numpy.zeros(nrepeat, dtype=self.out_dtype())
        for irepeat in xrange(nrepeat):
            
            if self['verbose']:
                stderr.write('-'*70 + '\n')
            # we always write this, although slower when not verbose
            if (nrepeat > 1) and (( (irepeat+1) % 10) == 0 or irepeat == 0):
                stderr.write("  %s/%s repeat done\n" % ((irepeat+1),nrepeat))

            iter=0
            while iter < self['itmax']:

                ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                         s2n_method=s2n_method,
                                         fluxfrac=s2ncalc_fluxfrac)
                if self['verbose']:
                    wlog("s2n_uw:",ci['s2n_uw'],"s2n_uw_psf:",ci['s2n_uw_psf'])
                    if iter == 0: stderr.write("%s " % str(ci.psf.shape))
                res = self.run(ci)

                if res['flags'] == 0:
                    st = self.copy_output(s2, ellip, s2n, ci, res)
                    out[irepeat] = st
                    break
                else:
                    iter += 1

            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        if self['verbose']:
            stderr.write("niter: %d\n" % (iter+1))

        if dowrite:
            write_output(self['run'], is2, is2n, out, itrial=itheta,
                         fs=self.fs)
        return out

    def process_trials_by_s2n(self, is2, is2n):
        """
        This only works for ring tests, processing
        all in the ring.  You can call the single one too
        """

        s2n = get_s2n(self, is2n)
        s2n_fac = self['s2n_fac']

        nring = self.simc['nring']
        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        ntot = nring*nrepeat
        out = numpy.zeros(ntot, dtype=self.out_dtype())

        ii = 0
        for i in xrange(nring):
            itheta=i

            if self['verbose']:
                stderr.write('-'*70)
                stderr.write('\n')
            # we always write this
            stderr.write("%d/%d %d%% done\n" % ((i+1),nring,
                                                100.*(i+1)/float(nring)))

            dolog=False
            if i==0:
                dolog=True
            st = self.process_trial_by_s2n(is2, is2n, itheta, dolog=dolog)
            out[ii:ii+nrepeat] = st
            ii += nrepeat

        write_output(self['run'], is2, is2n, out, fs=self.fs)
        return out

    def process_trial_by_e(self, is2, ie, itheta,
                           dowrite=False, 
                           dolog=False):

        """
        Process a singe element in the ring, with nrepeat
        possible noise realizations
        """
        s2,ellip = get_s2_e(self.simc, is2, ie)

        s2n = self['s2n']
        s2n_psf = self['s2n_psf']
        s2n_fac = self['s2n_fac']
        s2n_method = self['s2n_method']
        s2ncalc_fluxfrac =self['s2ncalc_fluxfrac']

        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        # do this once and get noise realizations.
        ci_nonoise = self.get_a_trial(self.shapesim, is2, ie, 
                                      itheta=itheta)
        if dolog:
            wlog("ring theta: %s/%s" % (itheta+1,self.simc['nring']))
            wlog('s2n:',s2n, 's2n_psf:',s2n_psf,'s2n_method:',s2n_method)

        out = numpy.zeros(nrepeat, dtype=self.out_dtype())
        for irepeat in xrange(nrepeat):
            if self['verbose']:
                stderr.write('-'*70 + '\n')
            # we always write this, although slower when not verbose
            if (nrepeat > 1) and (( (irepeat+1) % 10) == 0 or irepeat == 0):
                stderr.write("  %s/%s repeat\n" % ((irepeat+1),nrepeat))

            iter=0
            while iter < self['itmax']:

                ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf,
                                         s2n_method=s2n_method,
                                         fluxfrac=s2ncalc_fluxfrac)
                if self['verbose']:
                    wlog("s2n_uw:",ci['s2n_uw'],"s2n_uw_psf:",ci['s2n_uw_psf'])
                    if iter == 0: stderr.write("%s " % str(ci.psf.shape))
                res = self.run(ci)

                if res['flags'] == 0:
                    st = self.copy_output(s2, ellip, s2n, ci, res)
                    out[irepeat] = st
                    break
                else:
                    iter += 1

            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        if self['verbose']:
            stderr.write("niter: %d\n" % (iter+1))

        if dowrite:
            write_output(self['run'], is2, ie, out, itrial=itheta,
                         fs=self.fs)
        return out

           
 
    def process_trials_by_e(self, is2, ie):
        """
        Generate random realizations of a particular element in the s2 and
        ellip sequences.

        parameters
        ----------
        is2: integer
            A number between 0 and self['nums2']-1
        ie: integer
            A number is a number between 0 and self['nume']-1
        """
        import images 

        s2n = self['s2n']
        s2n_psf = self['s2n_psf']
        s2n_fac = self['s2n_fac']

        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        nring = self.simc['nring']
        ntot = nring*nrepeat
        out = numpy.zeros(ntot, dtype=self.out_dtype())

        s2,ellip = get_s2_e(self.simc, is2, ie)
        wlog('ellip:',ellip)

        ii = 0
        for i in xrange(nring):
            itheta=i

            if self['verbose']:
                stderr.write('-'*70)
                stderr.write('\n')

            stderr.write("%d/%d %d%% done\n" % ((i+1),nring,
                                                100.*(i+1)/float(nring)))
            dolog=False
            if i==0:
                dolog=True
            st = self.process_trial_by_e(is2, ie, itheta, dolog=dolog)
            out[ii:ii+nrepeat] = st
            ii += nrepeat

        write_output(self['run'], is2, ie, out, fs=self.fs)
        return out

    

    def get_a_trial(self, ss, is2, ie, itheta=None):
        """
        Get a trial.

        If not using the cache, a new image is created. If adding to the cache,
        the fits file is also written.
        """
        orient=self.simc['orient']
        if self['use_cache']:
            if itheta is not None:
                ci=ss.read_theta_cache(is2,ie,itheta)
            else:
                ci=ss.read_random_cache(is2,ie)
        else:
            #if orient == 'ring':
            #    raise ValueError("use a cache for ring tests")
            if self['add_to_cache']:
                ci=ss.write_trial(is2, ie, itheta=itheta)
            else:
                s2,ellip = get_s2_e(self.simc, is2, ie)
                theta = get_theta(self.simc, itheta=itheta)
                ci=ss.get_trial(s2,ellip,theta)

        retrim = self['retrim']
        if retrim:
            if 'retrim_fluxfrac' not in self:
                raise ValueError("you must set fluxfrac for a retrim")
            retrim_fluxfrac = self['retrim_fluxfrac']
            ci_full = ci
            ci = fimage.convolved.TrimmedConvolvedImage(ci_full, fluxfrac=retrim_fluxfrac)

            if self['verbose']:
                wlog("re-trimming with fluxfrac: %.12g" % retrim_fluxfrac)
                wlog("old dims:",str(ci_full.image.shape),"new dims:",str(ci.image.shape))
        return ci


    def get_a_ring_trial(self, ss, is2, ie):
        """
        Get a trial.

        If not using the cache, a new image is created. If adding to the cache,
        the fits file is also written.
        """
        orient=self.simc['orient']
        ci1,ci2=ss.read_random_ring_cache(is2,ie)
        return ci1, ci2




def get_s2_e(conf, is2, ie):
    """
    Extract the s2 and e corresponding to the input indices
    """
    check_is2_ie(conf, is2, ie)
    s2 = linspace(conf['mins2'],conf['maxs2'], conf['nums2'])[is2]
    ellip = linspace(conf['mine'],conf['maxe'], conf['nume'])[ie]

    return s2, ellip

def get_nums2n(conf):
    if 's2nvals' in conf:
        nums2n = len(conf['s2nvals'])
    else:
        nums2n = conf['nums2n']

    return nums2n
def get_s2n(conf, is2n):
    """
    Extract the s2n corresponding to index
    """

    if 's2nvals' in conf:
        s2n = conf['s2nvals'][is2n]
    else:
        s2n = linspace(conf['mins2n'],conf['maxs2n'], conf['nums2n'])[is2n]
    return s2n

def get_s2n_nrepeat(s2n, fac=0.4):
    """
    Number of repeats.  This is not enough now that I'm using the
    matched s/n

    The 0.4 gives rather noisy results for exp but can run less than a day.
    It gives *very* noisy results for dev, need to increase.
    """
    nrep = round( (fac/( s2n/100. )**2) )
    if nrep < 1:
        nrep = 1
    nrep = int(nrep)
    return nrep

def check_is2_ie(conf, is2, ie):
    """
    Verify the is2 and ie are within range
    """
    max_is2 = conf['nums2']-1
    max_ie  = conf['nume']-1
    if (is2 < 0) or (is2 > max_is2):
        raise ValueError("is2 must be within [0,%d], "
                         "got %d" % (max_is2,is2))
    if (ie < 0) or (ie > max_ie):
        raise ValueError("ie must be within [0,%d], "
                         "got %d" % (max_ie,ie))


def get_theta(conf, itheta=None):
    orient=conf['orient']
    if itheta is None:
        if orient == 'ring':
            raise ValueError("you must send itheta for ring orientation")
        theta = 180.0*numpy.random.random()
    else:
        if orient != 'ring':
            raise ValueError("itheta only makes sense for ring test "
                             "simulation types")
        if (conf['nring'] % 2) != 0:
            raise ValueError("ring sims must have nring even")

        thetas = zeros(conf['nring'])
        nhalf = conf['nring']/2
        thetas[0:nhalf] = linspace(0,90,nhalf)
        thetas[nhalf:] = thetas[0:nhalf] + 90

        theta = thetas[itheta]
    return theta



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

def get_simdir(fs=None):
    if fs=='hdfs':
        dir=os.environ.get('LENSDIR_HDFS')
    else:
        dir=os.environ.get('LENSDIR')
    return path_join(dir, 'shapesim')

def get_run_dir(run, fs=None):
    dir=get_simdir(fs=fs)
    return path_join(dir,run)


def get_wq_dir(run, bytrial=False, combine=False):
    dir=get_run_dir(run)
    dir=path_join(dir, 'wq')
    if bytrial:
        dir = path_join(dir, 'bytrial')
    elif combine:
        dir = path_join(dir, 'combine')
    return dir

def get_wq_url(run, is2, ie, itrial=None, combine=False):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    """
    dir=get_wq_dir(run, bytrial = (itrial is not None), combine=combine)
    f='%s' % run
    if combine:
        f += '-combine'

    f += '-%03i-%03i' % (is2,ie)
    if itrial is not None and not combine:
        f += '-%05d' % itrial

    f+='.yaml'
    return path_join(dir, f)

def get_plot_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'plots')
    return dir



def get_plot_file(run, type, s2min=None, yrng=None):
    d=get_plot_dir(run)
    f='%s' % run

    if s2min:
        f += '-s2min%0.2f' % s2min

    if yrng is not None:
        f += '-yr%0.3f-%0.3f' % tuple(yrng)
    f += '-%s.eps' % type
    f = path_join(d, f)
    return f

def get_output_dir(run, sub=None, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'outputs')
    if sub:
        dir = path_join(dir, sub)
    return dir

def get_output_url(run, is2, ie, itrial=None, fs=None):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    Note ie might actually be is2n

    """
    sub=None
    if itrial is not None:
        sub='bytrial'
    dir=get_output_dir(run, sub=sub, fs=fs)
    f='%s-%03i-%03i' % (run,is2,ie)
    if itrial is not None:
        if itrial == '*':
            f += '-*'
        else:
            f += '-%05i' % itrial
    f += '.rec'
    return path_join(dir, f)


def get_cache_output_dir(simname, is2, ie, fs=None):
    d=get_simdir(fs=fs)
    subd='%03i-%03i' % (is2,ie)
    d=os.path.join(d, 'cache', simname,'outputs',subd)
    return d

def get_cache_wq_dir(simname):
    d=get_simdir()
    d=os.path.join(d, 'cache', simname,'wq')
    return d

def get_cache_wq_url(simname, is2, ie):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    """
    dir=get_cache_wq_dir(simname)
    f='%s-%03i-%03i.yaml' % (simname,is2,ie)

    return path_join(dir, f)

def random_cache_url_exists(url):
    if url[0:4] == 'hdfs':
        return eu.hdfs.exists(url)
    else:
        return os.path.exists(url)

def get_theta_cache_url(simname, is2, ie, itheta, fs=None):
    import tempfile
    d=get_cache_output_dir(simname, is2, ie, fs=fs)
    f='%s-%03i-%03i-%05i.fits' % (simname,is2,ie, itheta)
    url=os.path.join(d,f)
    return url

def get_random_cache_url(simname, is2, ie, fs=None):
    import tempfile
    d=get_cache_output_dir(simname, is2, ie, fs=fs)
    f='%s-%03i-%03i-' % (simname,is2,ie)
    url=tempfile.mktemp(dir=d, prefix=f,suffix='.fits')
    while random_cache_url_exists(url):
        url=tempfile.mktemp(dir=d, prefix=f,suffix='.fits')
    return url

def get_cache_pattern(simname, is2, ie, fs=None):
    import tempfile
    d=get_cache_output_dir(simname, is2, ie, fs=fs)
    pattern='%s-%03i-%03i-*.fits' % (simname,is2,ie)
    return os.path.join(d,pattern)


def write_output(run, is2, ie, data, itrial=None, fs=None):
    f=get_output_url(run, is2, ie, itrial=itrial, fs=fs)
    wlog("Writing output:",f)
    eu.io.write(f, data, clobber=True)

def read_output(run, is2, ie, itrial=None, verbose=False, fs=None):
    f=get_output_url(run, is2, ie, itrial=itrial, fs=fs)
    if verbose:
        wlog("reading output:",f)
    return eu.io.read(f)


def get_averaged_url(run, is2, fs=None, docum=False):
    """
    All the trials are averaged in a given s2 bin, and all
    ellip/s2n bins are in a single struct.
    """
    dir=get_output_dir(run, fs=fs)
    f='%s-%03i-avg' % (run,is2)
    if docum:
        f += '-cum'
    f += '.rec'
    return path_join(dir, f)

def write_averaged_outputs(run, data, docum=False, skip1=[], fs=None):
    if fs != 'hdfs':
        d = get_output_dir(run, fs=fs)
        if not os.path.exists(d):
            os.makedirs(d)

    c=read_config(run)
    cs=read_config(c['sim'])
    numi1 = cs['nums2']
    idata=0
    for i1 in xrange(numi1):
        if i1 not in skip1:
            f=get_averaged_url(run, i1, fs=fs, docum=docum)
            wlog(f)
            eu.io.write(f, data[idata], clobber=True)
            idata += 1

def read_averaged_outputs(run, docum=False, skip1=[], fs=None):
    c=read_config(run)
    cs=read_config(c['sim'])
    numi1 = cs['nums2']
    data=[]
    for i1 in xrange(numi1):
        if i1 not in skip1:
            f=get_averaged_url(run, i1, fs=fs, docum=docum)
            wlog(f)
            d=eu.io.read(f)
            data.append(d)
    return data

def make_averaged_outputs(run, docum=True, 
                          skip1=[], 
                          skip2=[]):
    """
    We write to both local and hdfs, for speed
    """
    data=[]
    c=read_config(run)
    cs=read_config(c['sim'])
    runtype=c['runtype']

    straight_avg=False
    if run[0:5] == 'deswl':
        straight_avg=True
        wlog("Doing straight average for deswl")

    numi1 = cs['nums2']
    if runtype == 'byellip':
        numi2 = cs['nume']
    else:
        numi2 = get_nums2n(c)

    orient=cs['orient']
    for i1 in xrange(numi1):
        if i1 in skip1:
            continue
        s2data=[]
        for i2 in xrange(numi2):
            if i2 in skip2:
                continue
            try:
                edata = read_output(run, i1, i2, verbose=True, fs='hdfs')
                s2data.append(edata)
            except:
                pass
        s2data = average_outputs(s2data, straight_avg=straight_avg)
        data.append(s2data)

    wlog("writing averaged outputs")
    write_averaged_outputs(run, data, skip1=skip1)
    write_averaged_outputs(run, data, skip1=skip1,fs='hdfs')

    if docum:
        cumdata = accumulate_outputs(data)
        wlog("writing acummulated averaged outputs")
        write_averaged_outputs(run, cumdata, docum=docum, skip1=skip1)
        write_averaged_outputs(run, cumdata, docum=docum, skip1=skip1,fs='hdfs')
    return data

def read_all_outputs(run, 
                     verbose=False, 
                     skip1=[], 
                     skip2=[], 
                     fs=None):
    """
    Data are grouped as a list by is2 and then sublists by ie/is2n

    If docum=True, we accumulate the sums for s2 < s2i.
    """
    data=[]
    c=read_config(run)
    cs=read_config(c['sim'])
    runtype=c['runtype']

    numi1 = cs['nums2']
    if runtype == 'byellip':
        numi2 = cs['nume']
    else:
        numi2 = get_nums2n(c)

    orient=cs['orient']
    for i1 in xrange(numi1):
        if i1 in skip1:
            continue
        s2data=[]
        for i2 in xrange(numi2):
            if i2 in skip2:
                continue
            try:
                edata = read_output(run, i1, i2, verbose=verbose,fs=fs)
                s2data.append(edata)
            except:
                pass
        data.append(s2data)

    return data

def accumulate_outputs(data):
    """
    Should already be summed over realizations using
    average_outputs

    Need to implement deswl averaging
    """

    if 'gamma1sum' in data[0].dtype.names:
        straight_avg=True
    else:
        straight_avg=False

    dt =data[0].dtype.descr
    dt += [('shear1cum','f8'), # sums so we can do cumulative
           ('shear2cum','f8'),
           ('shear1cum_err','f8'),
           ('shear2cum_err','f8')]
    if not straight_avg:
        dt += [('Rshearcum','f8')]
    out=[]

    num=len(data[0])
    for i,s2data in enumerate(data):
        nd = len(s2data)
        if nd != num:
            raise ValueError("when accumulating, there can be no missing values "
                             "use skip1/skip2 to help")

        d = zeros(num, dtype=dt)
        eu.numpy_util.copy_fields(s2data, d)

        if i > 0:
            # add previous
            if straight_avg:
                d['gamma1sum'] += dold['gamma1sum']
                d['gamma2sum'] += dold['gamma2sum']
            else:
                d['e1sum']  += dold['e1sum']
                d['e2sum']  += dold['e2sum']
                d['esqsum'] += dold['esqsum']

                #print 'this:',d['e1err2invsum']
                #print 'adding:',dold['e1err2invsum']
                d['e1err2invsum'] += dold['e1err2invsum']
                d['e2err2invsum'] += dold['e2err2invsum']


            d['nsum']   += dold['nsum']

        if straight_avg:
            d['shear1cum'] = d['gamma1sum']/d['nsum']
            d['shear2cum'] = d['gamma2sum']/d['nsum']

        else:
            d['Rshearcum'] = 1-0.5*d['esqsum']/d['nsum']
            d['shear1cum'] = 0.5*d['e1sum']/d['nsum']/d['Rshearcum']
            d['shear2cum'] = 0.5*d['e2sum']/d['nsum']/d['Rshearcum']

        dold = d.copy()
        if not straight_avg:
            # this only works for ring test with no shape noise!
            # might want 1/Rshear here
            d['shear1cum_err'] = 0.5*sqrt(1/d['e1err2invsum'])
            d['shear2cum_err'] = 0.5*sqrt(1/d['e2err2invsum'])

        out.append(d)

    return out

def average_outputs(data, straight_avg=False):
    """
    Input should be a list of arrays.  The output will be
    an array with length of list, values averaged over the
    elements
    """
    dt = data[0].dtype.descr

    if 'e1_meas' in data[0].dtype.names:
        if straight_avg:
            # deswl we can just do straight_avg sums on gamma
            dt_extra = [('gamma1sum','f8'),
                   ('gamma2sum','f8'),
                   ('nsum','i8')]
        else:
            dt_extra = [('e1sum','f8'), # sums so we can do cumulative
                   ('e2sum','f8'),
                   ('e1err2invsum','f8'),
                   ('e2err2invsum','f8'),
                   ('esqsum','f8'),
                   ('nsum','i8'),
                   ('shear1','f8'),
                   ('shear1err','f8'),  # average in particular bin
                   ('shear2','f8'),
                   ('shear2err','f8'),
                   ('Rshear_true','f8'),('Rshear','f8')]
    else:
        raise ValueError('DEAL WITH e1meas missing')

    dt += dt_extra
    name_extra = [dd[0] for dd in dt_extra]

    d=zeros(len(data),dtype=dt)
    for i,edata in enumerate(data): # over different ellipticities
        #shear1,shear2,R,shear1err,shear2err \
        #    = lensing.util.average_shear(edata['e1_meas'],edata['e2_meas'],
        #                                 doerr=True)

        if straight_avg:
            g1 = edata['gamma1_meas']
            g2 = edata['gamma2_meas']
            num = g1.size

            d['shear1'][i] = g1
            d['shear2'][i] = g2
            d['gamma1sum'][i] = g1.sum()
            d['gamma2sum'][i] = g2.sum()
            d['nsum'][i] = num
        else:
            e1 = edata['e1_meas']
            e2 = edata['e2_meas']
            esq = e1**2 + e2**2

            num=e1.size

            e1sum = e1.sum()
            e2sum = e2.sum()
            esqsum = esq.sum()

            atype='uw'
            if atype=='weighted':
                wtstot = 1/(2*0.3**2 + edata['pars_err'][:,2]**2 + edata['pars_err'][:,3]**2)
                wts1 = 1/(0.3**2 + edata['pars_err'][:,2]**2)
                wts2 = 1/(0.3**2 + edata['pars_err'][:,3]**2)
                me1,e1err0 = wmom(e1, wts1, calcerr=True)
                me2,e2err0 = wmom(e2, wts2, calcerr=True)

                mcrap,e1err = wmom(e1-edata['e1true'], wts1, calcerr=True)
                mcrap,e2err = wmom(e2-edata['e2true'], wts2, calcerr=True)

                mesq,mesq_err = wmom(esq, wtstot, calcerr=True)

                e1err2inv = 1/e1err**2
                e2err2inv = 1/e2err**2

                print e1err0,e1err
                print e2err0,e2err
            else:

                mesq = esqsum/num
                me1 = e1sum/num
                me2 = e2sum/num

                # we use the scatter from true, becuase with ring tests
                # we don't have shape noise
                e1scatt = (e1-edata['e1true']).var()
                e2scatt = (e2-edata['e2true']).var()

                e1err = sqrt(e1scatt/num)
                e2err = sqrt(e2scatt/num)
                e1err2inv = num/e1scatt
                e2err2inv = num/e2scatt

                #mesq_err = esq.std()/sqrt(num)
                #e1err = sqrt(1/e1err2inv)
                #e2err = sqrt(1/e2err2inv)
                
                #g1err = g1*sqrt( (mesq_err/mesq)**2 + (e1err/me1)**2 )
                #g2err = g2*sqrt( (mesq_err/mesq)**2 + (e1err/me2)**2 )

            R = 1-.5*mesq
            g1 = 0.5*me1/R
            g2 = 0.5*me2/R

            g1err = 0.5*e1err/R
            g2err = 0.5*e2err/R


            d['Rshear'][i] = R
            d['Rshear_true'][i] = (1-0.5*edata['etrue']**2).mean()
            d['shear1'][i] = g1
            d['shear2'][i] = g2
            d['shear1err'][i] = g1err
            d['shear2err'][i] = g2err

            d['e1sum'][i] = e1sum
            d['e2sum'][i] = e2sum
            d['e1err2invsum'][i] = e1err2inv
            d['e2err2invsum'][i] = e2err2inv
            d['esqsum'][i] = esqsum
            d['nsum'][i] = num

        for n in d.dtype.names:
            if n not in name_extra:
                if edata[n].dtype.names is None and len(edata[n].shape) == 1:
                    #d[n][i] = median(edata[n])
                    d[n][i] = edata[n].mean()

    return d


def average_outputs_old(data):
    """
    data is a list of lists of arrays

    Take the results from read_all_outputs and average the trials for each
    ellip value.  The result will be a list of arrays, one for each s2.  Each
    array will have one entry for each ellipticity
    """
    out=[]
    dt = data[0][0].dtype.descr

    if 'e1_meas' in data[0][0].dtype.names:
        dt += [('shear1','f8'),('shear1err','f8'),
               ('shear2','f8'),('shear2err','f8'),
               ('Rshear_true','f8'),('Rshear','f8')]
    else:
        raise ValueError('DEAL WITH e1meas missing')

    for s2data in data: # over different values of s2
        d=zeros(len(s2data),dtype=dt)
        for i,edata in enumerate(s2data): # over different ellipticities
            shear1,shear2,R,shear1err,shear2err \
                = lensing.util.average_shear(edata['e1_meas'],edata['e2_meas'],
                                             doerr=True)
            for n in d.dtype.names:
                if n == 'Rshear':
                    d['Rshear'][i] = R
                elif n == 'Rshear_true':
                    d['Rshear_true'][i] = (1-0.5*edata['etrue']**2).mean()
                elif n == 'shear1':
                    d['shear1'][i] = shear1
                elif n == 'shear2':
                    d['shear2'][i] = shear2
                elif n == 'shear1err':
                    d['shear1err'][i] = shear1err
                elif n == 'shear2err':
                    d['shear2err'][i] = shear2err
                else:
                    if edata[n].dtype.names is None and len(edata[n].shape) == 1:
                        #d[n][i] = median(edata[n])
                        d[n][i] = edata[n].mean()

        out.append(d)
    return out


def plot_signal_vs_rad(im, cen):
    import biggles
    row,col=ogrid[0:im.shape[0], 0:im.shape[1]]
    rm = array(row - cen[0], dtype='f8')
    cm = array(col - cen[1], dtype='f8')
    radm = sqrt(rm**2 + cm**2)

    radii = numpy.arange(0,im.shape[0]/2)
    cnts=numpy.zeros(radii.size)
    for ir,r in enumerate(radii):
        w=where(radm <= r)
        if w[0].size > 0:
            cnts[ir] = im[w].sum()

    xlog=True
    ylog=True
    cnts /= cnts.max()
    plt=eu.plotting.bscatter(radii, cnts,
                             xlabel=r'$r [pix]$',
                             ylabel='counts/max',
                             xlog=xlog,ylog=ylog,
                             show=False)

    w,=where(cnts > 0.999937)
    if w.size > 0:
        plt.add(biggles.Point(radii[w[0]], cnts[w[0]], 
                              type='filled circle', color='red'))

    plt.show()

def combine_trials(run, is2, ie):
    pattern=get_output_url(run, is2, ie, itrial='*', fs='hdfs')
    outfile=get_output_url(run, is2, ie, fs='hdfs')

    flist = eu.hdfs.ls(pattern, full=True)
    flist.sort()
    datalist=[]
    for f in flist:
        print f
        t=eu.io.read(f)
        datalist.append(t)

    data = eu.numpy_util.combine_arrlist(datalist)
    print 'data.size:',data.size
    print 'writing:',outfile
    eu.io.write(outfile, data, clobber=True)
