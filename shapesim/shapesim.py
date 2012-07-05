from sys import stderr
import os
from os.path import join as path_join
import esutil as eu
from esutil.misc import wlog
import numpy
from numpy import ogrid, array, sqrt, where, linspace, median, zeros
from numpy.random import standard_normal
import fimage
from fimage import add_noise_admom, add_noise_dev, add_noise_uw
from fimage.conversions import mom2sigma, cov2sigma, etheta2e1e2, ellip2mom
from fimage.convolved import NoisyConvolvedImage
import lensing

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

    def write_ring_trial(self, is2, ie):
        """
        Write a simulated image pair, plus psf to the cache, to be used later.

        A random angle is chosen, and paired with one at theta+90

        parameters
        ----------
        is2: integer
            A number between 0 and self['nums2']-1
        ie: integer
            A number is a number between 0 and self['nume']-1

        All keys are written
        """

        s2,ellip = get_s2_e(self, is2, ie)
        theta1 = get_theta(self)
        theta2 = theta1+90

        ci1=self.get_trial(s2,ellip,theta1)
        ci2=self.get_trial(s2,ellip,theta2)

        wlog("dims: [%s,%s]" % ci1.image.shape)
        wlog("theta1:",theta1)
        wlog("theta2:",theta2)

        fits_file = get_random_cache_url(self['name'], is2, ie, fs=self.fs)

        wlog("writing cache file:",fits_file)
        shear=self.get_shear()
        if shear:
            h={'delta1':shear.e1,'delta2':shear.e2,
               'shear1':shear.g1,'shear2':shear.g2}
        else:
            h=None

        with eu.hdfs.HDFSFile(fits_file,verbose=True) as hdfs_file:
            self.write_ring_fits(hdfs_file.localfile, ci1, ci2, extra_keys=h)
            hdfs_file.put()

    def write_ring_fits(self, fits_file, ci1, ci2, extra_keys=None):
        """
        Write the images and metadata to a fits file

        The images are in separate extensions 
            'image1','image1_0'
            'image2','image2_0',
            'psf'
        the metadata are in a binary tables 'table1','table2'

        parameters
        ----------
        fits_file: string
            Name of the file to write
        ci1: child of ConvolverBase
        ci2: child of ConvolverBase
        """
        import fitsio

        shear=self.get_shear()
        cilist = [ci1,ci2]
        with fitsio.FITS(fits_file,mode='rw',clobber=True) as fitsobj:
            for i in [1,2]:
                ci = cilist[i-1]
                dt=[]
                if shear is not None:
                    dt+=[('delta1','f8'),('delta2','f8'),
                         ('shear1','f8'),('shear2','f8')]
                for k,v in ci.iteritems():
                    if isinstance(v,int) or isinstance(v,long):
                        dt.append( (k, 'i8') )
                    elif isinstance(v,float):
                        dt.append( (k, 'f8') )
                    elif isinstance(v,numpy.ndarray):
                        this_t = v.dtype.descr[0][1]
                        this_n = v.size
                        if this_n > 1:
                            this_dt = (k,this_t,this_n)
                        else:
                            this_dt = (k,this_t)
                        dt.append(this_dt)
                    else:
                        raise ValueError("unsupported type: %s" % type(v))
                table = numpy.zeros(1, dtype=dt)
                for k,v in ci.iteritems():
                    table[k][0] = v
                if shear is not None:
                    table['delta1'],table['delta2'] = shear.e1,shear.e2
                    table['shear1'],table['shear2'] = shear.g1,shear.g2

                h={}
                # note not all items will be written, only basic types,
                # so this is not for feeding to the sim code.  The full
                # metadata are in the table
                for k,v in ci.iteritems():
                    h[k] = v
                if extra_keys:
                    for k,v in extra_keys.iteritems():
                        h[k] = v

                fitsobj.write(ci.image, header=h, extname='image%d' % i)
                fitsobj.write(ci.psf, extname='psf%d' % i)
                fitsobj.write(ci.image0, extname='image%d_0' % i)
                fitsobj.write(table, extname='table%d' % i)




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

        dotrim=self.get('dotrim',True)

        ci_full = self.new_convolved_image(s2, ellip, theta)
        if dotrim:
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
        wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=True) as hdfs_file:
            hdfs_file.stage()
            ci = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile)
        return ci

    def read_random_ring_cache(self, is2, ie):
        """
        Read data from a random cache file
        """
        f=self.get_random_cache_file(is2,ie)
        wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=True) as hdfs_file:
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
        wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=True) as hdfs_file:
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
        sh = self.get('shear',None)
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
        psf_cov=fimage.ellip2mom(2*self['psf_sigma']**2,
                                 e=self['psf_ellip'],theta=0)
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


    def process_trials(self, is2, ie_or_is2n):
        runtype=self.get('runtype','byellip')
        if runtype == 'byellip':
            self.process_trials_by_e(is2, ie_or_is2n)
        else:
            self.process_trials_by_s2n(is2, ie_or_is2n)

    def process_trials_by_s2n(self, is2, is2n):
        """
        This only works for ring tests
        """
        orient=self.simc.get('orient','rand')
        if orient != 'ring':
            raise ValueError("by s2n only works for ring tests")

        # fixed ie
        ie = self['ie']
        s2,ellip = get_s2_e(self.simc, is2, ie)
        s2n = get_s2n(self, is2n)


        nring = self.simc['nring']
        nrepeat = get_s2n_nrepeat(s2n)

        ntot = nring*nrepeat
        out = numpy.zeros(ntot, dtype=self.out_dtype())

        simpars=self.get('simpars',{})
        ss = ShapeSim(self['sim'], **simpars)

        s2n_psf = self['s2n_psf']

        print 'ellip:',ellip

        ii = 0
        for i in xrange(nring):
            iter=0

            itheta=i

            ci_nonoise = self.get_a_trial(ss, is2, ie, itheta=itheta)

            for irepeat in xrange(nrepeat):
                
                stderr.write('-'*70)
                stderr.write("\n%d/%d %d%% done\n" % (ii+1,ntot,100.*(ii+1)/float(ntot)))
                while iter < self['itmax']:

                    ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf)

                    if iter == 0: stderr.write("%s " % str(ci.psf.shape))
                    res = self.run(ci)

                    if res['flags'] == 0:
                        st = self.copy_output(s2, ellip, s2n, ci, res)
                        out[ii] = st
                        ii += 1
                        break
                    else:
                        iter += 1

                if iter == self['itmax']:
                    raise ValueError("itmax %d reached" % self['itmax'])
            stderr.write("niter: %d\n" % (iter+1))
        write_output(self['run'], is2, is2n, out, fs=self.fs)
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

        orient=self.simc.get('orient','rand')
        if orient == 'ring':
            ntrial = self.simc['nring']
            nrepeat = self['nrepeat']
        else:
            ntrial = self['ntrial']
            nrepeat=1

        out = numpy.zeros(ntrial*nrepeat, dtype=self.out_dtype())

        simpars=self.get('simpars',{})
        ss = ShapeSim(self['sim'], **simpars)

        s2n = self['s2n']
        s2n_psf = self['s2n_psf']

        s2,ellip = get_s2_e(self.simc, is2, ie)
        print 'ellip:',ellip

        ii = 0
        for i in xrange(ntrial):
            stderr.write('-'*70)
            stderr.write("\n%d/%d " % (i+1,ntrial))
            iter=0

            if orient == 'ring':
                itheta=i
            else:
                itheta=None

            # for ring, we can do multiple noise realizations
            # of each image, so only get nonoise image once
            if orient == 'ring':
                ci_nonoise = self.get_a_trial(ss, is2, ie, itheta=itheta)

            for irepeat in xrange(nrepeat):
                while iter < self['itmax']:

                    if orient != 'ring':
                        # for not ring, we always grab a new angle/image
                        ci_nonoise = self.get_a_trial(ss, is2, ie)

                    if s2n > 0 or s2n_psf > 0:
                        ci = NoisyConvolvedImage(ci_nonoise, s2n, s2n_psf)
                    else:
                        ci = ci_nonoise

                    if iter == 0: stderr.write("%s " % str(ci.psf.shape))
                    res = self.run(ci)

                    if res['flags'] == 0:
                        st = self.copy_output(s2, ellip, s2n, ci, res)
                        out[ii] = st
                        ii += 1
                        break
                    else:
                        iter += 1

                if iter == self['itmax']:
                    raise ValueError("itmax %d reached" % self['itmax'])
            stderr.write("niter: %d\n" % (iter+1))
        write_output(self['run'], is2, ie, out, fs=self.fs)
        return out

    def process_ring_trials(self, is2, ie):
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

        orient=self.simc.get('orient','rand')
        ntrial = self['ntrial']

        out = numpy.zeros(ntrial*2, dtype=self.out_dtype())

        simpars=self.get('simpars',{})
        ss = ShapeSim(self['sim'], **simpars)

        s2n = self['s2n']
        s2n_psf = self['s2n_psf']

        s2,ellip = get_s2_e(self.simc, is2, ie)
        print 'ellip:',ellip

        ii = 0
        for i in xrange(ntrial):
            stderr.write('-'*70)
            stderr.write("\n%d/%d " % (i+1,ntrial))
            iter=0


            # for ring, we can do multiple noise realizations
            # of each image
            while iter < self['itmax']:

                ci1,ci2 = self.get_a_ring_trial(ss, is2, ie)
                if s2n > 0 or s2n_psf > 0:
                    ci1 = NoisyConvolvedImage(ci1, s2n, s2n_psf)
                    ci2 = NoisyConvolvedImage(ci2, s2n, s2n_psf)

                if iter == 0: stderr.write("%s " % str(ci1.psf.shape))
                res1 = self.run(ci1)
                if res1['flags'] == 0:
                    res2 = self.run(ci2)
                    if res2['flags'] == 0:
                        st1 = self.copy_output(s2, ellip, s2n, ci1, res1)
                        out[ii] = st1
                        ii += 1
                        st2 = self.copy_output(s2, ellip, s2n, ci2, res2)
                        out[ii] = st2
                        ii += 1
                        break
                    else:
                        iter+=1
                        continue
                else:
                    iter+=1
                    continue

            if iter == self['itmax']:
                raise ValueError("itmax %d reached" % self['itmax'])
        stderr.write("niter: %d\n" % (iter+1))
        write_output(self['run'], is2, ie, out, fs=self.fs)
        return out


    def get_a_trial(self, ss, is2, ie, itheta=None):
        """
        Get a trial.

        If not using the cache, a new image is created. If adding to the cache,
        the fits file is also written.
        """
        orient=self.simc.get('orient','rand')
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
        return ci


    def get_a_ring_trial(self, ss, is2, ie):
        """
        Get a trial.

        If not using the cache, a new image is created. If adding to the cache,
        the fits file is also written.
        """
        orient=self.simc.get('orient','rand')
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

def get_s2n(conf, is2n):
    """
    Extract the s2n corresponding to index
    """

    s2n = linspace(conf['mins2n'],conf['maxs2n'], conf['nums2n'])[is2n]
    return s2n

def get_s2n_nrepeat(s2n):
    """
    Number of repeats
    """
    ntrial = round( (0.4/( s2n/100. )**2) )
    if ntrial < 1:
        ntrial = 1
    ntrial = int(ntrial)
    return ntrial

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
    orient=conf.get('orient','rand')
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



def get_plot_file(run, type, s2max=None, s2meas=False, yrange=None):
    d=get_plot_dir(run)
    f='%s' % run

    if s2meas:
        f += '-s2meas'
    if s2max is not None:
        f += '-s2max%0.3f' % s2max

    if yrange is not None:
        f += '-yr%0.3f-%0.3f' % tuple(yrange)
    f += '-%s.eps' % type
    f = path_join(d, f)
    return f

def get_output_dir(run, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'outputs')
    return dir

def get_output_url(run, is2, ie, fs=None):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    Note ie might actually be is2n

    """
    dir=get_output_dir(run, fs=fs)
    f='%s-%03i-%03i.rec' % (run,is2,ie)
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


def write_output(run, is2, ie, data, fs=None):
    f=get_output_url(run, is2, ie, fs=fs)
    wlog("Writing output:",f)
    eu.io.write(f, data, clobber=True)

def read_output(run, is2, ie, verbose=False, fs=None):
    f=get_output_url(run, is2, ie, fs=fs)
    if verbose:
        wlog("reading output:",f)
    return eu.io.read(f)

def read_all_outputs(run, average=False, verbose=False, skip1=[], skip2=[], fs=None):
    """
    Data are grouped as a list by is2 and then sublists by ie/is2n
    """
    data=[]
    c=read_config(run)
    cs=read_config(c['sim'])
    runtype=c.get('runtype','byellip')

    numi1 = cs['nums2']
    if runtype == 'byellip':
        numi2 = cs['nume']
    else:
        numi2 = c['nums2n']

    orient=cs.get('orient','rand')
    for i1 in xrange(numi1):
        if i1 in skip1:
            continue
        s2data=[]
        for i2 in xrange(numi2):
            if i2 in skip2:
                continue
            if orient != 'ring':
                try:
                    edata = read_output(run, i1, i2, verbose=verbose,fs=fs)
                    s2data.append(edata)
                except:
                    pass
            else:
                # we require all for ring for cancellation
                edata = read_output(run, i1, i2, verbose=verbose,fs=fs)
                s2data.append(edata)
        data.append(s2data)

    if average:
        return average_outputs(data)
    else:
        return data

def average_outputs(data):
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
               ('Rshear','f8')]
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
