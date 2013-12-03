from sys import stderr
import os
from os.path import join as path_join
import esutil as eu
from esutil.misc import wlog
from esutil.stat import wmom
import numpy
from numpy import ogrid, array, sqrt, where, linspace, median, zeros, pi
from numpy import sin, cos, arctan2
from numpy.random import randn
import fimage
from fimage import add_noise_admom, add_noise_dev, add_noise_uw
from fimage.conversions import mom2sigma, cov2sigma, etheta2e1e2, ellip2mom
from fimage.convolved import NoisyConvolvedImage
import lensing
from lensing import Shear

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
        self.update(conf)
        self.update(keys)

        self._set_Tobj()
        self._set_counts()

        self.fs=get_default_fs()

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


    def get_trial(self, Tobj, ellip, theta, counts=1.0, center_offset=None):
        """
        Genereate a realization of the input size ratio squared and total
        ellipticity.

            - Unless dotrim is false, The image is trimmed to where 0.999937 is
            contained that is 4-sigma for a gaussian, but will have a different
            meaning for other profiles.

        parameters
        ----------
        Tobj: T value for object
        ellip: ellipticity value
        s2n: S/N ratio for object after convolution
        s2n_psf: S/N ratio for psf
        """

        ci_full = self.new_convolved_image(Tobj, ellip, theta, counts=counts, 
                                           center_offset=center_offset)
        if self['dotrim']:
            fluxfrac=self.get('fluxfrac',0.999937)
            if self['verbose']:
                wlog("trimming",fluxfrac)
            ci = fimage.convolved.TrimmedConvolvedImage(ci_full,fluxfrac=fluxfrac)
        else:
            ci=ci_full

        return ci

    def read_random_cache(self, iT, ie):
        """
        Read data from a random cache file
        """
        f=self.get_random_cache_file(iT,ie)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile)
        return ci

    def read_random_ring_cache(self, iT, ie):
        """
        Read data from a random cache file
        """
        f=self.get_random_cache_file(iT,ie)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci1 = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile,
                                                          hid=1)
            ci2 = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile,
                                                          hid=2)
        return ci1, ci2



    def read_theta_cache(self, iT, ie, itheta):
        """
        Read data from a random cache file
        """
        f=get_theta_cache_url(self['name'],iT,ie,itheta,fs=self.fs)
        if self['verbose']:
            wlog("\nreading from cache:",f)

        with eu.hdfs.HDFSFile(f,verbose=self['verbose']) as hdfs_file:
            hdfs_file.stage()
            ci = fimage.convolved.ConvolvedImageFromFits(hdfs_file.localfile)
        return ci


    def get_random_cache_file(self, iT, ie):
        """
        Get random file from cache, with replacement.
        """
        key = '%s-%03i-%03i' % (self['name'],iT,ie)
        flist = self.cache_list.get(key,None)

        pattern=get_cache_pattern(self['name'],iT,ie,fs=self.fs)

        if flist is None:
            if self.fs == 'hdfs':
                flist=eu.hdfs.ls(pattern)
            else:
                import glob
                flist=glob.glob(pattern)
        self.cache_list[key] = flist

        if len(flist) < 100:
            raise ValueError("less than 100 files in cache for "
                             "%s %s %s" % (self['name'],iT,ie))
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

    def new_convolved_image(self, Tobj, obj_ellip, obj_theta, counts=1.0,
                            center_offset=None):
        """
        Generate a convolved image with the input parameters and the psf and
        object models listed in the config.
        """
        # new thing using the gmix_image code for gaussian objects
        if self['objmodel'] in ['gexp','gdev','gauss','gbd']:
            return self.new_gmix_convolved_image(Tobj, obj_ellip, obj_theta, 
                                                 counts=counts,
                                                 center_offset=center_offset)

        psfmodel = self['psfmodel']
        objmodel = self['objmodel']

        if psfmodel in ['gauss','dgauss']:
            psfpars, psf_sigma_tot = self._get_gauss_psf_pars()
        elif psfmodel == 'turb':
            psfpars = {'model':'turb','psf_fwhm':self['psf_fwhm']}
            psf_sigma_tot = self['psf_fwhm']/fimage.convolved.TURB_SIGMA_FAC
        else:
            raise ValueError("unknown psf model: '%s'" % psfmodel)

        sigma = sqrt(Tobj)/2.

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
        if shear is not None:
            ci['shear1'] = shear.g1
            ci['shear2'] = shear.g2
        else:
            ci['shear1']=0.
            ci['shear2']=0.
        return ci

    def new_gmix_convolved_image(self, Tobj, obj_ellip, obj_theta, counts=1.0,
                                 center_offset=None):
        """
        Generate a convolved image with the input parameters and the psf and
        object models listed in the config.
        """
        import gmix_image
        psfmodel = self['psfmodel']
        objmodel = self['objmodel']
            
        g1psf = self['psf_g1']
        g2psf = self['psf_g2']

        shape_psf=Shear(g1=g1psf, g2=g2psf)

        if 'Tpsf' in self:
            Tpsf  = self['Tpsf']
        else:
            Tpsf  = 2*self['psf_sigma']**2
        if psfmodel in ['gauss','gturb']:
            psfpars=[-9., -9., shape_psf.g1, shape_psf.g2, Tpsf, 1.0]
            if psfmodel=='gauss':
                psf_gmix=gmix_image.GMixCoellip(psfpars)
            else:
                psf_gmix=gmix_image.GMixTurb(psfpars)
        else:
            raise ValueError("unsupported gmix psf type: '%s'" % psfmodel)

        e1,e2 = etheta2e1e2(obj_ellip, obj_theta)

        shape0 = Shear(e1=e1,e2=e2)
        shear=self.get_shear()

        if shear is not None:
            shape = shape0 + shear
            #shape = shear + shape0
        else:
            shape=shape0

        if objmodel in ['gbd','bd']:
            # we always set Texp to T
            frac_dev=self['frac_dev']
            Tfrac_dev=self['Tfrac_dev']
            
            Flux_exp = counts*(1.0-frac_dev)
            Flux_dev = counts*frac_dev

            # need to adjust a bit to get T
            T0 = (1-frac_dev)*(1-Tfrac_dev)*Tobj + frac_dev*Tfrac_dev*Tobj
            fac = Tobj/T0

            T_exp = (1-Tfrac_dev)*Tobj*fac
            T_dev = Tfrac_dev*Tobj*fac

            objpars=[-9., -9., shape.g1, shape.g2, T_exp, T_dev, Flux_exp, Flux_dev]
            obj_gmix=gmix_image.GMix(objpars,type='bd')
        else:
            objpars=[-9., -9., shape.g1, shape.g2, Tobj, counts]
            if objmodel=='gexp':
                obj_gmix=gmix_image.GMixExp(objpars)
            elif objmodel=='gdev':
                obj_gmix=gmix_image.GMixDev(objpars)
            elif objmodel=='gauss':
                obj_gmix=gmix_image.GMixCoellip(objpars)
            else:
                raise ValueError("unsupported gmix object type: '%s'" % objmodel)

        ci=fimage.convolved.ConvolverGMix(obj_gmix, psf_gmix, center_offset=center_offset, **self)

        ci['obj_theta'] = obj_theta
        if shear is not None:
            ci['shear1'] = shear.g1
            ci['shear2'] = shear.g2
        else:
            ci['shear1']=0.
            ci['shear2']=0.
        return ci


    def get_shear(self):
        if 'shearmag' in self:
            numpy.random.seed(None)

            shearmag=self['shearmag']+self['shearwidth']*randn()

            theta=numpy.random.random()*pi
            g1 = shearmag*cos(2*theta)
            g2 = shearmag*sin(2*theta)
            shear=Shear(g1=g1,g2=g2)
        elif 'shear' in self:
            sh = self['shear']
            if len(sh) != 2:
                raise ValueError("shear in config should have the "
                                 "form [g1,g2]")
            shear=Shear(g1=sh[0],g2=sh[1])
        else:
            shear=None
        return shear

    def get_cov(self, sigma, e, theta, shear=None):
        if shear:
            e1,e2 = etheta2e1e2(e, theta)
            shape=Shear(e1=e1,e2=e2)
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
        g1 = self['psf_g1']
        g2 = self['psf_g2']
        shape_psf=Shear(g1=g1, g2=g2)
        psf_cov=fimage.ellip2mom(2*self['psf_sigma']**2,
                                 e1=shape_psf.e1, e2=shape.e2)
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

    def _set_Tobj(self):
        self._T_dists=None
        self._Tvals=numpy.array(self['Tobj'])

        T_dist_type=self.get('T_dist',None)
        if T_dist_type is not None:
            self._T_dists=[]
            T_width_frac = self['T_width_frac']

            for T in self._Tvals:
                T_w = T*T_width_frac
                dist=eu.random.get_dist(T_dist_type,[T,T_w])
                self._T_dists.append(dist)

    def _get_Tobj(self, iT):
        if self._T_dists is not None:
            return self._T_dists[iT].sample()
        else:
            return self._Tvals[iT]
    
    def _set_counts(self):
        self._counts_mean=1.0
        self._counts_dist=None

        counts_dist=self.get('counts_dist',None)

        if counts_dist is not None:
            counts_width_frac=self['counts_width_frac']

            dist=eu.random.get_dist(counts_dist,
                                    [self._counts_mean, 
                                    self._counts_mean*counts_width_frac])
            self._counts_dist=dist

    def _get_counts(self):
        if self._counts_dist is not None:
            return self._counts_dist.sample()
        else:
            return self._counts_mean
    
class BaseSim(dict):
    def __init__(self, run):
        conf=read_config(run)
        for k,v in conf.iteritems():
            self[k] = v

        self.simc = read_config(self['sim'])

        self.fs=get_default_fs()

        self['verbose'] = self.get('verbose',False)

        simpars=self.get('simpars',{})
        simpars['verbose'] = self['verbose']
        self.shapesim = ShapeSim(self['sim'], **simpars)

        wlog("run self:")
        pprint.pprint(self, stream=stderr)

        orient=self.simc['orient']
        #if orient != 'ring':
        #    raise ValueError("no longer support anything but ring")

   
    def wlog(self, *args):
        if self['verbose']:
            wlog(*args)

    def run(self, ci):
        """
        Process the input convolved image.

        over-ride this
        """
        raise RuntimeError("Override the .run() method")

    def copy_output(self, Tobj, ellip, s2n, ci, res):
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


    def process_trials(self, iT, ie_or_is2n, itrial=None):
        runtype=self['runtype']
        if runtype == 'byellip':
            if itrial is not None:
                self.process_trial_by_e(iT, ie_or_is2n, itrial,
                                        dowrite=True, dolog=True)
            else:
                self.process_trials_by_e(iT, ie_or_is2n)
        else:
            if itrial is not None:
                self.process_trial_by_s2n(iT, ie_or_is2n, itrial,
                                          dowrite=True, dolog=True)
            else:
                self.process_trials_by_s2n(iT, ie_or_is2n)

    def process_trial_by_s2n(self, iT, is2n, itheta, 
                             dowrite=False, 
                             dolog=False):
        """
        Process a singe element in the ring, with nrepeat
        possible noise realizations
        """
        raise ValueError("fix this!")
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
            wlog("ring theta: %s/%s" % (itheta+1,self.simc['nsplit']))
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
                    wlog("s2n_admom:",ci['s2n_admom'],"s2n_uw:",ci['s2n_uw'],"s2n_matched:",ci['s2n_matched'])
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

        nsplit = self.simc['nsplit']
        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        ntot = nsplit*nrepeat
        out = numpy.zeros(ntot, dtype=self.out_dtype())

        ii = 0
        for i in xrange(nsplit):
            itheta=i

            if self['verbose']:
                stderr.write('-'*70)
                stderr.write('\n')
            # we always write this
            stderr.write("%d/%d %d%% done\n" % ((i+1),nsplit,
                                                100.*(i+1)/float(nsplit)))

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
            wlog("ring theta: %s/%s" % (itheta+1,self.simc['nsplit']))
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
        ie: integer
            A number is a number between 0 and self['nume']-1
        """
        import images 

        s2n = self['s2n']
        s2n_psf = self['s2n_psf']
        s2n_fac = self['s2n_fac']

        nrepeat = get_s2n_nrepeat(s2n, fac=s2n_fac)

        nsplit = self.simc['nsplit']
        ntot = nsplit*nrepeat
        out = numpy.zeros(ntot, dtype=self.out_dtype())

        s2,ellip = get_s2_e(self.simc, is2, ie)
        wlog('ellip:',ellip)

        ii = 0
        for i in xrange(nsplit):
            itheta=i

            if self['verbose']:
                stderr.write('-'*70)
                stderr.write('\n')

            stderr.write("%d/%d %d%% done\n" % ((i+1),nsplit,
                                                100.*(i+1)/float(nsplit)))
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


# from a BA13 prior run, exp galaxy
# /astro/u/esheldon/lensing/shapesim/cbafit-geg02r07/outputs/cbafit-geg02r07-000-avg.rec
s2n_ref_geg=[15,  20,  25,  30,  40,  50,  60,  70,  80,  90,  100,  120,  140,  160,  180,  200,  250,  300,  350,  400]

err_ref_geg=[6.04966948076, 6.27720009086, 6.39230479966, 6.46027631805, 6.31490331065, 5.07534788684, 4.24058192767, 3.63903754734, 3.18832441297, 2.83577943389, 2.55346528496, 2.12927717929, 1.82605936882, 1.5983608504, 1.42135553436, 1.27964356694, 1.02561949518, 0.857061929443, 0.737296676262, 0.648170473541]
npair_ref_geg=[1333000,  750000,  480000,  333000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000]

s2n_ref_geg=numpy.array(s2n_ref_geg,dtype='f8')
err_ref_geg=numpy.array(err_ref_geg)*1.e-5
npair_ref_geg=numpy.array(npair_ref_geg)


s2n_ref_deg=numpy.array([15,  20,  25,  30,  40,  50,  60,  70,  80,  90,  100,  120,  140,  160,  180,  200,  250,  300,  350,  400],dtype='f8')


err_ref_deg=numpy.array([  5.72424798e-05,   6.10264401e-05,   6.32783893e-05,
                         6.47026643e-05,   6.63254272e-05,   6.71194744e-05,
                         6.75744094e-05,   6.78835605e-05,   6.80540579e-05,
                         6.81622436e-05,   6.82656513e-05,   6.84032611e-05,
                         6.84544518e-05,   6.84826729e-05,   6.85100035e-05,
                         6.85641846e-05,   6.84825446e-05,   6.83578669e-05,
                         6.82696058e-05,   6.80481557e-05])
npair_ref_deg=numpy.array([3768000, 2280984, 1514069, 1072740,  615414,  397620,  277595,
                           204484,  156980,  124212,  100687,   69975,   51480,   39438,
                           31178,   25272,   16236,   11336,    8397,    6489])

def get_npair_by_noise(s2n, desired_err, run):
    """
    given the desired final error, determine the required number of pairs
    """

    if 'geg' in run:
        npairii = numpy.interp([s2n], s2n_ref_geg, npair_ref_geg)
        errii = numpy.interp([s2n], s2n_ref_geg, err_ref_geg)
    elif 'deg' in run:
        npairii = numpy.interp([s2n], s2n_ref_deg, npair_ref_deg)
        errii = numpy.interp([s2n], s2n_ref_deg, err_ref_deg)

    # desired_err = errii*sqrt(npairii/npair)
    # thus npair = npairii*(errii/desired_err)^2
    npair = npairii*(errii/desired_err)**2

    return npair[0]
    

def get_npair_nsplit_by_noise(c, is2n):
    from math import ceil
    s2n = c['s2n_vals'][is2n]
    npair_tot = get_npair_by_noise(s2n, c['desired_err'],c['run'])
    #print 'desired_err:',c['desired_err']
    #print 'npair_tot:',npair_tot

    # to keep equal time, normalize to zeroth
    nsplit0 = c['nsplit0']

    if is2n==0:
        nsplit=nsplit0
    else:
        npair_tot0 = get_npair_by_noise(c['s2n_vals'][0], c['desired_err'],c['run'])
        nsplit = int( ceil( nsplit0*float(npair_tot)/npair_tot0 ))
    
    npair_per = int(ceil(npair_tot/float(nsplit)))

    return npair_per, nsplit

def get_npair_nsplit(c, is2n):
    """
    Get number of pairs per split and number of splits

    For equal_time, we take number per split from is2n==0
    """
    if 'desired_err' in c:
        return get_npair_nsplit_by_noise(c, is2n)
    else:
        s2n = c['s2n_vals'][is2n]

        npair = get_s2n_nrepeat(s2n, fac=c['s2n_fac'])
        if npair < c['min_npair']:
            npair = c['min_npair']
        nsplit=c['nsplit']

        return npair,nsplit



def get_s2_e(conf, is2, ie):
    """
    Extract the s2 and e corresponding to the input indices
    """
    check_is2_ie(conf, is2, ie)
    s2 = linspace(conf['mins2'],conf['maxs2'], conf['nums2'])[is2]
    ellip = linspace(conf['mine'],conf['maxe'], conf['nume'])[ie]

    return s2, ellip

def get_numT(conf):
    return len(conf['Tobj'])
def get_nums2n(conf):
    return len(conf['s2nvals'])


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

        thetas = get_ring_thetas(conf['nsplit'])
        theta = thetas[itheta]
    return theta

def get_ring_thetas(nsplit):
    nsplit=int(nsplit)
    if (nsplit % 2) != 0:
        raise ValueError("nsplit must be even")

    thetas = zeros(nsplit)
    nhalf = nsplit/2
    thetas[0:nhalf] = linspace(0,90,nhalf)
    thetas[nhalf:] = thetas[0:nhalf] + 90

    return thetas




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

def get_default_fs():
    if os.environ.get('SHAPESIM_FS')=='hdfs':
        fs='hdfs'
    else:
        fs='nfs'
    return fs

def get_simdir(fs=None):
    if fs=='hdfs':
        dir=os.environ.get('LENSDIR_HDFS')
        dir=path_join(dir, 'shapesim')
    elif fs=='local':
        dir='/data/esheldon/lensing/shapesim'
    else:
        dir=os.environ['SHAPESIM_DIR']

    return dir

def get_run_dir(run, fs=None):
    dir=get_simdir(fs=fs)
    return path_join(dir,run)

def get_pbs_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'pbs')
    return dir

def get_condor_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'condor')
    return dir

def get_condor_job_url(run):
    d=get_condor_dir(run)
    return path_join(d,'%s.condor' % run)

def get_condor_master_url(run):
    d=get_condor_dir(run)
    return path_join(d,'%s.sh' % run)


def get_minions_url(run, i1):
    d=get_pbs_dir(run)
    return path_join(d,'%s-minions-%03d.pbs' % (run,i1))

def get_minions_script_url(run):
    d=get_pbs_dir(run)
    return path_join(d,'%s-minions.sh' % run)


def get_commands_url(run,i1):
    d=get_pbs_dir(run)
    return path_join(d,'%s-commands-%03d.txt' % (run,i1))


def get_wq_dir(run, bytrial=False, combine=False, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'wq')
    if bytrial:
        dir = path_join(dir, 'bytrial')
    elif combine:
        dir = path_join(dir, 'combine')
    return dir

def get_wq_url(run, is2, ie, itrial=None, combine=False, fs=None):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    """
    dir=get_wq_dir(run, bytrial = (itrial is not None), combine=combine, fs=fs)
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


def get_bias_file(run, type):
    d=get_plot_dir(run)
    f='%s' % run

    f += '-%s.fits' % type
    f = path_join(d, f)
    return f


def get_plot_file(run, type, s2n_name=None, s2min=None, yrng=None, use_pqr=False):
    d=get_plot_dir(run)
    f='%s' % run

    if s2min:
        f += '-s2min%0.2f' % s2min

    if yrng is not None:
        f += '-yr%0.3f-%0.3f' % tuple(yrng)
    f += '-%s' % type

    if use_pqr:
        f += '-pqr'

    if s2n_name is not None:
        f += '-%s' % s2n_name

    f += '.eps'
    f = path_join(d, f)
    return f

def get_output_dir(run, sub=None, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'outputs')
    if sub:
        dir = path_join(dir, sub)
    return dir

def get_output_url(run, is2, ie, itrial=None, fs=None, ext=None):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    Note ie might actually be is2n

    """
    if ext is None:
        if 'ngmix' in run:
            ext='fits'
        else:
            ext='rec'
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
    f += '.%s' % ext
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
    if 'hdfs' not in f:
        d=os.path.dirname(f)
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass
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
    numi1 = get_numT(cs)
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
    numi1 = get_numT(cs)
    data=[]
    for i1 in xrange(numi1):
        if i1 not in skip1:
            f=get_averaged_url(run, i1, fs=fs, docum=docum)
            wlog(f)
            d=eu.io.read(f)
            data.append(d)
    return data

def make_averaged_outputs(run, docum=False, 
                          skip1=[], 
                          skip2=[]):
    """
    if fs is hdfs, we write to both local and hdfs, for speed
    """

    fs=get_default_fs()

    data=[]
    c=read_config(run)
    cs=read_config(c['sim'])
    runtype=c['runtype']

    straight_avg=False
    bayes=False
    if 'bayes' in run or 'mixmc' in run or 'bafit' in run:
        wlog("doing bayes averaging")
        bayes=True
    elif 'deswl' in run:
        straight_avg=True
        wlog("Doing straight average for deswl")

    numi1 = get_numT(cs)
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
            edata = read_output(run, i1, i2, verbose=True, fs=fs)
            edata['shear_true'][:,0] = cs['shear'][0]
            edata['shear_true'][:,1] = cs['shear'][1]
            edata['gtrue'][:,0] = cs['shear'][0]
            edata['gtrue'][:,1] = cs['shear'][1]
            s2data.append(edata)

        if 'shearmag' in cs:
            s2data = average_randshear_outputs(s2data)
        else:
            s2data = average_outputs(s2data, straight_avg=straight_avg, bayes=bayes, orient=cs['orient'])
        data.append(s2data)

    wlog("writing averaged outputs")
    write_averaged_outputs(run, data, skip1=skip1, fs=fs)
    if fs=='hdfs':
        write_averaged_outputs(run, data, skip1=skip1,fs='nfs')

    if docum:
        cumdata = accumulate_outputs(data)
        wlog("writing acummulated averaged outputs")
        write_averaged_outputs(run, cumdata, docum=docum, skip1=skip1,fs=fs)
        if fs=='hdfs':
            write_averaged_outputs(run, cumdata, docum=docum, skip1=skip1,fs='nfs')
    return data


def read_all_outputs(run, 
                     verbose=False, 
                     skip1=[], 
                     skip2=[], 
                     fs=None):
    """
    Data are grouped as a list by is2 and then sublists by ie/is2n

    """
    data=[]
    c=read_config(run)
    cs=read_config(c['sim'])
    runtype=c['runtype']

    numi1 = get_numT(cs)
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
            edata = read_output(run, i1, i2, verbose=verbose,fs=fs)
            s2data.append(edata)
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

def average_outputs(data, straight_avg=False, bayes=False, orient='ring'):
    """
    Input should be a list of arrays.  The output will be
    an array with length of list, values averaged over the
    elements
    """
    dt = data[0].dtype.descr

    if straight_avg:
        # deswl we can just do straight_avg sums on gamma
        dt_extra = [('gamma1sum','f8'),
                    ('gamma2sum','f8'),
                    ('nsum','i8')]
    elif bayes:
        npars=data[0]['pars'].shape[1]
        dt_extra = [('g1sum','f8'),
                    ('g2sum','f8'),
                    ('g1sensum','f8'),
                    ('g2sensum','f8'),
                    ('g1err2invsum','f8'),
                    ('g2err2invsum','f8'),
                    ('pmeans','f8',npars),
                    ('pvarmeans','f8',npars),  # average of square error for each par
                    ('psums','f8',npars),
                    ('pvarsums','f8',npars),   # sum of square error for each par
                    ('nsum','i8'),
                    ('shear1','f8'),
                    ('shear1err','f8'),
                    ('shear2','f8'),
                    ('shear2err','f8')]
        if 'Ts2n' not in data[0].dtype.names:
            dt_extra+=[('Ts2n','f8')]
        dt_extra += [('Ts2n_sum','f8')]

        dt_extra += [('flux_s2n_sum','f8')]

        if 'flux' in data[0].dtype.names:
            dt_extra += [('flux_sum','f8'),
                         ('flux_err2invsum','f8')]


        if 'gcov0' in data[0].dtype.names:
            dt_extra += [('g1err0sum2','f8'),
                         ('g2err0sum2','f8'),
                         ('g1err0_mean','f8'), # this is the rms: sqrt(sum(sigma^2)/n)
                         ('g2err0_mean','f8')]

        if 'Q' in data[0].dtype.names:
            dt_extra+=[('Q_sum','f8',2),
                       ('Cinv_sum','f8',(2,2)),
                       ('bashear','f8',2),
                       ('bashear_cov','f8',(2,2))]
    else:
        dt_extra = [('g1sum','f8'), # sums so we can do cumulative
                    ('g2sum','f8'),
                    ('g1err2invsum','f8'),
                    ('g2err2invsum','f8'),
                    ('nsum','i8'),
                    ('shear1','f8'),
                    ('shear1err','f8'),  # average in particular bin
                    ('shear2','f8'),
                    ('shear2err','f8')]

    dt += dt_extra
    name_extra = [dd[0] for dd in dt_extra]

    # other names not to avarage
    name_extra += ['flux_err']

    d=zeros(len(data),dtype=dt)
    for i,edata in enumerate(data): # over different ellipticities

        if straight_avg:
            g1 = edata['gamma1_meas']
            g2 = edata['gamma2_meas']
            num = g1.size

            d['shear1'][i] = g1
            d['shear2'][i] = g2
            d['gamma1sum'][i] = g1.sum()
            d['gamma2sum'][i] = g2.sum()
            d['nsum'][i] = num
        elif bayes:
            g1 = edata['g'][:,0]
            g2 = edata['g'][:,1]
            g1sens = edata['gsens'][:,0]
            g2sens = edata['gsens'][:,1]

            g1sens_mean=g1sens.mean()
            g2sens_mean=g2sens.mean()

            g1sum = g1.sum()
            g2sum = g2.sum()
            g1sensum = g1sens.sum()
            g2sensum = g2sens.sum()

            # we use the scatter from true, becuase with ring tests
            # we don't have shape noise

            num = g1.size
            #g1err = sqrt(g1scatt/num)
            #g2err = sqrt(g2scatt/num)
            if orient == 'ring':
                if True:
                    SN=0.0
                    g1err2invsum = ( 1/(SN**2 + edata['gcov'][:,0,0]) ).sum()
                    g2err2invsum = ( 1/(SN**2 + edata['gcov'][:,1,1]) ).sum()
                    g1err = sqrt(1/g1err2invsum)/g1sens_mean
                    g2err = sqrt(1/g2err2invsum)/g2sens_mean
                else:
                    g1var = (g1-edata['gtrue'][:,0]).var()
                    g2var = (g2-edata['gtrue'][:,1]).var()
                    g1err2invsum = num/g1var
                    g2err2invsum = num/g2var
                    g1err = sqrt(g1var/num)
                    g2err = sqrt(g2var/num)
            else:
                g1corr = g1/g1sens.mean()
                g2corr = g2/g2sens.mean()
                g1err = g1corr.std()/sqrt(num)/g1sens_mean
                g2err = g2corr.std()/sqrt(num)/g2sens_mean

            d['nsum'][i] = g1.size
            d['shear1'][i] = g1sum/g1sensum
            d['shear2'][i] = g2sum/g2sensum
            d['shear1err'][i] = g1err
            d['shear2err'][i] = g2err
            print 'shear1: %.16g +/- %.16g' % (d['shear1'][i],d['shear1err'][i])

            d['g1sum'][i] = g1sum
            d['g2sum'][i] = g2sum
            d['g1sensum'][i] = g1sensum
            d['g2sensum'][i] = g2sensum

            d['g1err2invsum'][i] = g1err2invsum
            d['g2err2invsum'][i] = g2err2invsum


            if 'pars' in edata.dtype.names:
                for pi in xrange(npars):
                    d['psums'][i,pi] = edata['pars'][:,pi].sum()
                    d['pvarsums'][i,pi] = edata['pcov'][:,pi,pi].sum()
                    d['pmeans'][i,pi] = d['psums'][i,pi]/num
                    d['pvarmeans'][i,pi] = d['pvarsums'][i,pi]/num

                if 'Ts2n' in edata.dtype.names:
                    Ts2n_vals=edata['Ts2n']
                else:
                    # ack this won't work for complex models!
                    Ts2n_vals=edata['pars'][:,4]/sqrt(edata['pcov'][:,4,4])
                d['Ts2n_sum'][i] = Ts2n_vals.sum()
                d['Ts2n'][i] = d['Ts2n_sum'][i]/num

                if 'flux_s2n' in edata.dtype.names:
                    flux_s2n_vals=edata['flux_s2n']
                    d['flux_s2n_sum'][i] = flux_s2n_vals.sum()
                    d['flux_s2n'][i] = d['flux_s2n_sum'][i]/num

            if 'flux' in edata.dtype.names:
                flux_sum=edata['flux'].sum()
                d['flux_sum'][i] = flux_sum
                d['flux'][i] = flux_sum/num
                d['flux_err2invsum'][i] = (1./edata['flux_err']**2).sum()
                d['flux_err'][i] = sqrt(1./d['flux_err2invsum'][i])

            if 'gcov0' in data[0].dtype.names:
                d['g1err0sum2'][i] = edata['gcov0'][:,0,0].sum()
                d['g2err0sum2'][i] = edata['gcov0'][:,1,1].sum()
                d['g1err0_mean'][i] = sqrt(d['g1err0sum2'][i]/num)
                d['g2err0_mean'][i] = sqrt(d['g2err0sum2'][i]/num)

            if 'Q' in data[0].dtype.names:
                # NOTE!  can't use errors because this is a ring
                # test.  Using errors from above
                P = edata['P']
                Q = edata['Q']
                R = edata['R']
                g1g2, C, Q_sum, Cinv_sum = \
                        lensing.shear.get_shear_pqr(P,Q,R,get_sums=True)
                d['Q_sum'][i] = Q_sum
                d['Cinv_sum'][i] = Cinv_sum
                d['bashear'][i] = g1g2
                d['bashear_cov'][i,0,0] = d['shear1err'][i]**2
                d['bashear_cov'][i,1,1] = d['shear2err'][i]**2
                print 'bashear1: %.16g +/- %.16g' % (g1g2[0],sqrt(d['bashear_cov'][i,0,0]))
        else:

            g1 = edata['g'][:,0]
            g2 = edata['g'][:,1]

            num=g1.size

            # use median because sometimes the error is very small (e.g. wrong)
            med_invcov1 = numpy.median(1./edata['gcov'][:,0,0])
            med_invcov2 = numpy.median(1./edata['gcov'][:,1,1])
            g1err2invsum = num*med_invcov1
            g2err2invsum = num*med_invcov2

            #g1err2invsum = ( 1/(SN1**2 + edata['gcov'][:,0,0]) ).sum()
            #g2err2invsum = ( 1/(SN1**2 + edata['gcov'][:,1,1]) ).sum()
            g1err = sqrt(1/g1err2invsum)
            g2err = sqrt(1/g2err2invsum)

            g1sum = g1.sum()
            g2sum = g2.sum()

            mg1 = g1sum/num
            mg2 = g2sum/num

            d['shear1'][i] = mg1
            d['shear2'][i] = mg2
            d['shear1err'][i] = g1err
            d['shear2err'][i] = g2err

            d['g1sum'][i] = g1sum
            d['g2sum'][i] = g2sum
            d['g1err2invsum'][i] = g1err2invsum
            d['g2err2invsum'][i] = g2err2invsum
            d['nsum'][i] = num
            print 'shear1: %.16g +/- %.16g' % (d['shear1'][i],d['shear1err'][i])

        for n in d.dtype.names:
            if n in edata.dtype.names and n != 'model':
                if n not in name_extra:
                    if edata[n].dtype.names is None and len(edata[n].shape) == 1:
                        #d[n][i] = median(edata[n])
                        d[n][i] = edata[n].mean()

    return d


def average_runs(runlist, new_run_name, skip1=[]):
    """
    The runs must already be averaged and have the same size
    in all indices of relevance
    """

    fs=get_default_fs()

    if 'bayes' in runlist[0] or 'bafit' in runlist[0]:
        bayes=True
        sumlist=['g1sum',
                 'g2sum',
                 'g1sensum',
                 'g2sensum',
                 'g1err2invsum',
                 'g2err2invsum',
                 'nsum']

        if 'bafit' in runlist[0]:
            sumlist += ['Q_sum', 'Cinv_sum']

        f=get_averaged_url(runlist[0], 0, fs=fs)
        t=eu.io.read(f)
        if 'g1err0sum2' in t.dtype.names:
            sumlist+=['g1err0sum2','g2err0sum2']
        if 'Ts2n_sum' in t.dtype.names:
            print 'doing Ts2n'
            sumlist +=['Ts2n_sum']

        if 'flux_sum' in t.dtype.names:
            sumlist += ['flux_sum','flux_err2invsum','flux_s2n_sum']
    else:
        sumlist=['g1sum',
                 'g2sum',
                 'g1err2invsum',
                 'g2err2invsum',
                 'nsum']
        bayes=False

    dir = get_output_dir(new_run_name)
    if not os.path.exists(dir):
        wlog("making dir:",dir)
        os.makedirs(dir)
 
    c=read_config(runlist[0])
    cs0=read_config(c['sim'])

    numi1 = get_numT(cs0)
    for i1 in xrange(numi1):
        if i1 in skip1:
            continue

        for irun,run in enumerate(runlist):
            f=get_averaged_url(run, i1, fs=fs)
            wlog("Reading:",f)
            d=eu.io.read(f)

            if irun == 0:
                data = d.copy()
            else:
                for field in sumlist:
                    data[field] += d[field]
            
        # this is an array with mean for each in second index
        if bayes:
            g1sens_mean = data['g1sensum']/data['nsum']
            g2sens_mean = data['g2sensum']/data['nsum']
            data['shear1'] = data['g1sum']/data['g1sensum']
            data['shear2'] = data['g2sum']/data['g2sensum']
            data['shear1err'] = sqrt(1/data['g1err2invsum'])/g1sens_mean
            data['shear2err'] = sqrt(1/data['g2err2invsum'])/g2sens_mean
            for d in data:
                print 'shear1: %.16g +/- %.16g' % (d['shear1'],d['shear1err'])

            if 'g1err0sum2' in sumlist:
                data['g1err0_mean'] = sqrt(data['g1err0sum2']/data['nsum'])
                data['g2err0_mean'] = sqrt(data['g2err0sum2']/data['nsum'])

            if 'Ts2n_sum' in t.dtype.names:
                data['Ts2n'] = data['Ts2n_sum']/data['nsum']

            if 'flux_sum' in sumlist:
                data['flux'] = data['flux_sum']/data['nsum']
                data['flux_err'] = sqrt(1./data['flux_err2invsum'])

            if 'Q_sum' in sumlist:

                for i in xrange(data.size):
                    Q_sum = data['Q_sum'][i]
                    Cinv_sum = data['Cinv_sum'][i]
                    C = numpy.linalg.inv(Cinv_sum)
                    g1g2 = numpy.dot(C,Q_sum)

                    data['bashear'][i] = g1g2
                    data['bashear_cov'][i,0,0] = data['shear1err'][i]**2
                    data['bashear_cov'][i,1,1] = data['shear2err'][i]**2
                    print 'bashear1: %.16g +/- %.16g' % (g1g2[0],sqrt(data['bashear_cov'][i,0,0]))


        else:
            data['shear1'] = data['g1sum']/data['nsum']
            data['shear2'] = data['g2sum']/data['nsum']
            data['shear1err'] = sqrt(1/data['g1err2invsum'])
            data['shear2err'] = sqrt(1/data['g2err2invsum'])


        fout=get_averaged_url(new_run_name, i1, fs=fs)
        wlog("    writing:",fout)
        eu.io.write(fout, data, clobber=True)
        if fs=='hdfs':
            fout=get_averaged_url(new_run_name, i1, fs='nfs')
            wlog("    writing:",fout)
            eu.io.write(fout, data, clobber=True)


def average_randshear_outputs(data):
    """
    New version where we specify a shear magnitude but use random orientations
    for each realization.  In this case we average the shear difference from
    truth rather than the shapes

    Set up for bayes only right now
    """
    dt = data[0].dtype.descr
    dt_extra = [('sheardiff','f8'),
                ('osheardiff','f8'),
                ('sheardifferr','f8'),
                ('osheardifferr','f8')]

    dt += dt_extra
    name_extra = [dd[0] for dd in dt_extra]

    d=zeros(len(data),dtype=dt)
    for i,edata in enumerate(data): # over different ellipticities

        twotheta = -arctan2( edata['shear_true'][:,1],edata['shear_true'][:,0] )
        shearmag=sqrt(edata['shear_true'][:,0]**2 + edata['shear_true'][:,1]**2)

        g1 = edata['g'][:,0]
        g2 = edata['g'][:,1]
        shear  = g1*cos(twotheta) - g2*sin(twotheta)
        oshear = g1*sin(twotheta) + g2*cos(twotheta)

        g1sens = edata['gsens'][:,0]
        g2sens = edata['gsens'][:,1]
        g1sensum = g1sens.sum()
        g2sensum = g2sens.sum()

        g1sens_mean=g1sens.mean()
        g2sens_mean=g2sens.mean()

        shearsum = shear.sum()
        oshearsum = oshear.sum()

        # using g1sensum...should be same either way?
        shear = shearsum/g1sensum
        oshear = oshearsum/g1sensum

        sheardiff = shear-shearmag.mean()

        SN=0.0
        g1err2inv = ( 1/(SN**2 + edata['gcov'][:,0,0]) ).sum()
        g2err2inv = ( 1/(SN**2 + edata['gcov'][:,1,1]) ).sum()
        shearerr = sqrt(1/g1err2inv)
        oshearerr = sqrt(1/g1err2inv)

        d['sheardiff'][i] = sheardiff
        d['osheardiff'][i] = oshear
        d['sheardifferr'][i] = shearerr/g1sens_mean
        d['osheardifferr'][i] = oshearerr/g1sens_mean
        print 'shear diff:',d['sheardiff'][i],'+/-',shearerr

        for n in d.dtype.names:
            if n not in name_extra:
                if edata[n].dtype.names is None and len(edata[n].shape) == 1:
                    #d[n][i] = median(edata[n])
                    d[n][i] = edata[n].mean()

    return d



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

def combine_trials(run, is2, ie, allow_missing=True):
    fs=get_default_fs()
    c = read_config(run)
    cs = read_config(c['sim'])

    orient=cs.get('orient','rand')
    if orient == 'ring':
        ntrial = cs['nsplit']
    else:
        ntrial = cs['ntrial']


    outfile=get_output_url(run, is2, ie, fs=fs)

    datalist=[]
    for itrial in xrange(ntrial):
        f=get_output_url(run, is2, ie, itrial=itrial, fs=fs)
        print f
        if allow_missing and not os.path.exists(f):
            continue
        t=eu.io.read(f)
        datalist.append(t)

    data = eu.numpy_util.combine_arrlist(datalist)
    print 'data.size:',data.size
    print 'writing:',outfile
    eu.io.write(outfile, data, clobber=True)

def combine_ctrials(run, is2n, allow_missing=True):
    fs=get_default_fs()
    c = read_config(run)

    npair,nsplit = get_npair_nsplit(c, is2n)
    ntot=npair*2

    outfile=get_output_url(run, 0, is2n, fs=fs, ext='fits')
    print 'writing to:',outfile

    with fitsio.FITS(outfile,mode="rw",clobber=True) as output:

        for isplit in xrange(nsplit):
            f=get_output_url(run, 0, is2n, itrial=isplit, fs=fs)
            if allow_missing and not os.path.exists(f):
                continue
            print f
            t=eu.io.read(f)

            if t.size != ntot:
                raise ValueError("expected %d, got %d" % (npair,t.size))

            if isplit == 0:
                output.write(t)
            else:
                output[-1].append(t)

    
    print 'wrote:',outfile
