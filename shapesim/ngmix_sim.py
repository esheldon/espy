"""
Bernstein & Armstrong using the ngmix code
"""

import os
from sys import stderr
import ngmix
import numpy
import time

from .shapesim import read_config

NSIGMA_RENDER=5.0

class TryAgainError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

class NGMixSim(dict):
    def __init__(self, run, s2n, npairs, **keys):
        """
        Simulate and fit the requested number of pairs at
        the specified s/n
        """

        self.run=run
        self.s2n=s2n
        self.npairs=npairs

        self.update(keys)
        self.conf = read_config(run)
        self.update(self.conf)
        self.simc = read_config(self.conf['sim'])
        self.shear=self.simc['shear']
        self.nsub=self.simc['nsub']

        self.obj_model=self.simc['obj_model']

        self.make_struct()
        self.set_priors()
        self.make_psf()
        self.set_noise()

    def run_sim(self):
        """
        Run the simulation, fitting psf and all pairs
        """
        self.fit_psf()

        tm=time.time()
        i=0
        npairs=self.npairs
        for ipair in xrange(npairs):
            print >>stderr,'%s/%s' % (ipair+1,npairs)
            while True:
                try:
                    reslist=self.process_pair()
                    break
                except TryAgainError as err:
                    print >>stderr,str(err)

            self.copy_to_output(reslist[0], i)
            i += 1
            self.copy_to_output(reslist[1], i)
            i += 1

        tm=time.time()-tm
        print >>stderr,'time per image:',tm/(2*npairs)

    def process_pair(self):
        """
        Create a simulated image pair and perform the fit
        """

        imdicts = self.get_noisy_image_pair()
        reslist=[]
        for key in imdicts:
            res=self.fit_galaxy(imdicts[key])
            if res['flags'] != 0:
                raise TryAgainError("failed at %s" % key)

            reslist.append(res)
            self.print_res(res)

        return reslist

    def fit_galaxy(self, imdict):
        """
        Fit the model to the galaxy
        """

        gm=imdict['gm_pre']
        T_guess      = gm.get_T()
        counts_guess = gm.get_psum()

        fitter=ngmix.fitting.MCMCSimple(imdict['image'],
                                        imdict['wt'],
                                        imdict['jacobian'],
                                        self.obj_model,

                                        cen_prior=self.cen_prior,
                                        g_prior=self.g_prior,
                                        T_prior=self.T_prior,
                                        counts_prior=self.counts_prior,

                                        T_guess=T_guess,
                                        counts_guess=counts_guess,

                                        psf=self.psf_gmix_fit,
                                        nwalkers=self['nwalkers'],
                                        nstep=self['nstep'],
                                        burnin=self['burnin'],
                                        mca_a=self['mca_a'],
                                        do_pqr=True,
                                        do_lensfit=True)
        fitter.go()
        #fitter.make_plots(show=True)
        return fitter.get_result()

    def print_res(self,res):
        """
        print some stats
        """
        print >>stderr,'    arate:',res['arate']
        ngmix.fitting.print_pars(res['pars'],front='    pars: ',stream=stderr)
        ngmix.fitting.print_pars(res['perr'],front='    perr: ',stream=stderr)

    def fit_psf(self):
        """
        Fit the pixelized psf to a model
        """
        from ngmix.gexceptions import GMixRangeError, GMixMaxIterEM

        print >>stderr,'fitting psf'
        imsky,sky=ngmix.em.prep_image(self.psf_image)

        em=ngmix.em.GMixEM(imsky)
        guess=self.psf_gmix_true.copy()
        print 'psf guess:'
        print guess
        em.go(guess, sky, tol=1.e-5)

        self.psf_gmix_fit=em.get_gmix()
        print 'psf fit:'
        print self.psf_gmix_fit

    def set_priors(self):
        """
        Set all the priors
        """

        print >>stderr,"setting priors"
        T=self.simc['obj_T_mean']
        T_sigma = self.simc['obj_T_sigma_frac']*T
        counts=self.simc['obj_counts_mean']
        counts_sigma = self.simc['obj_counts_sigma_frac']*counts

        self.g_prior=ngmix.priors.GPriorBA(0.3)
        self.cen_prior=ngmix.priors.CenPrior(0.0, 0.0, 0.1, 0.1)
        self.T_prior=ngmix.priors.LogNormal(T, T_sigma)
        self.counts_prior=ngmix.priors.LogNormal(counts, counts_sigma)

    def make_psf(self):
        """
        make the psf gaussian mixture model
        """

        print >>stderr,"making psf"

        self.psf_dims, self.psf_cen=self.get_dims_cen(self.simc['psf_T'])

        pars=[self.psf_cen[0],
              self.psf_cen[1],
              self.simc['psf_shape'][0],
              self.simc['psf_shape'][1],
              self.simc['psf_T'],
              1.0]
        self.psf_gmix_true=ngmix.gmix.GMixModel(pars, self.simc['psf_model'])
        
        self.psf_image=self.psf_gmix_true.make_image(self.psf_dims,
                                                     nsub=self.nsub)
    
    def set_noise(self):
        """
        Find gaussian noise that when added to the image 
        produces the requested s/n.  Use a matched filter.

         sum(pix^2)
        ------------ = S/N^2
          skysig^2

        thus
            
        sum(pix^2)
        ---------- = skysig^2
          (S/N)^2
        """
        
        from numpy.random import randn

        print >>stderr,"setting noise"

        imdict=self.get_image_pair(random=False)
        im=imdict['im1']['image']
        skysig2 = (im**2).sum()/self.s2n**2
        skysig = numpy.sqrt(skysig2)

        noise_image = skysig*randn(im.size).reshape(im.shape)
        new_im = im + noise_image

        s2n_check = numpy.sqrt( (im**2).sum()/skysig**2 )
        print >>stderr,"S/N goal:",self.s2n,"found:",s2n_check

        self.skysig=skysig
        self.ivar=1.0/skysig**2


    def get_noisy_image_pair(self, random=True):
        """
        Get an image pair, with noise added
        """
        imdict=self.get_image_pair(random=random)
        self.add_noise(imdict['im1']['image'])
        self.add_noise(imdict['im2']['image'])

        wt=numpy.zeros(imdict['im1']['image'].shape) + self.ivar
        imdict['im1']['wt']=wt
        imdict['im2']['wt']=wt
        return imdict

    def add_noise(self, im):
        """
        Add gaussian random noise
        """

        from numpy.random import randn
        im[:,:] += self.skysig*randn(im.size).reshape(im.shape)

    def get_image_pair(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors
        """

        cen_offset, shape1, shape2, T, counts=self.get_pair_pars(random=random)

        # center is just placeholder for now
        pars1=[0.0, 0.0, shape1.g1, shape1.g2, T, counts]
        pars2=[0.0, 0.0, shape2.g1, shape2.g2, T, counts]

        gm1_pre=ngmix.gmix.GMixModel(pars1, self.obj_model)
        gm2_pre=ngmix.gmix.GMixModel(pars2, self.obj_model)

        gm1  = gm1_pre.convolve(self.psf_gmix_true)
        gm2  = gm2_pre.convolve(self.psf_gmix_true)

        T = gm1.get_T()
        dims, cen = self.get_dims_cen(T)

        # jacobian is at center before offset
        j=ngmix.jacobian.UnitJacobian(cen[0], cen[1])

        cen[0] += cen_offset[0]
        cen[1] += cen_offset[1]

        gm1.set_cen(cen[0], cen[1])
        gm2.set_cen(cen[0], cen[1])

        nsub = self.nsub
        im1=gm1.make_image(dims, nsub=nsub)
        im2=gm2.make_image(dims, nsub=nsub)

        out={'im1':{'gm_pre':gm1_pre,'gm':gm1,'image':im1,'jacobian':j},
             'im2':{'gm_pre':gm2_pre,'gm':gm2,'image':im2,'jacobian':j}}
        return out

    def get_pair_pars(self, random=False):
        """
        Get pair parameters
        """
        from numpy.random import random as randu

        if random:
            cen_offset=self.cen_prior.sample()
            g = self.g_prior.sample1d(1)
            g=g[0]
            rangle1 = randu()*2*numpy.pi
            rangle2 = rangle1 + numpy.pi/2.0
            g1_1 = g*numpy.cos(rangle1)
            g2_1 = g*numpy.sin(rangle1)
            g1_2 = g*numpy.cos(rangle2)
            g2_2 = g*numpy.sin(rangle2)

            T=self.T_prior.sample()
            counts=self.counts_prior.sample()
        else:
            cen_offset=[0.0, 0.0]
            g1_1=0.0
            g2_1=0.0
            g1_2=0.0
            g2_2=0.0
            T=self.T_prior.mean
            counts=self.counts_prior.mean

        shape1=ngmix.shape.Shape(g1_1, g2_1)
        shape2=ngmix.shape.Shape(g1_2, g2_2)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        return cen_offset, shape1, shape2, T, counts

    def get_dims_cen(self, T):
        """
        Based on T, get the required dimensions and a center
        """
        sigma=numpy.sqrt(T/2.)
        dims = [2.*sigma*NSIGMA_RENDER]*2
        cen = [(dims[0]-1.)/2.]*2

        return dims, cen

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data

    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        d=self.data
        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']
        d['P'][i] = res['P']
        d['Q'][i,:] = res['Q']
        d['R'][i,:,:] = res['R']
        d['g'][i,:] = res['g']
        d['gsens'][i,:] = res['g_sens']

    def make_struct(self):
        """
        Make the output array
        """
        dt=[('pars','f8',6),
            ('pcov','f8',(6,6)),
            ('P','f8'),
            ('Q','f8',2),
            ('R','f8',(2,2)),
            ('g','f8',2),
            ('gsens','f8',2)]
        self.data=numpy.zeros(self.npairs*2, dtype=dt)

def srandu(num=None):
    """
    Generate random numbers in the symmetric distribution [-1,1]
    """
    return 2*(numpy.random.random(num)-0.5)


