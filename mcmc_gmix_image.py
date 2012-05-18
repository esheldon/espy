import os
import numpy
from numpy import zeros, array, where, isfinite, sqrt, pi, log, linspace
from numpy.random import randn
from numpy.random import random

import fimage
from fimage import model_image
import gmix_image
from gmix_image import ogrid_image
import mcmc
import esutil as eu

class MCMCCoellip:
    """
    Fit a gaussian mixture model to an image using MCMC.

    Force all gaussians to be same center and ellipticity, for now

    A gaussian mixture model is represented by a list of dictionaries with the
    gaussian parameters.

    construction parameters
    -----------------------
    image: numpy array
        A two dimensional, 64-bit float numpy array.  If the image
        is a different type it will be converted.
    sky: 
        sky level
    skysig: 
        sigma of sky in image

    prior: array
        prior mean parameters

        row: center of the whole model
        col: center of the whole model
        irr: Covariance matrix element for row*row for first gauss
        irc: Covariance matrix element for row*col for first gauss
        icc: Covariance matrix element for col*col for first gauss
        pi: A normalization for each gaussian
        fi: relative size of covariance matrix for each gauss

        For three gaussians
        [row,col,irr,irc,icc,p1,p2,p3,f1,f2,f3]

        For N gaussians, there are 5+2*N parameters

        The normalizations will be recalculated internally such that 
            sum(p) == 1
    widths: array
        Widths of priors by type

            [cen,cov,p,f]

    stepsize: array
        step widths for each parameter.

        A sigma value for each type of parameter.

            [censig,covsig,psig,fsig]

    """
    def __init__(self, image, sky, skysig, prior, width, stepsize):
        self._image = image
        self._model = zeros(image.shape)
        self._counts = (image-sky).sum()

        self._sky = sky
        self._skysig = skysig

        self._prior = array(prior,copy=True,dtype='f8')
        self._npars = self._prior.size
        ngauss = (len(self._prior) - 5)/2
        self._ngauss = ngauss
        self._nsub=1

        self._width_short = array(width,copy=True,dtype='f8')

        self._width=zeros(len(self._prior))
        self._width[0:2] = width[0]
        self._width[2:5] = width[1]
        self._width[5:5+ngauss] = width[2]
        self._width[5+ngauss:5+ngauss*2] = width[3]


        self._stepsize_short = array(stepsize,copy=True,dtype='f8')

        self._stepsize=zeros(len(self._prior))
        self._stepsize[0:2] = stepsize[0]
        self._stepsize[2:5] = stepsize[1]
        self._stepsize[5:5+ngauss] = stepsize[2]
        self._stepsize[5+ngauss:5+ngauss*2] = stepsize[3]

    def step(self,pars):

        newpars = pars + self._stepsize*randn(self._npars)
        if self._ngauss == 1:
            # for one gaussian, f is 1, avoid degeneracy
            newpars[6] = 1
        # keep normalization
        newpars[5:5+self._ngauss] /= newpars[5:5+self._ngauss].sum()

        return newpars

    def likelihood(self, pars):
        """
        log likelihood
        """
        det = pars[2]*pars[4]-pars[3]**2
        if (det < 0 
                or pars[0] < 0 or pars[0] > (self._image.shape[0]-1)
                or pars[1] < 0 or pars[1] > (self._image.shape[1]-1)):
            chi2=9.999e9
        else:
            model = self.make_model(pars)
            chi2 = ( (self._image - model)**2 ).sum()/self._skysig**2
        loglike = -0.5*chi2
        priors = self.priors(pars)
        return loglike + priors

    def priors(self, pars):
        """
        1-d priors
            exp(-0.5*(par-prior)**2/2/sigma**2)/sqrt(2*pi*sigma**2)
        """
        sigma = self._width
        priors = -0.5*(pars-self._prior)**2/sigma**2 - 0.5*log(2*pi*sigma**2)
        priors = priors.sum()
        return priors

    def make_model(self, pars):
        """
        pars = [row,col,irr,irc,icc,p1,p2,p3,f1,f2,f3]
        """
        self._model[:] = self._sky
        psum=0.0
        cen = pars[0:2]
        for i in xrange(self._ngauss):
            pi = pars[5+i]
            fi = pars[5+self._ngauss+i]
            covari = pars[2:5]*fi
            self._model += model_image('gauss',
                                       self._image.shape,
                                       cen,covari,
                                       counts=pi*self._counts, # counts over sky
                                       nsub=self._nsub)
            psum += pi

        self._model /= psum
        return self._model

    def guess(self):
        return self._guess


def test_mcmc1(nstep=10000, burnin=1000, ntrial=1, s2n=35.0):
    """
    S/N is adaptive weighted S/N
    """
    sigma = sqrt(2)
    covar=[sigma**2,0.0,sigma**2]
    ngauss=1
    dim=int(2*4*sigma)
    if (dim % 2) == 0:
        dim += 1
    dims=array([dim,dim])
    cen=(dims-1)/2.
    det = covar[0]*covar[2] - covar[1]**2


    counts=1.0
    sky=0.0
    #skysig = counts/s2n/sqrt(4*pi*det)
    skysig = counts/sqrt(dims[0]*dims[1])/s2n

    print 'dims:  ',dims
    print 'cen:   ',cen
    print 'skysig:',skysig

    im0 = model_image('gauss',
                      dims,
                      cen,covar,
                      counts=counts,
                      nsub=1)
    allmeans = zeros( (ntrial, 2+3+2*ngauss) )
    allerrs  = zeros( (ntrial, 2+3+2*ngauss) )
    for j in xrange(ntrial):
        print '-'*70
        print '%d/%d' % ((j+1),ntrial)
        im = im0 + skysig*randn(dims[0]*dims[1]).reshape(dims)

        # guess is prior
        guess=zeros(2+3+2*ngauss)
        guess[0:2] = cen + 0.1*(random(2)-0.5)
        guess[2:5] = covar + 0.1*(random(3)-0.5)
        guess[5] = 1.0
        guess[6] = 1.0

        print 'guess:',guess


        # prior widths, generally broad
        width=array([0.1, # pixels
                     1.0, # pixels**2
                     1.0, 
                     1.0])

        # step widths
        #[censig,covsig,psig,fsig]
        stepsize=array([0.01,0.01,0.01,0.01])
        
        obj=MCMCCoellip(im, sky, skysig, guess, width, stepsize)
        m=mcmc.MCMC(obj)

        res = m.run(nstep, guess)

        means, errs = mcmc.extract_stats(res, burnin,sigma_clip=False)
        print 'means +/- err'
        for i in xrange(len(means)):
            print '  %.16g +/- %.16g' % (means[i],errs[i])

        allmeans[j,:] = means
        allerrs[j,:] = errs
     
    return allmeans, allerrs

def test_admom_multi():
    covar=[3.0,0.0,3.0]
    ntrial=10000
    n_s2n=20
    s2n_array = linspace(20.0,100.0,n_s2n)
    dt=[('s2n','f8'),
        ('Irr','f8'),('Irc','f8'),('Icc','f8'),
        ('Irr_mean','f8'),('Irc_mean','f8'),('Icc_mean','f8')]
    data=zeros(n_s2n,dtype=dt)
    for i,s2n in enumerate(s2n_array):
        data['s2n'][i] = s2n
        data['Irr'][i] = covar[0]
        data['Irc'][i] = covar[1]
        data['Icc'][i] = covar[2]

        res=test_admom(covar=covar, ntrial=ntrial, s2n=s2n)
        data['Irr_mean'][i] = res[:,0].mean()
        data['Irc_mean'][i] = res[:,1].mean()
        data['Icc_mean'][i] = res[:,2].mean()
    eu.io.write('~/tmp/test-admom-bias.fits',data,clobber=True)
    return data

def test_admom(covar, ntrial=1, s2n=35.0):
    print 'covar:',covar,'s2n:',s2n
    import admom
    nsub=1
    ngauss=1
    dims=array([21,21])
    cen=(dims-1)/2.
    det = covar[0]*covar[2] - covar[1]**2

    counts=1.0
    sky=0.0
    # weighted admom s2n of about 35
    #skysig=0.004
    skysig = counts/s2n/sqrt(4*pi*det)

    s2n_exp = counts/skysig/sqrt(4*pi*det)

    im0 = model_image('gauss',
                     dims,
                     cen,covar,
                     counts=counts,
                     nsub=nsub)
    allmeans = zeros( (ntrial, 3) )
    for j in xrange(ntrial):
        im = im0 + skysig*randn(dims[0]*dims[1]).reshape(dims)

        guess = (covar[1]+covar[2])/2 + 0.1*(random(2)-0.5)
        res = admom.admom(im, cen[0], cen[1], sigsky=skysig, nsub=nsub)

        #print s2n_exp,res['s2n']
        allmeans[j,:] = [res['Irr'],res['Irc'],res['Icc']]
     
    return allmeans

def test_coellip2(nstep=10000, burnin=1000, s2n=None):
    ngauss=2
    dims=array([21,21])
    cen=(dims-1)/2.
    covar1=[2.0,0.0,2.0]
    covar2=[2.5,0.0,2.5]
    counts1=1.0
    counts2=0.3
    sky=0.0
    skysig=0.001

    im1 = model_image('gauss',
                     dims,
                     cen,covar1,
                     counts=counts1,
                     nsub=1)
    im2 = model_image('gauss',
                     dims,
                     cen,covar2,
                     counts=counts2,
                     nsub=1)

    im = im1 + im2 + skysig*randn(dims[0]*dims[1]).reshape(dims)

    guess=zeros(2+3+2*ngauss)
    guess[0:2] = cen + 0.1*(random(2)-0.5)
    guess[2:5] = covar + 0.1*(random(3)-0.5)
    guess[5:7] = array([0.5,0.5]) # p vals
    guess[7:9] = array([1.0,1.0]) # f vals

    print 'guess:',guess

    # step widths
    #[censig,covsig,psig,fsig]
    stepsize=array([0.01,0.01,0.01,0.01])

    # prior widths
    width=array([0.2,0.2,0.2,0.2])

    obj=MCMCCoellip(im, sky, skysig, guess, width, stepsize)
    m=mcmc.MCMC(obj)

    res = m.run(nstep, guess)

    means, errs = mcmc.extract_stats(res, burnin,sigma_clip=False)
    print 'means +/- err'
    for i in xrange(len(means)):
        print '  %.16g +/- %.16g' % (means[i],errs[i])

     
    return res
