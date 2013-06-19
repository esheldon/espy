"""
Importance sampling
"""

import numpy
from numpy import exp, sqrt
import esutil as eu

def test_dist(dist_true, dist_guess,
              nsample=200,
              ntrial=5000):

    dt=[('mean','f8'),('sigma','f8')]
    data=numpy.zeros(ntrial, dtype=dt)
    
    for i in xrange(ntrial):

        # first draw from the guess
        samp=dist_guess.sample(nsample)

        # evaluate prob using true distribution
        ptrue = dist_true.prob(samp)

        # evaluate prob using guess distribution
        pguess = dist_guess.prob(samp)

        weights = ptrue/pguess

        mn,err,sdev=eu.stat.wmom(samp, weights, sdev=True) 

        data['mean'][i] = mn
        data['sigma'][i] = sdev


    trial_mean_mean = data['mean'].mean()
    trial_mean_err  = data['mean'].std()/sqrt(ntrial)

    trial_sigma_mean = data['sigma'].mean()
    trial_sigma_err  = data['sigma'].std()/sqrt(ntrial)

    return trial_mean_mean, trial_mean_err, trial_sigma_mean, trial_sigma_err

def test_normal(mean_true=10.0,
                sigma_true=1.0, 
                mean_guess=12.0,
                sigma_guess=2.0,
                nsample=200,
                ntrial=5000):
    """
    tester using gaussians
    """
    from esutil.random import Normal

    dist_true=Normal(mean_true, sigma_true)
    dist_guess=Normal(mean_guess, sigma_guess)

    res=test_dist(dist_true, dist_guess, nsample=nsample, ntrial=ntrial)

    (trial_mean_mean, 
     trial_mean_err, 
     trial_sigma_mean, 
     trial_sigma_err) = res

    print 'trial mean mean:  %g +/- %g' % (trial_mean_mean,trial_mean_err)
    print 'trial sigma mean: %g +/- %g' % (trial_sigma_mean,trial_sigma_err)
    print 'mean frac err:       %g +/- %g' % (trial_mean_mean/mean_true-1,trial_mean_err/mean_true)
    print 'sigma mean frac err: %g +/- %g' % (trial_sigma_mean/sigma_true-1,trial_sigma_err/sigma_true)

def test_lognormal(mean_true=10.0,
                   sigma_true=1.0, 
                   mean_guess=12.0,
                   sigma_guess=2.0,
                   nsample=200,
                   ntrial=5000):
    """
    tester using gaussians
    """
    from esutil.random import LogNormal

    dist_true=LogNormal(mean_true, sigma_true)
    dist_guess=LogNormal(mean_guess, sigma_guess)

    res=test_dist(dist_true, dist_guess, nsample=nsample, ntrial=ntrial)

    (trial_mean_mean, 
     trial_mean_err, 
     trial_sigma_mean, 
     trial_sigma_err) = res

    print 'trial mean mean:  %g +/- %g' % (trial_mean_mean,trial_mean_err)
    print 'trial sigma mean: %g +/- %g' % (trial_sigma_mean,trial_sigma_err)
    print 'mean frac err:       %g +/- %g' % (trial_mean_mean/mean_true-1,trial_mean_err/mean_true)
    print 'sigma mean frac err: %g +/- %g' % (trial_sigma_mean/sigma_true-1,trial_sigma_err/sigma_true)


