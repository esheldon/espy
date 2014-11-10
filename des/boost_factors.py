"""
calculate boost factors
"""
from __future__ import print_function
import numpy
from numpy import sqrt, diag, where, zeros

def calc_boost_factors(data, rand):
    """
    calculate the boost factor and its error

    parameters
    ----------
    data: array
        An array with fields 'wsum_mean' and 'wsum_mean_err'
    rand: array
        An array with fields 'wsum_mean' and 'wsum_mean_err'

    returns
    -------
    corr, corr_err

    corr: array
        The correction. An array with the same shape as the input wsum_mean
    corr_err: array
        The error, an array with the same shape as the input wsum_mean
    """

    corr = data['wsum_mean']/rand['wsum_mean']

    derr = data['wsum_mean_err']
    rerr = rand['wsum_mean_err']

    corr_err = zeros(corr.shape) + 9999

    w=where( (data['wsum_mean'] > 0) & (rand['wsum_mean'] > 0) )
    corr_err[w] = corr[w]*sqrt(   (derr[w]/data['wsum_mean'][w])**2
                                + (rerr[w]/rand['wsum_mean'][w])**2 )

    return corr, corr_err
