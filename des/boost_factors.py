"""
calculate boost factors, corrections for clustering of lenses
with sources
"""
from __future__ import print_function
import numpy
from numpy import sqrt, diag, where, zeros

def add_boost_factors(data, rand):
    """
    calculate the boost factor and its error, multiply
    dsig etc. by the factor.

    parameters
    ----------
    data: array
        An array with fields 'wsum_mean' and 'wsum_mean_err'
    rand: array
        An array with fields 'wsum_mean' and 'wsum_mean_err'

    returns
    -------
    corrected struct
    """
    import jackknife
    import esutil as eu

    nrad=data['dsig'].shape[1]


    add_dt=[('boost','f8',nrad),
            ('boost_err','f8',nrad)]
    odata=eu.numpy_util.add_fields(data, add_dt)

    corr, corr_err = calc_boost_factors(data,rand)

    odata['boost']     = corr
    odata['boost_err'] = corr_err

    # apply clipped version
    cclip=corr.clip(min=1.0)
    cclip2 = cclip**2

    mcols=['dsum','osum','dsig','dsigerr','osig']
    mcols2=['dsigcov']

    for col in mcols:
        odata[col] *= cclip

    for col in mcols2:
        odata[col] *= cclip2[:,:,numpy.newaxis]

    nbin=odata.size
    for i in xrange(nbin):
        odata['dsigcor'][i,:,:] = jackknife.covar2corr(odata['dsigcov'][i,:,:])

    return odata

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

    if w[0].size > 0:
        dfrac = derr[w]/data['wsum_mean'][w]
        rfrac = rerr[w]/rand['wsum_mean'][w]
        corr_err[w] = corr[w]*sqrt( dfrac**2 + rfrac**2 )

    return corr, corr_err
