"""
Some operations on output files
"""
from __future__ import print_function
import os
import sys
import numpy
from numpy import log10,sqrt,linspace,where
import lensing
import esutil as eu
from esutil.numpy_util import where1


def add_lensums(l1, l2):
    """

    Add the sums from l2 to l1.

    The rows of l1 must correspond to those of l2; zindex must match

    """

    w=where1(l1['zindex'] != l2['zindex'])
    if w.size > 0:
        raise ValueError("zindex do not line up")
    names = ['weight','sshsum','npair','rsum','wsum','dsum','osum']
    if 'totpairs' in l1.dtype.names:
        names.append('totpairs')
    for n in names:
        l1[n] += l2[n]




#
# Codes for combining the lensout "lensum" data and
# getting averages
#

def average_lensums(lout, weights=None):
    """

    combine the lens-by-lens lensums by summing over
    all the individual sums and producing averages

    This uses the averaged_dtype

    This is used by the binner routines

    """

    if weights is not None:
        return average_lensums_weighted(lout,weights)

    nlens = lout.size
    nbin = lout['rsum'][0].size

    comb = averaged_struct(nbin)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses
    comb['weightsum'][0] = lout['weight'].sum()
    comb['sshsum'][0] = lout['sshsum'].sum()
    comb['totpairs'][0] = lout['totpairs'].sum()

    comb['ssh'] = comb['sshsum']/comb['weightsum']

    for i in xrange(nbin):
        npair = lout['npair'][:,i].sum()
        rsum  = lout['rsum'][:,i].sum()


        wsum = lout['wsum'][:,i].sum()
        dsum = lout['dsum'][:,i].sum()
        osum = lout['osum'][:,i].sum()

        comb['npair'][0,i] = npair
        comb['rsum'][0,i] = rsum
        comb['wsum'][0,i] = wsum
        comb['dsum'][0,i] = dsum
        comb['osum'][0,i] = osum

        # averages
        comb['r'][0,i] = rsum/npair
        comb['dsig'][0,i] = dsum/wsum
        comb['osig'][0,i] = osum/wsum

        comb['dsigerr'][0,i] = numpy.sqrt(1.0/wsum)

        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        comb['wsum_mean'][0,i] = wsum/nlens

    return comb

def average_lensums_weighted(lout, weights):
    """

    Reduce the lens-by-lens lensums by summing over
    all the individual sums and producing averages

    """

    nlens = lout.size
    nbin = lout['rsum'][0].size

    if weights.size != nlens:
        raise ValueError("weights not same size as lensout, "
                         "%d instead of %d" % (weights.size,nlens))

    totweights = weights.sum()

    comb = averaged_struct(nbin)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses

    comb['weightsum'][0] = (lout['weight']*weights).sum()
    comb['sshsum'][0] = (lout['sshsum']*weights).sum()

    # should not use this for anything since weights
    # make it non-integer
    comb['totpairs'][0] = lout['totpairs'].sum()

    comb['ssh'] = comb['sshsum']/comb['weightsum']

    for i in xrange(nbin):

        npair = lout['npair'][:,i].sum()
        w_npair = (lout['npair'][:,i]*weights).sum()

        # not weighting?
        rsum  = lout['rsum'][:,i].sum()

        wsum = lout['wsum'][:,i].sum()

        w_wsum = (lout['wsum'][:,i]*weights).sum()
        w_dsum = (lout['dsum'][:,i]*weights).sum()
        w_osum = (lout['osum'][:,i]*weights).sum()


        comb['npair'][0,i] = npair
        comb['rsum'][0,i] = rsum
        comb['wsum'][0,i] = w_wsum
        comb['dsum'][0,i] = w_dsum
        comb['osum'][0,i] = w_osum

        # averages
        comb['r'][0,i] = rsum/npair

        comb['dsig'][0,i] = w_dsum/w_wsum
        comb['osig'][0,i] = w_osum/w_wsum

        # 1/err**2 = sum(1/sigma**2) = n*<1/sigma**2>
        # so for errors by using nlens*(average wsum over lenses)
        wsum_for_error = nlens*w_wsum/totweights
        comb['dsigerr'][0,i] = numpy.sqrt(1.0/wsum_for_error)


        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        comb['wsum_mean'][0,i] = w_wsum/totweights

    return comb



def averaged_struct(nbin, n=1):
    dt = averaged_dtype(nbin)
    return numpy.zeros(n, dtype=dt)

def averaged_dtype(nbin):
    dt=[('ssh','f8'),          # sshsum/weight
        ('weightsum','f8'),       # this is total of weight for each lens
        ('sshsum','f8'),       # this is total of sshsum for each lens
        ('totpairs','i8'),
        ('r','f8',nbin),
        ('dsig','f8',nbin),
        ('dsigerr','f8',nbin),
        ('osig','f8',nbin),
        ('wsum_mean','f8',nbin),
        ('npair','i8',nbin),
        ('rsum','f8',nbin),
        ('wsum','f8',nbin),
        ('dsum','f8',nbin),
        ('osum','f8',nbin)]
    return dt



def averaged_dtype_old():
    dt=[('r','f8'),
        ('dsig','f8'),
        ('dsigerr','f8'),
        ('osig','f8'),
        ('npair','i8'),
        ('weight','f8'),
        ('rsum','f8'),
        ('wsum','f8'),
        ('dsum','f8'),
        ('osum','f8')]
    return numpy.dtype(dt)

