"""
Some operations on output files
"""
from __future__ import print_function
import os
import sys
import numpy
from numpy import log10,sqrt,linspace,where,diag,zeros,ones
import lensing
import esutil as eu
from esutil.numpy_util import where1

import jackknife

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
    if 'sshsum' in lout.dtype.names:
        comb['sshsum'][0] = lout['sshsum'].sum()
        comb['ssh'] = comb['sshsum']/comb['weightsum']
    else:
        # this makes ssh unity
        comb['sshsum'][0] = comb['weightsum'][0]
        comb['ssh'] = comb['sshsum']/comb['weightsum']

    comb['totpairs'][0] = lout['totpairs'].sum()


    for i in xrange(nbin):
        npair = lout['npair'][:,i].sum()
        rsum  = lout['rsum'][:,i].sum()

        wsum = lout['wsum'][:,i].sum()
        wsum2 = (lout['wsum'][:,i]**2).sum()
        dsum = lout['dsum'][:,i].sum()
        osum = lout['osum'][:,i].sum()

        comb['npair'][0,i] = npair
        comb['rsum'][0,i] = rsum
        comb['wsum'][0,i] = wsum
        comb['dsum'][0,i] = dsum
        comb['osum'][0,i] = osum

        # averages
        comb['r'][0,i] = rsum/wsum
        comb['dsig'][0,i] = dsum/wsum
        comb['osig'][0,i] = osum/wsum

        comb['dsigerr_simple'][0,i] = numpy.sqrt(1.0/wsum)

        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        wsum_mean = wsum/nlens
        comb['wsum_mean'][0,i] = wsum_mean
        comb['wsum_mean_err_simple'][0,i] = sqrt(wsum2/nlens - wsum_mean**2)/sqrt(nlens)

    m,cov=jackknife.wjackknife(vsum=lout['dsum'], wsum=lout['wsum'])
    comb['dsigcov'][0,:,:] = cov
    comb['dsigcor'][0,:,:] = jackknife.covar2corr(cov)
    comb['dsigerr'][0,:] = sqrt(diag(cov))

    # turns out this agrees with the above one
    w=ones(lout['wsum'].shape)
    m,cov = jackknife.wjackknife(vsum=lout['wsum'], wsum=w)
    comb['wsum_mean_err'][0,:] = sqrt(diag(cov))

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

    if 'sshsum' in lout.dtype.names:
        comb['sshsum'][0] = (lout['sshsum']*weights).sum()
        comb['ssh'] = comb['sshsum']/comb['weightsum']
    else:
        # this makes ssh unity
        comb['sshsum'][0] = comb['weightsum'].sum()
        comb['ssh'] = comb['sshsum']/comb['weightsum']

    # should not use this for anything since weights
    # make it non-integer
    comb['totpairs'][0] = lout['totpairs'].sum()


    # these will get the extra weight
    # WE MUST MAKE A COPY!! ARRGGG THIS BIT ME!!!
    jwsum = lout['wsum'].copy()
    jdsum = lout['dsum'].copy()
    for i in xrange(nbin):

        npair = lout['npair'][:,i].sum()

        rsum  = lout['rsum'][:,i].sum()

        w_wsum = (lout['wsum'][:,i]*weights).sum()
        w_dsum = (lout['dsum'][:,i]*weights).sum()
        w_osum = (lout['osum'][:,i]*weights).sum()


        comb['npair'][0,i] = npair
        comb['rsum'][0,i]  = rsum
        comb['wsum'][0,i]  = w_wsum
        comb['dsum'][0,i]  = w_dsum
        comb['osum'][0,i]  = w_osum

        # averages
        comb['r'][0,i] = rsum/w_wsum

        comb['dsig'][0,i] = w_dsum/w_wsum
        comb['osig'][0,i] = w_osum/w_wsum


        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        comb['wsum_mean'][0,i] = w_wsum/totweights

        jwsum[:,i] *= weights
        jdsum[:,i] *= weights

    m,cov=jackknife.wjackknife(vsum=jdsum, wsum=jwsum)
    comb['dsigcov'][0,:,:] = cov
    comb['dsigcor'][0,:,:] = jackknife.covar2corr(cov)
    comb['dsigerr'][0,:] = sqrt(diag(cov))

 
    # make weights in shape of wsum
    w      = zeros(lout['wsum'].shape)
    w_wsum = zeros(lout['wsum'].shape)
    for i in xrange(lout.size):
        w_wsum[i,:] = lout['wsum'][i,:]*weights[i]
        w[i,:] = weights[i]
    m,cov = jackknife.wjackknife(vsum=w_wsum, wsum=w)
    comb['wsum_mean_err'][0,:] = sqrt(diag(cov))


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
        ('dsigerr_simple','f8',nbin),
        ('dsigerr','f8',nbin),
        ('dsigcov','f8',(nbin,nbin)),
        ('dsigcor','f8',(nbin,nbin)),
        ('osig','f8',nbin),
        ('wsum_mean','f8',nbin),
        ('wsum_mean_err','f8',nbin),
        ('wsum_mean_err_simple','f8',nbin),
        ('npair','i8',nbin),
        ('rsum','f8',nbin),
        ('wsum','f8',nbin),
        ('dsum','f8',nbin),
        ('osum','f8',nbin)]
    return dt

def add_lensums(l1, l2):
    """

    Add the sums from l2 to l1.

    The rows of l1 must correspond to those of l2; zindex must match

    Note we now do the combination from the splits using the c program
    redshear, but this could still be useful

    """

    w=where1(l1['zindex'] != l2['zindex'])
    if w.size > 0:
        raise ValueError("zindex do not line up")
    names = ['weight','sshsum','npair','rsum','wsum','dsum','osum']
    if 'totpairs' in l1.dtype.names:
        names.append('totpairs')
    for n in names:
        if n in l1.dtype.names and n in l2.dtype.names:
            l1[n] += l2[n]





def compare_hist_match(binsize, binnum=9, l=None, r=None):
    """
    Compare hist_match with weights and with remove method.
    """
    import weighting
    import lensing

    if l is None or r is None:
        l,r = load_test_data()

    binner=lensing.binning.LambdaBinner(12)
    print("selecting ",binner.bin_label(binnum))
    w=binner.select_bin(l, binnum)
    print("    kept %d/%d" % (w.size,l.size))

    print("using remove method")
    keep = weighting.hist_match_remove(r['z'], l['z'][w], binsize)
    perc = keep.size/(1.*r.size)
    print("    used number: %d/%d = %0.2f" % (keep.size,r.size, perc))

    print("    combining")
    comb_rm = average_lensums(r[keep])


    for i in xrange(comb_rm['r'].size):
        print("%0.3f %15.12f %15.12f" % \
              (comb_rm['r'][0,i], comb_rm['dsig'][0,i], comb_rm['wsum_mean'][0,i]))

    print("using weights method")
    weights = weighting.hist_match(r['z'], l['z'][w], binsize)
    effnum = weights.sum()
    effperc = effnum/r.size
    print("    effective number: %d/%d = %0.2f" % (effnum,r.size, effperc))
    print("    combining")
    comb = average_lensums(r, weights=weights)

    for i in xrange(comb_rm['r'].size):
        print("%0.3f %15.12f %15.12f %15.12f %15.12f" % \
              (comb['r'][0,i], 
              comb['dsig'][0,i], comb_rm['dsig'][0,i],
              comb['wsum_mean'][0,i], comb_rm['wsum_mean'][0,i]))

    
def load_test_data():
    import lensing
    run='08'
    rrun='r01'
    l=lensing.files.sample_read(type='collated',sample=run)
    r=lensing.files.sample_read(type='collated',sample=rrun)

    return l,r


