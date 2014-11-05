"""
code to average lens outputs and certain tags
"""
from __future__ import print_function
import numpy
from numpy import sqrt, diag, where, zeros

_DEFAULT_SHEAR_STYLE='reduced'

def lens_wmom(data, tag, ind=None, sdev=False):
    """
    average a tag from a lensum struct using the lensing weights
    """
    import esutil as eu
    if ind is None:
        wts = data['weight']
        tdata = data[tag]
    else:
        wts = data['weight'][ind]
        tdata = data[tag][ind]

    return eu.stat.wmom(tdata, wts, calcerr=True, sdev=sdev)

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
    import jackknife

    if weights is not None:
        return average_lensums_weighted(lout,weights)

    nlens = lout.size
    nbin = lout['rsum'][0].size

    if 'dsensum' in lout.dtype.names:
        shear_style='lensfit'
    else:
        shear_style='reduced'

    comb = averaged_struct(nbin, shear_style=shear_style)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses
    comb['weightsum'][0] = lout['weight'].sum()

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

        if shear_style=='lensfit':
            comb['dsensum'][0,i] = lout['dsensum'][:,i].sum()
            comb['osensum'][0,i] = lout['osensum'][:,i].sum()
            comb['dsig'][0,i] = dsum/comb['dsensum'][0,i]
            comb['osig'][0,i] = osum/comb['osensum'][0,i]
        else:

            comb['dsig'][0,i] = dsum/wsum
            comb['osig'][0,i] = osum/wsum

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
    w=numpy.ones(lout['wsum'].shape)
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

    if 'dsensum' in lout.dtype.names:
        shear_style='lensfit'
    else:
        shear_style='reduced'

    comb = averaged_struct(nbin, shear_style=shear_style)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses

    comb['weightsum'][0] = (lout['weight']*weights).sum()


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

        if shear_style=='lensfit':
            comb['dsensum'][0,i] = lout['dsensum'][:,i].sum()
            comb['osensum'][0,i] = lout['osensum'][:,i].sum()
            comb['dsig'][0,i] = dsum/comb['dsensum'][0,i]
            comb['osig'][0,i] = osum/comb['osensum'][0,i]
        else:

            comb['dsig'][0,i] = dsum/wsum
            comb['osig'][0,i] = osum/wsum

        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        comb['wsum_mean'][0,i] = w_wsum/totweights

        jwsum[:,i] *= weights
        jdsum[:,i] *= weights

    m,cov=jackknife.wjackknife(vsum=jdsum, wsum=jwsum)

    if shear_style=='lensfit':
        dsensum=comb['dsensum']
        for i in xrange(nbin):
            for j in xrange(nbin):
                cov[i,j] *= 1.0/(dsensum[0,i]*dsensum[0,j])

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

def average_ratio(l1, l2):
    """
    jackknife the ratio of dsig
    """
    import jackknife

    comb1=average_lensums(l1)
    comb2=average_lensums(l2)

    nlens, nbin = l1['wsum'].shape

    dsum1=l1['dsum']
    dsum_tot1=comb1['dsum']
    dsum2=l2['dsum']
    dsum_tot2=comb2['dsum']

    if 'dsensum' in l1.dtype.names:
        wsum1=l1['dsensum']
        wsum_tot1=comb1['dsensum']
        wsum2=l2['dsensum']
        wsum_tot2=comb2['dsensum']
    else:
        wsum1=l1['wsum']
        wsum_tot1=comb1['wsum']
        wsum2=l2['wsum']
        wsum_tot2=comb2['wsum']

    r      = comb1['r'][0,:]
    ratio  = comb2['dsig'][0,:]/comb1['dsig'][0,:]
    jmean1 = zeros(nbin)
    jmean2 = zeros(nbin)
    jratio = zeros(nbin)
    jdiff  = zeros(nbin)
    jsum   = zeros((nbin, nbin))

    for i in xrange(nlens):

        jmean1[:] = (dsum_tot1[:]-dsum1[i,:])/(wsum_tot1[:]-wsum1[i,:])
        jmean2[:] = (dsum_tot2[:]-dsum2[i,:])/(wsum_tot2[:]-wsum2[i,:])

        jratio[:] = jmean2/jmean1

        jdiff[:] = jratio[:] - ratio[:]

        # now grab all the cross terms and add to the sum for
        # the covariance matrix
        for ix in xrange(nbin):
            for iy in xrange(ix,nbin):
                val = jdiff[ix]*jdiff[iy]

                jsum[ix,iy] += val
                if ix != iy:
                    jsum[iy,ix] += val

    covar = float(nlens-1)/nlens*jsum

    return r, ratio, covar






def add_lensums(l1, l2):
    """

    Add the sums from l2 to l1.

    The rows of l1 must correspond to those of l2; zindex must match

    Note we now do the combination from the splits using the c program
    redshear, but this could still be useful

    """

    w,=where(l1['index'] != l2['index'])
    if w.size > 0:
        raise ValueError("index do not line up")

    names = ['totpairs','weight','npair','rsum','wsum','dsum','osum']
    if 'dsensum' in l1.dtype.names:
        names+=['dsensum','osensum']

    for n in names:
        if n in l1.dtype.names and n in l2.dtype.names:
            l1[n] += l2[n]

def averaged_struct(nbin, n=1, shear_style=_DEFAULT_SHEAR_STYLE):
    """
    struct to hold the lens averages
    """
    dt = averaged_dtype(nbin, shear_style=shear_style)
    return numpy.zeros(n, dtype=dt)

def averaged_dtype(nbin, shear_style=_DEFAULT_SHEAR_STYLE):
    dt=[('weightsum','f8'),       # this is total of weight for each lens
        ('totpairs','i8'),
        ('r','f8',nbin),
        ('dsig','f8',nbin),
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

    if shear_style=='lensfit':
        dt+=[('dsensum','f8',nbin),
             ('osensum','f8',nbin)]
    return dt


def lensbin_struct(nrbin, shear_style=_DEFAULT_SHEAR_STYLE, bintags=None, n=1):
    """
    a combo of averaged_dtype with extra tags for averaged bin tags
    """
    dt = lensbin_dtype(nrbin, shear_style=shear_style, bintags=bintags)
    return numpy.zeros(n, dtype=dt)

def lensbin_dtype(nrbin, shear_style=_DEFAULT_SHEAR_STYLE, bintags=None):
    """
    This is the same as averaged_dtype but with the averages added
    """
    dt=[('nlenses','i8')]
    if bintags is not None:
        if not isinstance(bintags,list):
            bintags = [bintags]

        for bt in bintags:
            tn = bt+'_range'
            dt.append( (tn,'f8',2) )

            tn = bt+'_minmax'
            dt.append( (tn,'f8',2) )

            tn = bt+'_mean'
            dt.append( (tn,'f8') )

            tn = bt+'_err'
            dt.append( (tn,'f8') )

            tn = bt+'_sdev'
            dt.append( (tn,'f8') )

    dt += averaged_dtype(nrbin, shear_style=shear_style)

    return numpy.dtype(dt)


