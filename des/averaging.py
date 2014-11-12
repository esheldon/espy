"""
code to average lens outputs and certain tags
"""
from __future__ import print_function
import numpy
from numpy import sqrt, diag, where, zeros, ones, newaxis
from .files import get_shear_style

_DEFAULT_SHEAR_STYLE='reduced'

#
# Codes for combining the lensout "lensum" data and
# getting averages
#

def average_lensums(lout, weights=None, jackreg_col=None):
    """
    average over all the individual lensums

    The covariance matrix is estimated from jackknifing. If regions are sent,
    use them for jackknifing, otherwise jackknife one object at a time

    parameters
    ----------
    data: array
        Array containing the outputs from xshear
    weights: array, optional
        Optional weights
    jackreg_col: string, optional
        column name holding the jackknife region ids
    """
    import jackknife

    if weights is not None:
        return average_lensums_weighted(lout,weights,jackreg_col=jackreg_col)

    nlens = lout.size
    nrad = lout['rsum'][0].size

    shear_style=get_shear_style(lout)

    comb = averaged_struct(nrad, shear_style=shear_style)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses
    comb['weightsum'] = lout['weight'].sum()

    comb['totpairs'] = lout['totpairs'].sum()
    
    comb['npair'] = lout['npair'].sum(axis=0)
    comb['rsum']  = lout['rsum'].sum(axis=0)
    comb['wsum']  = lout['wsum'].sum(axis=0)
    comb['dsum']  = lout['dsum'].sum(axis=0)
    comb['osum']  = lout['osum'].sum(axis=0)

    # averages
    comb['r'] = comb['rsum']/comb['wsum']

    if shear_style=='lensfit':
        comb['dsensum'] = lout['dsensum'].sum(axis=0)
        comb['osensum'] = lout['osensum'].sum(axis=0)
        comb['dsig'] = comb['dsum']/comb['dsensum']
        comb['osig'] = comb['osum']/comb['osensum']
    else:
        comb['dsig'] = comb['dsum']/comb['wsum']
        comb['osig'] = comb['osum']/comb['wsum']

    # this is average wsum over lenses
    # we calculate boost factors from this, wsum_mean/wsum_mean_random
    comb['wsum_mean'] = comb['wsum']/nlens

    # jackknife to get the covariance matrix
    m,cov=jackknife_lensums(lout, jackreg_col=jackreg_col)

    comb['dsigcov'][0] = cov
    comb['dsigcor'][0] = jackknife.covar2corr(cov)
    comb['dsigerr'][0] = sqrt(diag(cov))

    # also jackknife the wsums, used for errors on boost factors
    w=ones(lout['wsum'].shape)
    m,cov = jackknife.wjackknife(vsum=lout['wsum'], wsum=w)

    comb['wsum_mean_err'][0] = sqrt(diag(cov))
    return comb

def average_lensums_weighted(lout, weights_in, jackreg_col=None):
    """
    average over all the individual lensums with additional weights

    The covariance matrix is estimated from jackknifing. If regions are sent,
    use them for jackknifing, otherwise jackknife one object at a time

    parameters
    ----------
    data: array
        Array containing the outputs from xshear
    weights: array
        Additional weights
    jackreg_col: string, optional
        column name holding the jackknife region ids
    """

    import jackknife

    nlens = lout.size
    nrad = lout['rsum'][0].size

    weights=weights_in.copy()
    weights *= (1.0/weights.max())

    if weights.size != nlens:
        raise ValueError("weights not same size as lensout, "
                         "%d instead of %d" % (weights.size,nlens))

    totweights = weights.sum()

    shear_style=get_shear_style(lout)

    # broadcast it
    wa=weights[:,newaxis]
    comb = averaged_struct(nrad, shear_style=shear_style)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses
    comb['weightsum'] = (lout['weight']*weights).sum()

    comb['totpairs'] = lout['totpairs'].sum()

    comb['wsum']  = (lout['wsum']*wa).sum(axis=0)
    comb['npair'] = lout['npair'].sum(axis=0)
    comb['rsum']  = (lout['rsum']*wa).sum(axis=0)
    comb['dsum']  = (lout['dsum']*wa).sum(axis=0)
    comb['osum']  = (lout['osum']*wa).sum(axis=0)

    # averages
    comb['r'] = comb['rsum']/comb['wsum']

    if shear_style=='lensfit':
        comb['dsensum'] = (lout['dsensum']*wa).sum(axis=0)
        comb['osensum'] = (lout['osensum']*wa).sum(axis=0)
        comb['dsig'] = comb['dsum']/comb['dsensum']
        comb['osig'] = comb['osum']/comb['osensum']
    else:
        comb['dsig'] = comb['dsum']/comb['wsum']
        comb['osig'] = comb['osum']/comb['wsum']

    # we calculate boost factors from this, wsum_mean/wsum_mean_random
    comb['wsum_mean'] = comb['wsum']/totweights

    # jackknife to get the covariance matrix
    m,cov=jackknife_lensums(lout, weights=weights, jackreg_col=jackreg_col)

    comb['dsigcov'][0] = cov
    comb['dsigcor'][0] = jackknife.covar2corr(cov)
    comb['dsigerr'][0] = sqrt(diag(cov))
 

    #
    # also jackknife the wsums, used for errors on boost factors
    #

    # this will broadcase
    w_wsum_all = lout['wsum']*wa

    # but here we need fully expanded version
    weights_big=ones( (lout.size, nrad) )*wa

    m,cov = jackknife.wjackknife(vsum=w_wsum_all, wsum=weights_big)

    comb['wsum_mean_err'][0] = sqrt(diag(cov))
    return comb

def jackknife_lensums(data, jackreg_col=None, weights=None):
    """
    jackknife the data. If regions are sent, use them for jackknifing,
    otherwise jackknife one object at a time

    parameters
    ----------
    data: array
        An array with fields 'dsum' and 'wsum'. If shear style
        is lensfit, dsensum is needed rather than wsum.
    jackreg_col: string, optional
        column name holding the jackknife region ids
    weights: array, optional
        Additional weights per lens

    returns
    -------
    dsig, dsig_cov

    dsig: array
        The delta sigma in radial bins [nrad]
    dsig_cov: array
        The covariance matrix of delta sigma [nrad,nrad]
    """
    import jackknife

    jdsum, jwsum = get_jackknife_sums(data,
                                      weights=weights,
                                      jackreg_col=jackreg_col)
    dsig,dsig_cov=jackknife.wjackknife(vsum=jdsum, wsum=jwsum)

    return dsig, dsig_cov

def get_jackknife_sums(data, jackreg_col=None, weights=None):
    """
    the sums for jackknifing.  If regions are sent, use them for jackknifing,
    otherwise jackknife one object at a time

    parameters
    ----------
    data: array
        An array with fields 'dsum' and 'wsum'. If shear style
        is lensfit, dsensum is needed rather than wsum.
    jackreg_col: string, optional
        column name holding the jackknife region ids
    weights: array, optional
        Additional weights
    """
    from esutil.stat import histogram

    if weights is not None:
        return get_jackknife_sums_weighted(data, weights,
                                           jackreg_col=jackreg_col)

    shear_style=get_shear_style(data)

    dcol='dsum'
    if shear_style=='lensfit':
        wcol = 'wsum'
    else:
        wcol='dsensum'

    if jackreg_col is None:
        jdsum = data[dcol]
        jwsum = data[wcol]
    else:
        print("using jackreg_col:",jackreg_col)
        regions=data[jackreg_col]

        h,rev=histogram(regions, rev=True)

        nbin=h.size
        nrad=data[dcol].shape[1]
        jdsum=zeros( (nbin, nrad) )
        jwsum=zeros( (nbin, nrad) )

        for i in xrange(nbin):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                jdsum[i] = data[dcol][w].sum(axis=0)
                jwsum[i] = data[wcol][w].sum(axis=0)

        w,=where(h > 0)
        jdsum=jdsum[w,:]
        jwsum=jwsum[w,:]

    return jdsum, jwsum

def get_jackknife_sums_weighted(data, weights, jackreg_col=None):
    """
    the sums for jackknifing.  If regions are sent, use them for jackknifing,
    otherwise jackknife one object at a time

    parameters
    ----------
    data: array
        An array with fields 'dsum' and 'wsum'. If shear style
        is lensfit, dsensum is needed rather than wsum.
    weights: array
        Additional weights
    jackreg_col: string, optional
        column name holding the jackknife region ids
    """
    from esutil.stat import histogram


    shear_style=get_shear_style(data)

    dcol='dsum'
    if shear_style=='lensfit':
        wcol = 'wsum'
    else:
        wcol='dsensum'

    if jackreg_col is None:
        # broadcast it
        wa=weights[:,newaxis]
        jdsum = data[dcol]*wa
        jwsum = data[wcol]*wa
    else:
        print("using jackreg_col:",jackreg_col)
        regions=data[jackreg_col]

        h,rev=histogram(regions, rev=True)

        nbin=h.size
        nrad=data[dcol].shape[1]
        jdsum=zeros( (nbin, nrad) )
        jwsum=zeros( (nbin, nrad) )

        for i in xrange(nbin):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                # broadcast it
                wa=weights[w]
                wa=wa[:,newaxis]

                # note leaving off trailing axis in subscripts
                jdsum[i] = (data[dcol][w]*wa).sum(axis=0)
                jwsum[i] = (data[wcol][w]*wa).sum(axis=0)

        w,=where(h > 0)
        jdsum=jdsum[w,:]
        jwsum=jwsum[w,:]

    return jdsum, jwsum




def average_ratio(l1, l2):
    """
    jackknife (singly) the ratio of dsig1/dsig2
    """

    comb1=average_lensums(l1)
    comb2=average_lensums(l2)

    nlens, nrad = l1['wsum'].shape

    dsum1=l1['dsum']
    dsum_tot1=comb1['dsum']
    dsum2=l2['dsum']
    dsum_tot2=comb2['dsum']

    wsum1, wsum_tot1=_get_ratio_wsums(l1, comb1)
    wsum2, wsum_tot2=_get_ratio_wsums(l2, comb2)

    r      = comb1['r'][0,:]
    ratio  = comb1['dsig'][0,:]/comb2['dsig'][0,:]
    jmean1 = zeros(nrad)
    jmean2 = zeros(nrad)
    jratio = zeros(nrad)
    jdiff  = zeros(nrad)
    jsum   = zeros((nrad, nrad))

    for i in xrange(nlens):

        jmean1[:] = (dsum_tot1[:]-dsum1[i,:])/(wsum_tot1[:]-wsum1[i,:])
        jmean2[:] = (dsum_tot2[:]-dsum2[i,:])/(wsum_tot2[:]-wsum2[i,:])

        jratio[:] = jmean1/jmean2

        jdiff[:] = jratio[:] - ratio[:]

        # now grab all the cross terms and add to the sum for
        # the covariance matrix
        for ix in xrange(nrad):
            for iy in xrange(ix,nrad):
                val = jdiff[ix]*jdiff[iy]

                jsum[ix,iy] += val
                if ix != iy:
                    jsum[iy,ix] += val

    covar = float(nlens-1)/nlens*jsum

    return r, ratio, covar



def _get_ratio_wsums(data, comb):
    """
    helper routine for average_ratio
    get the appropriate denominator for dsum.  Will be dsensum for lensfit,
    which is sum(weight*sensitivity), otherwise wsum which is sum(weight)
    """
    if 'dsensum' in data.dtype.names:
        wsum=data['dsensum']
        wsum_tot=comb['dsensum']
    else:
        wsum=data['wsum']
        wsum_tot=comb['wsum']

    return wsum, wsum_tot


def average_ratio_allrad(l1, l2):
    """
    jackknife summed over all radii
    """

    nlens,nrad = l1['dsum'].shape

    dsum1     = l1['dsum'][:,0:nrad-1].sum(axis=1)
    dsum_tot1 = l1['dsum'][:,0:nrad-1].sum()
    dsum2     = l2['dsum'][:,0:nrad-1].sum(axis=1)
    dsum_tot2 = l2['dsum'][:,0:nrad-1].sum()

    wsum1, wsum_tot1 = _get_ratio_wsums_allrad(l1)
    wsum2, wsum_tot2 = _get_ratio_wsums_allrad(l2)

    dsig1 = dsum_tot1/wsum_tot1
    dsig2 = dsum_tot2/wsum_tot2

    ratio = dsig1/dsig2
    jsum  = 0.0

    for i in xrange(nlens):

        jmean1 = (dsum_tot1-dsum1[i])/(wsum_tot1-wsum1[i])
        jmean2 = (dsum_tot2-dsum2[i])/(wsum_tot2-wsum2[i])

        jratio = jmean1/jmean2

        jdiff = jratio - ratio

        jsum += jdiff*jdiff

    var = float(nlens-1)/nlens*jsum
    err = sqrt(var)

    return ratio, err

def _get_ratio_wsums_allrad(data):
    """
    helper routine for average_ratio
    get the appropriate denominator for dsum.
    """
    nlens,nrad = data['dsum'].shape
    if 'dsensum' in data.dtype.names:
        wsum=data['dsensum'][:,0:nrad-1].sum(axis=1)
        wsum_tot=data['dsensum'][:,0:nrad-1].sum()
    else:
        wsum=data['wsum'][:,0:nrad-1].sum(axis=1)
        wsum_tot=data['wsum'][:,0:nrad-1].sum()

    return wsum, wsum_tot





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

def averaged_struct(nrad, n=1, shear_style=_DEFAULT_SHEAR_STYLE):
    """
    struct to hold the lens averages
    """
    dt = averaged_dtype(nrad, shear_style=shear_style)
    return numpy.zeros(n, dtype=dt)

def averaged_dtype(nrad, shear_style=_DEFAULT_SHEAR_STYLE):
    dt=[('weightsum','f8'),       # this is total of weight for each lens
        ('totpairs','i8'),
        ('r','f8',nrad),
        ('dsig','f8',nrad),
        ('dsigerr','f8',nrad),
        ('dsigcov','f8',(nrad,nrad)),
        ('dsigcor','f8',(nrad,nrad)),
        ('osig','f8',nrad),
        ('wsum_mean','f8',nrad),
        ('wsum_mean_err','f8',nrad),
        ('npair','i8',nrad),
        ('rsum','f8',nrad),
        ('wsum','f8',nrad),
        ('dsum','f8',nrad),
        ('osum','f8',nrad)]

    if shear_style=='lensfit':
        dt+=[('dsensum','f8',nrad),
             ('osensum','f8',nrad)]
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

def average_lensums_slow(lout, weights=None):
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
    nrad = lout['rsum'][0].size

    if 'dsensum' in lout.dtype.names:
        shear_style='lensfit'
    else:
        shear_style='reduced'

    comb = averaged_struct(nrad, shear_style=shear_style)

    # weight is the weight for the lens.  Call this weightsum to
    # indicate a sum over multiple lenses
    comb['weightsum'] = lout['weight'].sum(axis=0)

    comb['totpairs'] = lout['totpairs'].sum(axis=0)


    for i in xrange(nrad):
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

        comb['wsum_mean'][0,i] = wsum/nlens

    if shear_style=='lensfit':
        m,cov=jackknife.wjackknife(vsum=lout['dsum'], wsum=lout['dsensum'])
    else:
        m,cov=jackknife.wjackknife(vsum=lout['dsum'], wsum=lout['wsum'])

    comb['dsigcov'][0,:,:] = cov
    comb['dsigcor'][0,:,:] = jackknife.covar2corr(cov)
    comb['dsigerr'][0,:] = sqrt(diag(cov))

    w=numpy.ones(lout['wsum'].shape)
    m,cov = jackknife.wjackknife(vsum=lout['wsum'], wsum=w)
    comb['wsum_mean_err'][0,:] = sqrt(diag(cov))

    return comb



def average_lensums_weighted_slow(lout, weights):
    """

    Reduce the lens-by-lens lensums by summing over
    all the individual sums and producing averages

    """
    import jackknife

    nlens = lout.size
    nrad = lout['rsum'][0].size

    if weights.size != nlens:
        raise ValueError("weights not same size as lensout, "
                         "%d instead of %d" % (weights.size,nlens))

    totweights = weights.sum()

    if 'dsensum' in lout.dtype.names:
        shear_style='lensfit'
    else:
        shear_style='reduced'

    comb = averaged_struct(nrad, shear_style=shear_style)

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
    for i in xrange(nrad):

        npair = lout['npair'][:,i].sum()

        w_rsum = (lout['rsum'][:,i]*weights).sum()

        w_wsum = (lout['wsum'][:,i]*weights).sum()
        w_dsum = (lout['dsum'][:,i]*weights).sum()
        w_osum = (lout['osum'][:,i]*weights).sum()


        comb['npair'][0,i] = npair
        comb['rsum'][0,i]  = w_rsum
        comb['wsum'][0,i]  = w_wsum
        comb['dsum'][0,i]  = w_dsum
        comb['osum'][0,i]  = w_osum

        # averages
        comb['r'][0,i] = w_rsum/w_wsum

        jdsum[:,i] *= weights
        if shear_style=='lensfit':
            comb['dsensum'][0,i] = lout['dsensum'][:,i].sum()
            comb['osensum'][0,i] = lout['osensum'][:,i].sum()
            comb['dsig'][0,i] = w_dsum/comb['dsensum'][0,i]
            comb['osig'][0,i] = w_osum/comb['osensum'][0,i]
            jwsum[:,i] = comb['dsensum']*weights
        else:

            comb['dsig'][0,i] = w_dsum/w_wsum
            comb['osig'][0,i] = w_osum/w_wsum
            jwsum[:,i] = comb['wsum']*weights


        # this is average wsum over lenses
        # we calculate clustering correction from this, wsum_mean/wsum_mean_random
        comb['wsum_mean'][0,i] = w_wsum/totweights

    # jwsum will be wsum*weights or dsensum*weights
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



