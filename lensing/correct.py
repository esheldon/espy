from __future__ import print_function

import numpy
from numpy import where, sqrt, zeros, ones
import esutil as eu
from esutil.numpy_util import where1

import lensing


def correct(data, rand, subtract_rand=False, minrad=None):
    """

    Apply various corrections.

    1) if requested, subtract random points within the specified
       radius.
    2) 1/ssh
    3) clustering correction

    parameters
    ----------
    data: ndarray
        This should be the result of binning the data
    rand: ndarray
        This should be the result of matching to the binned data

    """

    nrad = data['r'][0].size
    nbin = data.size
    if rand.size != data.size:
        raise ValueError("random size (%d) not same size as "
                         "data (%d)" % (rand.size, data.size))

    if subtract_rand:
        if minrad is None:
            raise ValueError("send minrad if you want subtraction of randoms")
        dsig_rand, dsigerr_rand = calc_rand_to_subtract(data, rand, minrad)
    else:
        dsig_rand= zeros((nbin,nrad))
        dsigerr_rand = zeros((nbin,nrad))


    ssh = calc_ssh(data)
    clust_corr, clust_corr_err = calc_clustering_correction(data, rand)

    out = add_corrections(data, ssh, 
                          clust_corr, clust_corr_err, 
                          rand['wsum_mean'], rand['wsum_mean_err'],
                          dsig_rand, dsigerr_rand)

    return out

def add_corrections(datain, 
                    ssh, 
                    clust_corr, clust_corr_err, 
                    wsum_mean_rand, wsum_mean_rand_err,
                    dsig_rand, dsigerr_rand):

    nrad = datain['r'][0].size
    data = eu.numpy_util.add_fields(datain,
                                    [('dsig_orig','f8',nrad),
                                     ('dsigerr_orig','f8',nrad),
                                     ('dsigcov_orig','f8',(nrad,nrad)),
                                     ('wsum_mean_rand','f8',nrad),
                                     ('wsum_mean_rand_err','f8',nrad),
                                     ('dsig_rand','f8',nrad),
                                     ('dsigerr_rand','f8',nrad),
                                     ('clust_corr','f8',nrad),
                                     ('clust_corr_err','f8',nrad)])

    data['dsig_orig'] = datain['dsig']
    data['dsigerr_orig'] = datain['dsigerr']
    data['dsigcov_orig'] = datain['dsigcov']

    # note copying over ssh that is there already
    data['ssh']            = ssh

    data['wsum_mean_rand'] = wsum_mean_rand
    data['wsum_mean_rand_err'] = wsum_mean_rand_err
    data['dsig_rand']      = dsig_rand
    data['dsigerr_rand']   = dsigerr_rand
    data['clust_corr']     = clust_corr
    data['clust_corr_err'] = clust_corr_err

    # The subtraction should always happen first
    err = sqrt(data['dsigerr']**2 + data['dsigerr_rand']**2)
    data['dsig'] -= dsig_rand
    data['dsigerr'] = err

    # these are multiplications, can be in any order
    data['dsig'] /= ssh
    data['dsigerr'] /= ssh

    corr = data['clust_corr']
    corr_clipped = corr.clip(1.0, corr.max())
    data['dsig'] *= corr_clipped
    data['dsigerr'] *= corr_clipped

    # for now, just rescale the correlation matrix.
    for i in xrange(data.size):
        data['dsigcov'][i] = eu.stat.cor2cov(data['dsigcor'][i],data['dsigerr'][i])

    return data


def calc_rand_to_subtract(datain, rand, minrad):
    nbin = datain.size
    if rand.size != nbin:
        raise ValueError("random size (%d) not same size as "
                         "data (%d)" % (rand.size, data.size))




    w,=where(data['r'] > minrad)
    wr,=where(rand['r'] > minrad)
    if w.size != wr.size:
        raise ValueError("Found > minrad %d from data but %d "
                         "from rand" % (w.size,wr.size))
    if w.size == 0:
        raise ValueError("no radii greater than minrad %0.2f" % minrad)


    nrad = datain['r'][0,:].size
    rand = zeros((nbin, nrad))
    rand_err = zeros((nbin, nrad))

    for binnum in xrange(nbin):
        rand[binnum,w] = rand['dsig'][binnum,w]
        rand_err[binnum,w] = rand['dsigerr'][binnum,w]

    return rand, rand_err


def calc_ssh(data):
    """

    Correct for the "shear responsivity"

    """

    # just use the mean ssh, should be more stable when
    # we only have a few lenses
    ssh = data['sshsum'].sum()/data['weightsum'].sum()

    return ssh


def calc_clustering_correction(data, rand):
    """

    Calculate the clustering correction, add corresponding fields to the
    structure, and return

    """

    corr = data['wsum_mean']/rand['wsum_mean']

    if 'wsum_mean_err' in data.dtype.names:
        derr = data['wsum_mean_err']
    else:
        derr = data['wsum_err']

    if 'wsum_mean_err' in rand.dtype.names:
        rerr = rand['wsum_mean_err']
    else:
        rerr = rand['wsum_err']

    corr_err = zeros(data['wsum_mean'].shape) + 9999
    w=where( (data['wsum_mean'] > 0) & (rand['wsum_mean'] > 0) )
    corr_err[w] = corr[w]*sqrt(   (derr[w]/data['wsum_mean'][w])**2
                                + (rerr[w]/rand['wsum_mean'][w])**2 )

    return corr, corr_err

