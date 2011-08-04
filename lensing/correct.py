from __future__ import print_function

import numpy
from numpy import sqrt
import esutil as eu
from esutil.numpy_util import where1

import lensing

def correct_ssh(data):
    """

    Correct for the "shear responsivity"

    """
    out=data.copy()

    # just use the mean ssh, should be more stable when
    # we only have a few lenses
    ssh = out['sshsum'].sum()/out['weightsum'].sum()

    out['dsig'] /= ssh

    return out

def correct_clustered(data, rand):
    """
    Correct for clustering of sources with lenses.
    """
    nrad = data['r'][0].size

    add_dt=[('wsum_mean_rand','f8',nrad),
            ('wsum_mean_rand_err','f8',nrad),
            ('clust_corr','f8',nrad),
            ('clust_corr_err','f8',nrad)]

    out = eu.numpy_util.add_fields(data, add_dt)
    out['wsum_mean_rand'] = rand['wsum_mean']
    out['wsum_mean_rand_err'] = rand['wsum_err']


    corr = data['wsum_mean']/rand['wsum_mean']

    corr_err = corr*sqrt(   (data['wsum_err']/data['wsum_mean'])**2
                          + (rand['wsum_err']/rand['wsum_mean'])**2 )

    out['clust_corr'] = corr
    out['clust_corr_err'] = corr_err

    # we only apply the clipped correction
    corr_clipped = corr.clip(1.0, corr.max())
    out['dsig'] *= corr_clipped

    return out


