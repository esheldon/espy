"""

jackknife the lensums.

"""

from __future__ import print_function
import os
import sys
import tempfile

import numpy
from numpy import log10,sqrt,linspace,where,zeros

import lensing
import esutil as eu
from esutil.numpy_util import where1
import recfile

def wjackknife(vsum=None, wsum=None, cleanup=True, slow=False):
    """
    jackknife the weighted data

    parameters
    ----------
    vsum: array
        An (nsample, nvar) array.  Each element should be interpreted as
        weight*var or sum(weights*var) such that the means will be
            vsum/wsum
    wsum: array
        An (nsample, nvar) array.  Each element should be interpreted as
        weight or sum(weights)

    cleanup: bool, optional
    """

    if vsum is None or wsum is None:
        raise ValueError("send vsum and wsum")

    if slow:
        return jackknife_purepy(vsum=vsum, wsum=wsum)

    f=tempfile.mktemp(prefix='jackknife-input-', suffix='.bin')
    fout = f.replace('.bin','-out.dat')
    try:

        with open(f,'w') as fobj:
            nsample, nvar = wsum.shape

            nsample=numpy.array([nsample], dtype='i8')
            nvar=numpy.array([nvar], dtype='i8')
            vsum=numpy.array(vsum, ndmin=1, dtype='f8', copy=False)
            wsum=numpy.array(wsum, ndmin=1, dtype='f8', copy=False)

            nsample.tofile(fobj)
            nvar.tofile(fobj)
            vsum.tofile(fobj)
            wsum.tofile(fobj)
        
        command='jackknife '+f+' '+fout
        if os.system(command) != 0:
            raise RuntimeError("command '%s' failed" % command)

        nvar=int(nvar[0])

        dtype=[('mean','f8'),('err','f8')]
        fobj = recfile.Open(fout,dtype=dtype,delim=' ',nrows=nvar)

        d=fobj.read()
        mean = d['mean'].copy()

        covar=numpy.fromfile(fobj.fobj, sep=' ', dtype='f8', count=nvar*nvar)
        covar=covar.reshape((nvar,nvar))

    finally:
        if cleanup:
            if os.path.exists(f):
                os.remove(f)
            if os.path.exists(fout):
                os.remove(fout)

    return mean, covar

def jackknife_purepy(vsum=None, wsum=None):
    """

    Simple jackknife removing one at a time instead of spacially

    """

    nsample, nvar = wsum.shape

    # sum over samples
    wsum_tot = wsum.sum(axis=0)
    vsum_tot = vsum.sum(axis=0)

    # mean in the radial bins
    mean = vsum_tot/wsum_tot

    jmean  = zeros(nvar)
    jdiff = zeros(nvar)
    jsum = zeros((nvar, nvar))

    for i in xrange(nsample):
        # calculate mean minus the ith sample

        jmean[:] = (vsum_tot[:]-vsum[i,:])/(wsum_tot[:]-wsum[i,:])

        # diff in each radial bin
        jdiff[:] = jmean[:]-mean[:]

        # now grab all the cross terms and add to the sum for
        # the covariance matrix
        for ix in xrange(nvar):
            for iy in xrange(ix,nvar):
                val = jdiff[ix]*jdiff[iy]

                jsum[ix,iy] += val
                if ix != iy:
                    jsum[iy,ix] += val



    covar = float(nsample-1)/nsample*jsum

    return mean, covar


def covar2corr(cov):
    corr = zeros(cov.shape)

    for ix in xrange(cov.shape[0]):
        for iy in xrange(cov.shape[1]):
            corr[ix,iy] = cov[ix,iy]/sqrt(cov[ix,ix]*cov[iy,iy])

    return corr
