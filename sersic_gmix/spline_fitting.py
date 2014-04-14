from __future__ import print_function

import os
import sys
import numpy
from numpy import linspace, zeros

from numba import autojit, jit, float64

from . import fitting
from .fitting import get_plot_fname

#@jit(argtypes=[float64[:], float64])
@autojit
def binary_search(a, x):
    """
    Index of closest value from a smaller than x

    however, defaults to edges when out of bounds
    """

    up=a.size
    down=-1

    if x < a[0]:
        return 0
    if x > a[up-1]:
        return up-1

    while up-down > 1:
        mid = down + (up-down)//2
        val=a[mid]

        if x >= val:
            down=mid
        else:
            up=mid
        
    return down

#@autojit
@jit(argtypes=[float64[:], float64[:,:], float64, float64[:]])
def interp_multi_scalar(xref, yref, x, output):
    """
    parameters
    ----------
    xref: array
        shape (n,)
    yref: array
        shape (n,ndim)
    x: scalar
        point at which to interpolate
    output: array
        shape (ndim,)
    """

    np=xref.size
    ndim=output.size

    ilo = binary_search(xref, x)
    if (ilo >= (np-1)):
        ilo = np-2
    ihi = ilo + 1

    for i in xrange(ndim):
        output[i] = (x-xref[ilo])*(yref[ihi,i] - yref[ilo,i])/(xref[ihi]-xref[ilo]) + yref[ilo,i]



def interp_multi_array(xref, yref, x):
    """
    parameters
    ----------
    xref: array
        shape (n,)
    yref: array
        shape (n,ndim)
    x: array
        points at which to interpolate
    """

    ndim=yref.shape[1]
    npoints=x.size
    output=zeros( (npoints, ndim) )
    res=zeros(ndim)
    for i in xrange(npoints):
        interp_multi_scalar(xref, yref, x[i], res)
        output[i,:] = res

    return output

def convert_hogg(pars):
    """
    Hogg's pars are [f1,f2,f2,f3,...,T1,T2,T3,T4,...]

    Also normalize

    They are already sorted
    """

    ngauss=pars.size/2
    Fvals = pars[0:ngauss].copy()
    Tvals = pars[ngauss:].copy()

    # Fvals now sum to unity
    Fvals /= Fvals.sum()

    Tmean = (Tvals*Fvals).sum()
    Tvals /= Tmean

    newpars=pars.copy()
    newpars[0:ngauss] = Tvals
    newpars[ngauss:] = Fvals

    return newpars

def fit_spline(nvals, pars):
    import biggles
    import pcolors

    ngauss=pars.shape[1]/2
    colors=pcolors.rainbow(ngauss)

    ninterp=10000
    nvals_interp=linspace(nvals[0]-0.001, nvals[-1]+0.001, ninterp)

    vals_interp_all = interp_multi_array(nvals, pars, nvals_interp)

    for type in ['T','flux']:

        if type=='T':
            start=0
        else:
            start=ngauss

        plt=biggles.FramedPlot()
        plt.xlabel='Sersic n'
        plt.ylabel=type

        for i in xrange(ngauss):

            color=colors[i]
            vals=pars[:,start+i]
            vals_interp=vals_interp_all[:,start+i]

            pts=biggles.Points(nvals, vals,
                               type='filled circle',color=color,size=1.0)
            crv=biggles.Curve(nvals_interp, vals_interp, type='solid',color=color)

            plt.add( pts, crv )

        plt.ylog=True

        minval=pars.min()
        maxval=pars.max()
        plt.xrange = [0.2,6.5]
        plt.yrange = [0.5*minval, 1.5*maxval]

        dir=fitting.get_dir()

        eps=get_plot_fname(ngauss, type)
        print(eps)
        plt.write_eps(eps)


def fit_spline_old(pars, nvals, type, order=3):
    import biggles
    import pcolors

    from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

    if order != 3 and order != 1:
        raise ValueError("order 1 or 3")

    ngauss=pars.shape[1]/2
    colors=pcolors.rainbow(ngauss)

    if type=='T':
        start=0
    else:
        start=ngauss

    plt=biggles.FramedPlot()
    plt.xlabel='Sersic n'
    plt.ylabel=type

    nvals_interp=linspace(nvals[0], nvals[-1],10000)


    for i in xrange(ngauss):

        color=colors[i]

        vals=pars[:,start+i]

        if order==3:
            interpolator=InterpolatedUnivariateSpline(nvals,vals,k=order)
            vals_interp=interpolator(nvals_interp)
        else:
            #                                                    vals)
            interpolator=interp1d(nvals,vals)
            vals_interp=interpolator(nvals_interp)

        pts=biggles.Points(nvals, vals, type='filled circle',color=color,size=1.0)
        crv=biggles.Curve(nvals_interp, vals_interp, type='solid',color=color)

        plt.add( pts, crv )

    plt.ylog=True

    minval=pars.min()
    maxval=pars.max()
    plt.xrange = [0.2,6.5]
    plt.yrange = [0.5*minval, 1.5*maxval]

    #dir=fitting.get_dir()

    #eps=get_plot_fname(ngauss, order, type)
    #print(eps)
    #plt.write_eps(eps)
    plt.show()


def sort_by_T(data):
    ndata=data.size
    ngauss=(data['pars'].shape[1]-4)/2

    for name in ['pars','pars_norm']:
        for i in xrange(ndata):
            Tvals=data[name][i,4:4+ngauss].copy()
            fvals=data[name][i,4+ngauss:].copy()

            s=Tvals.argsort()

            Tvals=Tvals[s]
            fvals=fvals[s]

            data[name][i,4:4+ngauss] = Tvals
            data[name][i,4+ngauss:] = fvals

    return data



