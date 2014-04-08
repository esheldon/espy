from __future__ import print_function

import os
import sys
import numpy
from numpy import linspace

from . import fitting
from .fitting import get_plot_fname


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

def fit_spline(pars, nvals, type, order=3):
    import biggles
    import pcolors

    from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

    if order != 3 and order != 1:
        raise ValueError("order 1 or 3")

    ngauss=pars.shape[1]
    colors=pcolors.rainbow(ngauss)

    plt=biggles.FramedPlot()
    plt.xlabel='Sersic n'
    plt.ylabel=type

    nvals_interp=linspace(nvals[0], nvals[-1],10000)


    for i in xrange(ngauss):

        color=colors[i]

        vals=pars[:,i]

        #interpolator=InterpolatedUnivariateSpline(nvals,vals,k=order)
        #vals_interp=interpolator(nvals_interp)
        if order==3:
            interpolator=InterpolatedUnivariateSpline(nvals,vals,k=order)
            vals_interp=interpolator(nvals_interp)
        else:
            #vals_interp=numpy.lib.function_base.compiled_interp(nvals_interp,
            #                                                    nvals,
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

    dir=fitting.get_dir()

    eps=get_plot_fname(ngauss, order, type)
    print(eps)
    plt.write_eps(eps)


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



