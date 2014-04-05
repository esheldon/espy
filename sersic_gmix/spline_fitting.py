from __future__ import print_function
import sys
import numpy
from numpy import linspace

def renorm(data):
    n=data.size

    ngauss=(data['pars'].shape[1]-4)/2

    for i in xrange(n):
        Tvals=data['pars'][i,4:4+ngauss].copy()
        fvals=data['pars'][i,4+ngauss:].copy()

        # fvals now sum to unity
        fvals /= fvals.sum()

        Tmean = (Tvals*fvals).sum()
        Tvals /= Tmean

        data['pars'][i,4:4+ngauss] = Tvals
        data['pars'][i,4+ngauss:] = fvals

    return data

def fit_spline(data, nvals, type):
    import biggles
    import pcolors
    from scipy.interpolate import InterpolatedUnivariateSpline

    ngauss=(data['pars'].shape[1]-4)/2
    colors=pcolors.rainbow(ngauss)

    plt=biggles.FramedPlot()
    plt.xlabel='Sersic n'
    plt.ylabel=type

    if type=='T':
        start=4
    else:
        start=4+ngauss

    nvals_interp=linspace(nvals[0], nvals[-1])
    for i in xrange(ngauss):

        color=colors[i]

        vals=data['pars'][:,start+i]
        pts=biggles.Points(nvals, vals, type='filled circle',color=color)

        interpolator=InterpolatedUnivariateSpline(nvals,vals,k=3)
        vals_interp=interpolator(nvals_interp)
        #crv=biggles.Curve(nvals, vals, type='solid',color=color)
        crv=biggles.Curve(nvals_interp, vals_interp, type='solid',color=color)
        plt.add( pts, crv )

    plt.ylog=True

    allvals = data['pars'][:,start:start+ngauss]
    minval=allvals.min()
    maxval=allvals.max()
    plt.xrange = [0.2,6.0]
    plt.yrange = [0.5*minval, 1.5*maxval]
    #plt.show()

    eps='%s-vs-n.eps' % type
    print(eps)
    plt.write_eps(eps)

