from __future__ import print_function, division
import os
import numpy

def fit_prior(show=False, gmax=0.975):
    import esutil as eu
    import ngmix

    g=cache_data(gmax)
    hd=eu.stat.histogram(g, nbin=100, min=0, max=gmax,more=True)

    xdata=hd['center']
    ydata=hd['hist'].astype('f8')

    #pg3=ngmix.priors.GPriorGreat3Exp()
    #pg3=ngmix.priors.GPriorGreatDES(gmax=gmax)
    pg3=ngmix.priors.GPriorGreatDES2(gmax=gmax)

    pg3.dofit(hd['center'], hd['hist'].astype('f8'), show=show)

def cache_data(gmax):
    import fitsio
    tmpdir=os.environ['TMPDIR']
    cachename=os.path.join(tmpdir,'sfit-noisefree-02-g00-gvals.fits')

    if os.path.exists(cachename):
        print("reading:",cachename)
        data=fitsio.read(cachename)
        g=data['g']
    else:
        fname='/astro/u/esheldon/lensing/great-des/sfit-noisefree-02/collated/sfit-noisefree-02-g00.fits'
        print("reading:",fname)
        t=fitsio.read(fname)

        g = numpy.sqrt( t['g'][:,0]**2 + t['g'][:,1]**2 )
        w,=numpy.where(t['flags']==0)

        g=g[w]

        out=numpy.zeros(g.size, dtype=[('g','f4')])
        out['g'] = g
        print("writing:",cachename)
        fitsio.write(cachename, out, clobber=True)

    w,=numpy.where(g < 0.99)
    g=g[w]
    return g
