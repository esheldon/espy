from __future__ import print_function, division
import os
import numpy

def fit_prior(type='great-des', show=False, gmax_data=0.975, gmax=1.0):
    import esutil as eu
    import ngmix

    g=cache_data(gmax_data)
    hd=eu.stat.histogram(g, nbin=100, min=0, max=gmax,more=True)

    xdata=hd['center']
    ydata=hd['hist'].astype('f8')

    if type=='no-exp':
        p=ngmix.priors.GPriorGreatDESNoExp(gmax=gmax)
    elif type=='great-des':
        p=ngmix.priors.GPriorGreatDES(gmax=gmax)
    elif type=='great-des2':
        p=ngmix.priors.GPriorGreatDES2()
    elif type=='ba':
        p=ngmix.priors.GPriorBA(0.2)
    elif type=='merf':
        p=ngmix.priors.GPriorMErf()
    else:
        raise ValueError("bad type '%s'" % type)

    p.dofit(hd['center'], hd['hist'].astype('f8'), show=show)

def cache_data(gmax_data):
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

    w,=numpy.where(g < gmax_data)
    g=g[w]
    return g
