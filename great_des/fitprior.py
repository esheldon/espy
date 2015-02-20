from __future__ import print_function, division
import numpy

def fit_prior(show=False):
    import esutil as eu
    import ngmix
    import fitsio

    fname='/astro/u/esheldon/lensing/great-des/sfit-noisefree-02/collated/sfit-noisefree-02-g00.fits'
    print("reading:",fname)
    t=fitsio.read(fname)

    g2 = numpy.sqrt( t['g'][:,0]**2 + t['g'][:,1]**2 )
    w,=numpy.where(t['flags']==0)

    hd=eu.stat.histogram(g2[w], nbin=100, min=0, max=0.985,more=True)

    xdata=hd['center']
    ydata=hd['hist'].astype('f8')

    #pg3=ngmix.priors.GPriorGreat3Exp()
    pg3=ngmix.priors.GPriorGreatDES(gmax=0.985)

    pg3.dofit(hd['center'], hd['hist'].astype('f8'), show=show)
