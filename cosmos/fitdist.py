"""
Fit distributions from the cosmos field
"""

import numpy
from . import files


def many_sigma2rhalf(model, sigma_vals, nsig=5, smooth=0.001):
    rhalf_vals=numpy.zeros(sigma_vals.size)

    for i in xrange(sigma_vals.size):
        rhalf_vals[i] = sigma2rhalf(model,sigma_vals[i],nsig=nsig,smooth=smooth)

    return rhalf_vals

def sigma2rhalf(model, sigma, nsig=5, smooth=0.001, doplot=False):
    """
    map sigma ( sqrt(T/2) ) to half life radii for the input model

    for exp we get
        sigma/rhalf = 2.36
        used smooth=0.001
        nsig=1
    for dev
        used smooth=0.001
        nsig=0.01
        sigma/rhalf = 558.

    Note this is for the gaussian mixtures 
    """
    import ngmix

    T = 2*sigma**2

    dim = 2*nsig*sigma
    cen=dim/2.0

    dims=[dim]*2

    print '-'*70
    print 'dims:',dims

    flux=1.0

    pars=[cen,cen, 0.0, 0.0, T, flux]
    gm=ngmix.gmix.GMixModel(pars, model)

    image=gm.make_image(dims)
    image /= image.max()

    fwhm=measure_image_width_erf(image, 0.5,smooth=smooth)
    rhalf=fwhm/2.0

    print 'sigma:',sigma
    print 'rhalf:',rhalf
    print 'T:',T
    # sigma = fac*rhalf

    fac = sigma/rhalf
    print 'sigma/rhalf:',fac

    if doplot:
        import esutil as eu
        vals=image[:, cen]
        ivals=numpy.arange(vals.size)
        plt=eu.plotting.bscatter(ivals, vals, show=False)
        eu.plotting.bscatter(ivals, [0.5]*vals.size, color='blue',
                             type='solid', plt=plt)

    return rhalf

def measure_image_width_erf(image, thresh_vals, smooth=0.1):
    """
    Measure width at the given threshold using an erf to smooth the contour.
    
    e.g. 0.5 would be the FWHM

    parameters
    ----------
    image: 2-d darray
        The image to measure
    thresh_vals: scalar or sequence
        threshold is, e.g. 0.5 to get a Full Width at Half max
    smooth: float
        The smoothing scale for the erf.  This should be between 0 and 1. If
        you have noisy data, you might set this to the noise value or greater,
        scaled by the max value in the images.  Otherwise just make sure it
        smooths enough to avoid pixelization effects.

    output
    ------
    widths: scalar or ndarray
        sqrt(Area)/pi where Area is,

            nim=image.image.max()
            arg =  (nim-thresh)/smooth
            vals = 0.5*( 1 + erf(arg) )
            area = vals.sum()
            width = 2*sqrt(area/pi)
    """
    from numpy import array, sqrt, zeros, pi, where
    from scipy.special import erf

    if isinstance(thresh_vals, (list,tuple,numpy.ndarray)):
        is_seq=True
    else:
        is_seq=False

    thresh_vals=array(thresh_vals,ndmin=1,dtype='f8')

    nim = image.copy()
    maxval=image.max()
    nim *= (1./maxval)

    widths=zeros(len(thresh_vals))
    for i,thresh in enumerate(thresh_vals):
        arg = (nim-thresh)/smooth

        vals = 0.5*( 1+erf(arg) )
        area = vals.sum()
        widths[i] = 2*sqrt(area/pi)

    if is_seq:
        return widths
    else:
        return widths[0]


