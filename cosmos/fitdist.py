"""
Fit distributions from the cosmos field
"""

import numpy
from . import files



def fit_gprior_exp_mcmc(a=0.25, g0=0.1, gmax=0.87, gmax_min=None, Awidth=1.0,
                        binsize=0.02, doplot=False):
    """
    This works much better than an lm fitter

    for all cosmos galaxies I get
        [840.0, 1.05, 0.087, 0.810]
    """
    import mcmc
    import emcee
    import esutil as eu
    from esutil.random import srandu

    fits_cat=files.read_fits_cat()

    # b/a
    r = fits_cat['sersicfit'][:, 3]
    g = (1-r)/(1+r)


    bs=eu.stat.Binner(g)
    bs.dohist(binsize=binsize)
    bs.calc_stats()
    xdata=bs['center']
    ydata=bs['hist']

    nwalkers=200
    burnin=100
    nstep=100

    print 'fitting exp'

    A=ydata.sum()*(xdata[1]-xdata[0])

    pcen=[A,a,g0,gmax]
    npars=4
    guess=numpy.zeros( (nwalkers,npars) )
    guess[:,0] = pcen[0]*(1.+0.1*srandu(nwalkers))
    guess[:,1] = pcen[1]*(1.+0.1*srandu(nwalkers))
    guess[:,2] = pcen[2]*(1.+0.1*srandu(nwalkers))
    guess[:,3] = pcen[3]*(1.+0.1*srandu(nwalkers))

    ivar = numpy.ones(xdata.size)
    w,=numpy.where(ydata > 0)
    ivar[w] = 1./ydata[w]
    gfitter=GPriorExpFitter(xdata, ydata, ivar, Aprior=A, Awidth=Awidth, gmax_min=gmax_min)

    print 'pcen:',pcen

    sampler = emcee.EnsembleSampler(nwalkers, 
                                    npars,
                                    gfitter.get_lnprob,
                                    a=2)

    pos, prob, state = sampler.run_mcmc(guess, burnin)
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, nstep)

    trials  = sampler.flatchain

    pars,pcov=mcmc.extract_stats(trials)

    d=numpy.diag(pcov)
    perr = numpy.sqrt(d)

    res={'A':pars[0],
         'A_err':perr[0],
         'a':pars[1],
         'a_err':perr[1],
         'g0':pars[2],
         'g0_err':perr[2],
         'gmax': pars[3],
         'gmax_err':perr[3],
         'pars':pars,
         'pcov':pcov,
         'perr':perr}


    fmt="""
A:    %(A).6g +/- %(A_err).6g
a:    %(a).6g +/- %(a_err).6g
g0:   %(g0).6g +/- %(g0_err).6g
gmax: %(gmax).6g +/- %(gmax_err).6g
    """.strip()

    print fmt % res

    if doplot:
        import ngmix
        p=ngmix.priors.GPriorExp(pars)
        gsamp=p.sample1d(g.size)
        plt=eu.plotting.bhist(g, binsize=binsize, show=False)
        eu.plotting.bhist(gsamp, binsize=binsize, plt=plt, color='blue')

    return res

class GPriorExpFitter:
    def __init__(self, xvals, yvals, ivar, Aprior=None, Awidth=None, gmax_min=None):
        """
        Fit with gmax free
        Input is the histogram data
        """
        self.xvals=xvals
        self.yvals=yvals
        self.ivar=ivar

        self.Aprior=Aprior
        self.Awidth=Awidth

        self.gmax_min=gmax_min


    def get_lnprob(self, pars):
        w,=numpy.where(pars < 0)
        if w.size > 0:
            return -9.999e20
        if pars[3] > 1:
            return -9.999e20

        if self.gmax_min is not None:
            if pars[-1] < self.gmax_min:
                return -9.999e20

        model=gprior1d_exp_vec(pars, self.xvals)

        chi2 = (model - self.yvals)**2
        chi2 *= self.ivar

        lnprob = -0.5*chi2.sum()

        if self.Aprior is not None and self.Awidth is not None:
            aprior=-0.5*( (self.Aprior-pars[0])/self.Awidth )**2
            lnprob += aprior

        return lnprob

def gprior1d_exp_vec(pars, g):
    from numpy import pi
    return 2*pi*g*gprior2d_exp_vec(pars, g)
def gprior2d_exp_vec(pars, g):
    from numpy import exp, sqrt
    A=pars[0]
    a=pars[1]
    g0=pars[2]
    gmax=pars[3]

    prior=numpy.zeros(g.size)

    w,=numpy.where(g < gmax)
    if w.size > 0:
        numer = A*(1-exp( (g-gmax)/a ))
        denom = (1+g)*sqrt(g**2 + g0**2)

        prior[w]=numer/denom

    return prior



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


