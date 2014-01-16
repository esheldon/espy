"""
Fit distributions from the cosmos field
"""

import numpy

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


def T_rat_from_fix_rhalf_rat(sigma_psf=1.414, rhalf_exp=4.0, show=False):
    """
    Generate models with r_{1/2}^{dev} = 0.6 r_{1/2}^{exp} and
    measure the T ratio

    parameters
    ----------
    rhalf: float, optional
        Half light radius of the exponential component

    dependencies
    ------------
    galsim
    """
    import galsim 
    import ngmix
    from numpy.random import randn
    from esutil.random import srandu
    import pprint

    flux=100.0
    s2n=10000.0

    shape=(32,32)
    size=shape[0]*shape[1]
    pixel_scale=1.0
    j=ngmix.jacobian.UnitJacobian(0.5*shape[0], 0.5*shape[1])

    rhalf_dev = 0.6*rhalf_exp

    pix = galsim.Pixel(pixel_scale)
    psf0=galsim.Gaussian(flux=1.0, sigma=sigma_psf)

    exp0=galsim.Exponential(flux=1.0, half_light_radius=rhalf_exp)
    dev0=galsim.DeVaucouleurs(flux=1.0, half_light_radius=rhalf_dev)
    gal0 = galsim.Add([exp0, dev0])
    gal0.setFlux(flux)

    gal = galsim.Convolve([gal0, psf0, pix])
    psf = galsim.Convolve([psf0, pix])

    galsim_image     = galsim.ImageD(shape[0], shape[1])
    galsim_psf_image = galsim.ImageD(shape[0], shape[1])

    gal.draw(galsim_image, dx=pixel_scale)
    psf.draw(galsim_psf_image, dx=pixel_scale)

    image_nonoise=galsim_image.array
    psf_image=galsim_psf_image.array

    skysig2 = (image_nonoise**2).sum()/s2n**2
    skysig = numpy.sqrt(skysig2)

    image = image_nonoise + skysig*randn(size).reshape(shape)

    if show:
        import images
        images.multiview(psf_image,title='psf')
        images.multiview(image,title='image')

    imsky,sky=ngmix.em.prep_image(psf_image)
    em=ngmix.em.GMixEM(imsky, jacobian=j)
    Tpsf_guess=2*sigma_psf**2
    psf_guess=ngmix.gmix.GMixModel([0.0, 0.0, 0.0, 0.0, Tpsf_guess, 1.0],'gauss')
    em.go(psf_guess, sky)

    psf_gmix = em.get_gmix()

    nwalkers=200
    guess=numpy.zeros( (nwalkers, 8) )
    guess[:,0] = 0.1*srandu(nwalkers)
    guess[:,1] = 0.1*srandu(nwalkers)
    guess[:,2] = 0.1*srandu(nwalkers)
    guess[:,3] = 0.1*srandu(nwalkers)
    guess[:,4] = 40*(1.0 + 0.1*srandu(nwalkers))
    guess[:,5] = 30*(1.0 + 0.1*srandu(nwalkers))
    guess[:,6] = flux*(1.0 + 0.1*srandu(nwalkers))
    guess[:,7] = flux*(1.0 + 0.1*srandu(nwalkers))

    T_prior = ngmix.priors.FlatPrior(0.001, 200)
    counts_prior = ngmix.priors.FlatPrior(0.001, 1000)

    wt=0*image + 1.0/skysig2
    fitter=ngmix.fitting.MCMCBDC(image, wt, j, 'bdc',
                                 psf=psf_gmix,
                                 full_guess=guess,
                                 T_b_prior=T_prior,
                                 T_d_prior=T_prior,
                                 counts_b_prior=counts_prior,
                                 counts_d_prior=counts_prior,
                                 nwalkers=nwalkers,
                                 burnin=800,
                                 nstep=100,
                                 mca_a=3.0,
                                 min_arate=0.3)
    fitter.go()
    if show:
        fitter.make_plots(show=show)
    res=fitter.get_result()

    pars=res['pars']
    perr=res['pars_err']
    Tb = pars[4] 
    Td = pars[5]
    Tb_err = perr[4]
    Td_err = perr[5]

    Trat = Tb/Td
    Trat_err = Trat*( (Tb/Tb_err)**2 + (Td/Td_err)**2 )
    pprint.pprint(res)
    print
    print '%s +/- %s' % (Trat, Trat_err)
