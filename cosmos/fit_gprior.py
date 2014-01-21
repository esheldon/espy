import numpy
from . import files
from . import analysis

def get_shapes(cat_type, version=None):
    if cat_type=="galfit":
        fits_cat=files.read_fits_cat()
        # b/a
        r = fits_cat['sersicfit'][:, 3]
        g = (1-r)/(1+r)

    elif cat_type in ["ngmix-exp","ngmix-dev","ngmix-bdf"]:
        data=files.read_output(version)

        if cat_type=="ngmix-exp":
            model="exp"
        elif cat_type=="ngmix-dev":
            model="dev"
        elif cat_type=="ngmix-bdf":
            model="bdf"

        w=analysis.select_by_s2n_flux(data, model)

        pname="%s_pars" % model
        g1 = data[pname][w,2]
        g2 = data[pname][w,3]

        g=numpy.sqrt(g1**2 + g2**2)

    else:
        raise ValueError("bad cat type: '%s'" % (cat_type))

    return g
                            
def fit_gprior_m_style(cat_type, version=None,
                       a=0.25, g0=0.1, gmax=0.87, gmax_min=None, Awidth=1.0,
                       binsize=0.02, doplot=False):
    """
    cat_type should be "galfit" or "ngmix-exp" "ngmix-dev" "ngmix-bdf"

    If cat_type=="galfit" then fit to the shapes from the sersic fits.
    
    If cat=="ngmix-exp" use my fits, same for dev.  Must send version= as well

    This works much better than an lm fitter

    for all cosmos galaxies I get
        [840.0, 1.05, 0.087, 0.810]
    """
    import mcmc
    import emcee
    import esutil as eu
    from esutil.random import srandu

    g=get_shapes(cat_type, version=version)

    bs=eu.stat.Binner(g)
    bs.dohist(binsize=binsize)
    bs.calc_stats()
    xdata=bs['center']
    ydata=bs['hist']

    nwalkers=200
    burnin=500
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
    gfitter=GPriorMFitter(xdata, ydata, ivar, Aprior=A, Awidth=Awidth, gmax_min=gmax_min)

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
        import mcmc
        import ngmix
        mcmc.plot_results(trials,names=['A','a','g0','gmax'],
                          title=cat_type)
        p=ngmix.priors.GPriorM(pars)
        gsamp=p.sample1d(g.size)
        plt=eu.plotting.bhist(g, binsize=binsize, show=False)
        eu.plotting.bhist(gsamp, binsize=binsize,
                          plt=plt, color='blue',
                          xlabel='|g|',
                          xrange=[0.,1.],
                          title=cat_type)

    return res

class GPriorMFitter:
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



