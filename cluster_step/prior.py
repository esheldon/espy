import numpy
from numpy import diag, where, cos, sin, exp, sqrt, zeros, random
from math import pi

from esutil.random import srandu

class GPrior(object):
    """
    This is the base class.  You need to over-ride a few of
    the functions, see below
    """
    def __init__(self, pars):
        self.pars=pars
        self.maxval = self(0., 0.)

        # sub-class may want to over-ride this, see GPriorExp
        self.gmax=1.0

    def __call__(self, g1, g2):
        """
        Get the 2d prior
        """
        g = sqrt(g1**2 + g2**2)
        return self.prior2d_gabs(g)

    def prior2d_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        raise RuntimeError("over-ride")

    def prior2d_gabs_scalar(self, g):
        """
        Get the 2d prior for the input |g| scalar value
        """
        raise RuntimeError("over-ride")


    def prior1d(self, g):
        """
        Get the 1d prior for an input |g| value(s).
        """
        return 2*pi*g*self.prior2d_gabs(g)


    def dbyg1(self, g1, g2, h=1.e-6):
        """
        Derivative with respect to g1 at the input g1,g2 location

        Uses central difference and a small enough step size
        to use just two points
        """
        ff = self(g1+h/2, g2)
        fb = self(g1-h/2, g2)

        return (ff - fb)/h

    def dbyg2(self, g1, g2, h=1.e-6):
        """
        Derivative with respect to g2 at the input g1,g2 location

        Uses central difference and a small enough step size
        to use just two points
        """
        ff = self(g1, g2+h/2)
        fb = self(g1, g2-h/2)
        return (ff - fb)/h



    def sample1d(self, nrand):
        """
        Get random |g| from the 1d distribution

        Set self.gmax appropriately

        parameters
        ----------
        nrand: int
            Number to generate
        """

        if not hasattr(self,'maxval1d'):
            self.set_maxval1d()

        g = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g in [0,1)
            grand = self.gmax*random.random(nleft)

            # now the height from [0,maxval)
            h = self.maxval1d*random.random(nleft)

            pvals = self.prior1d(grand)

            w,=where(h < pvals)
            if w.size > 0:
                g[ngood:ngood+w.size] = grand[w]
                ngood += w.size
                nleft -= w.size
   
        return g


    def sample2d(self, nrand):
        """
        Get random g1,g2 values

        parameters
        ----------
        nrand: int
            Number to generate
        """

        grand=self.sample1d(nrand)
        rangle = random.random(nrand)*2*pi
        g1rand = grand*cos(rangle)
        g2rand = grand*sin(rangle)
        return g1rand, g2rand

    def set_maxval1d(self):
        """
        Use a simple minimizer to find the max value of the 1d 
        distribution
        """
        import scipy.optimize
        
        (minvalx, fval, iterations, fcalls, warnflag) \
                = scipy.optimize.fmin(self.prior1dneg, 0.1, full_output=True, 
                                      disp=False)
        if warnflag != 0:
            raise ValueError("failed to find min: warnflag %d" % warnflag)
        self.maxval1d = -fval

    def prior1dneg(self, g, *args):
        """
        So we can use the minimizer
        """
        return -self.prior1d(g)




class GPriorExp(GPrior):
    def __init__(self, pars):
        """
        [A, a, g0, gmax]
        """
        super(GPriorExp,self).__init__(pars)
        self.gmax=pars[-1]

    def prior2d_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        return gprior2d_exp_vec(self.pars, g)

    def prior2d_gabs_scalar(self, g):
        """
        Get the 2d prior for the input |g| scalar value
        """
        return gprior2d_exp_scalar(self.pars, g)


def gprior2d_exp_vec(pars, g):
    A=pars[0]
    a=pars[1]
    g0=pars[2]
    gmax=pars[3]

    prior=zeros(g.size)

    w,=where(g < gmax)
    if w.size > 0:
        numer = A*(1-exp( (g-gmax)/a ))
        denom = (1+g)*sqrt(g**2 + g0**2)

        prior[w]=numer/denom

    return prior

def gprior2d_exp_scalar(pars, g):
    from math import exp, sqrt
    A=pars[0]
    a=pars[1]
    g0=pars[2]
    gmax=pars[3]

    if g > gmax:
        return 0.0

    numer = A*(1-exp( (g-gmax)/a ))
    denom = (1+g)*sqrt(g**2 + g0**2)

    prior=numer/denom

    return prior


def gprior1d_exp_vec(pars, g):
    return 2*pi*g*gprior2d_exp_vec(pars, g)


class GPriorExpFitterFixedGMax:
    def __init__(self, xvals, yvals, gmax=0.87):
        """
        Input is the histogram data
        """
        self.xvals=xvals
        self.yvals=yvals
        self.gmax=gmax

    def __call__(self, pars):
        w,=where(pars < 0)
        if w.size > 0:
            return zeros(self.xvals.size) + numpy.inf

        send_pars=list(pars) + [self.gmax]

        model=gprior1d_exp_vec(send_pars, self.xvals)
        return model-self.yvals

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
        w,=where(pars < 0)
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


def fit_gprior_exp_mcmc(xdata, ydata, ivar, a=0.25, g0=0.1, gmax=0.87, gmax_min=None, Awidth=1.0):
    """
    This works much better than the lm fitter
    Input is the histogram data.
    """
    import mcmc
    import emcee

    nwalkers=200
    burnin=100
    nstep=100

    print 'fitting exp'

    A=ydata.sum()*(xdata[1]-xdata[0])

    pcen=[A,a,g0,gmax]
    npars=4
    guess=zeros( (nwalkers,npars) )
    guess[:,0] = pcen[0]*(1.+0.1*srandu(nwalkers))
    guess[:,1] = pcen[1]*(1.+0.1*srandu(nwalkers))
    guess[:,2] = pcen[2]*(1.+0.1*srandu(nwalkers))
    guess[:,3] = pcen[3]*(1.+0.1*srandu(nwalkers))


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

    d=diag(pcov)
    perr = sqrt(d)

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

    return res


class GPriorDev(GPrior):
    def __init__(self, pars):
        """
        I'm not using this for anything right now.  In clusterstep
        the exp model works better in all cases pretty much.
        [A,b,c]
        """
        super(GPriorDev,self).__init__(pars)

    def prior2d_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        return gprior2d_dev_vec(self.pars, g)

    def prior2d_gabs_scalar(self, g):
        """
        Get the 2d prior for the input |g| scalar value
        """
        return gprior2d_dev_scalar(self.pars, g)




def gprior2d_dev_vec(pars, g):
    A=pars[0]
    b=pars[1]
    c=pars[2]
    return A*exp( -b*g - c*g**2 )

def gprior1d_dev_vec(pars, g):
    return 2*pi*g*gprior2d_dev_vec(pars, g)

def gprior2d_dev_scalar(pars, g):
    from math import exp
    A=pars[0]
    b=pars[1]
    c=pars[2]
    return A*exp( -b*g - c*g**2 )



class GPriorDevFitter:
    def __init__(self, xvals, yvals):
        """
        Input is the histogram data
        """
        self.xvals=xvals
        self.yvals=yvals

    def __call__(self, pars):
        w,=where(pars < 0)
        if w.size > 0:
            return zeros(self.xvals.size) + numpy.inf

        model=gprior1d_dev_vec(pars, self.xvals)
        return model-self.yvals



def fit_gprior_dev(xdata, ydata):
    """
    Input is the histogram data, should be close to
    normalized
    """
    from scipy.optimize import leastsq


    A=ydata.sum()*(xdata[1]-xdata[0])
    b=2.3
    c=6.7

    pstart=[A,b,c]
    print 'fitting dev'
    print 'pstart:',pstart
    gfitter=GPriorDevFitter(xdata, ydata)
    res = leastsq(gfitter, pstart, full_output=1)

    pars, pcov0, infodict, errmsg, ier = res

    if ier == 0:
        raise ValueError("bad args")

    if ier > 4:
        raise ValueError("fitting failed with\n    %s" % errmsg)

    pcov=None
    perr=None
    if pcov0 is None:
        raise ValueError("pcov0 is None")

    dof=xdata.size-pars.size

    ydiff=gfitter(pars)
    s_sq = (ydiff**2).sum()/dof
    pcov = pcov0 * s_sq 

    d=diag(pcov)
    w,=where(d < 0)

    if w.size > 0:
        raise ValueError("negative diag: %s" % d[w])

    perr = sqrt(d)

    print """    A:  %.6g +/- %.6g
    b:  %.6g +/- %.6g
    c:  %.6g +/- %.6g
    """ % (pars[0],perr[0],
           pars[1],perr[1],
           pars[2],perr[2])

    return {'A':pars[0],
            'b':pars[1],
            'c':pars[2],
            'pars':pars,
            'pcov':pcov,
            'perr':perr}

def fit_gprior_2gauss_cut(xdata, ydata, ivar):
    """
    This works much better than the lm fitter
    Input is the histogram data.
    """
    import mcmc
    import emcee

    nwalkers=800
    burnin=1000
    nstep=100

    A=ydata.sum()#*(xdata[1]-xdata[0])

    A1 = 0.6*A
    A2 = 0.4*A

    sigma1 = 0.02
    sigma2 = 0.3

    pcen = numpy.array([A1,sigma1,A2,sigma2])

    npars=pcen .size
    guess=zeros( (nwalkers,npars) )
    guess[:,0] = pcen[0]*(1.+0.2*srandu(nwalkers))
    guess[:,1] = pcen[1]*(1.+0.2*srandu(nwalkers))
    guess[:,2] = pcen[2]*(1.+0.2*srandu(nwalkers))
    guess[:,3] = pcen[3]*(1.+0.2*srandu(nwalkers))

    gfitter=GPrior2GaussCutFitter(xdata, ydata, ivar)

    print 'pcen:',pcen

    sampler = emcee.EnsembleSampler(nwalkers, 
                                    npars,
                                    gfitter.get_lnprob,
                                    a=2)

    pos, prob, state = sampler.run_mcmc(guess, burnin)
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, nstep)

    arate = sampler.acceptance_fraction.mean()
    print 'arate:',arate
    trials  = sampler.flatchain
    mcmc.plot_results(trials, ptypes=['log','linear','log','linear'])
    

    pars,pcov=mcmc.extract_stats(trials)

    d=diag(pcov)
    perr = sqrt(d)

    gprior=GPrior2GaussCut(pars)

    res={'A1':pars[0],
         'A1_err':perr[0],
         'sigma1':pars[1],
         'sigma1_err':perr[1],
         'A2':pars[2],
         'A2_err':perr[2],
         'sigma2':pars[3],
         'sigma2_err':perr[3],

         'pars':pars,
         'pcov':pcov,
         'perr':perr}


    fmt="""
A1:        %(A1).6g +/- %(A1_err).6g
sigma1:    %(sigma1).6g +/- %(sigma1_err).6g
A2:        %(A2).6g +/- %(A2_err).6g
sigma2:    %(sigma2).6g +/- %(sigma2_err).6g
    """.strip()

    print fmt % res

    return gprior,res


class GPrior2GaussCutFitter(object):
    def __init__(self, xvals, yvals, ivar):
        """
        Input is the histogram data in 1d
        """
        self.xvals=xvals
        self.yvals=yvals
        self.ivar=ivar

    def get_lnprob(self, pars):
        w,=where(pars < 0)
        if w.size > 0:
            return -9.999e20

        model=gprior1d_2gauss_cut(pars, self.xvals)

        lnprob = model
        lnprob -= self.yvals
        lnprob *= lnprob
        lnprob *= self.ivar

        lnprob = lnprob.sum()
        lnprob *= (-0.5)

        return lnprob

class GPrior2GaussCut(GPrior):
    def __init__(self, pars):
        """
        Input is the histogram data in 1d
        """
        self.pars=pars
        self.gmax=1.0

    def prior2d_gabs(self, g):
        """
        Get the 2d prior for the input |g| value(s)
        """
        return gprior2d_2gauss_cut(self.pars,g)

    def prior2d_gabs_scalar(self, g):
        """
        Get the 2d prior for the input |g| value(s)

        using same for now
        """
        return gprior2d_2gauss_cut(self.pars,g)

def gprior2d_2gauss_cut(pars, g):
    gsq = g**2

    amp1=pars[0]
    ivar1 = 1./pars[1]**2
    amp2=pars[2]
    ivar2 = 1./pars[3]**2

    n1=ivar1/(2*pi)
    n2=ivar2/(2*pi)
    p = (amp1*n1*exp(-0.5*gsq*ivar1) + amp2*n2*exp(-0.5*gsq*ivar2))*(1-gsq)**2 
    return p

def gprior1d_2gauss_cut(pars, g):
    return 2*pi*g*gprior2d_2gauss_cut(pars,g)




def fit_gprior_gmix_em(g1, g2, ngauss, n_iter=4000, min_covar=1.e-6, n_init=10):
    import esutil as eu
    #from scikits.learn import mixture
    from sklearn import mixture
    #gmm = mixture.GMM(n_states=ngauss)
    #gmm = mixture.gmm.GMMCenZero(n_components=ngauss,
    gmm = mixture.gmm.GMM(n_components=ngauss,
                                 n_iter=n_iter,
                                 n_init=n_init,
                                 min_covar=min_covar,
                                 covariance_type='spherical')

    #gmm.means_ = numpy.zeros( (ngauss,2) )
    vals = zeros( (g1.size, 2) )
    vals[:,0] = g1
    vals[:,1] = g2

    gmm.fit(vals)

    if not gmm.converged_:
        raise RuntimeError("not converged")
    return gmm

class TPrior(object):
    """
    Prior on T.  
    
    The actual underlying distribution is a lognormal on 

        sigma=sqrt(T/2)

    And it is the mean sigma and the width that are passed to the constructor

    """
    def __init__(self, sigma_mean, sigma_width):
        self.sigma_mean=sigma_mean
        self.sigma_width=sigma_width
        
        self._set_prior()

    def lnprob(self, T):
        sigma=sqrt(T/2)
        return self.ln.lnprob(sigma)

    def _set_prior(self):
        from esutil.random import LogNormal

        self.ln=LogNormal(self.sigma_mean, self.sigma_width)

class MultiGauss:
    def __init__(self, gmm):
        """
        Takes a gmm object

        eval(x) returns normalized evaluation of gaussians
        """
        self.gmm = gmm

    def __repr__(self):
        mess=[]
        gmm = self.gmm
        for i in xrange(gmm.n_states):
            mean = gmm.means[i,0]
            var = gmm.covars[i][0][0]
            weight = gmm.weights[i]
            mess0='p: %.6g x0: %.6g s: %.6g' 
            mess0=mess0 % (weight, mean, sqrt(var))
            mess.append(mess0)

        return '\n'.join(mess)

    def eval(self, x):
        """
        Actually this can just be exp(gmm.eval(x))
        """
        model = numpy.zeros(x.size, dtype='f8')

        gmm = self.gmm
        for i in xrange(gmm.n_states):
            mean = gmm.means[i,0]
            var = gmm.covars[i][0][0]
            weight = gmm.weights[i]
            g = self.gauss(x, mean, var)
            model[:] += weight*g


        return model

    def evalone(self, x, i):
        """
        Just evaluate one of the gaussians
        """

        gmm = self.gmm
        mean = gmm.means[i,0]
        var = gmm.covars[i][0][0]
        weight = gmm.weights[i]
        return weight*self.gauss(x, mean, var)


    def gauss(self, x, mean, var):
        siginv2 = 1./var

        g = exp(-0.5*(x-mean)**2*siginv2)
        norm = sqrt(siginv2/2./numpy.pi)

        g *= norm
        return g



