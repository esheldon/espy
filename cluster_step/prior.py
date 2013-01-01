import numpy
from numpy import diag, where, cos, sin, exp, sqrt, zeros, random
from math import pi

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

class GPriorDev(GPrior):
    def __init__(self, pars):
        """
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


class GPriorExpFitter:
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




def fit_gprior_exp(xdata, ydata):
    """
    Input is the histogram data, should be close to
    normalized
    """
    from scipy.optimize import leastsq


    A=ydata.sum()*(xdata[1]-xdata[0])
    a=0.25
    g0=0.1

    pstart=[A,a,g0]
    print 'pstart:',pstart
    gfitter=GPriorExpFitter(xdata, ydata)
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

    print """
    A:    %.6g +/- %.6g
    a:    %.6g +/- %.6g
    g0:   %.6g +/- %.6g
    """ % (pars[0],perr[0],
           pars[1],perr[1],
           pars[2],perr[2])

    return {'A':pars[0],
            'a':pars[1],
            'g0':pars[2],
            #'gmax':pars[3],
            'gmax':gfitter.gmax,
            'pars':pars,
            'pcov':pcov,
            'perr':perr}


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

    print """
    A:  %.6g +/- %.6g
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



def fit_gprior_exp_gmix(g1, g2):
    import esutil as eu
    from scikits.learn import mixture
    ngauss=5
    gmm = mixture.GMM(n_states=ngauss)

    vals = zeros(g1.size, 2)
    vals[:,0] = g1
    vals[:,1] = g2

    gmm.fit(vals, n_iter=400)#, min_covar=1.e-6)

    mg = MultiGauss(gmm)
    print mg

    return mg

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



