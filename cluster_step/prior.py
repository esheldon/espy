import numpy
from numpy import diag, where, exp, sqrt, zeros
from math import pi

class GPriorExp:
    def __init__(self, A, a, g0, gmax):
        self.A=A
        self.a=a
        self.g0=g0
        self.gmax=gmax

        self.pars = [A, a, g0, gmax]

        self.maxval = self(0., 0.)

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
        return gprior2d_exp_vec(self.pars, g)

    def prior1d(self, g):
        """
        Get the 1d prior for an input |g| value(s).
        """
        return 2*pi*g*self.prior2d_gabs(g)

    def sample1d(self, nrand):
        """
        Get random |g| from the 1d distribution

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

            # generate total g**2 in [0,1)
            grand = random.random(nleft)

            # now finally the height from [0,maxval)
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
        g1 = zeros(nrand)
        g2 = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand2 = random.random(nleft)
            grand = sqrt(grand2)
            # now uniform angles
            rangle = random.random(nleft)*2*pi

            # now get cartesion locations in g1,g2 plane
            g1rand = grand*cos(rangle)
            g2rand = grand*sin(rangle)

            # now finally the height from [0,maxval)
            h = self.maxval*random.random(nleft)

            pvals = self(g1rand, g2rand)

            w,=where(h < pvals)
            if w.size > 0:
                g1[ngood:ngood+w.size] = g1rand[w]
                g2[ngood:ngood+w.size] = g2rand[w]
                ngood += w.size
                nleft -= w.size

        return g1, g2


    def set_maxval1d(self):
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

def gprior1d_exp_vec(pars, g):
    return 2*pi*g*gprior2d_exp_vec(pars, g)


class GPriorFitter:
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



def fit_gprior_exp(xdata, ydata):
    """
    Input is the histogram data, should be close to
    normalized
    """
    from scipy.optimize import leastsq


    A=ydata.sum()*(xdata[1]-xdata[0])
    a=0.25
    g0=0.1
    gmax=0.87

    #pstart=[A,a,g0,gmax]
    pstart=[A,a,g0]
    print 'pstart:',pstart
    gfitter=GPriorFitter(xdata, ydata, gmax=gmax)
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
            'gmax':gmax,
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



