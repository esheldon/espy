import numpy
from numpy import diag, where, exp, sqrt, zeros

class GPriorExp:
    def __init__(self, A, a, g0, gmax=0.90):
        self.A=A
        self.a=a
        self.g0=g0
        self.gmax=gmax

        self.pars = [A, a, g0]

    def __call__(self, g1, g2):
        g = sqrt(g1**2 + g2**2)
        return self.prior_gabs(g)

    def prior_gabs(self, g):
        return  gprior_exp_vec(self.pars, g, gmax=self.gmax)


def gprior_exp_vec(pars, g, gmax=0.90):
    A=pars[0]
    a=pars[1]
    g0=pars[2]

    prior=zeros(g.size)

    w,=where(g < gmax)
    if w.size > 0:
        numer = A*g*(1-exp( (g-gmax)/a ))
        denom = (1+g)*sqrt(g**2 + g0**2)

        prior[w]=numer/denom

    return prior

class GPriorFitter:
    def __init__(self, xvals, yvals, gmax=0.90):
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

        model=gprior_exp_vec(pars, self.xvals, gmax=self.gmax)
        return model-self.yvals

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

def fit_gprior_exp(xdata, ydata):
    """
    Input is the histogram data, should be close to
    normalized
    """
    from scipy.optimize import leastsq

    #gmax=0.9
    gmax=0.85

    A=ydata.sum()*(xdata[1]-xdata[0])
    a=0.25
    g0=0.1

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
    A:  %.6g +/- %.6g
    a:  %.6g +/- %.6g
    g0: %.6g +/- %.6g
    """ % (pars[0],perr[0],
           pars[1],perr[1],
           pars[2],perr[2])

    return {'A':pars[0],
            'a':pars[1],
            'g0':pars[2],
            'gmax':gmax,
            'pars':pars,
            'pcov':pcov,
            'perr':perr}

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



