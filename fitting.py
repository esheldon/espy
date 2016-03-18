from __future__ import print_function, division
import numpy
import scipy.optimize
from pprint import pprint

def histogauss(data, guess=None, **keys):
    fitter=GaussFitter(data, **keys)

    if guess is None:
        mn=data.mean()
        sig=data.std()
        A=float(data.size)
        guess=[mn, sig, A]

        print("generated guess:",guess)

    fitter.dofit(guess)
    res=fitter.get_result()
    pprint(res)

    plt=fitter.make_plot(**keys)

    return plt


class GaussFitter(object):
    """
    Fit a 1-d gaussian
    """
    def __init__(self, data, **keys):
        """
        The data are histogrammed and fit according to the keywords
        """

        self.data=data
        self.conf=keys

        self.use_error = keys.get('use_error',False)


    def get_result(self):
        return {'pars':self.pars,
                'pcov':self.pcov,
                'perr':self.perr}

    def dofit(self, guess):
        """
        Run the lm fitter

        guess is [mean, sigma, amp]
        """

        self._make_hist()

        res=scipy.optimize.leastsq(self._errfunc,
                                   guess,
                                   maxfev=4000,
                                   full_output=1)

        self.pars, self.pcov0, self.infodict, self.errmsg, self.ier = res

        if self.ier == 0:
            # wrong args, this is a bug
            raise ValueError(self.errmsg)

        self.numiter = self.infodict['nfev']
        self.pcov=None
        self.perr=None

        if self.pcov0 is not None:
            self.pcov = self._scale_leastsq_cov(self.pars, self.pcov0)

            d=numpy.diag(self.pcov)
            w,=numpy.where(d < 0)

            if w.size == 0:
                # only do if non negative
                self.perr = numpy.sqrt(d)

    def eval_pars(self, pars):
        """
        [cen, sigma, amp]
        """
        from numpy import exp, sqrt, pi
        mean=pars[0]
        sigma=pars[1]
        amp=pars[2]

        norm = 1.0/sqrt(2*pi)/sigma

        modvals = amp*norm*exp( -0.5*(self.x-mean)**2/sigma**2 )

        return modvals

    def _errfunc(self, pars):
        model = self.eval_pars(pars)
        if not self.use_error is None:
            diff = model-self.y
        else:
            diff = (model-self.y)/self.yerr

        return diff


    def _scale_leastsq_cov(self, pars, pcov):
        """
        Scale the covariance matrix returned from leastsq; this will
        recover the covariance of the parameters in the right units.
        """
        dof = (self.x.size-len(pars))
        s_sq = (self._errfunc(pars)**2).sum()/dof
        return pcov * s_sq 

    def _make_hist(self):
        import esutil as eu

        if not hasattr(self, 'x'):
            print('histogramming')
            self.conf['more']=True
            h=eu.stat.histogram(self.data, **self.conf)

            self.x=h['center']
            self.y=h['hist']

            if self.use_error:
                self.yerr=numpy.sqrt(h['hist'])

    def make_plot(self, **keys):
        """
        compare fit to data
        """
        import biggles

        model=self.eval_pars(self.pars)

        plt=biggles.FramedPlot(**keys)

        x0=self.x[0]
        binsize=self.x[1]-self.x[0]
        h=biggles.Histogram(self.y, x0=x0, binsize=binsize, color='black')
        h.label='data'

        plt.add(h)
        if self.use_error:
            ep=biggles.SymmetricErrorBarsY(self.x, self.y, self.yerr)
            plt.add(ep)

        fh=biggles.Histogram(model, x0=x0, binsize=binsize, color='red')
        fh.label='fit'

        plt.add(fh)

        key=biggles.PlotKey(0.1, 0.9, [h, fh], halign='left')
        plt.add(key)

        mnstr=r'$\mu: %g +/- %g$' % (self.pars[0],self.perr[0])
        sigstr=r'$\sigma: %g +/- %g$' % (self.pars[1],self.perr[1])
        ampstr='amp: %g +/- %g' % (self.pars[2],self.perr[2])

        mnlab = biggles.PlotLabel(0.9,0.3,mnstr, halign='right')
        siglab = biggles.PlotLabel(0.9,0.2,sigstr, halign='right')
        amplab = biggles.PlotLabel(0.9,0.1,ampstr, halign='right')

        plt.add(mnlab, siglab, amplab)

        show=keys.get('show',True)
        if show:
            plt.show()

        return plt

class LogNormalFitter(GaussFitter):
    def eval_pars(self, pars):
        """
        [cen, sigma, amp]
        """
        import ngmix
        p=ngmix.priors.LogNormal(pars[0],pars[1])

        return pars[2]*p.get_prob_array(self.x)


def fit_line(x, y, yerr=None, **kw):
    lf=LineFitter(x, y, yerr=yerr, **kw)
    lf.dofit()
    return lf

class LineFitter(object):
    def __init__(self, x, y, yerr=None, method='max', **kw):
        self.x=x
        self.y=y
        self.yerr=yerr
        self.method=method

        self.npars=2

        if self.method=='mcmc':
            if yerr is None:
                raise RuntimeError("send yerr= for mcmc")
            self.nwalkers=kw['nwalkers']
            self.burnin=kw['burnin']
            self.nstep=kw['nstep']
            self.a = 2.0

        self._set_guess()

    def get_result(self):
        return self._result

    def dofit(self):
        if self.method=='max':
            self._dofit_max()
        elif self.method=='mcmc':
            self._dofit_mcmc()
        else:
            raise ValueError("bad method: '%s'" % self.method)

    def _dofit_mcmc(self):
        import emcee
        import mcmc
        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        self.npars, 
                                        self._get_lnprob,
                                        a=self.a)

        pos_burn, prob, state = sampler.run_mcmc(self.guess, self.burnin)
        pos, prob, state = sampler.run_mcmc(pos_burn, self.nstep)

        trials  = sampler.flatchain
        pars, pcov = mcmc.extract_stats(trials)

        self.sampler=sampler
        self.trials=trials
        perr=numpy.sqrt(numpy.diag(pcov))

        self._result={'pars':pars, 'pcov':pcov, 'perr':perr}

    def _dofit_max(self):
        res=scipy.optimize.leastsq(self._errfunc, self.guess,
                                   full_output=1)
        pars, pcov0, infodict, errmsg, ier = res
        if ier == 0:
            # wrong args, this is a bug
            raise ValueError(errmsg)

        numiter = infodict['nfev']
        pcov=None
        perr=None

        if pcov0 is not None:
            pcov = self._scale_leastsq_cov(pars, pcov0)

            d=numpy.diag(pcov)
            w,=numpy.where(d < 0)

            if w.size == 0:
                # only do if non negative
                perr = numpy.sqrt(d)
        
        self._result={'pars':pars, 'pcov':pcov, 'perr':perr}


    def _set_guess(self):
        best_fit=numpy.polyfit(self.x, self.y, 1)
        
        if self.method=='mcmc':
            guess = numpy.zeros( (self.nwalkers, self.npars) )

            rnums = 2.0*(numpy.random.random(self.nwalkers)-0.5)
            guess[:,0] += 0.01*rnums
            rnums = 2.0*(numpy.random.random(self.nwalkers)-0.5)
            guess[:,1] += 0.01*rnums

            self.guess=guess
        else:
            self.guess=best_fit

    def eval_pars(self, pars):
        return pars[0]*self.x + pars[1]

    def _get_lnprob(self, pars):
        model = self.eval_pars(pars)

        chi2 = ( (model-self.y)/self.yerr )**2

        return -0.5*chi2.sum()

    def _errfunc(self, pars):
        model = self.eval_pars(pars)
        if self.yerr is None:
            diff = model-self.y
        else:
            diff = (model-self.y)/self.yerr

        return diff

    def _scale_leastsq_cov(self, pars, pcov):
        """
        Scale the covariance matrix returned from leastsq; this will
        recover the covariance of the parameters in the right units.
        """
        dof = (self.x.size-len(pars))
        s_sq = (self._errfunc(pars)**2).sum()/dof
        return pcov * s_sq 


    def get_poly(self):
        return numpy.poly1d(self.pars)

    def __call__(self, x):
        """
        pars order same as for numpy.poly1d
        """
        return self.pars[0]*x + self.pars[1]

    def __repr__(self):
        if hasattr(self,'_result'):
            pars=self._result['pars']
            perr=self._result['perr']

            if perr is not None:
                rep = """y = p0*x + p1
    p0: %s +/- %s
    p1: %s +/- %s""" % (pars[0],perr[0],pars[1],perr[1])
            else:
                rep = """y = p0*x + p1
    p0: %s +/- None
    p1: %s +/- None""" % (pars[0],pars[1])

        else:
            rep=""
        return rep


def test_line(show=False):
    pars = [1,3]
    x = numpy.arange(20)
    y = pars[0]*x + pars[1]
    yerr = numpy.zeros(len(x)) + 2 + 0.1*numpy.random.random(len(x))
    y += yerr*numpy.random.randn(len(x))

    lf = LineFitter(x, y, yerr=yerr)
    print('guess:',lf.guess)
    print(lf)

    if show:
        import biggles
        plt=biggles.FramedPlot()
        plt.add(biggles.Points(x, y,type='filled circle'))
        plt.add(biggles.SymmetricErrorBarsY(x, y, yerr))
        plt.add(biggles.Curve(x, lf(x),color='blue'))
        plt.show()
