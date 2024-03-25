"""
Module:
    mcmc

Classes:
    MH: A class for running Monte Carlo Markov Chains using
        metropolis hastings
    MHTester: A class for testing the MH class.

stats:
    extract_stats: extract mean and covariance
    plot_results: Plot points and histograms for MCMC trials.

testing:
    test: Run the MHTester
    testmany:  Run mutliple realizations of the tester and plot a histogram
        of all the means for all realizations.

See the docs for these individual classes for more details.

Revision History:
    Created: 2010-04-02, Erin Sheldon, BNL
"""

import numpy as np
from numpy import array, zeros, sqrt, where, diag, log
from numpy.random import randn
from .mcmc_stats import extract_mcmc_stats


class GaussStepper(object):
    """
    A class to take gaussian steps.

    gs=GausssianStepper(sigmas)
    newpars = gs(oldpars)
    """

    def __init__(self, sigmas, rng):
        import numpy as np
        self._rng = rng

        self._sigmas = np.array(sigmas, dtype='f8')

    def __call__(self, pars):
        sigmas = self._sigmas
        return pars + sigmas * self._rng.normal(size=sigmas.size)


def print_stats(means, cov, names=None):
    npar = len(means)
    for i in range(npar):
        if names is not None:
            name = names[i]
        else:
            name = "%s" % i
        print("%s: %.16g +/- %.16g" % (name, means[i], sqrt(cov[i, i])))


def extract_maxlike_stats(data, burnin):
    nuse = data.size - burnin
    npar = data["pars"].shape[1]

    maxi = data["loglike"][burnin:].argmax()

    max_like = zeros(npar, dtype="f8")
    error = zeros(npar, dtype="f8")

    for i in range(npar):
        max_like[i] = data["pars"][burnin + maxi, i]

        # variance around this point
        vi = ((data["pars"][burnin:, i] - max_like[i]) ** 2).sum()
        vi /= nuse - 1.0

        error[i] = sqrt(vi)

    return max_like, error


def test_line(burnin=1000, nstep=10000, show=False):
    """
    run all steps at once so we can plot burnin phase
    """
    from .mh import run_mh

    rng = np.random.default_rng()

    pars = [2.0, 1.0]
    xmin = -1.0
    xmax = 1.0
    nx = 10
    yerr = 0.1
    x, y, yerr = noisy_line(pars, xmin, xmax, nx, yerr)

    lineprob = LineProb(x, y, yerr)

    stepper = GaussStepper(sigmas=[0.02, 0.02], rng=rng)
    res = run_mh(
        func=lineprob,
        stepper=stepper,
        rng=rng,
        start=pars,
        nstep=nstep,
    )

    if show:
        import matplotlib.pyplot as mplt
        import corner

        # plot the burnin
        trials = res['trials']
        steps = np.arange(res['loglike'].size, dtype="i4")

        fig, axs = mplt.subplots(nrows=2, ncols=2)
        offburn_plot = axs[0, 0]
        slopeburn_plot = axs[0, 1]

        offburn_plot.set(ylabel="offset")

        slopeburn_plot.set(ylabel="slope", xlabel="step number")

        offburn_plot.plot(steps[0:burnin], trials[:burnin, 0], color="red")
        slopeburn_plot.plot(steps[0:burnin], trials[:burnin, 1], color="red")
        offburn_plot.plot(steps[burnin:], trials[burnin:, 0], color="black")
        slopeburn_plot.plot(steps[burnin:], trials[burnin:, 1], color="black")

        # get status for chain
        parfit, cov = extract_mcmc_stats(trials[burnin:, :])
        yfit = parfit[0] * x + parfit[1]

        ax = axs[1, 0]
        ax.errorbar(x, y, yerr, color='black', label='data', ls='', marker='.')
        ax.plot(x, yfit, color='blue', label='fit')
        ax.legend()

        axs[1, 1].axis('off')

        mplt.show()

        _ = corner.corner(
            trials[burnin:], labels=['offset', 'slope'], show_titles=True,
            bins=20,
        )
        mplt.show()


class EmceeFitter(object):
    """
    Base class to fit the using emcee

    the user must over-ride these functions
        - get_guess() - return array of shape (nwalkers,ndims)
        - get_lnprob(pars) - return the log likelihood for the input pars
        - get_npoints() - return total number of data points
            this is for doing statistics such as the chisq probability

    """

    def __init__(self, nwalkers, burnin, nstep, a):

        self._nwalkers = nwalkers
        self._burnin = burnin
        self._nstep = nstep
        self._a = a

    def get_result(self):
        """
        A dictionary with the results

        pars,pcov,perr,lnprob,aic,bic,chi2,dof,chi2per,prob
        """
        return self._result

    def get_trials(self):
        """
        Return all points in the chain
        """
        return self._trials

    def get_lnprobs(self):
        """
        Return all log probabilities in the chain
        """
        return self._sampler.lnprobability.reshape(
            self._nwalkers * self._nstep
        )

    def get_guess(self):
        """
        array shape (nwalkers,ndims)
        """
        raise RuntimeError("over-ride")

    def get_npoints(self):
        """
        Total number of data points
        """
        raise RuntimeError("over-ride")

    def get_lnprob(self, pars):
        """
        scalar log likelihood or probability
        """
        raise RuntimeError("over-ride")

    def _run_trials(self):
        import emcee

        guess = self.get_guess()
        npars = guess.shape[1]

        sampler = emcee.EnsembleSampler(
            self._nwalkers, npars, self.get_lnprob, a=self._a
        )

        pos, prob, state = sampler.run_mcmc(guess, self._burnin)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self._nstep)

        self._trials = sampler.flatchain
        self._sampler = sampler

        self._calc_stats()

    def _calc_stats(self):
        pars, pcov = extract_mcmc_stats(self._trials)
        perr = sqrt(diag(pcov))

        npars = len(pars)

        lnprob = self.get_lnprob(pars)

        npts = self.get_npoints()

        chi2 = lnprob / (-0.5)
        dof = npts - npars
        chi2per = chi2 / dof

        aic = -2 * lnprob + 2 * npars
        bic = -2 * lnprob + npars * log(npts)

        res = {
            "pars": pars,
            "perr": perr,
            "pcov": pcov,
            "lnprob": lnprob,
            "aic": aic,
            "bic": bic,
            "chi2": chi2,
            "dof": dof,
            "chi2per": chi2per,
        }
        try:
            import scipy.stats

            prob = scipy.stats.chisqprob(chi2, dof)

            res["prob"] = prob
        except Exception:
            pass

        self._result = res

    def corner(self, burnin, **keys):
        import corner
        return corner.corner(
            self._trials,
            **keys
        )

    def __repr__(self):
        pars = self._result["pars"]
        perr = self._result["perr"]
        npars = len(pars)

        rep = []
        for i in range(npars):
            r = "  p%d: %.4g +/- %.4g" % (i, pars[i], perr[i])
            rep.append(r)

        rep = "\n".join(rep)
        return rep


class LogNormalFitter(EmceeFitter):
    def __init__(self, x, y, guess, nwalkers, burnin, nstep, a=2, **keys):

        super(LogNormalFitter, self).__init__(nwalkers, burnin, nstep, a)

        self._x = array(x, dtype="f8", ndmin=1, copy=False)
        self._y = array(y, dtype="f8", ndmin=1, copy=False)
        self._guess0 = guess

        self._xmax = self._x.max()

        if self._guess0 is not None:
            self._guess0 = array(self._guess0, dtype="f8", ndmin=1, copy=False)

        self._yerr = keys.get("yerr", None)
        self._ivar = keys.get("ivar", None)

        self._width = keys.get("width", False)
        if self._width is not None:
            self._width = array(self._width, dtype="f8", ndmin=1, copy=False)

        self._set_ivar()

        self._check_data()
        self._check_guess_and_width()

        self._set_full_guess()

        self._lowval = -9.999e20
        self._run_trials()

    def get_model(self):
        """
        Get a normalized lognormal model at the expectation value of the
        parameters.  This function integrates to unity.

            model = lnf.get_model()
            y = model(x)

        To get the scaled function with the right overall normalization

            y=model.scaled(x)
        """
        from esutil.random import LogNormal

        pars = self._result["pars"]
        A = pars[0]
        mean = pars[1]
        sigma = pars[2]

        return LogNormal(mean, sigma, norm=A)

    def get_guess(self):
        """
        over-ride for base class
        """
        return self._guess

    def get_npoints(self):
        """
        over-ride for base class
        """
        return self._x.size

    def get_lnprob(self, pars):
        (w,) = where(pars <= 1.0e-6)
        if w.size > 0:
            return self._lowval
        if pars[1] > self._xmax:
            return self._lowval

        model = self._get_model_at_pars(pars)

        chi2 = ((model - self._y) ** 2) * self._ivar
        lnprob = -0.5 * chi2.sum()

        if self._width is not None:
            w = self._width
            lnprior = ((pars - self._guess0) / self._width) ** 2
            lnprior = -0.5 * lnprior.sum()

            lnprob += lnprior

        return lnprob

    def _get_model_at_pars(self, pars):
        from esutil.random import LogNormal

        A = pars[0]
        mean = pars[1]
        sigma = pars[2]

        ln = LogNormal(mean, sigma, norm=A)
        return ln.scaled(self._x)

        # some optimizations here; could construct with norm=A and call
        # scaled(x).  This way we avoid some error checking that has already
        # been done

        # lnobj=LogNormal(mean, sigma)
        # lnprob=lnobj._lnprob(self._x)

        # prob=exp(lnprob)

        # return A*prob

    def corner(self, **keys):
        import corner
        labels = ["A", "mean", "sigma"]
        return corner.corner(
            self._trials,
            labels=labels,
            show_titles=True,
        )

    def _set_ivar(self):
        if self._yerr is not None:
            ivar = 1.0 / self._yerr**2
            self._ivar = array(ivar, dtype="f8", ndmin=1, copy=False)
        elif self._ivar is not None:
            self._ivar = array(self._ivar, dtype="f8", ndmin=1, copy=False)
        else:
            self._ivar = None

    def _check_data(self):
        (wbad,) = where(~np.isfinite(self._x))
        if wbad.size != 0:
            raise ValueError("%d x values are not finite" % wbad.size)
        (wbad,) = where(~np.isfinite(self._y))
        if wbad.size != 0:
            raise ValueError("%d y values are not finite" % wbad.size)

        if self._x.size != self._y.size:
            raise ValueError("x,y must be same size")

        if self._ivar is not None:
            (wbad,) = where(~np.isfinite(self._ivar))
            if wbad.size != 0:
                raise ValueError("%d ivar values are not finite" % wbad.size)

            if (self._x.size != self._ivar.size) and (self._ivar.size != 1):
                raise ValueError("x,y, and yerr/ivar must be same size")

        (wbad,) = where(self._x <= 0)
        if wbad.size != 0:
            raise ValueError("x values must all be > 0")

    def _check_guess_and_width(self):
        # make sure there is a guess for each dimension
        guess0 = self._guess0
        if guess0.size != 3:
            raise ValueError("guess should be length 3 for [norm,mean,sigma]")

        (wbad,) = where(~np.isfinite(guess0))
        if wbad.size != 0:
            mess = []
            mess.append("bad guess:")
            for i in wbad:
                mess.append("  %i %g" % (i, guess0[i]))
            mess = "\n".join(mess)
            raise ValueError(mess)

        if self._width is not None:
            if self._width.size != 3:
                raise ValueError(
                    "width should be length 3 for [norm,mean,sigma]"
                )
            (wbad,) = where(~np.isfinite(self._width))
            if wbad.size != 0:
                raise ValueError(
                    "Some width are not finite: "
                    "[%s, %s, %s]" % tuple(self._width)
                )

    def _set_full_guess(self):
        from esutil.random import srandu

        guess0 = self._guess0

        npars = len(guess0)
        nwalkers = self._nwalkers
        guess = zeros((nwalkers, npars))

        for i in range(npars):
            if guess0[i] == 0:
                guess[:, i] = guess0[i] + 0.1 * srandu(nwalkers)
            else:
                guess[:, i] = guess0[i] * (1 + 0.1 * srandu(nwalkers))

        self._guess = guess


def test_lognormal():
    import biggles
    import esutil as eu
    from esutil.random import LogNormal, srandu
    from esutil.stat import histogram

    n = 1000
    nwalkers = 100
    burnin = 100
    nstep = 100

    mean = 8
    sigma = 3
    ln = LogNormal(mean, sigma)
    vals = ln.sample(n)

    binsize = 0.5

    plt = eu.plotting.bhist(vals, binsize=binsize, show=False)

    h = histogram(vals, binsize=binsize, more=True)
    herr = sqrt(h["hist"])
    herr = herr.clip(1.0, herr.max())

    guess = [
        n * (1.0 + 0.1 * srandu()),
        mean * (1.0 + 0.1 * srandu()),
        sigma * (1.0 + 0.1 * srandu()),
    ]
    guess = [n * binsize, mean, sigma]

    print("guess:", guess)
    nlf = LogNormalFitter(
        h["center"], h["hist"], guess, nwalkers, burnin, nstep, yerr=herr
    )

    print(nlf)

    model = nlf.get_model()

    yvals = model.scaled(h["center"])
    plt.add(biggles.Curve(h["center"], yvals, color="blue"))
    plt.show()


class PolyFitter(object):
    """
    Fit a polygon to the input points using an affine invariant MCMC chain

    The emcee module is used for the MCMC chain.
    """

    def __init__(
        self,
        order,
        x,
        y,
        nwalkers,
        burnin,
        nstep,
        guess=None,
        a=2,
        yerr=None,
        ivar=None,
    ):

        self.order = order
        self.x = array(x, dtype="f8", ndmin=1, copy=False)
        self.y = array(y, dtype="f8", ndmin=1, copy=False)

        self.nwalkers = nwalkers
        self.burnin = burnin
        self.nstep = nstep
        self.a = a

        self.yerr = yerr
        self.ivar = ivar

        self._set_ivar()

        self._check_data()
        self._set_guess(guess)

        self._run_trials()

    def get_result(self):
        return self._result

    def get_poly(self):
        return self._ply

    def __call__(self, x):
        return self._ply(x)

    def _set_ivar(self):
        if self.yerr is not None:
            ivar = 1.0 / self.yerr**2
            self.ivar = array(ivar, dtype="f8", ndmin=1, copy=False)
        elif self.ivar is not None:
            self.ivar = array(self.ivar, dtype="f8", ndmin=1, copy=False)
        else:
            self.ivar = None

    def _check_data(self):
        (wbad,) = where(~np.isfinite(self.x))
        if wbad.size != 0:
            raise ValueError("%d x values are not finite" % wbad.size)
        (wbad,) = where(~np.isfinite(self.y))
        if wbad.size != 0:
            raise ValueError("%d y values are not finite" % wbad.size)

        if self.x.size != self.y.size:
            raise ValueError("x,y must be same size")

        if self.ivar is not None:
            (wbad,) = where(~np.isfinite(self.ivar))
            if wbad.size != 0:
                raise ValueError("%d ivar values are not finite" % wbad.size)

            if (self.x.size != self.ivar.size) and (self.ivar.size != 1):
                raise ValueError("x,y, and yerr/ivar must be same size")

    def _check_guess(self, guess):
        if guess.size != (self.order + 1):
            raise ValueError("guess should be length order+1")
        (wbad,) = where(~np.isfinite(guess))
        if wbad.size != 0:
            mess = []
            mess.append("bad guess:")
            for i in wbad:
                mess.append("  %i %g" % (i, guess[i]))
            mess = "\n".join(mess)
            raise ValueError(mess)

    def _set_guess(self, guess0):
        from esutil.random import srandu

        if guess0 is not None:
            guess0 = array(guess0, dtype="f8", ndmin=1, copy=True)
            self._check_guess(guess0)
        else:
            guess0 = np.polyfit(self.x, self.y, self.order)
            self.guess = guess0

        npars = len(guess0)
        nwalkers = self.nwalkers
        guess = zeros((nwalkers, npars))

        for i in range(npars):
            if guess0[i] == 0:
                guess[:, i] = guess0[i] + 0.1 * srandu(nwalkers)
            else:
                guess[:, i] = guess0[i] * (1 + 0.1 * srandu(nwalkers))

        self._guess = guess

    def _run_trials(self):
        import emcee

        guess = self._guess
        npars = guess.shape[1]
        sampler = emcee.EnsembleSampler(
            self.nwalkers, npars, self.get_lnprob, a=self.a
        )

        pos, prob, state = sampler.run_mcmc(guess, self.burnin)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self.nstep)

        self.trials = sampler.flatchain

        self._calc_stats()

    def _calc_stats(self):
        pars, pcov = extract_mcmc_stats(self.trials)
        perr = sqrt(diag(pcov))

        npars = len(pars)

        lnprob = self.get_lnprob(pars)

        chi2 = lnprob / (-0.5)
        dof = self.x.size - npars
        chi2per = chi2 / dof

        aic = -2 * lnprob + 2 * npars
        bic = -2 * lnprob + npars * log(self.x.size)

        res = {
            "pars": pars,
            "perr": perr,
            "pcov": pcov,
            "lnprob": lnprob,
            "aic": aic,
            "bic": bic,
            "chi2": chi2,
            "dof": dof,
            "chi2per": chi2per,
        }
        try:
            import scipy.stats

            prob = scipy.stats.chisqprob(chi2, dof)

            res["prob"] = prob
        except Exception:
            pass

        self._result = res
        try:
            self._ply = np.poly1d(res["pars"])
        except Exception:
            print("could not set poly")
            self._ply = None

    def get_lnprob(self, pars):
        from numpy import poly1d

        ply = poly1d(pars)
        model = ply(self.x)
        chi2 = (self.y - model) ** 2
        if self.ivar is not None:
            chi2 *= self.ivar

        return -0.5 * chi2.sum()

    def corner(self, **kw):
        import corner
        return corner.corner(self.trials, **kw)

    def __repr__(self):
        pars = self._result["pars"]
        perr = self._result["perr"]
        npars = len(pars)

        rep = []
        header = []
        for i in range(npars):
            h = "p%d" % i
            power = npars - i - 1
            if power > 0:
                h += " x^%d" % power
            header.append(h)

            r = "  p%d: %.4g +/- %.4g" % (i, pars[i], perr[i])
            rep.append(r)

        header = " + ".join(header)
        rep = [header] + rep
        rep = "\n".join(rep)
        return rep


class LineProb:
    def __init__(self, x, y, yerr):
        self.x = x
        self.y = y
        self.yerr = yerr
        self.ivar = 1.0 / yerr**2
        assert x.size == y.size == yerr.size

    def __call__(self, pars):
        yfunc = self.line_func(pars)
        chi2 = self.ivar * (self.y - yfunc) ** 2
        return -0.5 * chi2.sum()

    def line_func(self, pars):
        return pars[0] * self.x + pars[1]


def noisy_line(pars, xmin, xmax, nx, yerr):

    x = np.linspace(xmin, xmax, nx)
    y = pars[0] * x + pars[1]

    y += yerr * randn(nx)

    yerr_vals = np.array([yerr] * x.size, dtype="f8")

    return x, y, yerr_vals


def gaussfunc(mean, sigma, xvals):

    gnorm = 1.0 / np.sqrt(2.0 * np.pi * sigma**2)
    gauss = np.exp(-0.5 * (xvals - mean) ** 2 / sigma**2)

    return gauss * gnorm
