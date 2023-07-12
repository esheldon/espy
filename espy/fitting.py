import numpy as np


def histogauss(
    data, guess=None, nbin=None, binsize=None, min=None, max=None,
    **plot_keys,
):
    from pprint import pprint

    fitter = GaussFitter(
        data,
        nbin=nbin, binsize=binsize, min=min, max=max,
    )

    if guess is None:
        mn = data.mean()
        sig = data.std()
        A = float(data.size)
        guess = [mn, sig, A]

        print("generated guess:", guess)

    fitter.dofit(guess)
    res = fitter.get_result()
    pprint(res)

    plt = fitter.doplot(**plot_keys)

    return plt, res


class GaussFitter(object):
    """
    Fit a 1-d gaussian
    """
    def __init__(self, data, **keys):
        """
        The data are histogrammed and fit according to the keywords
        """

        self.data = data
        self.conf = keys

        self.use_error = keys.get('use_error', False)

    def get_result(self):
        return {
            'pars': self.pars,
            'pcov': self.pcov,
            'perr': self.perr,
        }

    def dofit(self, guess):
        """
        Run the lm fitter

        guess is [mean, sigma, amp]
        """
        import scipy.optimize

        self._make_hist()

        res = scipy.optimize.leastsq(self._errfunc,
                                     guess,
                                     maxfev=4000,
                                     full_output=1)

        self.pars, self.pcov0, self.infodict, self.errmsg, self.ier = res

        if self.ier == 0:
            # wrong args, this is a bug
            raise ValueError(self.errmsg)

        self.numiter = self.infodict['nfev']
        self.pcov = None
        self.perr = None

        if self.pcov0 is not None:
            self.pcov = self._scale_leastsq_cov(self.pars, self.pcov0)

            d = np.diag(self.pcov)
            w, = np.where(d < 0)

            if w.size == 0:
                # only do if non negative
                self.perr = np.sqrt(d)

    def eval_pars(self, pars):
        """
        [cen, sigma, amp]
        """
        from numpy import exp, sqrt, pi
        mean = pars[0]
        sigma = pars[1]
        amp = pars[2]

        norm = 1.0/sqrt(2*pi)/sigma

        modvals = amp*norm*exp(-0.5*(self.x-mean)**2/sigma**2)

        return modvals

    def _errfunc(self, pars):
        model = self.eval_pars(pars)
        if self.use_error is not None:
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
            self.conf['more'] = True
            h = eu.stat.histogram(self.data, **self.conf)

            self.x = h['center']
            self.y = h['hist']

            if self.use_error:
                self.yerr = np.sqrt(h['hist'])

    def doplot(
        self, aspect=1.618, alpha=0.5, show=True, file=None, dpi=100, **keys,
    ):
        """
        compare fit to data
        """
        import matplotlib.pyplot as pplt

        model = self.eval_pars(self.pars)

        fig, ax = pplt.subplots(**keys)
        ymax = self.y.max()
        ax.set(ylim=(0, 1.1*ymax))
        ax.step(
            self.x, self.y,
            label='data',
            alpha=alpha,
            where='mid',
        )

        if self.use_error:
            ax.errorbar(self.x, self.y, self.yerr)

        ax.plot(self.x, model, label='model', color='red')

        mnstr = r'$\mu: %.3g +/- %.3g$' % (self.pars[0], self.perr[0])
        sigstr = r'$\sigma: %.3g +/- %.3g$' % (self.pars[1], self.perr[1])
        ampstr = 'amp: %.3g +/- %.3g' % (self.pars[2], self.perr[2])

        ax.text(0.1, 0.9, mnstr, ha='left', transform=ax.transAxes)
        ax.text(0.1, 0.8, sigstr, ha='left', transform=ax.transAxes)
        ax.text(0.1, 0.7, ampstr, ha='left', transform=ax.transAxes)
        ax.legend()

        if show:
            pplt.show()

        if file is not None:
            fig.savefig(file, dpi=dpi)

        return fig, ax


class LogNormalFitter(GaussFitter):
    def eval_pars(self, pars):
        """
        [cen, sigma, amp]
        """
        import ngmix
        p = ngmix.priors.LogNormal(pars[0], pars[1])

        return pars[2]*p.get_prob_array(self.x)


def fit_line(x, y, yerr=None):
    """
    Fit a line to the input data

    Parameters
    ----------
    x: array like
        The independent variable
    y: array like
        The dependent variable
    yerr: array like, optional
        Uncertainties on y

    Returns
    -------
    dict with entries
        pars: array
            [slope, offset]
        perr: array
            errors on slope/offset
        slope: float
            Copy of pars[0]
        slope_err: float
            Copy of perr[0]
        offset: flot
            Copy of pars[1]
        offset_err: float
            Copy of perr[1]
        poly: np.poly1d(pars)
    """
    xsum, ysum, xysum, x2sum, wsum = _get_fit_line_quantities(x, y, yerr)

    slope = (xsum * ysum - xysum * wsum) / (xsum**2 - x2sum * wsum)
    offset = (xysum - slope * x2sum) / xsum

    slope_var = wsum / (x2sum * wsum - xsum**2)
    offset_var = x2sum / (x2sum * wsum - xsum**2)
    slope_err = np.sqrt(slope_var)
    offset_err = np.sqrt(offset_var)

    pars = np.array([slope, offset])
    perr = np.array([slope_err, offset_err])

    return {
        'pars': pars,
        'perr': perr,
        'slope': slope,
        'slope_err': slope_err,
        'offset': offset,
        'offset_err': offset_err,
        'poly': np.poly1d(pars),
    }


def _get_fit_line_quantities(x, y, yerr=None):
    x = np.array(x, ndmin=1)
    y = np.array(y, ndmin=1)

    if yerr is None:
        yerr = np.ones(x.size)
    else:
        yerr = np.array(yerr, ndmin=1)

    if x.size != y.size:
        raise ValueError(
            f'x[{x.size}] and y[{y.size}] must be same size for line fit'
        )

    if x.size != yerr.size:
        raise ValueError(
            f'x[{x.size}] and yerr[{yerr.size}] must be same size for line fit'
        )

    wts = 1/yerr**2
    xsum = (x * wts).sum()
    ysum = (y * wts).sum()
    xysum = (x * y * wts).sum()
    x2sum = (x * x * wts).sum()
    wsum = wts.sum()
    return xsum, ysum, xysum, x2sum, wsum


def llsq(X, y, W=None):
    """
    Perform generalized linear least squares fitting

    Parameters
    ----------
    X: array
        The independent variables, (nsamples, npars)
    y: array
        Data with size nsamples
    W: array, optional
        The weights for the data, either shape (nsamples, ) or
        (nsamples,nsamples).  This would usually be the inverse of the data
        covariance matrix.

        If not sent, the identity matrix is used.

    Returns
    -------
    A result dictionary.  The entries are
        pars: The parameter array
        pcov: The covariance of the parameter array
        pcor: The correlation matrix of the parameter array
        perr: The uncertainty on each parameter, the square root
          of the covariance matrix (pcov in the result dict)
        ypred: The prediction of the model, shape (nsamples, )
        chi2: chi^2 including weights
        dof: degrees of freedom
        chi2per: chi^2/dof

    Examples
    --------

    # fit a model y = a * x + b * y

    rng = np.random.RandomState(55)
    nsamples = 20
    npars = 2
    a = 1.0
    b = 2.0
    X = np.zeros((nsamples, npars))
    X[:, 0] = np.linspace(1, 20, nsamples)
    X[:, 1] = np.linspace(-5, 5, nsamples)

    noise = 5
    y = a * X[:, 0] + b * X[:, 1]
    y += rng.normal(scale=noise, size=nsamples)
    yerr = y * 0 + noise

    res = llsq(X, y)
    print('pars:', res['pars'])
    print('perr:', res['perr'])

    # fit a model with an offset and weights
    # y = a * x + b

    rng = np.random.RandomState(88)

    slope = 1.0
    offset = 3.0

    nsamples = 20
    X = np.zeros((nsamples, 2))
    X[:, 0] = np.arange(nsamples)
    X[:, 1] = 1.0  # all 1.0 for the constant offset term

    yerr = 1 + 0.5*rng.uniform(size=nsamples)
    y = slope * X[:, 0] + offset + yerr

    # diagonal covariance, and thus diagonal weights
    W = 1.0/yerr**2
    res = llsq(X, y, W=W)

    print('pars:', res['pars'])
    print('perr:', res['perr'])
    """

    W = _check_llsq_shapes(X=X, y=y, W=W)

    pars = np.linalg.inv(X.T @ W @ X) @ X.T @ W @ y

    ypred = X @ pars.T

    pcov, chi2, dof, chi2per = get_llsq_errors(
        X=X,
        y=y,
        W=W,
        pars=pars,
        ypred=ypred,
    )

    return {
        'pars': pars,
        'pcov': pcov,
        'pcor': cov2cor(pcov),
        'perr': np.sqrt(np.diag(pcov)),
        'ypred': ypred,
        'chi2': chi2,
        'dof': dof,
        'chi2per': chi2per,
    }


def get_llsq_errors(
    X,
    y,
    W,
    pars,
    ypred,
):
    nsamples = y.size
    npars = len(pars)

    residuals = (y - ypred)
    chi2 = residuals.T @ W @ residuals
    dof = nsamples - npars
    chi2per = chi2 / dof

    inv = np.linalg.inv(X.T @ W @ X)
    pcov = inv * chi2per

    return pcov, chi2, dof, chi2per


def _check_llsq_shapes(X, y, W):
    nsamples = y.size
    if X.shape[0] != nsamples:
        raise ValueError(
            f'X has shape {X.shape} but should be '
            f'({nsamples}, {X.shape[1]}) to be compabilty with '
            f'y shape {y.shape}'
        )

    if W is None:
        W = np.diag(np.ones(nsamples))
    elif W.ndim == 1:
        if W.size != nsamples:
            raise ValueError(
                f'W has shape {W.shape} but should be '
                f'({nsamples}, ) or ({nsamples}, {nsamples}) to '
                f'be compatible with y shape {y.shape}'
            )
        W = np.diag(W)
    else:
        if W.shape[0] != nsamples or W.shape[1] != nsamples:
            raise ValueError(
                f'W has shape {W.shape} but should be '
                f'({nsamples}, {nsamples}) to '
                f'be incompatible with y shape {y.shape}'
            )

    return W


def cov2cor(cov):
    """
    Convert the input covariance matrix to a correlation matrix

    corr[i,j] = cov[i,j]/sqrt(cov[i,i]*cov[j,j])

    Parameters
    ----------
    cov: square array
        An NxN covariance matrix

    Returns
    -------
    cor: square array
        The NxN correlation matrix
    """
    cor = np.zeros(cov.shape)

    for ix in range(cov.shape[0]):
        cxx = cov[ix, ix]
        if cxx <= 0.0:
            raise ValueError(
                "diagonal cov[%d,%d]=%e is not positive" % (ix, ix, cxx)
            )
        for iy in range(cov.shape[1]):
            cyy = cov[iy, iy]
            if cyy <= 0.0:
                raise ValueError(
                    "diagonal cov[%d,%d]=%e is not positive" % (iy, iy, cyy)
                )
            cor[ix, iy] = cov[ix, iy] / np.sqrt(cxx * cyy)

    return cor


class PowerLawFitter(object):
    def __init__(self, x, y, yerr=None, method='max', **kw):

        self.x = np.array(x, dtype='f8', ndmin=1)
        self.y = np.array(y, dtype='f8', ndmin=1)

        if yerr is not None:
            yerr = np.array(yerr, dtype='f8', ndmin=1)
        self.yerr = yerr

        self.method = method

        self.npars = 2

        if self.method == 'mcmc':
            if self.yerr is None:
                raise RuntimeError("send yerr= for mcmc")

            self.nwalkers = kw['nwalkers']
            self.burnin = kw['burnin']
            self.nstep = kw['nstep']
            self.a = kw.get("a", 2.0)

        self._set_guess(**kw)

    def get_result(self):
        return self._result

    def dofit(self):
        if self.method == 'max':
            self._dofit_max()
        elif self.method == 'mcmc':
            self._dofit_mcmc()
        else:
            raise ValueError("bad method: '%s'" % self.method)

    def _dofit_mcmc(self):
        import emcee
        import mcmc
        sampler = emcee.EnsembleSampler(
            self.nwalkers,
            self.npars,
            self._get_lnprob,
            a=self.a,
        )

        pos_burn, prob, state = sampler.run_mcmc(self.guess, self.burnin)
        pos, prob, state = sampler.run_mcmc(pos_burn, self.nstep)

        trials = sampler.flatchain
        pars, pcov = mcmc.extract_stats(trials)

        self.sampler = sampler
        self.trials = trials
        perr = np.sqrt(np.diag(pcov))

        self._result = {
            'pars': pars,
            'pcov': pcov,
            'perr': perr,
        }

    def _get_lnprob(self, pars):
        model = self.eval_pars(pars)

        chi2 = ((model-self.y) / self.yerr)**2

        return -0.5*chi2.sum()

    def _dofit_max(self):
        import scipy.optimize
        res = scipy.optimize.leastsq(
            self._errfunc, self.guess,
            full_output=1,
        )
        pars, pcov0, infodict, errmsg, ier = res
        if ier == 0:
            # wrong args, this is a bug
            raise ValueError(errmsg)

        pcov = None
        perr = None

        if pcov0 is not None:
            pcov = self._scale_leastsq_cov(pars, pcov0)

            d = np.diag(pcov)
            w, = np.where(d < 0)

            if w.size == 0:
                # only do if non negative
                perr = np.sqrt(d)

        self._result = {
            'pars': pars,
            'pcov': pcov,
            'perr': perr,
        }

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

    def _set_guess(self, **kw):
        if 'guess' not in kw:
            raise ValueError("send guess= for power law fit")

        self.guess = np.array(kw['guess'], dtype='f8')
        if self.guess.size != self.npars:
            raise ValueError("guess should have size "
                             "%d, got %d" % (self.npars, self.guess.size))

    def eval_pars(self, pars):
        return pars[0]*self.x**pars[1]

    def __call__(self, x):
        res = self.get_result()
        pars = res['pars']
        return pars[0] * x**pars[1]

    def __repr__(self):
        if hasattr(self, '_result'):
            pars = self._result['pars']
            perr = self._result['perr']

            if perr is not None:
                rep = """y = p0*x^p1
    p0: %s +/- %s
    p1: %s +/- %s""" % (pars[0], perr[0], pars[1], perr[1])
            else:
                rep = """y = p0*x^p1
    p0: %s +/- None
    p1: %s +/- None""" % (pars[0], pars[1])

        else:
            rep = ""
        return rep


def test_line(show=False):
    pars = [1.0, 3.0]

    npts = 20
    x = np.arange(npts)
    y = pars[0]*x + pars[1]
    yerr = 1 + 0.5*np.random.uniform(size=npts)
    y += yerr*np.random.normal(size=npts)

    res = fit_line(x, y, yerr=yerr)
    ply = res['poly']

    if show:
        from matplotlib import pyplot as mplt
        fig, ax = mplt.subplots()

        ax.errorbar(x, y, yerr, label='data')
        ax.plot(x, ply(x), label='fit')
        ax.legend()
        mplt.show()


def test_line_err(seed=5, ntrial=10000):
    rng = np.random.RandomState(seed)
    pars = [1.0, 3.0]

    npts = 20
    x = np.arange(npts)
    y0 = pars[0]*x + pars[1]

    slopes = np.zeros(ntrial)
    slope_errors = np.zeros(ntrial)
    offsets = np.zeros(ntrial)
    offset_errors = np.zeros(ntrial)

    for i in range(ntrial):
        yerr = 1 + 0.5*rng.uniform(size=npts)
        y = y0 + yerr*rng.normal(size=npts)

        res = fit_line(x, y, yerr=yerr)

        slopes[i] = res['slope']
        slope_errors[i] = res['slope_err']
        offsets[i] = res['offset']
        offset_errors[i] = res['offset_err']

    slope_mean = slopes.mean()
    offset_mean = offsets.mean()

    print(f'slope mean: {slope_mean} true: {pars[0]}')
    print(f'offset mean: {offset_mean} true: {pars[1]}')

    slope_std = slopes.std()
    offset_std = offsets.std()

    med_slope_err = np.median(slope_errors)
    med_offset_err = np.median(offset_errors)

    print(f'slope err: {slope_std} predicted: {med_slope_err}')
    print(f'offset err: {offset_std} predicted: {med_offset_err}')


def test_line_full(seed=5, ntrial=10000):
    rng = np.random.RandomState(seed)

    pars = [1.0, 3.0]

    npts = 20
    X = np.ones((npts, 2))
    X[:, 0] = np.arange(npts)
    y0 = pars[0]*X[:, 0] + pars[1]

    slopes = np.zeros(ntrial)
    slope_errors = np.zeros(ntrial)
    offsets = np.zeros(ntrial)
    offset_errors = np.zeros(ntrial)

    for i in range(ntrial):
        yerr = 1 + 0.5*rng.uniform(size=npts)
        y = y0 + yerr*rng.normal(size=npts)

        W = 1.0/yerr**2
        res = llsq(X, y, W=W)

        slopes[i] = res['pars'][0]
        slope_errors[i] = res['perr'][0]
        offsets[i] = res['pars'][1]
        offset_errors[i] = res['perr'][1]

    slope_mean = slopes.mean()
    offset_mean = offsets.mean()

    print(f'slope mean: {slope_mean} true: {pars[0]}')
    print(f'offset mean: {offset_mean} true: {pars[1]}')

    slope_std = slopes.std()
    offset_std = offsets.std()

    med_slope_err = np.median(slope_errors)
    med_offset_err = np.median(offset_errors)

    print(f'slope err: {slope_std} predicted: {med_slope_err}')
    print(f'offset err: {offset_std} predicted: {med_offset_err}')


def test_lin2(show=False):
    seed = 55
    rng = np.random.RandomState(seed)

    # fit a model y = a * x + b * y

    nsamples = 20
    npars = 2
    atrue = 1.0
    btrue = 2.0
    X = np.zeros((nsamples, npars))
    X[:, 0] = np.linspace(1, 20, nsamples)
    X[:, 1] = np.linspace(-5, 5, nsamples)

    noise = 5
    y = atrue * X[:, 0] + btrue * X[:, 1]
    y += rng.normal(scale=noise, size=nsamples)
    yerr = y * 0 + noise

    res = llsq(X, y)
    print('pars:', res['pars'])
    print('pcov:', res['pcov'])

    if show:
        import proplot as pplt
        fig, axs = pplt.subplots(nrows=2, share=False)

        kw = {'marker': 'o', 'markersize': 4, 'ls': ''}
        axs[0].errorbar(X[:, 0], y, yerr, **kw)
        axs[0].plot(X[:, 0], res['ypred'])

        s = X[:, 1].argsort()
        axs[1].errorbar(X[s, 1], y[s], yerr[s], **kw)
        axs[1].plot(X[s, 1], res['ypred'][s])
        pplt.show()
