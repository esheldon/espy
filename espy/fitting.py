import numpy
import scipy.optimize
from pprint import pprint


def histogauss(
    data, guess=None, nbin=None, binsize=None, min=None, max=None,
    **plot_keys,
):
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

            d = numpy.diag(self.pcov)
            w, = numpy.where(d < 0)

            if w.size == 0:
                # only do if non negative
                self.perr = numpy.sqrt(d)

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
                self.yerr = numpy.sqrt(h['hist'])

    def doplot(self, aspect=1.618, alpha=0.5, show=True, file=None, dpi=100, **keys):
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
        p=ngmix.priors.LogNormal(pars[0],pars[1])

        return pars[2]*p.get_prob_array(self.x)


def fit_line(x, y, yerr=None, **kw):
    lf=LineFitter(x, y, yerr=yerr, **kw)
    lf.dofit()
    return lf

class LineFitter(object):
    def __init__(self, x, y, yerr=None, method='max', **kw):

        self.x=numpy.array(x,dtype='f8',ndmin=1)
        self.y=numpy.array(y,dtype='f8',ndmin=1)
        if yerr is not None:
            yerr=numpy.array(yerr,dtype='f8',ndmin=1)
        self.yerr=yerr

        self.method=method

        self.npars=2

        if self.method=='mcmc':
            if self.yerr is None:
                raise RuntimeError("send yerr= for mcmc")
            self.nwalkers=kw['nwalkers']
            self.burnin=kw['burnin']
            self.nstep=kw['nstep']
            self.a = kw.get("a",2.0)

        self._set_guess(**kw)

    def doplot(self, **keys):
        import biggles

        show=keys.pop('show',False)
        file=keys.pop('file',None)

        res=self.get_result()
        ply = self.get_poly()

        plt=biggles.FramedPlot()

        pts = biggles.Points(
            self.x, self.y,
            type='filled circle',
            size=2,
        )
        pts.label='data'

        plt.add(pts)
        if self.yerr is not None:
            err = biggles.SymmetricErrorBarsY(
                self.x, self.y, self.yerr,
            )
            plt.add(err)

        line = biggles.Curve(
            self.x, ply(self.x),
            color='blue',
        )
        line.label='%.3g x + %.3g' % tuple(res['pars'])
        plt.add(line)

        key=biggles.PlotKey(
            0.1,0.9,[pts,line],
            halign='left',
        )
        plt.add(key)

        if show:
            plt.show(**keys)
        if file is not None:
            print("writing:",file)
            if '.eps' in file:
                plt.write_eps(file)
            elif '.png' in file:
                width=keys.pop('width',800)
                height=keys.pop('height',800)
                plt.write_img(width,height,file)

        return plt


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


    def _set_guess(self, **kw):
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
        res=self.get_result()
        return numpy.poly1d(res['pars'])

    def __call__(self, x):
        """
        pars order same as for numpy.poly1d
        """
        res=self.get_result()
        pars=res['pars']
        return pars[0]*x + pars[1]

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


class PowerLawFitter(LineFitter):
    def _set_guess(self, **kw):
        if 'guess' not in kw:
            raise ValueError("send guess= for power law fit")

        self.guess=numpy.array(kw['guess'], dtype='f8')
        if self.guess.size != self.npars:
            raise ValueError("guess should have size "
                             "%d, got %d" % (self.npars,self.guess.size))

    def eval_pars(self, pars):
        return pars[0]*self.x**pars[1]

    def __call__(self, x):
        res=self.get_result()
        pars=res['pars']
        return pars[0]*x**pars[1]

    def __repr__(self):
        if hasattr(self,'_result'):
            pars=self._result['pars']
            perr=self._result['perr']

            if perr is not None:
                rep = """y = p0*x^p1
    p0: %s +/- %s
    p1: %s +/- %s""" % (pars[0],perr[0],pars[1],perr[1])
            else:
                rep = """y = p0*x^p1
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
