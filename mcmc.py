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

import numpy
from numpy import array, zeros, ones, sqrt, arange, isfinite, \
        where, diag, exp, log
from numpy.random import randn
from sys import stdout
import os
from sys import stderr


def extract_stats(data, weights=None):
    """
    Extract mean and covariance from mcmc trials

    parameters
    ----------
    trials: array
        (nsteps, npars) array.
    weights: optional
        Optional weights for each trial
    """
    if weights is not None:
        return _extract_weighted_stats(data, weights)
    else:
        return _extract_stats(data)

def _extract_stats(data):
    ntrials=data.shape[0]
    npar = data.shape[1]

    means = zeros(npar,dtype='f8')
    cov = zeros( (npar,npar), dtype='f8')

    for i in xrange(npar):
        means[i] = data[:, i].mean()

    num=ntrials

    for i in xrange(npar):
        idiff = data[:,i]-means[i]
        for j in xrange(i,npar):
            if i == j:
                jdiff = idiff
            else:
                jdiff = data[:,j]-means[j]

            cov[i,j] = (idiff*jdiff).sum()/(num-1)

            if i != j:
                cov[j,i] = cov[i,j]

    return means, cov

def _extract_weighted_stats(data, weights):
    if weights.size != data.shape[0]:
        raise ValueError("weights not same size as data")

    npar = data.shape[1]

    wsum = weights.sum()
    wsum2 = wsum**2

    means = zeros(npar,dtype='f8')
    cov = zeros( (npar,npar), dtype='f8')

    for i in xrange(npar):
        dsum = (data[:, i]*weights).sum()
        means[i] = dsum/wsum

    for i in xrange(npar):
        idiff = data[:,i]-means[i]
        for j in xrange(i,npar):
            if i == j:
                jdiff = idiff
            else:
                jdiff = data[:,j]-means[j]

            wvar = ( weights*idiff*jdiff ).sum()/wsum
            cov[i,j] = wvar

            if i != j:
                cov[j,i] = cov[i,j]

    return means, cov

def plot_results(trials, **keys):
    """
    Plot the points and histograms of trials from an MCMC chain.
    """
    import biggles
    import esutil

    npars=trials.shape[1]

    fontsize_min=keys.get('fontsize_min',1)
    biggles.configure( 'default', 'fontsize_min', fontsize_min)
    weights=keys.get('weights',None)

    binfac=keys.get('binfac',0.2)
    names=keys.get('names',None)
    show=keys.get('show',True)
    ptypes=keys.get('ptypes',['linear']*npars)

    means,cov = extract_stats(trials,weights=weights)
    errs=sqrt(diag(cov)) 

    plt=biggles.Table(npars,2)

    ind = numpy.arange(trials.shape[0])

    for i in xrange(npars):
        if names is not None:
            name=names[i]
        else:
            name=r'$p_{%d}$' % i

        burn_plot_i = esutil.plotting.bscatter(ind,trials[:,i],
                                               type='solid',
                                               xlabel='step',
                                               ylabel=name,
                                               show=False)
        plt[i,0] = burn_plot_i

        if ptypes[i] == 'linear':
            vals=trials[:,i]
            bsize = binfac*errs[i]
            xlabel=name
        else:
            vals=numpy.log10(trials[:,i])
            bsize=0.2*vals.std()
            xlabel=r'$log_{10}(%s)$' % name

        hdict = esutil.stat.histogram(vals,
                                      binsize=bsize, 
                                      weights=weights,
                                      more=True)
        if weights is not None:
            hist=hdict['whist']
            hplot = biggles.Curve(hdict['center'],
                                  hdict['whist'])
        else:
            hist=hdict['hist']
            hplot = biggles.Histogram(hdict['hist'], 
                                      x0=hdict['low'][0], 
                                      binsize=bsize)
        plti=biggles.FramedPlot()

        plti.xlabel=xlabel

        hmax=hist.max()
        plti.yrange=[-0.05*hmax, 1.2*hmax]

        plti.add(hplot)
            
        lab = r'$<%s> = %0.4g \pm %0.4g$' % (name,means[i],errs[i])
        plab = biggles.PlotLabel(0.1,0.8,lab,
                                 halign='left',
                                 color='blue')

        plti.add(plab)

        plt[i,1]=plti
    
    plt.title=keys.get('title',None)

    if show:
        plt.show()

    return plt

def plot_results_separate(trials, **keys):
    """
    Plot the points and histograms of trials from an MCMC chain.
    """
    import biggles
    import esutil
    import images

    npars=trials.shape[1]

    fontsize_min=keys.get('fontsize_min',1)
    biggles.configure( 'default', 'fontsize_min', fontsize_min)
    weights=keys.get('weights',None)

    binfac=keys.get('binfac',0.2)
    names=keys.get('names',None)
    show=keys.get('show',True)
    ptypes=keys.get('ptypes',['linear']*npars)

    means,cov = extract_stats(trials,weights=weights)
    errs=sqrt(diag(cov)) 

    nrows, ncols =images.get_grid(npars)
    burn_plt=biggles.Table(nrows, ncols)
    hist_plt=biggles.Table(nrows, ncols)

    ind = numpy.arange(trials.shape[0])

    prow=0
    pcol=0
    for i in xrange(npars):

        prow=i/ncols
        pcol=i % ncols

        if names is not None:
            name=names[i]
        else:
            name=r'$p_{%d}$' % i

        lab = r'$<%s> = %0.4g \pm %0.4g$' % (name,means[i],errs[i])
        plab = biggles.PlotLabel(0.1,0.8,lab,
                                 halign='left',
                                 color='blue')

        # steps
        burn_plot_i = esutil.plotting.bscatter(ind,trials[:,i],
                                               type='solid',
                                               xlabel='step',
                                               ylabel=name,
                                               show=False)

        burn_plot_i.add(plab)
        burn_plt[prow,pcol] = burn_plot_i

        # hist
        vals=trials[:,i]
        bsize = binfac*errs[i]

        hdict = esutil.stat.histogram(vals,
                                      binsize=bsize, 
                                      weights=weights,
                                      more=True)
        if weights is not None:
            hist=hdict['whist']
            hplot = biggles.Curve(hdict['center'],
                                  hdict['whist'])
        else:
            hist=hdict['hist']
            hplot = biggles.Histogram(hdict['hist'], 
                                      x0=hdict['low'][0], 
                                      binsize=bsize)
        plti=biggles.FramedPlot()

        plti.xlabel=name

        hmax=hist.max()
        plti.yrange=[-0.05*hmax, 1.2*hmax]

        plti.add(hplot)
        plti.add(plab)

        hist_plt[prow,pcol]=plti
   
    burn_plt.title=keys.get('title',None)
    hist_plt.title=keys.get('title',None)

    if show:
        burn_plt.show()
        hist_plt.show()

    return burn_plt, hist_plt


class MH(object):
    """
    Run a Monte Carlo Markov Chain (MCMC) using metropolis hastings.
    
    parameters
    ----------
    lnprob_func: function or method
        A function to calculate the log proability given the input
        parameters.  Can be a method of a class.
            ln_prob = lnprob_func(pars)
            
    stepper: function or method 
        A function to take a step given the input parameters.
        Can be a method of a class.
            newpars = stepper(pars)

    seed: floating point, optional
        An optional seed for the random number generator.

    examples
    ---------
    m=mcmc.MCMC(lnprob_func, stepper, seed=34231)
    m.run(pars_start, nstep)
    trials = m.get_trials()
    loglike = m.get_loglike()
    arate = m.get_acceptance_rate()

    """
    def __init__(self, lnprob_func, stepper,
                 seed=None, random_state=None):
        self._lnprob_func=lnprob_func
        self._stepper=stepper

        self.set_random_state(seed=seed, state=random_state)
        self.reset(seed=seed)

    def reset(self, seed=None):
        """
        Clear all data
        """
        self._trials=None
        self._loglike=None
        self._accepted=None

    def set_random_state(self, seed=None, state=None):
        """
        set the random state

        parameters
        ----------
        seed: integer, optional
            If state= is not set, the random state is set to
            numpy.random.RandomState(seed=seed)
        state: optional
            A random number generator with method .uniform()
        """
        if state is not None:
            self._random_state=state
        else:
            self._random_state=numpy.random.RandomState(seed=seed)

    def run(self, pars_start, nstep):
        """
        Run the MCMC chain.  Append new steps if trials already
        exist in the chain.

        parameters
        ----------
        pars_start: sequence
            Starting point for the chain in the n-d parameter space.
        nstep: integer
            Number of steps in the chain.
        """
        
        self._init_data(pars_start, nstep)

        for i in xrange(1,nstep):
            self._step()

        self._arate=self._accepted.sum()/float(self._accepted.size)
        return self._trials[-1,:]
    run_mcmc=run

    def get_trials(self):
        """
        Get the trials array
        """
        return self._trials

    def get_loglike(self):
        """
        Get the trials array
        """
        return self._loglike

    def get_acceptance_rate(self):
        """
        Get the acceptance rate
        """
        return self._arate

    def get_accepted(self):
        """
        Get the accepted array
        """
        return self._accepted

    def _step(self):
        """
        Take the next step in the MCMC chain.  
        
        Calls the stepper lnprob_func methods sent during construction.  If the
        new loglike is not greater than the previous, or a uniformly generated
        random number is greater than the the ratio of new to old likelihoods,
        the new step is not used, and the new parameters are the same as the
        old.  Otherwise the new step is kept.

        This is an internal function that is called by the .run method.
        It is not intended for call by the user.
        """

        index=self._current

        oldpars=self._oldpars
        oldlike=self._oldlike

        # Take a step and evaluate the likelihood
        newpars = self._stepper(oldpars)
        newlike = self._lnprob_func(newpars)

        log_likeratio = newlike-oldlike

        randnum = self._random_state.uniform()
        log_randnum = numpy.log(randnum)

        # we allow use of -infinity as a sign we are out of bounds
        if (isfinite(newlike) 
                and ( (newlike > oldlike) | (log_randnum < log_likeratio)) ):

            self._accepted[index]  = 1
            self._loglike[index]   = newlike
            self._trials[index, :] = newpars

            self._oldpars = newpars
            self._oldlike = newlike

        else:
            self._accepted[index] = 0
            self._loglike[index]  = oldlike
            self._trials[index,:] = oldpars

        self._current += 1

    def _init_data(self, pars_start, nstep):
        """
        Set the trials and accept array.
        """

        pars_start=array(pars_start,dtype='f8',copy=False)
        npars = pars_start.size

        self._trials   = numpy.zeros( (nstep, npars) )
        self._loglike  = numpy.zeros(nstep)
        self._accepted = numpy.zeros(nstep, dtype='i1')
        self._current  = 1

        self._oldpars = pars_start.copy()
        self._oldlike = self._lnprob_func(pars_start)

        self._trials[0,:] = pars_start
        self._loglike[0]  = self._oldlike
        self._accepted[0] = 1



class GaussStepper(object):
    """
    A class to take gaussian steps.

    gs=GausssianStepper(sigmas)
    newpars = gs(oldpars)
    """
    def __init__(self, sigmas):
        sigmas=numpy.asanyarray(sigmas, dtype='f8', copy=False)

        self._sigmas=sigmas
        self._ndim=sigmas.size

    def __call__(self, pars):
        sigmas=self._sigmas
        return pars + sigmas*randn(sigmas.size)
    
def print_stats(means, cov, names=None):
    npar=len(means)
    for i in xrange(npar):
        if names is not None:
            name=names[i]
        else:
            name='%s' % i
        print '%s: %.16g +/- %.16g' % (name,means[i],sqrt(cov[i,i]))

def extract_maxlike_stats(data, burnin):
    nuse = data.size-burnin
    npar = data['pars'].shape[1]

    maxi = data['loglike'][burnin:].argmax()

    max_like = zeros(npar,dtype='f8')
    error    = zeros(npar,dtype='f8')

    for i in xrange(npar):
        max_like[i] = data['pars'][burnin+maxi, i]

        # variance around this point
        vi = ( (data['pars'][burnin:, i] - max_like[i])**2 ).sum()
        vi /= (nuse-1.)

        error[i] = sqrt(vi)

    return max_like, error



class MHTester:
    """

    A simple class for testing the MCMC code.  This is a good template for how
    to begin your own MCMC project.

        import mcmc
        tc = mcmc.MHTester(type="constant")
        tc.run_test()

    You can then use the plot_results method to see a histogram of the results.

    inputs to constructur:
      type: The type of model.  Default is "constant"
    """
    def __init__(self, type="constant", verbose=False):
        self.setup(type, verbose=verbose)


    def setup(self,type="constant", verbose=False):
        self.type = type
        self.verbose=verbose
        if self.type is "constant":

            # Set up the "data"
            ny = 10
            val,sigma = self.get_truepars()

            self.true_pars = val
            self.true_error = sigma
            self.npars = 1

            self.npoints = ny
            self.y = zeros(ny, dtype='f8')
            self.y[:] = val
            self.y[:] += sigma*randn(ny)

            self.yerr = zeros(ny, dtype='f8')
            self.yerr[:] = sigma

            self.ivar = 1.0/self.yerr**2

            # step sizes
            self.step_sizes = sigma/sqrt(ny)

            #self.parguess=val + 3*sigma*randn()
            self.parguess=0.

            if verbose:
                stdout.write('  type: "constant"\n')
                stdout.write("  true_pars: %s\n  npars: %s\n  step_sizes: %s\n" %
                             (self.true_pars,self.npars,self.step_sizes))


    def get_truepars(self):
        if self.type is "constant":
            val = 10.0
            sigma = 1.0
            return val, sigma
        else:
            raise ValueError("only support type='constant'")



    # pars must be an array
    def step(self,pars):
        return pars + self.step_sizes*randn(self.npars)

    def get_loglike(self,pars):
        if self.type == 'constant':
            chi2 = self.ivar*(self.y-pars)**2
            chi2 = -0.5*chi2.sum()
        else:
            raise ValueError("only support type='constant'")
        return chi2

    def run_test(self, nstep=10000):

        if self.verbose:
            stdout.write("  nstep: %s\n" % nstep)

        m=MH(self.get_loglike, self.step)
        m.run(self.parguess, nstep)
        self.trials = m.get_trials()
        self.loglike = m.get_loglike()
        self.arate = m.get_acceptance_rate()

    def run_multiple(self, ntrial, nstep=1000, burnin=100):
        """
        Run multiple tests and save in self.meanvals
        Note default is fewer steps per mcmc
        """
        self.meanvals = zeros(ntrial, dtype='f4')

        for i in xrange(ntrial):
            # create a new realization
            self.setup()
            self.run_test(nstep)
            self.meanvals[i] = self.trials[burnin:,0].mean()



    def compare_results(self, burnin=100, doplot=False, hardcopy=False):
        if self.type == 'constant':
            par0 = self.trials[burnin:, :]

            mn  = par0.mean()
            err = par0.std()

            stdout.write("Using burn in of: %s\n" % burnin)
            stdout.write("True value: %s\n" % self.true_pars)
            stdout.write("True error on each of %s points: %s\n" % (self.npoints,self.true_error))
            expected_error = self.true_error/numpy.sqrt(self.npoints)
            stdout.write("Expected error from combined points: %s\n" % expected_error)
            stdout.write("  Mean from trials: %s\n" % mn)
            stdout.write("  Error from trials: %s\n" % err)

            if doplot:
                self.plot_results(burnin, hardcopy=hardcopy)


        else:
            raise ValueError("only support type='constant'")

    def plot_multiple(self, hardcopy=False):
        """
        Plot the results from run_multiple, saved in meanvals.
        """
        import biggles
        from esutil import stat
        if self.type == 'constant':

            meanvals = self.meanvals

            xmin = meanvals.min()
            xmax = meanvals.max()
            xstd = meanvals.std()

            binsize=xstd*0.2
            bindata = stat.histogram(meanvals, binsize=binsize, more=True)

            plt = biggles.FramedPlot()

            x_range = [xmin,xmax]
            #plt.x1.range = x_range

            d=biggles.Histogram( bindata['hist'], x0=xmin, binsize=binsize)
            d.label = 'Trials Means'


            # get the expected gaussian
            expected_error = self.true_error/numpy.sqrt(self.npoints)

            xvals = numpy.arange(x_range[0], x_range[1], 0.02, dtype='f8')
            gpoints=self.gaussfunc(self.true_pars, expected_error, xvals)

            gpoints *= meanvals.size*binsize

            g = biggles.Curve(xvals, gpoints, color='blue')
            g.label = 'Expected Distribution'

            k = biggles.PlotKey( .1, .9, [d,g] )

            plt.add(d,k,g)

            plt.xlabel = 'Trial Means'
            plt.ylabel = 'count'

            if hardcopy:
                fname='mcmc-constant-multi.eps'
                stdout.write('Writing test file hardcopy: %s\n' % fname)
                plt.write_eps(fname)
            plt.show()
        else:
            raise ValueError("only support type='constant'")



    def plot_results(self, burnin=100, hardcopy=False):
        import biggles
        from esutil import stat
        if self.type == 'constant':
            par0 = self.trials[burnin:, 0]

            xmin = par0.min()
            binsize=0.05
            bindata = stat.histogram(par0, binsize=binsize, more=True)

            plt = biggles.FramedPlot()

            x_range = [8.5,11.5]
            plt.x1.range = x_range

            d=biggles.Histogram( bindata['hist'], x0=xmin, binsize=binsize)
            d.label = 'trials'


            # get the expected gaussian
            expected_error = self.true_error/numpy.sqrt(self.npoints)

            xvals = numpy.arange(x_range[0], x_range[1], 0.02, dtype='f8')
            gpoints=self.gaussfunc(self.true_pars, expected_error, xvals)

            gpoints *= par0.size*binsize

            g = biggles.Curve(xvals, gpoints, color='blue')
            g.label = 'Expected Distribution'

            k = biggles.PlotKey( .1, .9, [d,g] )

            plt.add(d,k,g)

            plt.xlabel = 'trial values'
            plt.ylabel = 'count'

            if hardcopy:
                fname='mcmc-constant.eps'
                stdout.write('Writing test file hardcopy: %s\n' % fname)
                plt.write_eps(fname)
            plt.show()
        else:
            raise ValueError("only support type='constant'")

    def gaussfunc(self,mean,sigma,xvals):

        gnorm = 1.0/numpy.sqrt(2.0*numpy.pi*sigma**2)
        gauss = numpy.exp(-0.5*(xvals - mean)**2/sigma**2 )

        return gauss*gnorm



def test(burnin=1000, nstep=10000, doplot=False, hardcopy=False):
    """
    Name:
        test
    Purpose:
        Run the MHTester.
    """
    tc = MHTester("constant")
    tc.run_test(nstep)
    print 'acceptance rate:',tc.arate
    tc.compare_results(doplot=doplot, hardcopy=hardcopy)

def testmany(ntrial, **keys):
    """
    Just do a bunch of realizations and make sure we get the right mean
    """

    tc = MHTester("constant")
    tc.run_multiple(ntrial, **keys)
    tc.plot_multiple()

def test_line(burnin=1000, nstep=10000, doplot=False):
    """
    run all steps at once so we can plot burnin phase
    """
    import esutil
    pars = [1.0,1.0]
    xmin = -1.0
    xmax =  1.0
    nx = 10
    yerr = 0.1
    x, y, yerr = noisy_line(pars, xmin, xmax, nx, yerr)

    LF = LinFitter(x, y, yerr)

    fitter = MH(LF.get_loglike, LF.step)

    # bad guess
    parguess = [ pars[0] + 0.2, pars[1]-0.2 ]

    ntot=nstep+burnin
    pos = fitter.run(parguess, ntot)

    data=fitter.get_trials()

    if doplot:
        import biggles
        from esutil import stat

        burnin = 1000

        # plot the burnin
        tab = biggles.Table(2,1)

        steps = numpy.arange(ntot, dtype='i4')

        offset_steps_plot = biggles.FramedPlot()
        offset_steps_plot.ylabel = 'offset'

        slope_steps_plot = biggles.FramedPlot()
        slope_steps_plot.ylabel = 'slope'
        slope_steps_plot.xlabel = 'step number'

        offset_burnin_curve = biggles.Curve(steps[0:burnin], data[0:burnin,0], color='red')
        slope_burnin_curve = biggles.Curve(steps[0:burnin], data[0:burnin,1], color='red')
        offset_rest_curve = biggles.Curve(steps[burnin:], data[burnin:,0])
        slope_rest_curve = biggles.Curve(steps[burnin:], data[burnin:,1])


        offset_steps_plot.add( offset_burnin_curve, offset_rest_curve )
        slope_steps_plot.add( slope_burnin_curve, slope_rest_curve )


        tab[0,0] = offset_steps_plot
        tab[1,0] = slope_steps_plot

        tab.show()

        # get status for chain
        parfit, cov = extract_stats(data[burnin:, :])
        errfit = sqrt(diag(cov))

        # plot the histograms and comparison plot

        tab = biggles.Table(2,2)
        offsets = data[burnin:,0]
        slopes  = data[burnin:,1]

        offset_binsize = offsets.std()*0.2
        slope_binsize = slopes.std()*0.2

        offset_hist = stat.histogram(offsets, binsize=offset_binsize, more=True)  
        slope_hist = stat.histogram(slopes, binsize=slope_binsize, more=True)  


        offset_phist = biggles.FramedPlot()
        offset_phist.add( biggles.Histogram(offset_hist['hist'], 
                                            x0=offset_hist['low'][0], 
                                            binsize=offset_binsize) )
        offset_phist.xlabel = 'Offsets'
        offset_phist.add( biggles.PlotLabel(0.1,0.9,
                                            "offset=%0.2f +/- %0.2f" % (parfit[0],errfit[0]),
                                            halign='left') )
        offset_phist.yrange = [0,offset_hist['hist'].max()*1.2]

        slope_phist = biggles.FramedPlot()
        slope_phist.add( biggles.Histogram(slope_hist['hist'], 
                                           x0=slope_hist['low'][0], 
                                           binsize=slope_binsize) )
        slope_phist.xlabel = 'slopes'
        slope_phist.add( biggles.PlotLabel(0.1,0.9,"slope=%0.2f +/- %0.2f" % (parfit[1],errfit[1]),
                                           halign='left'))
        slope_phist.yrange = [0,slope_hist['hist'].max()*1.2]

        
        tab[0,0] = offset_phist
        tab[0,1] = slope_phist

        # now plot original data and best fit par

        yfit = parfit[0]*x + parfit[1]

        fitplt = biggles.FramedPlot()
        data_errbar = biggles.SymmetricErrorBarsY(x,y,yerr)
        data_points = biggles.Points(x,y,type='filled circle')
        data_points.label = 'Data'

        yfit_curve = biggles.Curve(x,yfit,color='blue')

        key = biggles.PlotKey( .1, .9, [data_points, yfit_curve] )

        fitplt.add( data_errbar, data_points, yfit_curve, key )

        tab[1,0] = fitplt

        tab.show()




    return data
                        
   
def get_grid(ntot):
    sq=int(sqrt(ntot))
    if ntot==sq*sq:
        return (sq,sq)
    elif ntot <= sq*(sq+1):
        return (sq,sq+1)
    else:
        return (sq+1,sq+1)

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

        self._nwalkers=nwalkers
        self._burnin=burnin
        self._nstep=nstep
        self._a=a

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
        return self._sampler.lnprobability.reshape(self._nwalkers*self._nstep)

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

        guess=self.get_guess()
        npars=guess.shape[1]

        sampler = emcee.EnsembleSampler(self._nwalkers, 
                                        npars,
                                        self.get_lnprob,
                                        a=self._a)

        pos, prob, state = sampler.run_mcmc(guess, self._burnin)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self._nstep)

        self._trials  = sampler.flatchain
        self._sampler = sampler 

        self._calc_stats()

    def _calc_stats(self):
        pars,pcov=extract_stats(self._trials)
        perr=sqrt(diag(pcov))

        npars=len(pars)

        lnprob=self.get_lnprob(pars)

        npts=self.get_npoints()

        chi2=lnprob/(-0.5)
        dof=npts-npars
        chi2per = chi2/dof

        aic = -2*lnprob + 2*npars
        bic = -2*lnprob + npars*log(npts)

        res={'pars':pars,
             'perr':perr,
             'pcov':pcov,
             'lnprob':lnprob,
             'aic':aic,
             'bic':bic,
             'chi2': chi2,
             'dof':dof,
             'chi2per':chi2per}
        try:
            import scipy.stats
            prob = scipy.stats.chisqprob(chi2, dof)

            res['prob']=prob
        except:
            pass

        self._result=res

    def plot_trials(self, **keys):
        plt=plot_results(self._trials, **keys)
        return plt

    def __repr__(self):
        pars=self._result['pars']
        perr=self._result['perr']
        npars=len(pars)

        rep=[]
        for i in xrange(npars):
            r= "  p%d: %.4g +/- %.4g" % (i,pars[i],perr[i])
            rep.append(r)
        
        rep='\n'.join(rep)
        return rep


class LogNormalFitter(EmceeFitter):
    def __init__(self, x, y, guess, nwalkers, burnin, nstep, a=2, **keys):

        super(LogNormalFitter,self).__init__(nwalkers, burnin, nstep, a)

        self._x=array(x, dtype='f8', ndmin=1, copy=False)
        self._y=array(y, dtype='f8', ndmin=1, copy=False)
        self._guess0=guess

        self._xmax=self._x.max()

        if self._guess0 is not None:
            self._guess0=array(self._guess0, dtype='f8', ndmin=1, copy=False)

        self._yerr=keys.get('yerr',None)
        self._ivar=keys.get('ivar',None)

        self._width=keys.get('width',False)
        if self._width is not None:
            self._width=array(self._width, dtype='f8', ndmin=1, copy=False)

        self._set_ivar()

        self._check_data()
        self._check_guess_and_width()

        self._set_full_guess()

        self._lowval=-9.999e20
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

        pars=self._result['pars']
        A=pars[0]
        mean=pars[1]
        sigma=pars[2]

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
        w,=where(pars <= 1.e-6)
        if w.size > 0:
            return self._lowval
        if pars[1] > self._xmax:
            return self._lowval


        model=self._get_model_at_pars(pars)

        chi2=((model-self._y)**2)*self._ivar
        lnprob = -0.5*chi2.sum()

        if self._width is not None:
            w=self._width
            g=self._guess0
            lnprior=( (pars-self._guess0)/self._width )**2
            lnprior = -0.5*lnprior.sum()
            
            lnprob += lnprior

        return lnprob

    def _get_model_at_pars(self, pars):
        from esutil.random import LogNormal
        A=pars[0]
        mean=pars[1]
        sigma=pars[2]

        ln=LogNormal(mean, sigma, norm=A)
        return ln.scaled(self._x)

        # some optimizations here; could construct with norm=A and call
        # scaled(x).  This way we avoid some error checking that has already
        # been done

        #lnobj=LogNormal(mean, sigma)
        #lnprob=lnobj._lnprob(self._x)

        #prob=exp(lnprob)

        #return A*prob

    def plot_trials(self, **keys):
        names=['A','mean','sigma']
        keys['names']=names
        plt=plot_results(self._trials, **keys)
        return plt

    def _set_ivar(self):
        if self._yerr is not None:
            ivar=1./self._yerr**2
            self._ivar=array(ivar, dtype='f8', ndmin=1, copy=False)
        elif self._ivar is not None:
            self._ivar=array(self._ivar, dtype='f8', ndmin=1, copy=False)
        else:
            self._ivar=None

    def _check_data(self):
        wbad,=where(numpy.isfinite(self._x) == False)
        if wbad.size!=0:
            raise ValueError("%d x values are not finite" % wbad.size)
        wbad,=where(numpy.isfinite(self._y) == False)
        if wbad.size!=0:
            raise ValueError("%d y values are not finite" % wbad.size)

        if self._x.size != self._y.size:
            raise ValueError("x,y must be same size")

        if self._ivar is not None:
            wbad,=where(numpy.isfinite(self._ivar) == False)
            if wbad.size!=0:
                raise ValueError("%d ivar values are not finite" % wbad.size)

            if (self._x.size != self._ivar.size) and (self._ivar.size != 1):
                raise ValueError("x,y, and yerr/ivar must be same size")

        wbad,=where(self._x <= 0)
        if wbad.size != 0:
            raise ValueError("x values must all be > 0")

    def _check_guess_and_width(self):
        # make sure there is a guess for each dimension
        guess0=self._guess0
        if guess0.size != 3:
            raise ValueError("guess should be length 3 for [norm,mean,sigma]")

        wbad,=where(numpy.isfinite(guess0) == False)
        if wbad.size != 0:
            mess=[]
            mess.append("bad guess:")
            for i in wbad:
                mess.append("  %i %g" % (i,guess0[i]))
            mess='\n'.join(mess)
            raise ValueError(mess)

        if self._width is not None:
            if self._width.size != 3:
                raise ValueError("width should be length 3 for [norm,mean,sigma]")
            wbad,=where(numpy.isfinite(self._width) == False)
            if wbad.size != 0:
                raise ValueError("Some width are not finite: [%s, %s, %s]" % tuple(self._width))


    def _set_full_guess(self):
        from esutil.random import srandu

        guess0=self._guess0

        npars=len(guess0)
        nwalkers=self._nwalkers
        guess=zeros((nwalkers,npars))

        for i in xrange(npars):
            if guess0[i]==0:
                guess[:,i] = guess0[i]+0.1*srandu(nwalkers)
            else:
                guess[:,i] = guess0[i]*(1+0.1*srandu(nwalkers))

        self._guess=guess

def test_lognormal():
    import biggles
    import esutil as eu
    from esutil.random import LogNormal, srandu
    from esutil.stat import histogram

    n=1000
    nwalkers=100
    burnin=100
    nstep=100

    mean=8
    sigma=3
    ln=LogNormal(mean,sigma)
    vals=ln.sample(n)

    binsize=0.5

    plt=eu.plotting.bhist(vals, binsize=binsize,show=False)

    h=histogram(vals, binsize=binsize,more=True)
    herr=sqrt(h['hist'])
    herr=herr.clip(1.0, herr.max())

    guess=[n*(1. + .1*srandu()),
           mean*(1. + .1*srandu()),
           sigma*(1. + .1*srandu())]
    guess=[n*binsize,mean,sigma]

    print 'guess:',guess
    nlf=LogNormalFitter(h['center'], h['hist'], guess, nwalkers, burnin, nstep,
                        yerr=herr)

    print nlf

    res=nlf.get_result()
    
    model=nlf.get_model()

    yvals=model.scaled(h['center'])
    plt.add(biggles.Curve(h['center'], yvals, color='blue'))
    plt.show()


                        
class PolyFitter(object):
    """
    Fit a polygon to the input points using an affine invariant MCMC chain

    The emcee module is used for the MCMC chain.
    """
    def __init__(self, order, x, y, nwalkers, burnin, nstep, 
                 guess=None, a=2, yerr=None, ivar=None):

        self.order=order
        self.x=array(x, dtype='f8', ndmin=1, copy=False)
        self.y=array(y, dtype='f8', ndmin=1, copy=False)

        self.nwalkers=nwalkers
        self.burnin=burnin
        self.nstep=nstep
        self.a=a

        self.yerr=yerr
        self.ivar=ivar

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
            ivar=1./self.yerr**2
            self.ivar=array(ivar, dtype='f8', ndmin=1, copy=False)
        elif self.ivar is not None:
            self.ivar=array(self.ivar, dtype='f8', ndmin=1, copy=False)
        else:
            self.ivar=None

    def _check_data(self):
        wbad,=where(numpy.isfinite(self.x) == False)
        if wbad.size!=0:
            raise ValueError("%d x values are not finite" % wbad.size)
        wbad,=where(numpy.isfinite(self.y) == False)
        if wbad.size!=0:
            raise ValueError("%d y values are not finite" % wbad.size)

        if self.x.size != self.y.size:
            raise ValueError("x,y must be same size")


        if self.ivar is not None:
            wbad,=where(numpy.isfinite(self.ivar) == False)
            if wbad.size!=0:
                raise ValueError("%d ivar values are not finite" % wbad.size)

            if (self.x.size != self.ivar.size) and (self.ivar.size != 1):
                raise ValueError("x,y, and yerr/ivar must be same size")



    def _check_guess(self,guess):
        if guess.size != (self.order+1):
            raise ValueError("guess should be length order+1")
        wbad,=where(numpy.isfinite(guess) == False)
        if wbad.size != 0:
            mess=[]
            mess.append("bad guess:")
            for i in wbad:
                mess.append("  %i %g" % (i,guess[i]))
            mess='\n'.join(mess)
            raise ValueError(mess)

    def _set_guess(self, guess0):
        from esutil.random import srandu

        if guess0 is not None:
            guess0=array(guess0, dtype='f8', ndmin=1, copy=True)
            self._check_guess(guess0)
        else:
            guess0=numpy.polyfit(self.x, self.y, self.order)
            self.guess=guess0

        npars=len(guess0)
        nwalkers=self.nwalkers
        guess=zeros((nwalkers,npars))

        for i in xrange(npars):
            if guess0[i]==0:
                guess[:,i] = guess0[i]+0.1*srandu(nwalkers)
            else:
                guess[:,i] = guess0[i]*(1+0.1*srandu(nwalkers))

        self._guess=guess

    def _run_trials(self):
        import emcee

        guess=self._guess
        npars=guess.shape[1]
        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        npars,
                                        self.get_lnprob,
                                        a=self.a)

        pos, prob, state = sampler.run_mcmc(guess, self.burnin)
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self.nstep)

        self.trials  = sampler.flatchain
        
        self._calc_stats()

    def _calc_stats(self):
        pars,pcov=extract_stats(self.trials)
        perr=sqrt(diag(pcov))

        npars=len(pars)

        lnprob=self.get_lnprob(pars)

        chi2=lnprob/(-0.5)
        dof=self.x.size-npars
        chi2per = chi2/dof

        aic = -2*lnprob + 2*npars
        bic = -2*lnprob + npars*log(self.x.size)

        res={'pars':pars,
             'perr':perr,
             'pcov':pcov,
             'lnprob':lnprob,
             'aic':aic,
             'bic':bic,
             'chi2': chi2,
             'dof':dof,
             'chi2per':chi2per}
        try:
            import scipy.stats
            prob = scipy.stats.chisqprob(chi2, dof)

            res['prob']=prob
        except:
            pass

        self._result=res
        try:
            self._ply=numpy.poly1d(res['pars'])
        except:
            print 'could not set poly'
            self._ply=None


    def get_lnprob(self, pars):
        from numpy import poly1d

        ply=poly1d(pars)
        model=ply(self.x)
        chi2 = (self.y-model)**2
        if self.ivar is not None:
            chi2 *= self.ivar

        return -0.5*chi2.sum()

    def plot_trials(self, **kw):
        plot_results(self.trials, **kw)
        '''
        import biggles
        import esutil as eu

        npars=len(self._result['pars'])

        trials=self.trials

        for i in xrange(npars):
            row=i/ncol
            col=i % ncol

            std=trials[:,i].std()
            binsize=0.2*std

            plti,hi=eu.plotting.bhist(trials[:,i], binsize=binsize, show=False,
                                      gethist=True)
        '''


    def __repr__(self):
        pars=self._result['pars']
        perr=self._result['perr']
        npars=len(pars)

        rep=[]
        header=[]
        for i in xrange(npars):
            h="p%d" % i
            power=npars-i-1
            if power>0:
                h += " x^%d" % power
            header.append(h)

            r= "  p%d: %.4g +/- %.4g" % (i,pars[i],perr[i])
            rep.append(r)
        
        header=' + '.join(header)
        rep = [header] + rep
        rep='\n'.join(rep)
        return rep



class LinFitter:
    def __init__(self, x, y, yerr, err_guess=None):
        self.x=x
        self.y=y
        self.yerr=yerr
        self.ivar = 1.0/yerr**2
        assert(x.size == y.size == yerr.size)

        self.npars = 2

        if err_guess is None:
            yerr_med = numpy.median(yerr)

            # guess at the error for the step() function
            # assign half var to each

            offset_err_guess = 0.5*yerr_med/numpy.sqrt(x.size)

            # pick value at some x, var(y) = x^2*var(slope)
            #  thus var(slope) = var(y)/x^2
            x_std = self.x.std() 
            slope_err_guess = 0.5*yerr_med/x_std

            self.err_guess = numpy.array([offset_err_guess, slope_err_guess], dtype='f8')
        else:
            if len(err_guess) != 2:
                raise ValueError("Error guess must be length 2")
            self.err_guess = err_guess


    def step(self,pars):
        newpars = pars + self.err_guess*randn(self.npars)
        return newpars
    
    def get_loglike(self, pars):
        yfunc = self.line_func(pars)
        chi2 = self.ivar*(self.y-yfunc)**2
        return -0.5*chi2.sum()

    def line_func(self, pars):
        # line is pars[0]*x + self.pars[1]
        return pars[0]*self.x + pars[1]

def noisy_line(pars, xmin, xmax, nx, yerr):

    x = numpy.linspace(xmin, xmax, nx)
    y = pars[0]*x + pars[1]

    y += yerr*randn(nx)
    
    yerr_vals = numpy.array([yerr]*x.size, dtype='f8')

    return x,y,yerr_vals

def gaussfunc(mean,sigma,xvals):

    gnorm = 1.0/numpy.sqrt(2.0*numpy.pi*sigma**2)
    gauss = numpy.exp(-0.5*(xvals - mean)**2/sigma**2 )

    return gauss*gnorm


