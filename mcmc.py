"""
Module:
    mcmc

Classes:
    MCMC: A class for running Monte Carlo Markov Chains.
    MCMCTester: A class for testing the MCMC class.

testing:
    test: Run the MCMCTester
    testmany:  Run mutliple realizations of the tester and plot a histogram
        of all the means for all realizations.

See the docs for these individual classes for more details.

Revision History:
    Created: 2010-04-02, Erin Sheldon, BNL
"""

import numpy
from numpy import array, zeros, ones, sqrt, arange, isfinite, \
        where, diag, log
from numpy.random import randn
from sys import stdout
import os


def extract_stats(data, weights=None):
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

class MCMC:
    """
    Class:
        MCMC
    Purpose:
        Run a Monte Carlo Markov Chain (MCMC).  The user inputs an object that
        has the methods "step" and "loglike" that can be used to generate the
        chain. Note the loglike method should return the log likelihood

    Calling Sequence:
        m=mcmc.MCMC(obj)
        result = m.run(nstep, par_guess, seed=None)

    Construction:
        Inputs:
            obj: 
                An object to use for evaluating likelihoods and taking steps in
                the chain.  The object must have the methods "step" and
                "loglike".  These methods must have the following
                signatures:
                    newpars = obj.step(pars)
                        Take a new step and return the new parameters.  
                    loglike = obj.loglike(pars)
                        Return the log likelihood of the input parameters.

    See docs for the .run method for more details.

    Revision History:
        Created: 2010-04-02, Erin Sheldon, BNL
    """
    def __init__(self, obj):
         self.init(obj)

    def init(self, obj):
        self.obj = obj
        self.fobj = None
        self.parguess=None

    def run(self, nstep, parguess, seed=None):
        """
        Class:
            MCMC
        Method Name:
            run
        Purpose:
            Run the MCMC chain.
        Calling Sequence:
            m=mcmc.MCMC(obj)

            # run the chain
            chain_data = m.run(nstep, par_guess, seed=None)

        Inputs:
            nstep: Number of steps in the chain.
            parguess:  Starting point for the chain in the n-dimensional
                parameters space.

        Optional Inputs:
            seed: A seed for the random number generator. If not given, 
                a seed is automatically generated by numpy.random.seed()

        Outputs:
            A rec array holding the chain data.  The data type is

                [('accepted','i1'),('pars',('f8',npar)), ('loglike','f8')]

            So output['pars'][i] is the parameters in step i, and
            output['loglike'][i] is the log likelihood of step i.
            Accepted is 1 if the step was accepted 0 otherwise.

            For example:
                data = m.run(1000)
                data['pars']
                data['loglike']

        """

        self.parguess=array(parguess,dtype='f8')

        # If seed sent use it, else just generate one
        numpy.random.seed(seed)


        self.oldpars = self.parguess.copy()
        self.oldlike = self.obj.loglike(self.oldpars)
        self.npar = self.oldpars.size

        output = self.result_struct(nstep)

        for i in xrange(nstep):
            self.step()

            output['pars'][i] = self.newpars
            output['loglike'][i] = self.newlike
            output['accepted'][i] = self.accepted

            self.oldpars = self.newpars
            self.oldlike = self.newlike

        return output

    def result_struct(self, num):
        dtype = self.result_dtype(self.npar)
        st = zeros(num, dtype=dtype)
        return st

    def result_dtype(self, npar):
        return [('accepted','i1'),('pars',('f8',npar)), ('loglike','f8')]

    def step(self):
        """
        Class:
            MCMC
        Method Name:
            step
        Purpose:
            Take the next step in the MCMC chain.  Calls the .step and
            .loglike methods of the object send during construction.  If
            the new loglike is not greater than the previous, or a
            uniformly generated random number is greater than the the ratio of
            new to old likelihoods, the new step is not used, and the new
            parameters are the same as the old.  Otherwise the new step is
            kept.

            This is an internal function that is called by the .run method.
            It is not intended for call by the user.
        """
        # Take a step and evaluate the likelihood
        newpars = self.obj.step(self.oldpars)
        newlike = self.obj.loglike(newpars)

        likeratio = newlike-self.oldlike

        randnum = numpy.random.random()
        randnum = numpy.log(randnum)

        # we allow use of -infinity as a sign we are out of bounds
        if (isfinite(newlike) 
                and ( (newlike > self.oldlike) | (randnum < likeratio)) ):
            self.newpars=newpars
            self.newlike=newlike
            self.accepted=1
        else:
            self.newpars=self.oldpars
            self.newlike=self.oldlike
            self.accepted=0







class MCMCTester:
    """

    A simple class for testing the MCMC code.  This is a good template for how
    to begin your own MCMC project.

        import mcmc
        tc = mcmc.MCMCTester(type="constant")
        tc.RunTest()

    You can then use the PlotResults method to see a histogram of the results.

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
            val,sigma = self.truepars()

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
            self.psigma = self.yerr[0]

            #self.parguess=val + 3*sigma*randn()
            self.parguess=0.

            if verbose:
                stdout.write('  type: "constant"\n')
                stdout.write("  true_pars: %s\n  npars: %s\n  psigma: %s\n" %
                             (self.true_pars,self.npars,self.psigma))


    def truepars(self):
        if self.type is "constant":
            val = 10.0
            sigma = 1.0
            return val, sigma
        else:
            raise ValueError("only support type='constant'")



    # pars must be an array
    def step(self,pars):
        return pars + self.psigma*randn(self.npars)

    def loglike(self,pars):
        if self.type == 'constant':
            chi2 = self.ivar*(self.y-pars)**2
            chi2 = -0.5*chi2.sum()
        else:
            raise ValueError("only support type='constant'")
        return chi2

    def RunTest(self, nstep=10000):

        if self.verbose:
            stdout.write("  nstep: %s\n" % nstep)


        m=MCMC(self)
        self.trials = m.run(nstep, self.parguess)

    def RunMultipleTests(self, ntrial, nstep=1000, burnin=100):
        """
        Run multiple tests and save in self.meanvals
        Note default is fewer steps per mcmc
        """
        self.meanvals = zeros(ntrial, dtype='f4')

        for i in xrange(ntrial):
            # create a new realization
            self.setup()
            self.RunTest(nstep)
            self.meanvals[i] = self.trials['pars'][burnin:].mean()



    def CompareResults(self, burnin=100, doplot=False, hardcopy=False):
        if self.type == 'constant':
            par0 = self.trials['pars'][burnin:]

            mn = par0.mean()
            err = par0.std()

            stdout.write("Using burn in of: %s\n" % burnin)
            stdout.write("True value: %s\n" % self.true_pars)
            stdout.write("True error on each of %s points: %s\n" % (self.npoints,self.true_error))
            expected_error = self.true_error/numpy.sqrt(self.npoints)
            stdout.write("Expected error from combined points: %s\n" % expected_error)
            stdout.write("  Mean from trials: %s\n" % mn)
            stdout.write("  Error from trials: %s\n" % err)

            if doplot:
                self.PlotResults(burnin, hardcopy=hardcopy)


        else:
            raise ValueError("only support type='constant'")

    def PlotMultipleResults(self, hardcopy=False):
        """
        Plot the results from RunMultipleTests, saved in meanvals.
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



    def PlotResults(self, burnin=100, hardcopy=False):
        import biggles
        from esutil import stat
        if self.type == 'constant':
            par0 = self.trials['pars'][burnin:]

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



    def get_results(self):
        return self.trials, self.like



def plot_results6(res, burnin, sigma_clip=True):
    import biggles
    import esutil
    tab = biggles.Table(2,3)
    means,cov= extract_stats(res, burnin)
    errs = sqrt(diag(cov)) 
    for i in xrange(6):
        irow = i//3
        icol = i % 3
        binsize=errs[i]*0.2
        min = means[i]-4.0*errs[i]
        max = means[i]+4.0*errs[i]
        hdict = esutil.stat.histogram(res['pars'][burnin:, i], 
                                      binsize=binsize, 
                                      min=min,max=max,
                                      more=True)

        hplot = biggles.Histogram(hdict['hist'], x0=hdict['low'][0], binsize=binsize)
        plt=biggles.FramedPlot()
        plt.add(hplot)
        plt.xlabel = 'parameter %i' % i

        tab[irow,icol] = plt

    tab.show()


def plot_results(trials, binfac=0.2, names=None, show=True):
    import biggles
    import esutil

    biggles.configure( 'default', 'fontsize_min', 1)

    means,cov = extract_stats(trials)
    errs=sqrt(diag(cov)) 

    npars=len(means)

    nrow,ncol=get_grid(npars)
    plt=biggles.Table(nrow,ncol)
    plt.aspect_ratio=float(nrow)/ncol

    for i in xrange(npars):
        row=i/ncol
        col=i % ncol

        bsize = binfac*errs[i]

        hdict = esutil.stat.histogram(trials[:, i], 
                                      binsize=bsize, 
                                      more=True)
        hplot = biggles.Histogram(hdict['hist'], 
                                  x0=hdict['low'][0], 
                                  binsize=bsize)
        plti=biggles.FramedPlot()
        if names is not None:
            name=names[i]
        else:
            name=r'$p_{%d}$' % i

        plti.xlabel=name
        plti.add(hplot)
            
        lab = r'$<%s> = %0.4g \pm %0.4g$' % (name,means[i],errs[i])
        plab = biggles.PlotLabel(0.1,0.9,lab,halign='left')

        plti.add(plab)

        plt[row,col]=plti
    
    if show:
        plt.show()

    return plt

def test(burnin=100, nstep=10000, doplot=False, hardcopy=False):
    """
    Name:
        test
    Purpose:
        Run the MCMCTester.
    """
    tc = MCMCTester("constant")
    tc.RunTest(nstep)
    plot_burnin(tc.trials['pars'], tc.trials['loglike'], burnin)
    tc.CompareResults(doplot=doplot, hardcopy=hardcopy)

def testmany(ntrial):
    """
    Just do a bunch of realizations and make sure we get the right mean
    """

    tc = MCMCTester("constant")
    tc.RunMultipleTests(ntrial)
    tc.PlotMultipleResults()

def test_line(nstep=10000, doplot=False):
    import esutil
    pars = [1.0,1.0]
    xmin = -1.0
    xmax =  1.0
    nx = 10
    yerr = 0.1
    x, y, yerr = noisy_line(pars, xmin, xmax, nx, yerr)

    LF = LinFitter(x, y, yerr)

    chain = MCMC(LF)

    parguess = [ pars[0] + 0.1, pars[1]-0.1 ]
    data = chain.run(nstep, parguess)


    if doplot:
        import biggles
        from esutil import stat

        burnin = 1000

        # plot the burnin
        tab = biggles.Table(2,1)

        ntoplot = burnin*10
        if ntoplot > data.size:
            ntoplot = data.size
        steps = numpy.arange(ntoplot, dtype='i4')

        offset_burnin = biggles.FramedPlot()
        offset_burnin.ylabel = 'offset'
        offset_burnin.add( biggles.Curve(steps, data['pars'][0:ntoplot,0] ) )

        slope_burnin = biggles.FramedPlot()
        slope_burnin.ylabel = 'slope'
        slope_burnin.xlabel = 'step number'
        slope_burnin.add( biggles.Curve(steps, data['pars'][0:ntoplot,1] ))

        tab[0,0] = offset_burnin
        tab[1,0] = slope_burnin

        tab.show()

        # get status for chain
        parfit, cov = extract_stats(data, burnin)
        errfit = sqrt(diag(cov))

        # plot the histograms and comparison plot

        tab = biggles.Table(2,2)
        offsets = data['pars'][burnin:,0]
        slopes  = data['pars'][burnin:,1]

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


class PolyFitter(object):
    def __init__(self, order, x, y, nwalkers, burnin, nstep, 
                 guess=None, a=3, yerr=None, ivar=None):

        self.order=order
        self.x=array(x, dtype='f8', ndmin=1, copy=False)
        self.y=array(y, dtype='f8', ndmin=1, copy=False)

        self.nwalkers=nwalkers
        self.burnin=burnin
        self.nstep=nstep
        self.a=a

        self.guess=guess

        if yerr is not None:
            ivar=1./yerr**2
            self.ivar=array(ivar, dtype='f8', ndmin=1, copy=False)
        elif ivar is not None:
            self.ivar=array(ivar, dtype='f8', ndmin=1, copy=False)
        else:
            self.ivar=ones(self.x.size, dtype='f8')
        
        if ((self.x.size != self.y.size)
                or
                (self.x.size != self.ivar.size)):
            raise ValueError("x,y, and yerr/ivar must be same size")

        self._go()

    def get_result(self):
        return self._result

    def get_poly(self):
        return numpy.poly1d(self._result['pars'])

    def _get_guess(self):
        from esutil.random import srandu

        if self.guess is not None:
            guess0=array(self.guess, dtype='f8', ndmin=1, copy=True)
            if guess0.size != (self.order+1):
                raise ValueError("guess should be length order+1")
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

        return guess

    def _go(self):
        import emcee
        guess=self._get_guess()
        npars=guess.shape[1]
        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        npars,
                                        self.loglike,
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

        loglike=self.loglike(pars)

        chi2=loglike/(-0.5)
        dof=self.x.size-npars
        chi2per = chi2/dof

        aic = -2*loglike + 2*npars
        bic = -2*loglike + npars*log(self.x.size)

        res={'pars':pars,
             'perr':perr,
             'pcov':pcov,
             'loglike':loglike,
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


    def loglike(self, pars):
        from numpy import poly1d

        ply=poly1d(pars)
        model=ply(self.x)
        chi2 = self.ivar*(self.y-model)**2
        return -0.5*chi2.sum()

    def plot_trials(self):
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
    
    def loglike(self, pars):
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


