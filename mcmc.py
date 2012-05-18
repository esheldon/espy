"""
Module:
    mcmc

Classes:
    MCMC: A class for running Monte Carlo Markov Chains.
    MCMCTester: A class for testing the MCMC class.

functions:
    read_results(filename):
        Read a result structure output by the mcmc code.  

testing:
    test: Run the MCMCTester
    testmany:  Run mutliple realizations of the tester and plot a histogram
        of all the means for all realizations.

See the docs for these individual classes for more details.

Revision History:
    Created: 2010-04-02, Erin Sheldon, BNL
"""
import numpy
from sys import stdout
import os

def read_results(filename):
    """
    Name:
        read_results
    Purpose:
        Read a result structure output by the mcmc code.  

    Description:

        The format of the file is binary in the native endian of the machine.
        The first 4 bytes are a 4-byte integer holding the number of
        parameters:

            4 byte integer                number of parameters

        then a large block of bytes representing each step in the chain. Each
        step is represented by npar 8-byte floats for each parameter, followed
        by an 8-byte float for the likelihood of that step

            npar*8 + 8 bytes              number of steps for each trial

        Thus the number of steps can be computed from the size of the file:
            
            nsteps = (nbytes-4)/(npar*8 + 8)

        The result is an array with fields, with data type:

            [('pars',('f8',npar)), ('like','f8')]


    """
    stat = os.stat(filename)
    nbytes = stat.st_size

    fobj = open(filename, 'r')

    npar = numpy.zeros(1, dtype='i4')
    npar = numpy.fromfile(fobj, dtype='i4', count=1)

    # each row is par0,par1,par2,...,parN,like
    # where the par and like are 8-byte floats
    # thus the number of trials is nbytes/(npar*8 + 8)
    # don't forget to remove 4 bytes for npar at 
    # the beginning

    nstep = (nbytes-4)/(npar*8 + 8)

    dtype = result_dtype(npar)
    data = numpy.fromfile(fobj, dtype=dtype, count=nstep)

    fobj.close()
    return data

def result_dtype(npar):
    return [('pars',('f8',npar)), ('like','f8')]

def extract_stats(data, burnin, sigma_clip=True):
    import esutil
    npar = data['pars'].shape[1]

    means = numpy.zeros(npar,dtype='f8')
    errs  = numpy.zeros(npar,dtype='f8')

    for i in xrange(npar):
        if not sigma_clip:
            means[i] = data['pars'][burnin:, i].mean()
            errs[i] = data['pars'][burnin:, i].std()
        else:
            means[i], errs[i] = esutil.stat.sigma_clip(data['pars'][burnin:, i])

    return means, errs

class MCMC:
    """
    Class:
        MCMC
    Purpose:
        Run a Monte Carlo Markov Chain (MCMC).  The user inputs an object that
        has the methods "step" and "likelihood" that can be used to generate
        the chain

    Calling Sequence:
        m=mcmc.MCMC(obj, log=True)
        result = m.run(nstep, par_guess, seed=None)

    Construction:
        Inputs:
            obj: 
                An object to use for evaluating likelihoods and taking steps in
                the chain.  The object must have the methods "step" and
                "likelihood".  These methods must have the following
                signatures:
                    .newpars = obj.step(pars)
                        Take a new step and return the new parameters.  
                    .likelihood = obj.likelihood(pars)
                        Return the likelihood of the input parameters.

        Optional Inputs:
            log:  If True, the object returns log likelihoods.  Default
                is True.

    See docs for the .run method for more details.

    Revision History:
        Created: 2010-04-02, Erin Sheldon, BNL
    """
    def __init__(self, obj, log=True):
         self.init(obj, log=log)

    def init(self, obj, log=True):
        self.obj = obj
        self.log = log
        self.fobj = None
        self.parguess=None

    def open_output(self, filename):
        if self.parguess is None:
            raise ValueError("parguess must be set before opening")

        self.fobj = open(filename, 'w')
        try:
            npar = len(self.parguess)
        except:
            npar = 1
        npar = numpy.array(npar, dtype='i4')
        npar.tofile( self.fobj )
 
        
    def close_output(self):
        if isinstance(self.fobj, file):
            self.fobj.close()
            self.fobj=None

    def run(self, nstep, parguess, seed=None, file=None):
        """
        Class:
            MCMC
        Method Name:
            run
        Purpose:
            Run the MCMC chain.
        Calling Sequence:
            m=mcmc.MCMC(obj, log=True)

            # run the chain
            chain_data = m.run(nstep, par_guess, seed=None)

        Inputs:
            nstep: Number of steps in the chain.
            parguess:  Starting point for the chain in the n-dimensional
                parameters space.

        Optional Inputs:
            seed: A seed for the random number generator. If not given, 
                a seed is automatically generated by numpy.random.seed()

            file: 
                a file name to write the results.  The file will contain a
                4-byte integer for the number of parameters, followed by the
                data, 8 byte floats for each parameter and and 8-byte float for
                the likelihood.  Repeated for each step.  See the read_results
                function for more info.

                if file is sent, the return value of this method is None

        Outputs:
            A rec array holding the chain data.  The data type is
                [('pars',('f8',npar)), ('like','f8')]
            So output['pars'][i] is the parameters in step i, and
            output['likelihood'][i] is the likelihood of step i.

            For example:
                data = m.run(1000)
                data['pars']
                data['likelihood']

            If file is sent, then None is returned.


        """

        # If seed sent use it, else just generate one
        numpy.random.seed(seed)

        if not numpy.isscalar(parguess):
            self.parguess = numpy.array(parguess, dtype='f8')
            self.npar = self.parguess.size
        else:
            self.parguess = parguess
            self.npar=1

        # these are the results
        if file is not None:
            self.open_output(file)
            output = None
        else:
            output = self.result_struct(nstep)

        self.oldpars = self.parguess.copy()
        self.oldlike = self.obj.likelihood(self.oldpars)

        for i in xrange(nstep):
            self.step()

            if file is not None:
                self.write_step()
            else:
                output['pars'][i] = self.newpars
                output['like'][i] = self.newlike

            self.oldpars = self.newpars
            self.oldlike = self.newlike

        if file is not None:
            self.close_output()
        # when writing to a file, will be None
        return output

    def result_struct(self, num):
        dtype = result_dtype(self.npar)
        st = numpy.zeros(num, dtype=dtype)
        return st

    def write_step(self):
        parout = numpy.array(self.newpars, dtype='f8', copy=False)
        parout.tofile(self.fobj)
        likeout = numpy.array(self.newlike, dtype='f8')
        likeout.tofile(self.fobj)



    def get_results(self):
        """
        Class:
            MCMC
        Method Name:
            get_results
        Purpose:
            return the parameters and likelihood for each trial
        Calling Sequence:
            mcmc=mcmc.MCMC(obj, log=True)
            mcmc.run(nstep, par_guess, seed=None)
            trials, liklihoods = mcmc.get_results()
        Outputs:
            trials, likelihoods:  A tuple containing trials and likelihoods.
                trials: 
                    an (npars, nstep) array containing the parameters at each
                    step in the chain.
                likelihoods: 
                    is an (nstep) length array containing the likelihood at
                    each step in the chain.
        """
        return self.pars, self.like



    def step(self):
        """
        Class:
            MCMC
        Method Name:
            step
        Purpose:
            Take the next step in the MCMC chain.  Calls the .step and
            .likelihood methods of the object send during construction.  If
            the new likelihood is not greater than the previous, or a
            uniformly generated random number is greater than the the ratio of
            new to old likelihoods, the new step is not used, and the new
            parameters are the same as the old.  Otherwise the new step is
            kept.

            This is an internal function that is called by the .run method.
            It is not intended for call by the user.
        """
        # Take a step and evaluate the likelihood
        newpars = self.obj.step(self.oldpars)
        newlike = self.obj.likelihood(newpars)

        if self.log:
            likeratio = newlike-self.oldlike
        else:
            likeratio = newlike/self.oldlike

        randnum = numpy.random.random()
        if self.log:
            randnum = numpy.log(randnum)

        if (newlike > self.oldlike) | (randnum < likeratio):
            self.newpars=newpars
            self.newlike=newlike
        else:
            self.newpars=self.oldpars
            self.newlike=self.oldlike







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
            self.y = numpy.zeros(ny, dtype='f8')
            self.y[:] = val
            self.y[:] += sigma*numpy.random.standard_normal(ny)

            self.yerr = numpy.zeros(ny, dtype='f8')
            self.yerr[:] = sigma

            self.ivar = 1.0/self.yerr**2
            self.psigma = self.yerr[0]

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
        return pars + self.psigma*numpy.random.standard_normal(self.npars)

    def likelihood(self,pars):
        if self.type == 'constant':
            chi2 = self.ivar*(self.y-pars)**2
            chi2 = -0.5*chi2.sum()
        else:
            raise ValueError("only support type='constant'")
        return chi2

    def RunTest(self, nstep=10000, file=None):

        if self.verbose:
            stdout.write("  nstep: %s\n" % nstep)

        parguess = self.true_pars

        m=MCMC(self)
        self.trials = m.run(nstep, parguess, file=file)

    def RunMultipleTests(self, ntrial, nstep=1000, burnin=100):
        """
        Run multiple tests and save in self.meanvals
        Note default is fewer steps per mcmc
        """
        self.meanvals = numpy.zeros(ntrial, dtype='f4')

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
    means,errs = extract_stats(res, burnin, sigma_clip=sigma_clip)
    
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


def test(nstep=10000, doplot=False, hardcopy=False):
    """
    Name:
        test
    Purpose:
        Run the MCMCTester.
    """
    tc = MCMCTester("constant")
    tc.RunTest(nstep)
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
        parfit, errfit = extract_stats(data, burnin)


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
        newpars = pars + self.err_guess*numpy.random.standard_normal(self.npars)
        return newpars
    
    def likelihood(self, pars):
        yfunc = self.line_func(pars)
        chi2 = self.ivar*(self.y-yfunc)**2
        return -0.5*chi2.sum()

    def line_func(self, pars):
        # line is pars[0]*x + self.pars[1]
        return pars[0]*self.x + pars[1]

def noisy_line(pars, xmin, xmax, nx, yerr):

    x = numpy.linspace(xmin, xmax, nx)
    y = pars[0]*x + pars[1]

    y += yerr*numpy.random.standard_normal(nx)
    
    yerr_vals = numpy.array([yerr]*x.size, dtype='f8')

    return x,y,yerr_vals
