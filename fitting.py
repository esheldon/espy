import numpy
import scipy.optimize

class LineFitter:
    def __init__(self, x, y, yerr=None):
        self.x=x
        self.y=y
        self.yerr=yerr
        self.set_guess()
        self.dofit()

    def set_guess(self):
        self.guess=numpy.polyfit(self.x, self.y, 1)

    def eval_pars(self, pars):
        return pars[0]*self.x + pars[1]

    def errfunc(self, pars):
        model = self.eval_pars(pars)
        if self.yerr is None:
            diff = model-self.y
        else:
            diff = (model-self.y)/self.yerr

        return diff

    def dofit(self):
        res=scipy.optimize.leastsq(self.errfunc, self.guess,
                                       full_output=1)
        self.pars, self.pcov0, self.infodict, self.errmsg, self.ier = res
        if self.ier == 0:
            # wrong args, this is a bug
            raise ValueError(self.errmsg)

        self.numiter = self.infodict['nfev']
        self.pcov=None
        self.perr=None

        if self.pcov0 is not None:
            self.pcov = self.scale_leastsq_cov(self.pars, self.pcov0)

            d=numpy.diag(self.pcov)
            w,=numpy.where(d < 0)

            if w.size == 0:
                # only do if non negative
                self.perr = numpy.sqrt(d)

    def get_result(self):
        return {'pars':self.pars,
                'pcov':self.pcov,
                'perr':self.perr}
    def scale_leastsq_cov(self, pars, pcov):
        """
        Scale the covariance matrix returned from leastsq; this will
        recover the covariance of the parameters in the right units.
        """
        dof = (self.x.size-len(pars))
        s_sq = (self.errfunc(pars)**2).sum()/dof
        return pcov * s_sq 


    def get_poly(self):
        return numpy.poly1d(self.pars)

    def __call__(self, x):
        """
        pars order same as for numpy.poly1d
        """
        return self.pars[0]*self.x + self.pars[1]

    def __repr__(self):
        if self.perr is not None:
            return """y = p0*x + p1
p0: %s +/- %s
p1: %s +/- %s""" % (self.pars[0],self.perr[0],self.pars[1],self.perr[1])
        else:
            return """y = p0*x + p1
p0: %s +/- None
p1: %s +/- None""" % (self.pars[0],self.pars[1])

def test_line(show=False):
    pars = [1,3]
    x = numpy.arange(20)
    y = pars[0]*x + pars[1]
    yerr = numpy.zeros(len(x)) + 2 + 0.1*numpy.random.random(len(x))
    y += yerr*numpy.random.randn(len(x))

    lf = LineFitter(x, y, yerr=yerr)
    print 'guess:',lf.guess
    print lf

    if show:
        import biggles
        plt=biggles.FramedPlot()
        plt.add(biggles.Points(x, y,type='filled circle'))
        plt.add(biggles.SymmetricErrorBarsY(x, y, yerr))
        plt.add(biggles.Curve(x, lf(x),color='blue'))
        plt.show()
