import numpy

class ConstFunc:
    """

    We use a minimizer to find the minimum of the function. The class
    representing the data must hold the data internally, and 
    it's evaluation must take "x" values and the pars.

    """

    def __init__(self, value, err, ndata):
        self.data = value*numpy.ones(ndata,dtype='f8')
        self.data += err*numpy.random.standard_normal(ndata)

    def chi2(self, pars):
        chi2 = (self.data-pars[0])**2
        chi2 = chi2.sum()
        return chi2

def fitconst():
    import scipy.optimize
    value = 1.0
    err=0.1

    guess = value + numpy.random.standard_normal(1)[0]

    ndata=50

    cf = ConstFunc(value, err, ndata)

    print scipy.optimize.leastsq(cf.chi2, guess)

class LinFunc:
    def __init__(self, a1, a2, err, ndata):
        self.x = numpy.linspace(-1.0, 1.0, ndata)
        self.y = a1 + a2*x
        self.y += numpy.random.standard_normal(ndata)

    def chi2(self, pars):
        return (self.data-pars[0])**2
