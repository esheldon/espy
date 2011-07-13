import numpy
import esutil
from scikits.learn import mixture

class GaussMix(mixture.GMM):
    """

    Inherit from the scikits gaussian mixture model, and implement some
    additional functionality

    Simpler calls for 1-d data.

    """

    def fit1d(self, data1d, **keywords):
        self.fit(data1d.reshape(data1d.size, 1), **keywords)
    def score1d(self, xvals, state=None):
        import numpy as np
        from numpy import log, pi as PI
        if state is not None:
            if self._cvtype != 'diag':
                raise ValueError("not sure this 1d works for non diag,check")
            obs = numpy.asanyarray(xvals)
            mean = self.means[state,0]
            covar = self.covars[state][0][0]
            log_weight = self._log_weights[state]

            lpr = -0.5*( (obs-mean)**2/covar + log(2*PI*covar) ) + log_weight
            return lpr
        else:
            return self.score( xvals.reshape(xvals.size, 1) )

    def compare1d(self, data1d, binsize, each=False):
        import biggles
        import pcolors

        colors = pcolors.rainbow(self.n_states,'hex')

        plt=biggles.FramedPlot()

        hdict = esutil.stat.histogram(data1d, binsize=binsize, more=True)

        h = biggles.Histogram(hdict['hist'], x0=hdict['low'][0], binsize=binsize)
        renorm = data1d.size*(hdict['high'][0] - hdict['low'][0])

        plt.add(h)

        # now the models
        logprob = self.score1d(hdict['center'])
        prob = renorm*numpy.exp(logprob)
        pprob = biggles.Curve(hdict['center'], prob, type='shortdashed')
        plt.add(pprob)

        for i in xrange(self.n_states):
            logprob = self.score1d(hdict['center'], state=i)
            prob = renorm*numpy.exp(logprob)
            pprob = biggles.Curve(hdict['center'], prob, color=colors[i])
            plt.add(pprob)


        plt.show()

def test_compare():
    n1=4000
    n2=6000
    m1 = 1.0
    sig1 = 1.0
    m2 = 2.0
    sig2 = 0.5
    data1 = m1 + sig1*numpy.random.randn(n1)
    data2 = m2 + sig2*numpy.random.randn(n2)

    data = numpy.zeros(n1+n2)
    data[0:n1] = data1
    data[n1:n1+n2] = data2

    gm = GaussMix(n_states=2)
    gm.fit1d(data)

    gm.compare1d(data, 0.05)

