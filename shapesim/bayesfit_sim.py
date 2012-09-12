from numpy import sqrt, cos, sin, exp, pi, zeros, random, where, array
class EPrior:
    """
    This is in g1,g2 space

    Prob = A cos(|e| pi/2) exp( - [ 2 |e| / B / (1 + |e|^D) ]^C )
    """
    def __init__(self):
        # A actually depends on norm when doing full thing
        self.A = 1
        self.B = 0.03
        self.C = 0.45
        self.D = 13.

        self.maxval = self.prior(0., 0.)

    def prior(self, g1, g2):
        """
        Prior actually only depends on the total mag of g
        """
        g1 = array(g1, ndmin=1, copy=False)
        g2 = array(g2, ndmin=1, copy=False)
        g = sqrt(g1**2 + g2**2)

        prior = zeros(g1.size)

        w,=where(g < 1)
        if w.size > 0:
            prior[w] = self.A * cos(g[w]*pi/2)*exp( - ( 2*g[w] / self.B / (1 + g[w]**self.D) )**self.C )
        return prior

    def sample(self, nrand):
        g1 = zeros(nrand)
        g2 = zeros(nrand)

        ngood=0
        nleft=nrand
        while ngood < nrand:

            # generate total g**2 in [0,1)
            grand2 = random.random(nleft)
            grand = sqrt(grand2)
            # now uniform angles
            rangle = random.random(nleft)*2*pi

            # now get cartesion locations in g1,g2 plane
            g1rand = grand*cos(rangle)
            g2rand = grand*sin(rangle)

            # now finally the height from [0,maxval)
            h = self.maxval*random.random(nleft)

            pvals = self.prior(g1rand, g2rand)

            w,=where(h < pvals)
            if w.size > 0:
                g1[ngood:ngood+w.size] = g1rand[w]
                g2[ngood:ngood+w.size] = g2rand[w]
                ngood += w.size
                nleft -= w.size
        return g1, g2
