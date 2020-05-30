import numpy
from numpy import log, arctan, arctanh, sqrt, log10, pi as PI

import esutil as eu
from esutil.numpy_util import where1, arrscl

class NFW:
    """

    Evaluate DeltaSigma or rho for an nfw.

    This has fixed omega_m and z, entered on construction.

    """
    def __init__(self, omega_m, z):
        self.omega_m = omega_m
        self.z = z

        # this is E(z) for flat universe
        Ez = omega_m*(1+z)**3 + (1.0-omega_m)

        self.Ez=Ez
        self.rhocrit = 0.277520*Ez
        self.rho0 = self.omega_m*self.rhocrit

        # buffer region around rs, linearly interpolate through
        # region near rs
        self.ep=0.001

    def rho(self, r, r200, c):
        del_c = (200/3.0)*c**3 /(log(1.0+c)-c/(1.0+c))
        rs = r200/c
        x = r/rs

        rho = del_c*self.rhocrit/x/(1+x)**2
        return rho

    def plot_rho(self, r200, c):
        from biggles import FramedPlot, Curve
        n=1000
        r = numpy.linspace(0.01, 20.0,n)
        rho = self.rho(r, r200, c)

        plt=FramedPlot()
        plt.add( Curve(r,rho) )
        plt.xlog=True
        plt.ylog=True

        plt.show()

    def m(self, r, r200, c):
        """
        Mass less than radius r in solar masses
        r and r200 in Mpc.
        """
        del_c = (200/3.0)*c**3 /(log(1.0+c)-c/(1.0+c))
        rs = r200/c
        x = r/rs

        m = 4*PI*del_c*self.rhocrit*rs**3*( log(1+x) - x/(1+x) )
        return m*1.e12

    def rhomean(self, r, r200, c):
        """
        Mean density within radius r for the given NFW parameters
        """
        v = 4.0/3.0*PI*r**3
        return self.m(r, r200, c)/v

    def m200(self, r200):
        """

        Gives mass in solar masses for r200 in Mpc

        Same as puttin r=r200 in the .m() method
        Note independent of c

        """
        m200 = 200*self.rhocrit*(4.0/3.0)*PI*r200**3
        return m200*1.e12


    def plot_m(self, r200, c):
        from biggles import FramedPlot, Curve
        n=1000
        r = numpy.linspace(0.01, 20.0,n)
        m = self.m(r, r200, c)

        plt=FramedPlot()
        plt.add( Curve(r,m) )
        plt.xlog=True
        plt.ylog=True

        plt.show()



    def dsig(self, r, r200, c):
        """

        The results here are exactly the same as Dave's
        r in Mpc
        pars[0] = r200 in Mpc
        pars[1] = concentration

        """

        ep = self.ep
        omep = 1-ep
        opep = 1+ep

        rs = r200/c
        xx = r/rs # = c*r/r200
        del_c = (200/3.0)*c**3 /(log(1.0+c)-c/(1.0+c))

        fac=rs*del_c*self.rhocrit

        w1=where1( xx < (1-ep) )
        w2=where1( (xx >= (1-ep)) & (xx <= (1+ep)) )
        w3=where1( xx > (1+ep))

        dsig=numpy.zeros(r.size, dtype='f8')

        if w1.size > 0:
            x=xx[w1]
            x2=x**2
            A=arctanh(sqrt((1-x)/(1+x)))
            s=8*A/(x2*sqrt(1-x2)) + 4*log(x/2)/x2 -2/(x2-1) + 4*A/((x2-1)*sqrt(1-x2))
            
            dsig[w1] = s

        # interpolate between the two regions
        if w2.size > 0:

            e = 1-ep
            e2 = e**2
            A=arctanh(sqrt((1-e)/(1+e)))
            s1=8*A/(e2*sqrt(1-e2)) + 4*log(e/2)/e2 -2/(e2-1) + 4*A/((e2-1)*sqrt(1-e2))

            e  = 1+ep
            e2 = e**2
            A=arctan(sqrt((e-1)/(1+e)))
            s2=8*A/(e2*sqrt(e2-1)) + 4*log(e/2)/(e2) -2/(e2-1) + 4*A/((e2-1)**(3/2.0))

            
            e1=1-ep
            e2=1+ep
            x = xx[w2]
            s = (x-e1)*s2/(e2-e1) + (x-e2)*s1/(e1-e2)

            dsig[w2] = s

        if w3.size > 0:
            x=xx[w3]
            x2=x**2

            A=arctan(sqrt((x-1)/(1+x)))
            s=8*A/(x2*sqrt(x2-1)) + 4*log(x/2)/x2 -2/(x2-1) + 4*A/((x2-1)**(3/2.0))
            dsig[w3] = s
        dsig *= fac

        return dsig

    def plot_dsig(self, r200, c):
        from biggles import FramedPlot, Curve
        n=1000
        r = numpy.linspace(0.01, 20.0,n)
        ds = self.dsig(r, r200, c)

        plt=FramedPlot()
        plt.add( Curve(r,ds) )
        plt.xlog=True
        plt.ylog=True

        plt.show()

    def r_fmean(self, r200, c, fac):
        """
        Define the mean density within radius r as
            rhor = M(r)/V(r)

        return the radius where rhor = rho0*fac
        where rho0 = omega_m*rhocrit is the mean density

        Use a root finder to get the value of r where the equation
            rhomean - fac*rho0 = 0 
            
        """
        import scipy.optimize
        guess = r200
        r_fmean = scipy.optimize.fsolve(self.r_fmean_kernel,
                                        guess,
                                        args=(r200,c,fac))
        return r_fmean[0]

    def r_fmean_kernel(self, r, r200, c, fac):
        """
        Used for solving for mmean

        rhor = M(r)/V(r) the mean density within radius r
        for the NFW profile.
        """

        return self.rhomean(r, r200,  c) - fac*self.rho0

    def r_fcrit(self, r200, c, fac):
        """
        Define the mean density within radius r as
            rhor = M(r)/V(r)

        return the radius where rhomean = rhocrit*fac. This
        this should give us back r200 for fac=200

        Use a root finder to get the value of r where the equation
            rhomean - fac*rhocrit = 0 
            
        """
        import scipy.optimize
        guess = r200
        r_fmean = scipy.optimize.fsolve(self.r_fcrit_kernel,
                                        guess,
                                        args=(r200,c,fac))
        return r_fmean[0]

    def r_fcrit_kernel(self, r, r200, c, fac):
        """
        Used for solving for mmean

        rhor = M(r)/V(r) the mean density within radius r
        for the NFW profile.
        """

        return self.rhomean(r, r200,  c) - fac*self.rhocrit





