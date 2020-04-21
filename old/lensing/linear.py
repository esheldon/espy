import numpy
from numpy import log, log10, sqrt, exp, sin, cos,\
        where, ceil, isfinite, \
        int32, roll, pi as PI, linspace, arange, interp
import esutil as eu
from esutil.numpy_util import where1
from esutil.stat.util import interplin

from lensing import project

try:
    import biggles
    from biggles import FramedPlot,Points, Table, Curve, PlotKey
except:
    pass

class Linear:
    """
    Compute P(k), xi(r), delta_sigma(r)

        ported from Dave Johnston's IDL code.
    """
    def __init__(self, **keys):
        self.omega_m=keys.get('omega_m',0.25)
        self.omega_b=keys.get('omega_b',0.055)
        self.sigma_8=keys.get('sigma_8',0.8)
        self.h=keys.get('h',1.0)
        self.ns=keys.get('ns',0.98)

        self.ehu = EisenHu(**keys)

        self.qg1000 = eu.integrate.QGauss(npts=1000)

        # rhocrit at z=0
        self.rhocrit0 = 0.277520

        try:
            import biggles
            biggles.configure('screen','width', 800)
            biggles.configure('screen','height', 800)
        except:
            pass

    def interpsp2(self, y, x, xinterp):
        """
        the x,y must be sorted by x
        """
        from scipy import interpolate


        m = x.size
        s = x.searchsorted(xinterp) -1
        s = numpy.where(s < 1, 1, s)
        s = numpy.where(s > (m-3), (m-3), s)

        nx = xinterp.size
        yinterp = numpy.zeros(nx,dtype='f8')

        sold = -1
        for i in xrange(nx):
            s0 = s[i]-1
            if sold != s0:
                x0 = x[s0: s0+4]
                y0 = y[s0: s0+4]
                tck = interpolate.splrep(x0,y0,k=3,s=0)
                sold = s0
            yinterp[i] = interpolate.splev(xinterp[i],tck,der=0)

        return yinterp


    def interpsp1(self, y, x, xinterp):
        """
        the x,y must be sorted by x
        """
        from scipy import interpolate

        tck = interpolate.splrep(x,y,k=3,s=0)
        yinterp = interpolate.splev(xinterp,tck,der=0)
        return yinterp

        # cubic (k=3) spline representation
        # no smoothing (s=0)
        nx = xinterp.size
        yinterp = numpy.zeros(nx,dtype='f8')
        for i in xrange(xinterp.size):
            # get closest point that is greater
            ii=x.searchsorted(xinterp[i])

            minii = ii-2
            if minii < 0:
                minii=0
            maxii = ii+4

            tck = interpolate.splrep(x[minii:maxii],y[minii:maxii],k=3,s=0)

            # evaluate the spline. der=0 means don't bother
            # with derivatives
            yinterp[i] = interpolate.splev(xinterp[i],tck,der=0)
        return yinterp

    def pk(self, k, norm=5.e6):
        """

        Get a P(k)
            k^ns * T^2

        A fake normalization of 5.e6 is applied and can be controlled
        through the keyword norm=

        To get the p(k) multiplied by (sigma_8/sigma)**2 call xi
        wit the getpk=True keyword
        """

        T = self.ehu.T(k)
        Pk=k**self.ns * T**2
        return Pk*norm

    def pkgen(self, log_kmin, log_kmax, **keys):
        """
        Name:
            pkgen
        Calling Sequence:
            k,pk = pkgen(log_kmin,log_kmax,**keys)
            
        Purpose:
            Generate k values equally spaced in log
            space and return (k,Pk)

        Inputs:
            log_kmin,log_kmax: log10 min and max
        Optional Inputs:
            nplk=160: the number of points per unit log(k)
                Then npts=(log_kmax-log_kmin)*nplk
            npts: Number of points.  This will over-ride nplk. 
            norm=5.e6

        """
        norm = keys.get('norm',5.e6)
        npts = keys.get('npts',None)
        if npts is None:
            nplk = keys.get('nplk',160)
            npts=(log_kmax-log_kmin)*nplk

        k=10**linspace(log_kmin, log_kmax, npts) 

        return k, self.pk(k,norm=norm)

    def plot_pk(self,k,pk):
        from biggles import FramedPlot, Curve

        plt=FramedPlot()
        plt.add(Curve(k,pk))
        plt.xlog=True
        plt.ylog=True

        plt.xlabel = r'$k [h/Mpc]$'
        plt.ylabel = r'$P_{lin}(k)$'
        plt.aspect_ratio = 1
        plt.show()

    def plot_pkgen(self):
        k,pk = self.pkgen(-4.1,1.1, npts=1000)
        self.plot_pk(k,pk)


 


    def xi(self, rin, **keys):
        """
        Name:
            xi
        Calling Sequence:
            xi(r, nplk=)
        Purpose:
            Calculate the 3d linear correlation function by integrating
            over the power spectrum.

        Inputs:
            r: radius in Mpc
            nplk=: The number of points in P(k) per log k.  See
                pkgen() for the default.
        """

        r = numpy.array(rin, dtype='f8', copy=False, ndmin=1)
        
        # first get P(k) on a log spaced grid
        log_kmin = -5.0
        log_kmax =  6.0


        k,Pk = self.pkgen(log_kmin, log_kmax, **keys)

        rmax=2000.0
        if r.max() > rmax:
            raise ValueError("Some rmax greater than %s" % rmax)

        xi = numpy.zeros(r.size, dtype='f8')

        for i in xrange(r.size):
            #xi[i] = self.xi_int_pk(k,Pk,r[i])
            xi[i] = self.xi_int_pk_slow(k,Pk,r[i])

        w=where1(r < 25.0)
        rr,ss = self.xi2sigmavar(r[w], xi[w])

        s8=exp(interplin(log(ss),log(rr),log(8.0)))

        xi=xi*(self.sigma_8/s8)**2

        return xi

    def drho(self, r, **keys):
        """

        In principle this should be xi*rhobar. But we leave it as

            xi*rhocrit

        This is so fitting for the amplitude B would give you

            \Omega_m \sigma_8^2 D(z)^2 b(M,z)

        
        """
        xi = self.xi(r, **keys)
        # not putting in real norm, doing relative to rhocrit
        drho=xi*self.rhocrit0
        return drho

    def j3(self, r, **keys):
        """
        This is proportional to the mass enclosed.  If you
        do a fit for B, then you should be able to just multiply
        B*rhocrit0*j3 to get the mass.  rhocrit0 is the rhocrit
        at redshift 0.

        See the m() for that very thing.

        You should use a small radius to get the mass right
        """
        xi = self.xi(r, **keys)

        w=where1(xi <= 0)
        if w.size > 0:
            raise ValueError("all xi must be > 0 for power law "
                             "interpolation")
        lr  = log(r)
        lxi = log(xi)

        al=( lxi-roll(lxi,1) )/( lr-roll(lr,1) )
        al[0]=al[1]

        A=xi/r**(al)

        ex=3.0+al
        Rb=r
        Ra=roll(r,1)
        Ra[0]=0.0
        int0=A*(1.0/ex) *(Rb**ex -Ra**ex)

        j3=4*PI*int0.cumsum()
        return j3

    def m(self, r, B, **keys):
        """
        The mass enclosed.  B is

            \Omega_m \sigma_8^2 D(z)^2 b(M,z)

        And would be the result, for example, of a 

            nfw.dsig() + B*lin.dsig() fit

        """
        j3 = self.j3(r, **keys)
        return j3*self.rhocrit0*B*1.e12

    def xigen(self, log_rmin, log_rmax, npts=1000):
        """
        Generate some xi values on log spaced grid
        """
        r = 10.0**linspace(log_rmin, log_rmax, npts)
        return r, self.xi(r)

    def plot_xi(self, r, xi):
        from biggles import FramedPlot, Curve
        minval = 1.e-4

        xi = where(xi < minval, minval, xi)

        plt=FramedPlot()
        plt.add(Curve(r,xi))
        plt.xlog=True
        plt.ylog=True

        plt.xlabel = r'$r [Mpc/h]$'
        plt.ylabel = r'$\xi_{lin}(r)$'
        plt.aspect_ratio=1
        plt.show()


    def plot_xigen(self):
        r,xi = self.xigen()
        self.plot_xi(r,xi)

    def xi_int_pk_slow(self, k, Pk, r, debug=False, doplot=False):
        """
        Maybe try a slow but more accurate version here
        get xi at position r by integrating P(k). r should be a scalar

        the integral is done in two parts

        """
        import scipy.integrate

        # some of this could be pre-computed, but it might
        # not matter

        NumPerPer = 30.0
        # Number of samples per period, should always be > 5 at least
        npd = 100
        # depends on cosmology but this is ballpark                  
        L_wiggle=170.0

        Pref=1.0 /(2.0 * PI**2 * r)


        dk1=2.0 *PI/(NumPerPer*L_wiggle)
        dk2=2.0 *PI/(NumPerPer*r)
        dktry=min([dk1,dk2])


        # first section
        kmin=0.0
        kmax=2.0
        numk=int32((kmax-kmin)/dktry)
        dk=(kmax-kmin)/(1.0*numk-1)


        kk=arange(numk,dtype='f8')*dk+kmin

        # Dave used *local* cubic spline interpolation (the /spline)
        # keyword for interpol.  This seems to give radically different
        # results!  Using interpsp2 above gets us within about 6%
        # of dave's answer.  Which is right?  I think the right thing
        # to do is use more values for Pk,k and do linear interpolation

        Pkk=interplin(Pk,k,kk)
        #Pkk=interp(kk,k,Pk)
        tab=Pkk*kk*sin(kk*r)


        integ = scipy.integrate.simps(tab,x=kk)
        xi_1=Pref*integ

        if debug:
            print 'r=',r
            if dk1 < dk2:
                print 'Baryon wiggle limited' 
            else: 
                print 'Sin(k*r) limited'
            print 'Numk=',numk
            print 'dk=',dk

            integ_trap = scipy.integrate.trapz(tab,x=kk)
            integ_qg1000 = self.qg1000.integrate(kk,tab)
            print 'integ=',integ
            print 'integ_trap=',integ_trap
            print 'integ_qg1000=',integ_qg1000
            nprint=20
            
            if doplot:
                wk=where1(k < kmax)
                wkk=where1(kk > 1.e-8)

                ptab=Table(2,1)

                lplt=FramedPlot()
                lplt.add( Points(k[wk],Pk[wk],type='filled circle',size=1) )
                lplt.add( Points(kk[wkk],Pkk[wkk], color='red', 
                                 type='filled circle', size=0.7))
                lplt.add( Curve(kk[wkk],Pkk[wkk], color='blue'))
                lplt.xlog=True
                ptab[0,0] = lplt

                splt = FramedPlot()
                splt.add(Points(k[wk],Pk[wk]*k[wk]*sin(k[wk]*r),
                                type='filled circle',size=1))
                splt.add(Points(kk[wkk],tab[wkk],
                                 color='red',type='filled circle',size=0.7))
                splt.add(Curve(kk[wkk],tab[wkk], color='blue'))
                splt.xlog=True

                ptab[1,0] = splt
                ptab.show()


        # next  section/method

        # now need integral \int_xmax^infinity x P(x/r) sin(x)
        Pref2=Pref/(r**2)

        xmax=r*kmax
        num_dec=5
        num_per_dec=npd

        numx=num_dec*num_per_dec
        x=xmax*10.0**(arange(numx,dtype='f8')/num_per_dec)
        Px=interplin(Pk,k,x/r)
        y=x*Px
        ly=log(y)
        n=numx

        al=(roll(ly,-1)-ly)/(roll(x,-1)-x)
        al[n-1]=al[n-2]

        Amp=y/exp(al*x)

        # integral boundaries
        a=x
        b=roll(x,-1)

        norm=y/(1+al**2)
        Ta=al*sin(a)-cos(a)
        Tb=al*sin(b)-cos(b)
        dif=exp(al*(b-a))*Tb-Ta
        d=norm*dif
        d=d[0:n-2]

        integ=d.sum()

        xi_2=integ*Pref2
        xi=xi_1+xi_2
        if debug:
            print 'xi_1:',xi_1
            print 'xi_2:',xi_2

        return xi



    def xi_int_pk(self, k, Pk, r):
        """

        get xi at position r by integrating P(k). r should be a scalar

        the integral is done in two parts

        """
        import scipy.integrate

        # some of this could be pre-computed, but it might
        # not matter

        NumPerPer = 30.0
        # Number of samples per period, should always be > 5 at least
        npd = 100
        # depends on cosmology but this is ballpark                  
        L_wiggle=170.0

        Pref=1.0 /(2.0 * PI**2 * r)

        dk1=2.0 *PI/(NumPerPer*L_wiggle)
        dk2=2.0 *PI/(NumPerPer*r)
        dktry=min([dk1,dk2])

        # first section
        kmin=0.0
        kmax=2.0
        numk=int32((kmax-kmin)/dktry)
        dk=(kmax-kmin)/(1.0*numk-1)


        kk=arange(numk,dtype='f8')*dk+kmin

        # Dave used cubic spline interpolation
        Pkk=interplin(Pk,k,kk)
        tab=Pkk*kk*sin(kk*r)

        integ = scipy.integrate.simps(tab,kk)
        xi_1=Pref*integ

        # next  section/method

        # now need integral \int_xmax^infinity x P(x/r) sin(x)
        Pref2=Pref/(r**2)

        xmax=r*kmax
        num_dec=5
        num_per_dec=npd

        numx=num_dec*num_per_dec
        x=xmax*10.0**(arange(numx,dtype='f8')/num_per_dec)
        Px=interplin(Pk,k,x/r)
        y=x*Px
        ly=log(y)
        n=numx

        al=(roll(ly,-1)-ly)/(roll(x,-1)-x)
        al[n-1]=al[n-2]

        Amp=y/exp(al*x)

        # integral boundaries
        a=x
        b=roll(x,-1)

        norm=y/(1+al**2)
        Ta=al*sin(a)-cos(a)
        Tb=al*sin(b)-cos(b)
        dif=exp(al*(b-a))*Tb-Ta
        d=norm*dif
        d=d[0:n-2]

        integ=d.sum()

        xi_2=integ*Pref2
        xi=xi_1+xi_2

        return xi



    def xi2sigmavar(self,r,xi):
        """

        Compute sigma variance from xi where xi is the 3d correlation
        function. returns rr and sig, rr is not exactly the same as r
        
        You can calculate sigma8 from the results

            rr,ss = self.xi2sigmavar(r, xi)
            s8=exp(interplin(log(ss),log(rr),log(8.0)))
        """

        w = where1( xi <= 0 )
        if w.size > 0:
            raise ValueError("xi must be positive to perform powerlawe interp")

        # choose a new set of points, logarithmically spaced such that a
        # whole number of points (nper2) correspond to a factor of two

        nper2 = 100

        rmax = r.max()
        rmin = r.min()

        tmp = int32( ceil(log(rmax/rmin)/log(2)) )
        num = 1+nper2*tmp

        tmp = arange(num,dtype='f8')/nper2
        rr = rmin*2.0**tmp

        log_r   = log(r)
        log_xi  = log(xi)
        log_rr  = log(rr)

        lxx = interplin(log_xi,log_r,log_rr)
        xxi = exp(lxx)

        # roll is same as shift in IDL
        al=(roll(lxx,-1)-lxx)/(roll(log_rr,-1)-log_rr)
        al[num-1]=al[num-2]

        A=xxi/rr**(al)
        Ra=rr.copy()
        Ra[0]=0.0

        # first one , integrate to zero
        Rb=roll(rr,-1)

        # 3 integrals to do

        beta=0.0
        ex=3.0+al+beta
        int0=A*(1.0/ex) *(Rb**ex -Ra**ex)

        beta=1.0
        ex=3.0+al+beta
        int1=A*(1.0/ex) *(Rb**ex -Ra**ex)

        beta=3.0
        ex=3.0+al+beta
        int2=A*(1.0/ex) *(Rb**ex -Ra**ex)

        # array of sub integrals, last one to be ignored

        # Now add up the sub integrals with 0 to rmin integral
        T0=roll(int0.cumsum(),-nper2)
        T1=roll(int1.cumsum(),-nper2)
        T2=roll(int2.cumsum(),-nper2)

        sig2=3.0*T0/(rr**3) + (-9.0/4)*T1/(rr**4) + (3.0/16)*T2/(rr**6)

        # now trim off last bunch
        rr=rr[0:num-nper2-2]
        log_rr=log_rr[0:num-nper2-2]
        sig2=sig2[0:num-nper2-2]
        sig=sqrt(sig2)

        # now interpolate back onto original points

        ss=exp(interplin(log(sig),log_rr,log_r))
        wcut=where1(r*2 < rmax) 
        rr=r[wcut]
        sig=ss
        sig=sig[wcut]

        return rr, sig

    def dsig(self,r, **keys):
        """
        Calculate \Delta\Sigma

        NOTE: the normalization is simply

            delta_rho = rhocrit*xi

        instead of 
        
            rhomean*xi

        The redshift dependence is left out of this calculation. You
        should add it in as needed, or leave the normalization as part
        of a fit.

        """

        drho = self.drho(r, **keys)

        sig = project.project3d(r, drho)
        dsig = self.sigma2dsig(r, sig)

        return dsig

    def sigma2dsig(self, r, sig):
        """

        Compute \Delta\Sigma from 2D sigma. Assumes power law
        interpolation and powerlaw extrapolation

        """

        # slopes
        al=log( sig/roll(sig,1) )/log( r/roll(r,1) )
        al[0] = al[1]
        slope_min = -1.95
        if al[0] < slope_min:
            print 'Warning: profile too steep to converge: ',al[0]
            print 'truncating to inner slope of: ',slope_min
            al[0] = slope_min

        A=sig/(r**al)

        RA=roll(r,1)
        RA[0]=0.0
        RB=r


        Ints=A/(al+2.0) *(RB**(al+2.0)-RA**(al+2.0))
        Icum=Ints.cumsum()
        avg_sig=2.0*Icum/(r**2)

        dsig=avg_sig-sig

        return dsig

class EisenHu:
    def __init__(self,**keys):

        self.omega_m=keys.get('omega_m',0.25)
        self.omega_b=keys.get('omega_b',0.055)
        self.sigma_8=keys.get('sigma_8',0.8)
        self.h=keys.get('h',1.0)

 
        #Tcmb=2.728
        self.Tcmb=2.72511  # What Bullocks code uses 

        self.calc_constants()

    def calc_constants(self):

        bfrac=self.omega_b/self.omega_m
        self.bfrac=bfrac

        Theta=self.Tcmb/2.7

        Om=self.omega_m*self.h**2
        Ob=self.omega_b*self.h**2

        # Equation 2
        Zeq=2.50e4 * Om * Theta**(-4)

        # Equation 3
        Keq=7.46e-2 * Om * Theta**(-2)
        self.Keq = Keq

        # Equations 4
        b1=0.313 * Om**(-0.419) * (1+ 0.607 * Om**0.674)
        b2=0.238 * Om**0.223
        Zd= 1291 * Om**0.251 * (1 + b1 * Ob**b2) / (1+0.659 * Om**0.828)

        # Equation 5
        Rconst=31.5 * Ob * Theta**(-4)
        Rd  = Rconst * (Zd/1e3)**(-1)
        Req = Rconst * (Zeq/1e3)**(-1)
        self.Rd=Rd

        # Equation 6

        s=2.0/(3.0*Keq) * sqrt(6.0/Req) * log((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
        self.s=s

        # Equation 7
        self.Ksilk = 1.6 * Ob**0.52 * Om**0.73 * (1 + (10.4 * Om)**(-0.95))

        # Equations 11,12

        a1 = (46.9 * Om)**0.670 * (1 + (32.1 * Om)**(-0.532))
        a2 = (12.0 * Om)**0.424 * (1 + (45.0 * Om)**(-0.582))
        self.alc= a1**(-bfrac) * a2**(-bfrac**3) 

        bb1 = 0.944 * (1 + (458.0 * Om)**(-0.708))**(-1)
        bb2 = (0.395 * Om)**(-0.0266)
        self.betac = (1+ bb1*((1.0-bfrac)**bb2-1))**(-1)

        # baryon constants
        # equations 14,15
        y=(1.0+Zeq)/(1.0+Zd)
        sq=sqrt(1.0+y)
        Gy=y * (-6.0*sq + (2+3*y)*log((sq+1.0)/(sq-1.0)))
        self.alb = 2.07*Keq*s*(1+Rd)**(-0.75)*Gy

        # Equations 22,23
        self.betanode = 8.41 * Om**0.435

        # Equation 24
        self.betab=0.5 + bfrac + (3-2*bfrac) * sqrt(1.0+(17.2*Om)**2)

    def T(self, kinput, retall=False):
        """
        Name:
            T
        Purpose:
            Calculate the matter transfer function using the formulas
            in Eisenstein and Hu 1998 , ApJ 496,605

            Ported from Dave Johnston's IDL code
        Calling Sequence:
            T(k, retall=False)
        Inputs:
            k: wavenumber in h/Mpc
            retall=False: if true, return (T,Tc,Tb)
        """

        # convert to physical 1/Mpc wavenumbers
        k = numpy.array(kinput, dtype='f8', copy=True)
        k = k*self.h



        # Equation 10
        q=k/(13.41 * self.Keq)

        # Equations 18-20

        e=numpy.e
        s=self.s
        nat   = log(e + 1.8 * self.betac * q)
        Term1 = 386.0/(1.0+69.9 * q**1.08)
        Co    = 14.2/self.alc + Term1
        To    = nat/(nat+ Co * q**2)
        f     = 1.0/(1.0+(k*s/5.4)**4)
        Co1   = 14.2 + Term1
        To1   = nat/(nat+ Co1 * q**2)

        # Equation 17
        Tc= f * To1 + (1.0-f) * To

        # Now the baryons

        # Equations 22,23

        stil=s*(1.0 + (self.betanode/(k*s))**3)**(-0.3333333)

        # Equation 24
        nat2 = log(e+1.8*q) 
        To11=nat2/(nat2+Co1*q**2)


        x=k*stil
        jo=sin(x)/x

        Tb=(To11/(1.0+(k*s/5.2)**2) + \
            self.alb/(1.0+(self.betab/(k*s))**3)*exp(-(k/self.Ksilk)**1.4))*jo

        T=self.bfrac*Tb + (1.0-self.bfrac)*Tc

        if retall:
            return T,Tc,Tb
        else:
            return T

    def plot_T(self, k, T, Tc=None, Tb=None):
        from biggles import FramedPlot, Curve, PlotKey

        plt=FramedPlot()

        c=Curve(k, T**2)
        c.label = '$T^2$'
        plt.add(c)
        plist = [c]

        if Tc is not None:
            cc=Curve(k,Tc**2,color='blue')
            cc.label = '$Tc^2$'
            plt.add(cc)
            plist.append(cc)

        if Tb is not None:
            tmp = where(Tb < 1.e-5, 1.e-5, Tb)
            cb=Curve(k,tmp**2,color='red')
            cb.label = '$Tb^2$'
            plt.add(cb)
            plist.append(cb)

        plt.xlog=True
        plt.ylog=True
        plt.ylabel = '$T^2'
        plt.xlabel = 'k'
        plt.yrange = [1.e-8,1.0]
        plt.aspect_ratio=1

        if Tc is not None or Tb is not None:
            key=PlotKey(0.1,0.9,plist)
            plt.add(key)

        plt.show()

    def plot_Tgen(self):
        rmin = 0.1 # Mpc
        rmax = 1000 # Mpc
        kmin = 1.0/rmax
        kmax = 1.0/rmin

        logk = linspace(log10(kmin), log10(kmax), 1000)
        k=10.0**logk

        T,Tc,Tb = self.T(k,retall=True)
        self.plot_T(k, T, Tc, Tb)



# convergence test for xi with the number of p(k) sample
# points used.
def test_xi_converge_nplk(epsfile=None):
    """
    Test how xi converges with the number of k points per log10(k)
    Note we should test other convergence factors too!
    """

    tab=Table(2,1) 
    pltbig = FramedPlot()
    pltzoom = FramedPlot()

    pltbig.xlabel='r'
    pltbig.ylabel='xi(r)'
    pltbig.xlog=True
    pltbig.ylog=True
    pltzoom.xlabel='r'
    pltzoom.ylabel='xi(r)'

    l=Linear()
    r=10.0**linspace(0.0,2.3,1000) 
    nplk_vals  = [20,60,100,140,160]
    color_vals = ['blue','skyblue','green','orange','magenta','red']

    plist=[]
    lw=2.4
    for nplk,color in zip(nplk_vals,color_vals):
        print 'nplk:',nplk
        xi =  l.xi(r, nplk=nplk)

        limxi = where(xi < 1.e-5, 1.e-5, xi)
        climxi = Curve(r,limxi,color=color,linewidth=lw)
        climxi.label = 'nplk: %i' % nplk
        pltbig.add(climxi)

        plist.append(climxi)

        w=where1(r > 50.0)
        cxi = Curve(r[w], xi[w], color=color, linewidth=lw)
        pltzoom.add(cxi)

    key = PlotKey(0.7,0.8, plist)
    pltzoom.add(key)
    tab[0,0] = pltbig
    tab[1,0] = pltzoom
    if epsfile is not None:
        tab.write_eps(epsfile)
    else:
        tab.show()


# compare do Dave's code
def compare_idl(ratios=False, epsfile=None):
    from biggles import FramedPlot,Points, Table
    omega_m=0.25
    omega_b=0.055
    sigma_8=1.0
    h=1.0
    ns=0.98

    r=eu.recfile.Open('/home/esheldon/tmp/linear-compare/pkidl.dat',
                      delim=' ',dtype=[('k','f8'),('pk','f8')])
    pkidlst = r.read()
    r.close()

    k,pkidl = pkidlst['k'],pkidlst['pk']

    r=eu.recfile.Open('/home/esheldon/tmp/linear-compare/xiidl.dat',
                      delim=' ',dtype=[('r','f8'),('xi','f8')])
    xiidlst = r.read()
    r.close()

    r,xiidl = xiidlst['r'],xiidlst['xi']

    l = Linear(omega_m=omega_m,omega_b=omega_b,sigma_8=sigma_8,h=h,ns=ns)

    pk = l.pk(k)
    xi = l.xi(r, nplk=1000)
    #xi = l.xi(r)

    pkrat = pk/pkidl
    xirat = xi/xiidl

    tab=Table(2,1)

    symsize=1

    if ratios:

        pkplt = FramedPlot()
        pkplt.xlabel='k'
        pkplt.xlog=True
        pkplt.ylabel='P(k)/P(k) dave'
        pkplt.add( Points(k,pkrat,type='filled circle',size=symsize) )
        pkplt.add( Curve(k,pkrat,color='blue') )

        tab[0,0] = pkplt

 
    if ratios:

        xiplt = FramedPlot()
        xiplt.xlabel='r'
        xiplt.xlog=True

        xiplt.ylabel='xi(r)/xi(r) dave'
        xiplt.add( Points(r,xirat,type='filled circle',size=symsize) )
        xiplt.add( Curve(r,xirat,color='blue') )

        tab[1,0] = xiplt
    else:

        # big log-log plot
        limxi = where(xi < 1.e-5, 1.e-5, xi)
        #pxi = Points(r,limxi,type='filled circle',size=symsize)
        cxi = Curve(r,limxi,color='blue')

        limxiidl = where(xiidl < 1.e-5, 1.e-5, xiidl)
        #pxiidl = Points(r,limxiidl,type='filled circle',size=symsize)
        cxiidl = Curve(r,limxiidl,color='red')

        xipltall = FramedPlot()
        xipltall.xlabel='r'
        xipltall.xlog=True
        xipltall.ylog=True
        xipltall.ylabel='xi(r)'

        #xipltall.add(pxi,cxi,pxiidl,cxiidl)
        xipltall.add(cxi,cxiidl)
       
        tab[0,0] = xipltall

        # zoomed in plot
        xiplt = FramedPlot()
        xiplt.xlabel='r'
        xiplt.xlog=True
        xiplt.ylabel='xi(r)'

        pxi = Points(r,xi,type='filled circle',size=symsize) 
        cxi = Curve(r,xi,color='blue')
        cxi.label = 'Mine'

        pxiidl = Points(r,xiidl,type='filled circle',size=symsize)
        cxiidl = Curve(r,xiidl,color='red')
        cxiidl.label = 'Dave'

        key=PlotKey(0.8,0.9, [cxi,cxiidl])

        xiplt.add(pxi,cxi,pxiidl,cxiidl)
        xiplt.xrange = [50.0,r.max()]
        xiplt.yrange=[-2.e-3,1.e-2]
        xiplt.add(key)

        tab[1,0] = xiplt



    if epsfile is not None:
        tab.write_eps(epsfile)
    else:
        tab.show()

    

