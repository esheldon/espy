import esutil as eu
from math import sqrt, atan2, tanh, atanh
import copy
from . import util
from .util import ShapeRangeError


class Shear:
    """
    A class for shapes and shears

    Examples
    --------
        s1 = Shear()
        s1.set_e1e2(e1=0.2, e2=0.3)
        
        print s1.e1, s1.e2, s1.g1, s1.g2, s1.eta1, s2.eta2

        # This performs proper shear addition
        s3 = s1 + s2
        s4 = s1 - s2
    """
    def __init__(self, 
                 e1=None, e2=None, 
                 g1=None, g2=None,
                 eta1=None,eta2=None):
        if e1 is not None and e2 is not None:
            self.set_e1e2(e1,e2)
        elif g1 is not None and g2 is not None:
            self.set_g1g2(g1,g2)
        elif eta1 is not None and eta2 is not None:
            self.set_eta1eta2(eta1,eta2)
        else:
            raise ValueError("send both e1,e2 or g1,g2 or eta1,eta2")

    def set_e1e2(self, e1, e2):
        etot = sqrt(e1**2 + e2**2)

        if etot >= 1:
            mess="e values must be < 1, found %.16g (%.16g,%.16g)"
            mess = mess % (etot,e1,e2)
            raise ShapeRangeError(mess)

        if etot==0:
            self.eta1,self.eta2,self.g1,self.g2=(0.,0.,0.,0.)
            self.e1,self.e2=(0., 0.)
            return 

        eta = atanh(etot)
        gtot = tanh(eta/2)

        cos2theta = e1/etot
        sin2theta = e2/etot

        self.e1=e1
        self.e2=e2
        self.g1=gtot*cos2theta
        self.g2=gtot*sin2theta
        self.eta1=eta*cos2theta
        self.eta2=eta*sin2theta

    def set_g1g2(self, g1, g2):
        gtot = sqrt(g1**2 + g2**2)
        if gtot >= 1:
            mess="g values must be < 1, found %.16g (%.16g,%.16g)"
            mess = mess % (gtot,g1,g2)
            raise ShapeRangeError(mess)

        if gtot==0:
            self.eta1,self.eta2,self.e1,self.e2=(0.,0.,0.,0.)
            self.g1,self.g2=(0., 0.)
            return 

        eta = 2*atanh(gtot)
        etot = tanh(eta)

        cos2theta = g1/gtot
        sin2theta = g2/gtot

        self.g1=g1
        self.g2=g2
        self.e1=etot*cos2theta
        self.e2=etot*sin2theta
        self.eta1=eta*cos2theta
        self.eta2=eta*sin2theta

    def set_eta1eta2(self, eta1, eta2):
        eta=sqrt(eta1**2 + eta2**2)

        if eta==0.:
            self.e1,self.e2,self.g1,self.g2=(0.,0.,0.,0.)
            return

        etot=tanh(eta)
        gtot=tanh(eta/2.)

        if etot >= 1.0:
            mess="e values must be < 1, found %.16g" % etot
            raise ShapeRangeError(mess)
        if gtot >= 1.0:
            mess="g values must be < 1, found %.16g" % gtot
            raise ShapeRangeError(mess)

        cos2theta = eta1/eta
        sin2theta = eta2/eta

        e1=etot*cos2theta
        e2=etot*sin2theta

        g1=gtot*cos2theta
        g2=gtot*sin2theta

        self.eta1,self.eta2=eta1,eta2
        self.e1,self.e2=e1,e2
        self.g1,self.g2=g1,g2

    def __neg__(self):
        s = Shear()
        s.set_e1e2(-self.e1, -self.e2)
        return s

    def __add__(self, s):
        """
        Add a shear.  If this object is "shape", the convention is

            newshape = shape + shear

        I see in galsim CppShear.cpp that they have a += operator which
        uses a different convention

            shape += shear -> shear+shape
        """
        if not isinstance(s,Shear):
            raise ValueError("Only Shear objects can be added")

        if s.e1 == 0 and s.e2 == 0:
            return copy.deepcopy(self)

        oneplusedot = 1.0 + self.e1*s.e1 + self.e2*s.e2
        if oneplusedot == 0:
            return Shear(e1=0, e2=0)

        esq = s.e1**2 + s.e2**2

        fac = (1.0 - sqrt(1.0-esq))/esq

        e1 = self.e1 + s.e1 + s.e2*fac*(self.e2*s.e1 - self.e1*s.e2)
        e2 = self.e2 + s.e2 + s.e1*fac*(self.e1*s.e2 - self.e2*s.e1)

        e1 /= oneplusedot
        e2 /= oneplusedot

        sout = Shear(e1=e1,e2=e2)
        return sout

    def add_by_g(self, s):
        """
        equivalent to the distortion based addition formula
        """
        if not isinstance(s,Shear):
            raise ValueError("Only Shear objects can be added")

        if s.e1 == 0 and s.e2 == 0:
            return copy.deepcopy(self)

        A = 1 + self.g1*s.g1 + self.g2*s.g2
        B = self.g2*s.g1 - self.g1*s.g2
        denom_inv = 1./(A**2 + B**2)

        g1 = A*(self.g1 + s.g1) + B*(self.g2 + s.g2)
        g2 = A*(self.g2 + s.g2) - B*(self.g1 + s.g1)

        g1 *= denom_inv
        g2 *= denom_inv

        sout = Shear(g1=g1, g2=g2)
        return sout

    def add_complex(self, s):
        gc = complex(self.g1, self.g2)
        sc = complex(s.g1, s.g2)

        newg = (gc + sc)/(1 + sc.conjugate()*gc)

        return Shear(g1=newg.real, g2=newg.imag)

    def __sub__(self, s):
        return self.__add__(-s)
    def __repr__(self):
        return '(%.16g, %.16g)' % (self.e1,self.e2)


class Shape:
    """
    A class for shapes and shears.  This differs from the old
    Shear class in that it uses complex ellipticities which
    simplifies things greatly

    Need to make sure it reproduces exactly the old results

    Examples
    --------
        s1 = Shape(g=complex(0.1, 0.2))
        s1.set_e1e2(g=complex(0.1, 0.2))
        
        print s1.e
        print s1.g

        # This performs proper shear addition
        s3 = s1 + s2
        s4 = s1 - s2
    """
    def __init__(self, e=None, g=None):

        if e is not None:
            self.set_e(e)
        elif g is not None:
            self.set_g(g)
        else:
            raise ValueError("send e= or g=")

    def set_e(self, e):
        """
        Set e with internal consistency, updating g as well

        parameters
        ----------
        e: complex
            complex distortion
        """
        if not isinstance(e,complex):
            raise ValueError("e must be complex")
        self.e=e
        self.g = e/(1+sqrt(1 - abs(e)**2))

    def set_g(self, g):
        """
        Set g with internal consistency, updating g as well

        parameters
        ----------
        g: complex
            complex ellipticity
        """

        if not isinstance(g,complex):
            raise ValueError("g must be complex")
        self.g=g
        self.e = 2*g/(1+abs(g)**2)

    def __add__(self, s):
        """
        Add a shear.  If this object is "shape", the convention is

            newshape = shape + shear

        I see in galsim CppShear.cpp that they have a += operator which
        uses a different convention

            shape += shear -> shear+shape
        """
        if not isinstance(s,Shape):
            raise ValueError("Only Shape objects can be added")

        if abs(s.g) == 0:
            return copy.deepcopy(self)

        newg = (self.g + s.g)/(1 + s.g.conjugate()*self.g)

        return Shape(g=newg)

    def __neg__(self):
        s = Shape(g=-self.g)
        return s

    def __sub__(self, s):
        return self.__add__(-s)

    def __repr__(self):
        return '%s' % self.g


def test():
    """
    make sure all the transformations work
    """

    e1=-0.85
    e2=0.49

    s=Shear(e1=e1, e2=e2)
    print 'input e1,e2:',s.e1,s.e2
    print 'output g1,g2:',s.g1,s.g2
    print 'output eta1,eta2:',s.eta1,s.eta2

    print
    print 'feeding g1,g2 back in'
    s2=Shear(g1=s.g1,g2=s.g2)
    print 'output e1,e2:',s.e1,s.e2
    print 'output eta1,eta2:',s.eta1,s.eta2

    print
    print 'feeding eta1,eta2 back in'
    s2=Shear(eta1=s.eta1,eta2=s.eta2)
    print 'output e1,e2:',s.e1,s.e2
    print 'output g1,g2:',s.g1,s.g2
    print 'output eta1,eta2:',s.eta1,s.eta2


    for eta1,eta2 in zip([-2.3,2.57],[1.79,0.0]):
        print
        print 'input eta1,eta2:',eta1,eta2
        s=Shear(eta1=eta1, eta2=eta2)
        print 'output e1,e2:',s.e1,s.e2
        print 'output g1,g2:',s.g1,s.g2
        print 'output eta1,eta2:',s.eta1,s.eta2

def test_average_shear(shapenoise=0.16, n=1000000):
    import numpy
    from .util import average_shear
    shear=Shear(g1=0.02, g2=-0.04)

    g1=9999+numpy.zeros(n)
    g2=9999+numpy.zeros(n)
    e1obs=numpy.zeros(n)
    e2obs=numpy.zeros(n)

    while True:
        gtot=g1**2 + g2**2
        w,=numpy.where(gtot > 1)
        if w.size == 0:
            break

        g1[w] = shapenoise*numpy.random.randn(w.size)
        g2[w] = shapenoise*numpy.random.randn(w.size)

    for i in xrange(n):
        s = Shear(g1=g1[i], g2=g2[i])
        s = s + shear

        e1obs[i] = s.e1
        e2obs[i] = s.e2

    sh1,sh2,R = average_shear(e1obs, e2obs) 

    sh1err=e1obs.std()/sqrt(n)
    sh2err=e2obs.std()/sqrt(n)
    print 'sh1: %.5f +/- %.5f' % (sh1,sh1err)
    print 'sh2: %.5f +/- %.5f' % (sh2,sh2err)


def dgs_by_dgo_jacob(g1, g2, s1, s2):
    """
    jacobian of the transformation
        |dgs/dgo|_{-g}

    parameters
    ----------
    g1,g2: numbers or arrays
        shape pars for "observed" image
    s1,s2: numbers or arrays
        shape pars for shear, applied negative
    """
    # usable in both C and python
    A = 1 - g1*s1 - g2*s2
    B = - g2*s1 + g1*s2
    C = 1./(B*B + A*A)
    C2 = C*C

    g2ms2 = g2-s2
    g1ms1 = g1-s1

    D = (g2ms2*B + g1ms1*A)
    E = s1*g2ms2 + g1ms1*s2
    F = g1ms1*s1 - g2ms2*s2
    G = (2*s2*B - 2*s1*A)
    H = -g1ms1*B + g2ms2*A
    I = (-2*s1*B - 2*s2*A)

    g1s_by_g1o =  -D*G*C2 -  ( A - F)*C

    g1s_by_g2o =  -D*I*C2 -  ( B - E)*C

    g2s_by_g1o =  -H*G*C2 -  (-B - E)*C

    g2s_by_g2o =  -H*I*C2 -  ( A + F)*C

    return g1s_by_g1o*g2s_by_g2o - g1s_by_g2o*g2s_by_g1o

def dgs_by_dgo_jacob_full(g1, g2, s1, s2):
    """
    same as above but without simplification
    |dgs/dgo|_{-g}
    """

    g1s_by_g1o = -((((g2 - s2)*(-(g2*s1) + g1*s2) + (g1 - s1)*(1 - g1*s1 - g2*s2))*(2*s2*(-(g2*s1) + g1*s2) - 2*s1*(1 - g1*s1 - g2*s2)))/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)**2) -  (1 - g1*s1 - (g1 - s1)*s1 - g2*s2 + (g2 - s2)*s2)/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)

    g1s_by_g2o = -((((g2 - s2)*(-(g2*s1) + g1*s2) + (g1 - s1)*(1 - g1*s1 - g2*s2))*(-2*s1*(-(g2*s1) + g1*s2) - 2*s2*(1 - g1*s1 - g2*s2)))/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)**2) -  (-(g2*s1) - s1*(g2 - s2) + g1*s2 - (g1 - s1)*s2)/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)

    g2s_by_g1o = -(((2*s2*(-(g2*s1) + g1*s2) - 2*s1*(1 - g1*s1 - g2*s2))*(-((g1 - s1)*(-(g2*s1) + g1*s2)) + (g2 - s2)*(1 - g1*s1 - g2*s2)))/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)**2) -  (g2*s1 - s1*(g2 - s2) - g1*s2 - (g1 - s1)*s2)/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)

    g2s_by_g2o =  -(((-((g1 - s1)*(-(g2*s1) + g1*s2)) + (g2 - s2)*(1 - g1*s1 - g2*s2))*(-2*s1*(-(g2*s1) + g1*s2) - 2*s2*(1 - g1*s1 - g2*s2)))/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)**2) -  (1 - g1*s1 + (g1 - s1)*s1 - g2*s2 - (g2 - s2)*s2)/((-(g2*s1) + g1*s2)**2 + (1 - g1*s1 - g2*s2)**2)

    return g1s_by_g1o*g2s_by_g2o - g1s_by_g2o*g2s_by_g1o
