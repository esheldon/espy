import esutil as eu
import numpy
from numpy import tanh, arctanh, sqrt, where
from numpy.random import random as randu
import copy

class ShapeRangeError(Exception):
    pass


def average_shear(e1, e2, e1err=None, e2err=None):
    """
    Calculate the mean shear from shapes.

    parameters
    ----------
    e1: array
        e1 values
    e2: array
        e2 values
    e1err,e2err: 
        Calculate the error on the mean and return
            g1,g2,g1err,g2err

    notes
    -----
    The formula used takes the average of the e_i^2
    compoments for the responsivity
    """


    num=e1.size
    if e2.size != num:
        raise ValueError("e1 and e2 must be same size")

    esq = e1**2 + e2**2
    mesq = esq.mean()
    R = 1-.5*mesq
    me1 = e1.mean()
    me2 = e2.mean()

    g1 = 0.5*me1/R
    g2 = 0.5*me2/R

    if e1err is not None and e2err is not None:
        e1_ivar = 1/e1err**2
        e2_ivar = 1/e2err**2

        g1err = 0.5/sqrt(e1_ivar.sum())/R
        g2err = 0.5/sqrt(e2_ivar.sum())/R

        return g1,g2,R,g1err,g2err
    else:
        return g1,g2, R


def e2g(e):
    """
    Convert "delta" to g

        delta=(a^2-b^2)/(a^2+b^2) 
        g=(a-b)/(a+2)
    """
    if e >= 1:
        raise ShapeRangeError("e must be <= 1, found %.16g" % e)
    g = tanh(0.5*arctanh(e))
    if g >= 1:
        raise ShapeRangeError("converted g must be <= 1, found %.16g" % g)
    return g

def g2e(g):
    """
    Convert g to delta.  
    
        delta=(a^2-b^2)/(a^2+b^2) 
        g=(a-b)/(a+2)
    """
    if isinstance(g,numpy.ndarray):
        e=numpy.zeros(g.size)
        for i in xrange(g.size):
            e[i] = g2e(g[i])
            return e

    if g >= 1:
        raise ShapeRangeError("ellipticity must be < 1, "
                              "found %.16g" % g)
    e = tanh(2*arctanh(g))
    if e >= 1:
        raise ShapeRangeError("converted e must be < 1, found %.16g" % e)
    return e

e2gamma=e2g
gamma2e=g2e

def e1e2_to_g1g2(e1, e2):
    if isinstance(e1,numpy.ndarray):
        g1=numpy.zeros(e1.size)
        g2=numpy.zeros(e1.size)
        for i in xrange(e1.size):
            g1[i],g2[i] = e1e2_to_g1g2(e1[i],e2[i])
    else: 
        e = sqrt(e1**2 + e2**2)
        if e == 0:
            return 0.,0.
        g = e2g(e)
        fac = g/e
        g1, g2 = fac*e1, fac*e2
    return g1,g2

def g1g2_to_e1e2(g1, g2):
    if isinstance(g1,numpy.ndarray):
        e1=numpy.zeros(g1.size)
        e2=numpy.zeros(g1.size)
        for i in xrange(g1.size):
            e1[i],e2[i] = g1g2_to_e1e2(g1[i],g2[i])
    else: 
        g = sqrt(g1**2 + g2**2)
        if g == 0:
            return 0.,0.
        e = g2e(g)
        fac = e/g
        e1, e2 = fac*g1, fac*g2
    return e1,e2


def g1g2_to_eta1eta2(g1,g2):
    gtot=sqrt(g1**2 + g2**2)
    if isinstance(gtot,numpy.ndarray):
        eta1=numpy.zeros(gtot.size)
        eta2=numpy.zeros(gtot.size)

        w,=numpy.where(gtot != 0)
        if w.size > 0:
            cos2theta=g1[w]/gtot[w]
            sin2theta=g2[w]/gtot[w]

            eta = 2*arctanh(gtot[w])
            eta1[w]=eta*cos2theta
            eta2[w]=eta*sin2theta
    else:
        from math import atanh        
        if gtot == 0:
            eta1,eta2=0.0,0.0
        else:
            cos2theta=g1/etot
            sin2theta=g2/etot

            eta = 2*atanh(gtot)
            eta1=eta*cos2theta
            eta2=eta*sin2theta

    return eta1,eta2

def shear_fracdiff(e, em, deriv=1.0):
    """
    e=etrue
    em=emeasured

    Hirata & Seljak eq 27
    putting in 1 for derivative d(emeas)/d(etrue)

    e=etrue
    em=emeas
    deriv deriviative of measured e with true e

    """
    return ((1-e**2)*deriv + em/e)/(2-em**2) - 1.0


def lens_wmom(data, tag, ind=None, sdev=False):
    if ind is None:
        wts = data['weight']
        tdata = data[tag]
    else:
        wts = data['weight'][ind]
        tdata = data[tag][ind]

    return eu.stat.wmom(tdata, wts, calcerr=True, sdev=sdev)

def randomize_e1e2(e1start,e2start):
    if e1start == 0 and e2start==0:
        e1rand = 0.05*(randu()-0.5)
        e2rand = 0.05*(randu()-0.5)
    else:
        nmax=100
        ii=0
        while True:
            e1rand = e1start*(1 + 0.2*(randu()-0.5))
            e2rand = e2start*(1 + 0.2*(randu()-0.5))
            etot = sqrt(e1rand**2 + e2rand**2)
            if etot < 0.95:
                break
            ii += 1
            if ii==nmax:
                wlog("---- hit max try on randomize e1e2, setting zero and restart")
                e1start=0
                e2start=0
                ii=0

    return e1rand, e2rand


def rotate_shape(e1, e2, theta, units='degrees'):
    from math import cos, sin, radians

    if units=='degrees':
        twotheta=2*radians(theta)
    else:
        twotheta=2*theta

    cos2angle=cos(twotheta)
    sin2angle=sin(twotheta)
    e1rot =  e1*cos2angle + e2*sin2angle
    e2rot = -e1*sin2angle + e2*cos2angle

    return e1rot,e2rot
