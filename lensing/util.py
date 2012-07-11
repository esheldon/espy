from __future__ import print_function
import esutil as eu
from numpy import tanh, arctanh, sqrt
import copy

def average_shear(e1, e2, weights=None, doerr=False):
    """
    Calculate the mean shear from shapes.

    parameters
    ----------
    e1: array
        e1 values
    e2: array
        e2 values
    weights: array, optional
    doerr: bool, optional
        Calculate the error. Not meaningful in ring tests where the
        orientations are not random.

    notes
    -----
    The formula used takes the average of the e_i^2
    compoments for the responsivity
    """

    if weights is not None:
        raise ValueError("implemented weighted")

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

    if doerr:
        mesq_err = esq.std()/sqrt(num)
        e1err = e1.std()/sqrt(num)
        e2err = e2.std()/sqrt(num)
        
        g1err = g1*sqrt( (mesq_err/mesq)**2 + (e1err/me1)**2 )
        g2err = g2*sqrt( (mesq_err/mesq)**2 + (e1err/me2)**2 )
        return g1,g2,R,g1err,g2err
    else:
        return g1,g2, R


def e2gamma(e):
    """
    Convert "delta", e.g. (a^2-b^2)/(a^2+b^2) to gamma which is 0.5*(a-b)/(a+2)
    """
    return tanh(0.5*arctanh(e))

def gamma2e(gamma):
    return tanh(2*arctanh(gamma))

def e1e2_to_g1g2(e1, e2):
    e = sqrt(e1**2 + e2**2)
    g = e2gamma(e)
    fac = g/e
    g1, g2 = fac*e1, fac*e2
    return g1,g2

def g1g2_to_e1e2(g1, g2):
    g = sqrt(g1**2 + g2**2)
    if g == 0:
        return 0.,0.
    e = gamma2e(g)
    fac = e/g
    e1, e2 = fac*g1, fac*g2
    return e1,e2


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


