from __future__ import print_function
import esutil as eu
from numpy import tanh, arctanh, sqrt
import copy

def add_distortion(e1_1, e2_1, e1_2, e2_2):
    """
    Add distortions delta_1 + delta2.  The order matters.

    Note 1,2 reversed from BJ02
    """

    esq_1 = e1_1**2 + e2_1**2
    esq_2 = e1_2**2 + e2_2**2
    if esq_1 > 1 or esq_2 > 1:
        raise ValueError("ellipticities must be <= 1")
    if esq_2 == 0:
        # no distortion
        return copy.copy(e1_1), copy.copy(e2_1)
    if esq_1 == 0:
        # first one is round, just return second
        return copy.copy(e1_2), copy.copy(e2_2)

    edot = e1_1*e1_2 + e2_1*e2_2
    oneplusedot = 1. + edot

    if oneplusedot == 0:
        return 0., 0.

    fac = (1.-sqrt(1.-esq_1) )/esq_1
    
    e1o = e1_1 + e1_2 + (e1_1 * e2_2 - e2_1 * e1_2)*fac*e2_1
    e2o = e2_1 + e2_2 + (e2_1 * e1_2 - e1_1 * e2_2)*fac*e1_1
    e1o /= oneplusedot
    e2o /= oneplusedot

    return e1o, e2o

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


