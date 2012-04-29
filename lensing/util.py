from __future__ import print_function
import esutil as eu

def e2gamma(e):
    """
    Convert "delta", e.g. (a^2-b^2)/(a^2+b^2) to gamma which is 0.5*(a-b)/(a+2)
    """
    return tanh(0.5*arctanh(e))

def e1e2_to_g1g2(e1, e2):
    e = sqrt(e1**2 + e2**2)
    g = e2gamma(e)
    fac = g/e
    g1, g2 = fac*e1, fac*e2
    return g1,g2



def lens_wmom(data, tag, ind=None, sdev=False):
    if ind is None:
        wts = data['weight']
        tdata = data[tag]
    else:
        wts = data['weight'][ind]
        tdata = data[tag][ind]

    return eu.stat.wmom(tdata, wts, calcerr=True, sdev=sdev)


