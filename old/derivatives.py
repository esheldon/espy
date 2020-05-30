"""
currently just test codes
"""
from __future__ import print_function
import numpy
from numpy import array

offsets3_1 = array([2., 1.0, -1.0, -2.0])
coeffs3_1 = array([1.0, -2.0, 2.0, -1.0])

def get_third1(func, x0, h):

    xvals = x0 + offsets3_1*h 

    fvals = func(xvals)

    return (fvals*coeffs3_1).sum()*0.5*(1.0/h)**3

def third_func1(x):
    return x**3

def test_third1():
    print("expect 6")

    h=1.0e-3
    val = get_third1(third_func1, 1.5, h)
    print("got:",val)


def get_partial2x_partialy(func, x0, y0, hx, hy):
    """
    test partial2 x partial y
    """
    xoffsets = array([2.0, 0.0, -2.0,  2.0,  0.0, -2.0])
    yoffsets = array([1.0, 1.0,  1.0, -1.0, -1.0, -1.0])
    coeffs = array([1.0, -2.0, 1.0, -1.0, 2.0, -1.0])
    
    xvals = x0 + xoffsets*hx
    yvals = y0 + yoffsets*hy

    fvals = func(xvals,yvals)

    return (fvals*coeffs).sum()*(1.0/8.0)*(1.0/hx**2/hy)

def get_partialx_partial2y(func, x0, y0, hx, hy):
    """
    test partial2 x partial y
    """
    xoffsets = array([1.0, 1.0,  1.0, -1.0, -1.0, -1.0])
    yoffsets = array([2.0, 0.0, -2.0,  2.0,  0.0, -2.0])
    coeffs = array([1.0, -2.0, 1.0, -1.0, 2.0, -1.0])
    
    xvals = x0 + xoffsets*hx
    yvals = y0 + yoffsets*hy

    fvals = func(xvals,yvals)

    return (fvals*coeffs).sum()*(1.0/8.0)*(1.0/hx/hy**2)


def third_func2(x,y):
    return x**2 * y

def third_func3(x,y):
    return x * y**2


def test_partial2x_partialy(x0=1.0, y0=1.0):
    print("expect 2")

    h=1.0e-3
    val = get_partial2x_partialy(third_func2, x0, y0, h, h)
    print("got:",val)

def test_partialx_partial2y(x0=1.0, y0=1.0):
    print("expect 2")

    h=1.0e-3
    val = get_partialx_partial2y(third_func3, x0, y0, h, h)
    print("got:",val)
