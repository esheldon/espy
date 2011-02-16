"""
project3d:

    Compute the 2D projection of points prepresenting a function of the
    radius r=sqrt(x**2 + y**2 + z**2). E.g. for a function rho(r) where
    r is a three-d radius, compute the integral along a line of sight.

    Uses a three point neighborhood to interpolate the
    function as a quadratic. Then uses analytic formulas to compute
    the integral endpoint correction by extrapolating.

    The last point is computed using a powerlaw extrapolation. You might
    want to just remove the last point if you don't trust this.

    Based on original IDL code by Dave Johnston

lag_quad_poly(x,y):
    Computes the uniq quadratic poly through three x,y points

"""
import numpy
from numpy import int32, log, log10, arange, sqrt

def project3d(r, rho, extrapolate=True):
    """

    Compute the 2D projection of points prepresenting a function of the
    radius r=sqrt(x**2 + y**2 + z**2). E.g. for a function rho(r) where
    r is a three-d radius, compute the integral along a line of sight.

    Uses a three point neighborhood to interpolate the
    function as a quadratic. Then uses analytic formulas to compute
    the integral endpoint correction by extrapolating.

    The last point is computed using a powerlaw extrapolation. You might
    want to just remove the last point if you don't trust this.

    """

    nr = r.size
    if extrapolate:
        # extrapolate then recurse
        r_ext, rho_ext = extrapolate_powerlaw(r,rho)
        sig = project3d(r_ext, rho_ext, extrapolate=False)
        sig = sig[0:nr]
        return sig

    sig = numpy.zeros(nr, dtype='f8')
    for j in xrange(nr-2):
        RR=r[j]
        RR2=RR**2
        num=nr-j
        Int=numpy.zeros(num,dtype='f8')

        for i in xrange(num-2):
            x=r[i+j:i+j+3]
            y=rho[i+j:i+j+3]
            p=lag_quad_poly(x,y)
            A = x[0]
            B = x[1]
     
            I0,I1,I2=getI(RR,RR2,A,B,p)

            Int[i]=2*(I0+I1+I2)
        sig[j] = Int.sum()
    return sig

def getI(RR,RR2,A,B,p):
    Rad=B
    Rad2 = Rad**2
    S=sqrt(Rad2-RR2)
    I0_B = S
    I1_B = S*Rad/2.0 + RR2*log(Rad+S)/2.0
    I2_B = S*(2*RR2 + Rad2)/3.0

    Rad=A
    Rad2 = Rad**2
    S=sqrt(Rad2-RR2)
    I0_A = S
    I1_A = S*Rad/2.0 + RR2*log(Rad+S)/2.0
    I2_A = S*(2*RR2 + Rad2)/3.0

    I0=(I0_B-I0_A)*p[0]
    I1=(I1_B-I1_A)*p[1]
    I2=(I2_B-I2_A)*p[2]

    return I0,I1,I2

def extrapolate_powerlaw(r,rho):
    nr = r.size
    # extrapolate for the integral, but then trim back

    # the fraction of the log interval to add on, 0.5 is probably good
    fac = 1.0
    n_extf=nr*fac
    if n_extf < 10.0: n_extf = 10.0
    n_ext = int32(n_extf)

    rmin=r.min()
    rmax=r.max()
    Lext=( log10(rmax)-log10(rmin) )*fac

    grid = 1+arange(n_ext,dtype='f8')
    r_ext=rmax*10.0**(Lext*grid/(n_ext-1))

    
    al=log(rho[nr-1]/rho[nr-2])/log(r[nr-1]/r[nr-2])
    A=rho[nr-1]/(r[nr-1]**al)
    rho_ext=A*r_ext**al

    r_ext = numpy.concatenate( (r,r_ext) )
    rho_ext = numpy.concatenate( (rho,rho_ext) )

    return r_ext, rho_ext


def lag_quad_poly(x,y):
    """
    NAME: 
        lag_quad_poly

    PURPOSE:
        Computes the uniq quadratic poly through three x,y points

    CALLING SEQUENCE:
        p=lag_quad_poly(x,y)

    INPUTS:
        x,y -  set of three points

    OUTPUTS:
        polynomial -  three numbers , usual IDL notation

    METHOD:
        Uses Lagrange formula
    """

    if x.size !=3 or y.size != 3:
        raise ValueError("x,y must be length 3")

    q=numpy.zeros(3,dtype='f8')
    p=numpy.zeros(3,dtype='f8')

    q[0]=y[0]/((x[0]-x[1])*(x[0]-x[2]))
    q[1]=y[1]/((x[1]-x[0])*(x[1]-x[2]))
    q[2]=y[2]/((x[2]-x[0])*(x[2]-x[1]))

    p[0]=q[0]*x[1]*x[2] + q[1]*x[0]*x[2] + q[2]*x[0]*x[1]
    p[1]=-q[0]*(x[1]+x[2])-q[1]*(x[0]+x[2]) -q[2]*(x[0]+x[1])
    p[2]=q.sum()

    return p
