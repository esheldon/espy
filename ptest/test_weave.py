from scipy import weave
from numpy import *

def weave_passarr():
    x = arange(10)

    # gets and returns length of the array x
    # Nx is always defined.
    code = \
         """
         int n = Nx[0];
         return_val = n;
         """
    n = weave.inline(code, ['x'],
                     type_converters = weave.converters.blitz)
    print "n = ",n,x.size

def weave_sum(a):

    size = a.size
    
    code = \
         """
         double sum=0.0;
         for (int i=0; i<size; i++)
           sum += a(i);
         return_val = sum;
         """

    sum = weave.inline(code, ['a','size'],
                       type_converters = weave.converters.blitz)

    return sum



def gauleg(x1, x2, npts):

    # outputs
    x = zeros(npts, dtype=float64)
    w = zeros(npts, dtype=float64)

    code = \
         """
         int i, j, m;
         double xm, xl, z1, z, p1, p2, p3, pp, pi, EPS, abszdiff;
         
         EPS = 3.e-11;
         pi=3.1415927;

         m = (npts + 1)/2;

         xm = (x1 + x2)/2.0;
         xl = (x2 - x1)/2.0;
         z1 = 0.0;

         for (i=1; i<= m; ++i) 
         {
      
           z=cos( pi*(i-0.25)/(npts+.5) );

           abszdiff = fabs(z-z1);

           while (abszdiff > EPS) 
           {
             p1 = 1.0;
             p2 = 0.0;
             for (j=1; j <= npts;++j)
             {
                p3 = p2;
                p2 = p1;
                p1 = ( (2.0*j - 1.0)*z*p2 - (j-1.0)*p3 )/j;
             }
             pp = npts*(z*p1 - p2)/(z*z -1.);
             z1=z;
             z=z1 - p1/pp;

             abszdiff = fabs(z-z1);

           }
      
           x(i-1) = xm - xl*z;
           x(npts+1-i-1) = xm + xl*z;
           w(i-1) = 2.0*xl/( (1.-z*z)*pp*pp );
           w(npts+1-i-1) = w(i-1);


         }

         return_val = 0;
         

         """
    
    weave.inline(code, ['x1', 'x2', 'npts', 'x', 'w'],
                 type_converters = weave.converters.blitz)

    return x,w


def time_reduce():
    import timeit

    t = timeit.Timer("sum = numpy.arange(1000000).sum()","import numpy")
    
    print min(t.repeat(10,1))
