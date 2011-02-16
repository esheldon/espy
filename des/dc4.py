"""
According to Huan this is the setup for DC4
note ra=-25 is really ra=335.0

These are in pixel coordinates
    x -> DEC
    y -> -RA
But we measure in RA/DEC space.  

RA          Dec        gamma1  gamma2
-----------------------------------------
>= -25 deg  >= -25 deg     0.05    0
>= -25      <  -25        -0.025   0.025
<  -25      >= -25         0.05    0.05
<  -25      <  -25         0      -0.025

rotated ra point is 335.0
"""

import numpy
from numpy import where

rabound = 335.0
decbound = -25.0
shearvals = {}
shearvals[1] = [0.05, 0.0]
shearvals[2] = [-0.025, 0.025]
shearvals[3] = [0.05, 0.05]
shearvals[4] = [0.0, -0.025]

def region_select(self, ra, dec, region):

    if region == 1:
        w, = numpy.where(  (ra >= rabound) & (dec >= decbound) )
    elif region==2:
        w, = numpy.where(  (ra >= rabound) & (dec < decbound) )
    elif region==3:
        w, = numpy.where(  (ra < rabound) & (dec >= decbound) )
    elif region==4:
        w, = numpy.where(  (ra < rabound) & (dec < decbound) )
    else:
        raise ValueError("region must be 1,2,3,4")
    return w


