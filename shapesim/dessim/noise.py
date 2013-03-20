"""
Define a set of reasonable values for sky noise and zero points based on
commissioning data.

The skies may be brighter than we expect for main survey, but it is hard to
say.

The zero points look pretty stable though.
"""

from numpy import array
from .util import FILTERNUM, FILTERCHAR

# these are super rough from reasonably dark times
# but need to get official "good" data for comparison
#SKY_ELECTRONS_PER_SEC=array([5., 13.0, 37., 70.0, 45.])
SKY_ELECTRONS_PER_SEC=array([7.05106,
                             12.5659,
                             41.9028,
                             95.1041,
                             45.0])

#ZEROPOINTS_ADU = [30.15, 30.34, 30.25, 29.97, 28.1]
ZEROPOINTS_E = array( [26.5670,
                       26.3742,
                       26.3213,
                       26.2623,
                       26.3] ) # Y is made up
MAGLIM = array( [24.9560,
                 24.4530,
                 23.7510,
                 23.2490,
                 23.0] ) # Y is made up

def get_sky(filter, exptime):
    """
    Get a sky value in electrons
    """
    fnum=FILTERNUM[filter]

    e_per_sec = SKY_ELECTRONS_PER_SEC[fnum]

    return e_per_sec*exptime

def get_skyvar(filter, exptime):
    """
    Get the background variance: pure sky for now
    """

    sky=get_sky(filter, exptime)
    return sky

def get_flux(filter, mag, exptime):
    """
    Get flux in electrons
    """
    fnum=FILTERNUM[filter]
    zp=ZEROPOINTS_E[fnum]

    arg=0.4*(zp-mag)
    e_per_sec=10.**arg

    return e_per_sec*exptime
