"""
Define a set of reasonable values for sky noise and zero points based on
commissioning data.

The skies may be brighter than we expect for main survey, but it is hard to
say.

The zero points look pretty stable though.

I can't find gain values anywhere in the DB, but jiangang is finding ~4.5
"""

from numpy import array
from .util import FILTERNUM, FILTERCHAR

# these are super rough from reasonably dark times
# but need to get official "good" data for comparison
SKY_ELECTRONS_PER_SEC=array([5., 13.0, 37., 70.0, 45.])

# can't find these anywhere
GAINS     = array( [4.5,4.5,4.5,4.5,4.5] ) # electrons/adu

# m = -2.5 x log10(DN / EXPTIME) + ZEROPOINT
#   note these are DN(ADU) not electrons.
# from ~oh/tmp/some-zeropoints.txt, the 90 second exposures
# except for Y which is the 45 second

ZEROPOINTS_ADU = [30.15, 30.34, 30.25, 29.97, 28.1]

def get_sky(filter, exptime, units='e'):
    """
    Get a sky value in electrons
    """
    fnum=FILTERNUM[filter]

    e_per_sec = SKY_ELECTRONS_PER_SEC[fnum]
    if units.lower() =='e':
        flux=e_per_sec
    elif units.upper()=='ADU':
        flux=e_per_sec/GAINS[fnum]
    else:
        raise ValueError("units should be e or ADU")

    return flux*exptime

def get_skyvar(filter, exptime, units='e'):
    """
    Get a sky variance
    """

    sky=get_sky(filter, exptime, units=units)
    return sky

def get_flux(filter, mag, exptime, units='e'):
    """
    Get flux in electrons or ADU
    """
    fnum=FILTERNUM[filter]
    zp=ZEROPOINTS_ADU[fnum]

    arg=0.4*(zp-mag)
    adu_per_second=10.**arg

    if units.lower() =='e':
        flux=adu_per_second*GAINS[fnum]
    elif units.upper() =='ADU':
        flux=adu_per_second
    else:
        raise ValueError("units should be e or ADU")

    return flux*exptime
