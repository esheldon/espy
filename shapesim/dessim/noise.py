from numpy import array
from .util import FILTERNUM, FILTERCHAR

MAGLIMS    = array( [24.9,24.5,23.7,23.2,21.5] )
TEFFS      = array( [900.,900.,900.,900.,900.,900.] )
READNOISES = array( [0., 0., 0., 0., 0.] )
GAINS      = array( [1.,1.,1.,1.,1.] ) # electrons/adu

def get_sky(filter, teff):
    """
    Get the sky in electrons
    """
    fnum=FILTERNUM[filter]

    teff_lim=TEFFS[fnum]
    gain=GAINS[fnum]

    arg=(22.5-MAGLIMS[fnum])/2.5
    flim=10.**arg

    # sky in one second
    sky1 = teff*( (flim**2*teff)/100. - flim )

    # sky over teff
    sky = sky1*teff/gain
    
    return sky

def get_skyvar(filter, teff):
    fnum=FILTERNUM[filter]
    sky=get_sky(filter, teff)

    return sky + READNOISES[fnum]
