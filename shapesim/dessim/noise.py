from numpy import array
from .util import FILTERNUM, FILTERCHAR

MAGLIMS    = array( [24.9,24.5,23.7,23.2,21.5] )
TEFFS      = array( [900.,900.,900.,900.,900.,900.] )
READNOISES = array( [0., 0., 0., 0., 0.] )
IGAINS     = array( [1.,1.,4.286,1.,1.] ) # electrons/adu

def get_sky(filter, exptime):
    """
    Get the sky in electrons

    Conventions from the aardvark paper
    """
    fnum=FILTERNUM[filter]
    
    mlim = MAGLIMS[fnum]
    flim = get_flux(mlim,1.0)

    # sky in one second
    sky1 = (flim**2*exptime)/100. - flim

    # sky over exptime
    sky = sky1*exptime
    return sky

def get_skyvar(filter, exptime):
    fnum=FILTERNUM[filter]
    sky=get_sky(filter, exptime=exptime)

    return sky + READNOISES[fnum]

def get_flux(mag, exptime):
    """
    Get flux in electrons
    """
    arg=(22.5-mag)/2.5
    flux=10.**arg

    return flux*exptime
