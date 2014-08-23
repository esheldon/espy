import numpy
from numpy import array

def get_modest(mag_auto, mag_psf, class_star, spread_model, spread_model_err, flags=None):
    """
    get logic for selecting galaxies.

    generally done in i-band

    parameters
    ----------
    mag_auto: scalar or array-like
        auto mag from sextractor
    mag_psf: scalar or array-like
        psf mag from sextractor
    class_star: scalar or array-like
        sextractor star galaxy classifier
    spread_model: scalar or array-like
        new sextractor star-galaxy separation variable
    spread_model_err: scalar or array-like
        new sextractor star-galaxy separation variable, error

    https://cdcvs.fnal.gov/redmine/projects/des-sci-verification/wiki/A_Modest_Proposal_for_Preliminary_StarGalaxy_Separation
    12/20/2013

    Eli's IDL expression

    (FLAGS_I <=3)
        AND NOT 
            (    ((CLASS_STAR_I > 0.3) AND (MAG_AUTO_I < 18.0))
              OR ((SPREAD_MODEL_I + 3*SPREADERR_MODEL_I) < 0.003)
              OR ((MAG_PSF_I > 30.0) AND (MAG_AUTO_I < 21.0))
            )

    i.e. all these must be false
        (class_star > 0.3) and (mag_auto < 18.0)
        (spread_model + 3*spread_model_err) < 0.003
        (mag_psf > 30.0) and (mag_auto < 21.0)

    if flags are sent, and with this
        flags <=3
    

    """

    mag_auto=array(mag_auto, ndmin=1, copy=False)
    mag_psf=array(mag_psf, ndmin=1, copy=False)
    class_star=array(class_star, ndmin=1, copy=False)
    spread_model=array(spread_model, ndmin=1, copy=False)
    spread_model_err=array(spread_model_err, ndmin=1, copy=False)

    bright_test    = ((class_star > 0.3) & (mag_auto < 18.0)) == False
    locus_test     = ((spread_model + 3*spread_model_err) < 0.003) == False
    faint_psf_test = ((mag_psf > 30.0) & (mag_auto < 21.0)) == False

    sg_test = (bright_test & locus_test & faint_psf_test)

    if flags is not None:
        flags=array(flags, ndmin=1, copy=False)
        sg_test = sg_test & (flags <= 3)

    return sg_test


