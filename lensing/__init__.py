"""

process
-------

Note when not doing randoms we will go ahead and combine the collated lens
splits

    First you need a catalog class for each catalog type.
    See 

    lensing.lcat.RedMapper
    lensing.scat.DR8RegaussCatalog
    lensing.scat.DR8GMixCatalog

    And an associated yaml config file with a sample name.
    Also for actual runs you need another config.  There
    are also scinv config and cosmo config.

    # if the scat are not already generated
    /bin/make-objshear-input.py -t scat -s scat_sample

    #
    # now generate the lens samples and wq scripts
    #

    # note lcat are now always split
    /bin/make-objshear-input.py -t lcat -s rm05

    # creates config,shear,src_reduce,collate
    # version can be gsens, mom etc.
    /bin/make-objshear-proc.py rm05gmix01 gsens

    # in the $LENSDIR/proc/run directory
    incsub run-${run}-[0-9]*.yaml
    
    # we then need to reduce across sources at fixed lens split
    incsub run-${run}-src-reduce-*.yaml

    # and then collate each of the lens splits
    incsub run-${run}-collate-*.yaml

    # combine the collations
    python $ESPY_DIR/lensing/bin/combine-collated-chunks.py rm05gmix01


    # This would bin by lambda into 12 bins, must have defined this binning
    /bin/bin-lenses.py rm05gmix01 lambda 12

    # plot the binning (not corrected, jackknife yet)
    /bin/plot-dsig-byrun.py -t binned rm05gmix01 lambda 12

    #
    # now randoms with sample rmrand01.
    #

    # generation of randoms is slow, so do it in chunks.
    # make sure to put nsplit into the lcat yaml file.
    # for pre-generated ra/dec

    /bin/make-random-chunk-scripts.py rmrand01

    # things are different when generating the ra/dec
    /bin/make-random-chunk-scripts.py --gen-radec rmrand01

    # scripts are in $LENSDIR/proc/{sample}/genrand/

    # now make the associated run and the wq scripts
    # let's call the run simply run-rmrand01gmix01
    # we need to reduce across sources then concatenate
    # the lens splits together

    /bin/make-objshear-proc.py rmrand01gmix01 gsens

    # now
    # - submit the shear wq scripts 
    # - submit the src reduce wq scripts (reduce across sources at 
    #   fixed lens split) *reduce*.yaml
    # - collate the results in the splits by submitting the *collate*.yaml
    # - combine the collated splits
    #       combine-collated-chunks.py rmrand01gmix01

    #
    # match randoms to the lens bins this can eat tons of memory, so probably
    # want to get an interactive job on a high mem machine, maybe 48G machine
    #

    /bin/match-randoms.py -t lambda -n 12 rm05gmix01 rmrand01gmix01
    # if you have additional weights in the randoms collated file
    /bin/match-randoms.py -w rmweight -t lambda -n 12 rm05gmix01 rmrand01gmix01

    # correct from randoms
    /bin/correct-shear.py

    # 

You can then fit NFW or NFW+lin masses

    # with linear term
    fitter = lensing.fit.NFWBiasFitter(omega_m, z, r, **linkeywords)
    p,cov = fitter.fit(dsig, dsigerr, guess, more=False)

    # without linear term
    fitter = lensing.fit.NFWBiasFitter(omega_m, z, r, withlin=False)
    p,cov = fitter.fit(dsig, dsigerr, guess, more=False)

To fit a run/binning and write out a file with the fits:

    fit_nfw_lin_dsig_byrun(run, name, withlin=True, rmax_from_true=False,
                           rmin=None, rmax=None):

Above the name will be something like 'n200-12' or 'm12z3'
    fit_nfw_lin_dsig_byrun('05', 'n200-12', withlin=True, 
                           rmin=None, rmax=None):


Some plots that are on a grid are done currently from the BINNER:
    

    mzbin.plot_dsig_byrun(run, dops=False)
    mzbin.plot_invmass_byrun(run, type='',
                             residual=False,
                             yrange=None,
                             dops=False)
    mzbin.plot_nfwmass_byrun(run, 
                             withlin=True,
                             residual=False,
                             noerr=False,
                             fudge=None,tag='m200',
                             yrange=None,
                             dops=False)
 
Others done one-by-one are done in fit.py, but note this
is actually specific to mz binning right now

    lensing.fit.plot_nfw_lin_fits_byrun(run, name, prompt=False)

"""

from . import util
from . import shear
from .shear import Shear
from .util import ShapeRangeError

from . import pqr

from . import rotation

from . import sigmacrit
from . import lens_efficiency
from . import test_estimators
from . import lcat
from . import scat
from . import convert
from . import objshear_config
from . import scripts
from . import condor
from . import wqsubmit


try:
    from . import regauss
    from . import regauss_test
    from . import regauss_sim
except:
    pass

try:
    from . import gmix_sdss_old
except:
    pass

try:
    from . import princeton
except:
    pass

from . import files
from . import pbslens

from . import binning
from . import outputs
from . import correct


from . import fit
from . import nfw
from . import linear
from . import project
from . import invert

from . import plotting

from . import testing

