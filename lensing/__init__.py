"""

quick version
-------------

I'll go through each below, but for now assume you have a lens run called l07
and a randoms run called r01, and the config files are set up, and scinv
exists.

    # if the scat are not already generated
    /bin/make-objshear-input.py scat scat_sample

    #
    # now generate the lens samples and processing scripts, condor
    #
    /bin/make-objshear-input.py lcat l07
    /bin/make-objshear-proc.py l07

    # in the $LENSDIR/proc/run directory
    condor_submit run-l07.condor

    /bin/reduce-lensout.py l07
    /bin/collate-reduced.py l07

    # bin by lambda into 12 bins, must have defined this binnign
    /bin/bin-lenses.py l07 lambda 12
    # plot the binning (not corrected, jackknife yet)
    /bin/plot-dsig-byrun.py -t binned l07 lambda 12

    #
    # now randoms
    #

    # generation of randoms is slow, so do it in chunks, here 100
    /bin/make-random-chunk-scripts.py r01 100

    # scripts are in the $LENSDIR/lcat/r01/condor
    for f in *.condor; do echo "$f"; condor_submit "$f"; done

    # then combine them
    /bin/combine-random-chunks.py r01 100

    # now make the scripts, condor
    /bin/make-objshear-proc.py r01

    # in the $LENSDIR/proc/r01 directory
    condor_submit run-r01.condor

    # and combine same as lenses
    /bin/reduce-lensout.py r01
    /bin/collate-reduced.py r01


    # match randoms to those bins
    /bin/match-randoms.py -t lambda -n 12 l07 r01

    # correct from randoms
    /bin/correct-shear.py

Detailed version
----------------

First you need a catalog class for each catalog type.
See 

    lensing.lcat.RedMapper
    lensing.scat.DR8Catalog

And an associated yaml config file with a sample name.

Then you must make sure you have source scinv in the
sweeps_reduct/regauss/04.cols by running

    /bin/add-scinv.py

Then create input catalogs. 

    /bin/make-objshear-input.py scat scat_sample
    /bin/make-objshear-input.py lcat lcat_sample

For random lcat this can be *very* slow so be prepared.
Note this requires dealing with catalog names in these

    lensing.lcat.create_input(sample)
    lensing.scat.create_input(sample)

Then create a lensing run yaml in $ESPY_DIR/lensing/config/ and run
the config, script, and condor file creators using 

    /bin/make-objshear-proc.py

The files go under ~/lensing/proc/run

Make sure to install objshear under ~/exports/objshear-work
    python build.py --prefix=~/exports/objshear-work install


Then submit the scripts

    condor_submit run-04-035.condor

Since the results are split by source chunks, you then need
to "reduce" the lens outputs, and collate with the original
catalog

    /bin/reduce-lensout.py
    /bin/collate-reduced.py

You can then bin up the results. Youll need to use a binner
from binning.py or make a new one.

    /bin/bin-lenses.py  type run nbin

    # also plots
    binner = N200Binner(nbin)

    binner.plot_dsig_byrun(run)

Then you need to run some randoms for corrections.  This involves
makine a new lcat.  Currently we have

    lcat.SDSSRandom

Then run as usual, reduce, collate, etc.

You must match the randoms to lenses using

    /bin/match-randoms.py

currently only for binned samples.  Then corrections are made with

    /bin/correct-shear.py

although this needs work.


Make some plots

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


from . import regauss
from . import regauss_test
from . import regauss_sim

from . import princeton
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

