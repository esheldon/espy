"""
First you need a catalog class for each catalog type.
See 
    lensing.lcat.DESMockLensCatalog
    lensing.scat.DESMockSrcCatalog

Then add an entry in the if statement for these functions and run them:

    lensing.lcat.create_input(catalog, version, sample, nsplit=)
    lensing.scat.create_input(catalog, version, sample)

You can run those with the script

    /bin/make-objshear-input.py

Then create a lensing run json in $ESPY_DIR/lensing/config/ and run
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

    lensing.outputs.make_reduced_lensout(run)
    lensing.outputs.make_collated_lensred(run)


You can then bin up the results. e.g.
    binner = N200Binner(nbin)
    binner.bin_byrun(run)

    # add clustering/Ssh corrections here

    binner.plot_dsig_byrun(run)

Should make a more general binner like we had in IDL.

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
from . import test_estimators
from . import lcat
from . import scat
from . import convert
from . import objshear_config
from . import scripts
from . import condor


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

