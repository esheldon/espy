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


You should then bin up the results. e.g.
    mzbin = lensing.binning.MZBinner(nmass, nz)
    mzbin.bin_byrun(run)
which calls this under the hood
    mzbin.bin(data)

Should make a more general binner like we had in IDL.

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

Some plots that are on a grid are done currently from the binner:
    

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
 
Others done one-by-one are done in fit.py

    lensing.outputs.plot_nfw_lin_fits_byrun(run, name, prompt=False)

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


from . import fit
from . import nfw
from . import linear
from . import project
from . import invert

from . import plotting

from . import testing
