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


UPDATE THIS FOR NEW SOURCE SPLITS

You should then bin up the results. e.g.

An example is:
    lensing.outputs.mzbin_byrun

Should make a more general binner like we had in IDL.

Then fits and plots:
    
    lensing.outputs.plot_mzbin(run, nmass, nz)

    lensing.outputs.nfw_lin_fits_byrun(run, name, rrange=None)
    lensing.outputs.plot_nfw_lin_fits_byrun(run, name, prompt=False)
    lensing.outputs.plot_mzbin_mass_byrun(run, nmass, nz, fudge=1.0,tag='m200',
                                          epsfile=None)

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
