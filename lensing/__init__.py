"""
First you need a catalog class for each catalog type.
See 
    lensing.lcat.DESMockLensCatalog
    lensing.scat.DESMockSrcCatalog

Then add an entry in the if statement for these functions and run them:

    lensing.lcat.create_input(catalog, version, sample, nsplit=)
    lensing.scat.create_input(catalog, version, sample)

Then create a lensing run json in $ESPY_DIR/lensing/config/ and run
the config file and pbs file creators

    lensing.config.write_config(run)
    lensing.pbslens.write_pbs(run)

Make sure to install objshear under ~/exports/objshear-work
    scons install prefix=~/exports/objshear-work

Then you can submit the pbs files which are under
    ${LENSDIR}/pbslens/{run}/


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
import sigmacrit
import test_estimators
import lcat
import scat
import convert
import config
import regauss
import regauss_sim
import princeton
import files
import pbslens

import outputs


import fit
import nfw
import linear
import project
import invert

import plotting

import testing
