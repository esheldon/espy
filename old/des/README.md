process for lensing work
========================

cosmology
---------
before anything you need the cosmology defined, e.g.

    $ESPY_DIR/des/config/cosmo-01.yaml

p(z) and scinv
--------------

First create the scinv from the p(z). These are made in chunks.  The relevant
information is the cosmology, the pz_vers and pz_type, which currently refer to
christophers hdf5 files, and a chunksize.

Create the wq scripts

    $ESPY_DIR/des/bin/make-scinv-wq $cosmo_vers $pz_vers $pz_type $chunksize
    $ESPY_DIR/des/bin/make-scinv-wq cosmo-01 v0.1 skynet_mag_auto 100000

These run the make-scinv script.  Then combine the outputs

    $ESPY_DIR/des/bin/combine-scinv $cosmo_vers $pz_vers $pz_type $chunksize

source catalogs
---------------

First there are the original catalogs, currently only the Daniel Gruen matched
catalogs, e.g. ngmix009-dg.

### when using point z

In this case we just read the scat directly and write it, so skip to the
next step

### when using full p(z)

We want to match these to the scinv calculated above

match the source catalog to the scinv in *tiles* not the chunks above.  This is
assuming the source catalogs are broken up by tilename

    $ESPY_DIR/des/bin/match-scinv $scat_name $cosmo_vers $pz_vers $pz_type
    $ESPY_DIR/des/bin/match-scinv ngmix009 cosmo-01 v0.1.5 tpz

This is high memory, so I've been running it as a single job rather than
splitting things up.

input source catalogs for xshear
--------------------------------

you need to define a source sample, e.g.

    $ESPY_DIR/des/config/scat-ngmix004.yaml

which will refer to an scat_name that points under the catalogs directory,
as well as p(z) version and type, and cosmology if using p(z)

Then make the source input for xshear.  This is pretty fast, no need currently
to split it up

    $ESPY_DIR/des/bin/make-xshear-scat $scat_vers
    $ESPY_DIR/des/bin/make-xshear-scat scat-ngmix004
    $ESPY_DIR/des/bin/make-xshear-scat -t DES0440-4623 scat-001

lens catalogs
---------------

you need to define a lens sample, e.g.

    $ESPY_DIR/des/config/lcat-001.yaml

which will refer to an lcat_name that points under the catalogs directory.  It
also refers to the cosmology.

Create the lens xshear input file

    $ESPY_DIR/des/bin/make-xshear-lcat $lcat_vers
    $ESPY_DIR/des/bin/make-xshear-lcat lcat-001
    $ESPY_DIR/des/bin/make-xshear-lcat --chunk 35 lcat-001

For randoms and large catalogs you will want to process each chunk
separately.  Make wq files to make the chunks.

    $ESPY_DIR/des/bin/make-lcat-wq $lcat_vers


xshear runs
-----------

define a run, e.g.

    $ESPY_DIR/des/config/run-rm004.yaml

then write the xshear config and all wq files

    $ESPY_DIR/des/config/make-run-wq

then submit the xshear wq scripts

    $LENSDIR/des-lensing/run/{run-name}/wq-xshear/lens-{chunk}/*.yaml 
    e.g.
    incsub -n 200 $(find . -name "*.yaml")

the reduction script

    $LENSDIR/des-lensing/run/{run-name}/wq-redshear/*.yaml 

the combination of all lense splits

    $LENSDIR/des-lensing/run/{run-name}/wq-combine/*.yaml 

and finally collation

    $LENSDIR/des-lensing/run/{run-name}/wq-collat/*.yaml 

binning
-------

Define a bin scheme, e.g. 

    $ESPY_DIR/des/config/bin-lambda01-z01.yaml

Then perform the binning


    $ESPY_DIR/des/bin/bin-lenses run bin_scheme
    $ESPY_DIR/des/bin/bin-lenses --jack run bin_scheme

matching random point lens runs
-------------------------------

    $ESPY_DIR/des/bin/match-rand lens_run rand_run bin_scheme

calculating boost factors
-------------------------

    $ESPY_DIR/des/bin/boost-correct lens_run rand_run bin_scheme

plotting the boost factors
---------------------------

    $ESPY_DIR/des/bin/plot-boost lens_run rand_run bin_scheme

generating your own randoms
---------------------------

First make a random config.  Then generate nrand randoms
and write to catalog cat_name

    $ESPY_DIR/des/bin/make-randoms rand_conf nrand cat_name

plots
-----

Plot a binned set of lenses

    $ESPY_DIR/des/bin/plot-dsig run bin_scheme

Plot ratio between two runs

    $ESPY_DIR/des/bin/plot-ratio run1 run2 bin_scheme
