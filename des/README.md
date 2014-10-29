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

you need to define a source sample, e.g.

    $ESPY_DIR/des/config/scat-001.yaml

which will refer to an scat_name that points under the catalogs directory.  It
also refers to the pz_vers and pz_type and cosmology.

match the source catalog in *tiles* not the chunks above.  This is assuming
the source catalogs are broken up by tilename

    $ESPY_DIR/des/bin/match-scinv $scat_vers
    $ESPY_DIR/des/bin/match-scinv scat-001

This is high memory, so I've been running it as a single job rather than
splitting things up.
