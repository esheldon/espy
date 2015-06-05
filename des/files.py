"""
getting file paths and some reading/writing

We need to put these all in one place to avoid import
confusion
"""

from __future__ import print_function
import os
import numpy
import fitsio

#
# configuration files
#

def config_dir():
    """
    dir holding config files
    """
    d=os.environ['ESPY_DIR']
    d=os.path.join(d, 'des', 'config')
    return d

def get_config_file(name):
    """
    config file name
    """
    d=config_dir()
    name='%s.yaml' % name

    return os.path.join(d, name)

def read_config(name):
    """
    read the config file
    """
    import yaml
    fname=get_config_file(name)
    #print("reading:",fname)
    with open(fname) as fobj:
        data=yaml.load(fobj)
    return data

def cascade_config(run):
    """
    read the run config as well as the lens, source, and cosmo
    configs
    """
    conf=read_config(run)
    conf['run']=run

    conf['lens_conf']=read_config(conf['lcat_vers'])
    conf['source_conf']=read_config(conf['scat_vers'])
    conf['cosmo_conf']=read_config(conf['lens_conf']['cosmo_vers'])
    conf['mask_conf']=read_config(conf['lens_conf']['mask_vers'])

    lc=conf['lens_conf']['cosmo_vers']
    sc=conf['source_conf']['cosmo_vers']
    if lc != sc:
        raise ValueError("cosmo mismatch: '%s' '%s'" % (lc,sc))

    lc=conf['lens_conf']['mask_vers']
    sc=conf['source_conf']['mask_vers']
    if lc != sc:
        raise ValueError("mask version: '%s' '%s'" % (lc,sc))

    return conf

#
# basic lens directory structure
#

def get_lensdir():
    """
    This is the root dir under which all data are stored
    """
    if 'LENSDIR' not in os.environ:
        raise ValueError("LENSDIR is not set")
    return os.environ['LENSDIR']

def get_des_lensdir():
    """
    root for des-specific lensing dir
    """
    d=get_lensdir()
    return os.path.join(d, 'des-lensing')

#
# Catalog files
#

def get_cat_basedir():
    """
    base dir holding all catalog dirs
    """
    catdir = os.path.join(get_lensdir(), 'catalogs')
    return catdir

def get_cat_dir(catname):
    """
    dir holding a catalog type
    """
    d=get_cat_basedir()
    d=os.path.join(d, catname)
    return d

def get_lcat_original_file(lcat_name):
    """
    get the original input file

    We always use a symlink to enforce uniform names
    """
    dir=get_cat_dir(lcat_name)
    return os.path.join(dir, '%s.fits' % lcat_name)

def read_lcat_original(lcat_name):
    """
    read the original input file
    """
    fname=get_lcat_original_file(lcat_name)
    print("reading:",fname)
    return fitsio.read(fname, lower=True)

#
# photoz related files
#

def get_pz_base_dir():
    return '/astro/u/astrodat/data/DES/EXTRA/photoz/combined'

def get_pz_vers_dir(pz_vers):
    d=get_pz_base_dir()
    return os.path.join(d, pz_vers)

def get_pz_version_info_file():
    dir=get_pz_base_dir()
    name='version-info.yaml'
    return os.path.join(dir, name)

def read_pz_vers_info():
    import yaml
    fname=get_pz_version_info_file()
    with open(fname) as fobj:
        data=yaml.load(fobj)
    return data

def get_pz_h5_file(pz_vers):
    dir=get_pz_vers_dir(pz_vers)
    name='DES_photoz_PDFS_%s.h5' % pz_vers
    return os.path.join(dir, name)

#
# sigmacrit files
#

def get_scinv_dir(pz_vers, pz_type, cosmo_vers):
    """
    sigmacrit files for the give p(z) and cosmology
    """
    dir=get_pz_vers_dir(pz_vers)
    return os.path.join(dir, 'scinv-%s-%s' % (pz_type,cosmo_vers))


def get_scinv_file(pz_vers, pz_type, cosmo_vers, chunk=None):
    """
    sigmacrit files for the give p(z) and cosmology
    """
    dir=get_scinv_dir(pz_vers,pz_type,cosmo_vers)
    name='DES_scinv_%s_%s_%s' % (pz_vers, pz_type,cosmo_vers)

    if chunk is not None:
        dir=os.path.join(dir,'chunks')
        name='%s_%06d' % (name,chunk)
    name='%s.fits' % name
    return os.path.join(dir, name)

def read_scinv(pz_vers, pz_type, cosmo_vers, chunk=None, get_header=False):
    """
    read the sigmacrit file for the give p(z) and cosmology
    """
    fname=get_scinv_file(pz_vers, pz_type, cosmo_vers, chunk=None)

    print("reading:",fname)
    with fitsio.FITS(fname) as fits:
        data=fits['scinv'][:]
        zlvals=fits['zlvals'][:]
        print("    read:",data.size)

        if get_header:
            h=fits['scinv'].read_header()
            ret=data, zlvals, h
        else:
            ret=data, zlvals

    return ret

def get_scinv_wq_dir(pz_vers, pz_type, cosmo_vers):
    """
    wq submit files for calculating sigmacrit
    """
    d=get_scinv_dir(pz_vers, pz_type, cosmo_vers)
    return os.path.join(d, 'wq')

def get_scinv_wq_file(pz_vers, pz_type, cosmo_vers, chunk):
    """
    wq submit files for calculating sigmacrit
    """
    dir=get_scinv_wq_dir(pz_vers, pz_type, cosmo_vers)
    name='DES_scinv_%s_%s_%s_%06d.yaml' % (pz_vers, pz_type, cosmo_vers, chunk)
    return os.path.join(dir, name)


#
# scat files
#

# these are more generic
def get_orig_scat_file(scat_name, tilename):
    """
    original scat files, split by tile

    parameters
    ----------
    scat_name: string
        e.g. ngmix011-v14
    tilename: string
        the tilename
    """
    if '-dg' in scat_name:
        fname=get_dg_scat_file(scat_name, tilename)
    else:
        d=get_cat_dir(scat_name)
        fname='{scat_name}-{tilename}.fits'.format(scat_name=scat_name,
                                                   tilename=tilename)
        fname=os.path.join(d, fname)

    return fname

def read_orig_scat(scat_name, tilename):
    """
    read scat files, split by tile

    parameters
    ----------
    scat_name: string
        e.g. ngmix011-v14
    tilename: string
        the tilename
    """

    fname=get_orig_scat_file(scat_name, tilename)
    print("reading:",fname)
    return fitsio.read(fname)


dg_name={'ngmix009':'{tilename}_{scat_name}m.fits.gz',
         'ngmix010':'{tilename}_im3shapev72_ngmix_009_010.fits.gz'}

def get_dg_scat_file(scat_name, tilename):
    """
    daniel's matched files
    """
    d=get_cat_dir(scat_name+'-dg')
    pattern=dg_name[scat_name]

    fname=pattern.format(scat_name=scat_name, tilename=tilename)
    return os.path.join(d, fname)

def read_dg_scat(scat_name, tilename):
    """
    daniel's matched files
    """
    fname=get_dg_scat_file(scat_name, tilename)
    print("reading:",fname)
    return fitsio.read(fname, lower=True)

def get_orig_scat_file_full(scat_name):
    """
    original file, all tiles combined

    parameters
    ----------
    scat_name: string
        e.g. ngmix011-v14
    """
    d=get_cat_dir(scat_name)
    fname='{scat_name}.fits'.format(scat_name=scat_name)
    return os.path.join(d, fname)

def read_orig_scat_full(scat_name):
    """
    read original file, all tiles combined

    parameters
    ----------
    scat_name: string
        e.g. ngmix011-v14
    """

    fname=get_orig_scat_file_full(scat_name)
    print("reading:",fname)
    return fitsio.read(fname)



#
# original source files matched to the scinv outputs
#

def get_scinv_matched_dir(scat_name, pz_vers, pz_type, cosmo_vers):
    """
    dir to hold the scinv matched files
    """
    d=get_cat_basedir()
    pattern='{scat_name}-{pz_vers}-{pz_type}-{cosmo_vers}-match'
    sub_dir=pattern.format(scat_name=scat_name,
                           pz_vers=pz_vers,
                           pz_type=pz_type,
                           cosmo_vers=cosmo_vers)
    return os.path.join(d, sub_dir)

def get_scinv_matched_file(scat_name, pz_vers, pz_type, cosmo_vers, tilename):
    """
    source catalog matched to scinv file
    """
    d=get_scinv_matched_dir(scat_name, pz_vers, pz_type, cosmo_vers)

    pattern='{tilename}-{scat_name}-{pz_vers}-{pz_type}-{cosmo_vers}-match.fits'
    fname=pattern.format(tilename=tilename,
                         scat_name=scat_name,
                         pz_vers=pz_vers,
                         pz_type=pz_type,
                         cosmo_vers=cosmo_vers)

    return os.path.join(d, fname)

def read_scinv_matched(scat_name, pz_vers, pz_type, cosmo_vers, tilename):
    """
    read the matched file
    """
    pass

#
# source catalog files for input to xshear
#

def get_scat_dir(scat_vers):
    """
    directory to hold source catalog files
    """
    d=get_des_lensdir()
    return os.path.join(d, 'scat', scat_vers)

def get_scat_file(scat_vers, tilename):
    """
    source catalog ascii files 
    """
    d=get_scat_dir(scat_vers)
    fname='{scat_vers}-{tilename}.dat'
    fname=fname.format(scat_vers=scat_vers,
                       tilename=tilename)
    return os.path.join(d,fname)

def write_scat(fname, data):
    """
    Write a source ascii file
    """
    from esutil.recfile import Recfile
    if os.path.exists(fname):
        os.remove(fname)
    with Recfile(fname,'w',delim=' ') as robj:
        robj.write(data)

#
# lens catalogs for input to xshear
#

def get_lcat_dir(lcat_vers):
    """
    directory for a lens catalog
    """
    d=get_des_lensdir()
    return os.path.join(d, 'lcat', lcat_vers)

def get_lcat_file_dir(lcat_vers):
    """
    directory holding the data
    """
    d=get_lcat_dir(lcat_vers)
    return os.path.join(d, 'data')

def get_lcat_file(lcat_vers, chunk):
    """
    lens catalog ascii files 
    """
    d=get_lcat_file_dir(lcat_vers)
    fname='%(lcat_vers)s-%(chunk)06d.dat'
    fname=fname % {'lcat_vers':lcat_vers, 'chunk':chunk}
    return os.path.join(d,fname)

def write_lcat(fname, data):
    """
    Write a source ascii file
    """
    from esutil.recfile import Recfile

    d=os.path.dirname(fname)
    if not os.path.exists(d):
        print("making dir:",d)
        os.makedirs(d)

    if os.path.exists(fname):
        os.remove(fname)

    print("writing:",fname)
    with Recfile(fname,'w',delim=' ') as robj:
        robj.write(data)


#
# run configuration
#

def get_run_basedir():
    """
    lensdir/run
    """
    d=get_des_lensdir()
    return os.path.join(d, 'run')

def get_run_dir(run):
    """
    lensdir/run/run_name
    """
    d=get_run_basedir()
    return os.path.join(d, run)

def get_xshear_config_dir(run):
    """
    directory holding the config file
    """
    d=get_run_dir(run)
    return os.path.join(d, 'config')

def get_xshear_config_file(run):
    """
    lensdir/run/{run_name}/{run_name}.cfg

    parameters
    ----------
    run:
        e.g. run-001
    """
    d=get_xshear_config_dir(run)
    fname='%s.cfg' % run
    return os.path.join(d, fname)

#
# xshear output files
#

def get_output_dir(run, lens_chunk):
    """
    lensdir/run/run_name
    """
    d=get_run_dir(run)
    return os.path.join(d, 'output', 'lens-%06d' % lens_chunk)

def get_output_file(run, lens_chunk, source_tilename):
    """
    the xshear output file for a lens chunk and source tilename
    """
    d=get_output_dir(run, lens_chunk)
    fname="%(run)s-lens-%(lens_chunk)06d-src-%(source_tilename)s.dat"
    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk,
                   'source_tilename':source_tilename}

    return os.path.join(d, fname)

def get_shear_style(data):
    """
    determine which shear style is used from the column names
    """
    if 'dsensum' in data.dtype.names:
        shear_style='lensfit'
    else:
        shear_style='reduced'
    return shear_style


def get_lensum_dtype(nbin, shear_style):
    """
    get the dtype for the number of bins and shear type
    """
    dt=[('index','i8'),
        ('weight','f8'),
        ('totpairs','i8'),

        ('npair','i8',nbin),

        ('rsum','f8',nbin),
        ('wsum','f8',nbin),
        ('dsum','f8',nbin),
        ('osum','f8',nbin)]

    if shear_style=='lensfit':
        dt+=[('dsensum','f8',nbin),
             ('osensum','f8',nbin)]

    return dt

def read_lensum(fname, nbin, shear_style):
    """
    read a generic lensum file. Used by the other readers

    parameters
    ----------
    filneame:
        location of file
    nbin: int
        number of radial bins
    shear_style: string
        Determine the columns that must be read
    """
    from esutil.recfile import Recfile

    dt=get_lensum_dtype(nbin, shear_style)

    print("reading:",fname)
    with Recfile(fname, 'r', dtype=dt, delim=' ') as robj:
        data=robj.read()

    return data

def _read_run_lensum(fname, run):
    """
    inernal routine to read a lensum, determining the nbin
    and shear_style from the run
    """
    conf=cascade_config(run)

    nbin=conf['lens_conf']['nbin']
    shear_style=conf['source_conf']['shear_style']

    return read_lensum(fname, nbin, shear_style)


def read_output(run, lens_chunk, source_tilename):
    """
    read the single output files

    parameters
    ----------
    run: string
        the run identifier, e.g. run-001
    lens_chunk: int
        lens chunk number
    source_tilename:
        the source tilename
    """

    fname=get_output_file(run, lens_chunk, source_tilename)
    return _read_run_lensum(fname, run)


#
# xshear reduced files
#

def get_reduced_dir(run):
    """
    lensdir/run/run_name/rediced
    """
    d=get_run_dir(run)
    return os.path.join(d, 'reduced')


def get_reduced_file(run, lens_chunk):
    """
    File reduced across sources
    """
    d=get_reduced_dir(run)
    if lens_chunk=='*':
        fname="%(run)s-lens-*-reduced.dat"
    else:
        fname="%(run)s-lens-%(lens_chunk)06d-reduced.dat"

    fname=fname % {'run':run,
                   'lens_chunk':lens_chunk}

    return os.path.join(d, fname)

def read_reduced(run, lens_chunk):
    """
    read a file for outputs reduced over sources

    parameters
    ----------
    run: string
        the run identifier, e.g. run-001
    lens_chunk: int
        lens chunk number
    """

    fname=get_reduced_file(run, lens_chunk)
    return _read_run_lensum(fname, run)

#
# xshear combined files (combine all lens chunks)
#

def get_combined_dir(run):
    """
    lensdir/run/run_name/combined
    """
    d=get_run_dir(run)
    return os.path.join(d, 'combined')

def get_combined_file(run):
    """
    File combined for all lenses
    """
    d=get_combined_dir(run)
    fname="%s-combined.dat" % run
    return os.path.join(d, fname)

def read_combined(run):
    """
    read a file for all lens chunks combined

    parameters
    ----------
    run: string
        the run identifier, e.g. run-001
    """

    fname=get_combined_file(run)
    return _read_run_lensum(fname, run)

#
# the collated file
#

def get_collated_dir(run):
    """
    lensdir/run/run_name/collated
    """
    d=get_run_dir(run)
    return os.path.join(d, 'collated')

def get_collated_file(run):
    """
    File collated for all lenses
    """
    d=get_collated_dir(run)
    fname="%s-collated.fits" % run
    return os.path.join(d, fname)

def read_collated(run):
    """
    read a file for all lens chunks collated

    parameters
    ----------
    run: string
        the run identifier, e.g. run-001
    """
    fname=get_collated_file(run)
    print("reading:",fname)
    return fitsio.read(fname)

#
# binned files
#

def get_binned_dir(run, bin_scheme):
    """
    lensdir/run/run_name/binned/bin_name
    """
    d=get_run_dir(run)
    return os.path.join(d, 'binned/%s' % bin_scheme)

def get_binned_file(run, bin_scheme, ext='fits'):
    """
    get the file holding binned data, or the basic plot file

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    ext: string
        default fits, could be eps etc.
    """

    d=get_binned_dir(run, bin_scheme)

    fname="%s-%s.%s" % (run, bin_scheme, ext)
    return os.path.join(d, fname)

def read_binned(run, bin_scheme):
    """
    read the file holding binned data

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    """

    fname=get_binned_file(run,bin_scheme)
    print("reading:",fname)
    return fitsio.read(fname)

#
# weights for randoms
#

def get_match_dir(lens_run, rand_run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname
    """
    totrun='%s-%s' % (lens_run, rand_run)
    d=get_run_dir(totrun)


def get_match_weights_dir(lens_run, rand_run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname/weights
    """

    d=get_match_binned_dir(lens_run, rand_run, bin_scheme)
    return os.path.join(d, 'weights')

def get_match_weights_file(lens_run, rand_run, bin_scheme, binnum=None, ext='fits'):
    """
    get the file holding weights for the random-matched data

    parameters
    ----------
    lens_run: string
        the run identifier
    rand_run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    binnum: int
        bin number for weights file
    ext: string
        default fits, could be eps etc.
    """

    d=get_match_weights_dir(lens_run, rand_run, bin_scheme)

    totrun='%s-%s' % (lens_run, rand_run)
    if binnum is not None:
        fname="%s-%s-%02d-weights.%s" % (totrun, bin_scheme, binnum, ext)
    else:
        fname="%s-%s-weights.%s" % (totrun, bin_scheme, ext)
    return os.path.join(d, fname)


# these are just convenience to combine the run names
def get_match_binned_dir(lens_run, rand_run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname/binned
    """

    totrun='%s-%s' % (lens_run, rand_run)
    return get_binned_dir(totrun, bin_scheme)

def get_match_binned_file(lens_run, rand_run, bin_scheme, ext='fits'):
    """
    get the file holding binned data, or the basic plot file

    parameters
    ----------
    lens_run: string
        the run identifier
    rand_run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    ext: string
        default fits, could be eps etc.
    """

    totrun='%s-%s' % (lens_run, rand_run)
    return get_binned_file(totrun, bin_scheme, ext=ext)

def read_match_binned(lens_run, rand_run, bin_scheme):
    """
    read the file holding binned data, or the basic plot file

    parameters
    ----------
    lens_run: string
        the run identifier
    rand_run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    """

    fname=get_match_binned_file(lens_run, rand_run, bin_scheme)
    print("reading:",fname)
    return fitsio.read(fname)
#
# jackknifed here
#

def get_jack_dir(run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname/jack
    """

    d=get_run_dir(run)
    return os.path.join(d, 'jack/%s' % bin_scheme)

def get_jack_file(run, bin_scheme, ext='fits'):
    """
    get the file holding jackknifed data, or the basic plot file

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    ext: string
        default fits, could be eps etc.
    """

    d=get_jack_dir(run, bin_scheme)

    fname="%s-%s-jack.%s" % (run, bin_scheme, ext)
    return os.path.join(d, fname)

def read_jack(run, bin_scheme):
    """
    read file holding jackknifed data, or the basic plot file

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    """
    fname=get_jack_file(run, bin_scheme)
    print("reading:",fname)
    return fitsio.read(fname)

#
# corrected files here
#

def get_corr_binned_dir(lens_run, rand_run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname/corr-binned
    """

    totrun='%s-%s' % (lens_run, rand_run)
    d=get_run_dir(totrun)
    return os.path.join(d, 'corr-binned/%s' % bin_scheme)

def get_corr_binned_file(lens_run, rand_run, bin_scheme, ext='fits'):
    """
    get the file holding corrected binned data, or the basic plot file

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    ext: string
        default fits, could be eps etc.
    """

    d=get_corr_binned_dir(lens_run, rand_run, bin_scheme)

    totrun='%s-%s' % (lens_run, rand_run)
    fname="%s-%s-corr.%s" % (totrun, bin_scheme, ext)
    return os.path.join(d, fname)

def read_corr_binned(lens_run, rand_run, bin_scheme):
    """
    read the file holding binned data, or the basic plot file

    parameters
    ----------
    lens_run: string
        the run identifier
    rand_run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    """

    fname=get_corr_binned_file(lens_run, rand_run, bin_scheme)
    print("reading:",fname)
    return fitsio.read(fname)

def get_corr_jack_dir(lens_run, rand_run, bin_scheme):
    """
    lensdir/run/lrun_name-rrunname/corr-jack
    """

    totrun='%s-%s' % (lens_run, rand_run)
    d=get_run_dir(totrun)
    return os.path.join(d, 'corr-jack/%s' % bin_scheme)

def get_corr_jack_file(lens_run, rand_run, bin_scheme, ext='fits'):
    """
    get the file holding corrected jackknifed data, or the basic plot file

    parameters
    ----------
    run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    ext: string
        default fits, could be eps etc.
    """

    d=get_corr_jack_dir(lens_run, rand_run, bin_scheme)

    totrun='%s-%s' % (lens_run, rand_run)
    fname="%s-%s-jack-corr.%s" % (totrun, bin_scheme, ext)
    return os.path.join(d, fname)


def read_corr_jack(lens_run, rand_run, bin_scheme):
    """
    read the file holding corrected jackknifed data, or the basic plot file

    parameters
    ----------
    lens_run: string
        the run identifier
    rand_run: string
        the run identifier
    bin_scheme: string
        name for the binning scheme, e.g. bin-lambda08-z01
    """

    fname=get_corr_jack_file(lens_run, rand_run, bin_scheme)
    print("reading:",fname)
    return fitsio.read(fname)



#
# jackknife region center files
#

def get_jackknife_centers_dir(des_region):
    """
    directory holding center files
    """
    d=get_des_lensdir()
    return os.path.join(d, 'jackknife-regions', des_region)

def get_jackknife_centers_file(des_region, ncen):
    """
    file holding centers for the specified region
    """
    d=get_jackknife_centers_dir(des_region)

    fname='jackknife-%s-%06d.fits' % (des_region, ncen)
    return os.path.join(d, fname)

def get_jackknife_centers_epsfile(des_region, ncen, extra=None):
    """
    file holding centers for the specified region
    """
    d=get_jackknife_centers_dir(des_region)

    fname=['jackknife',des_region,'%06d' % ncen]
    if extra is not None:
        fname += [extra]

    fname='-'.join(fname)
    fname = '%s.eps' % fname

    return os.path.join(d, fname)

def read_jackknife_centers(des_region, ncen):
    """
    read the file holding centers for the specified region
    """
    fname=get_jackknife_centers_file(des_region, ncen)
    print("reading:",fname)
    return fitsio.read(fname)

#
# cacheing tilenames from the db, since they are not stored
# in the flat catalogs
#

def get_tilenames(tablename, full=False):
    """
    parameters
    ----------
    tablename: string
        table where data is read originally
    full: bool
        If True, return the full data with coadd_objects_id
        and tilename for each object
    """

    print("getting tile list for:",tablename)

    data=load_tilenames_cache(tablename)

    if full:
        return data
    else:
        tilenames=numpy.unique(data['tilename'])
        print("found",tilenames.size,"tiles")
        return tilenames

def load_tilenames_cache(tablename):
    """
    write the tilenames to a file

    parameters
    ----------
    tablename: string
        table where data is read originally
    full: bool
        If True, return the full data with coadd_objects_id
        and tilename for each object
    """
    import desdb
    fname=get_tilename_cache_file(tablename)

    if not os.path.exists(fname):
        print("    caching from database")
        dir=get_tilename_cache_dir()
        if not os.path.exists(dir):
            print("making dir:",dir)
            os.makedirs(dir)

        with desdb.Connection() as conn:
            q="""
    select
        coadd_objects_id, tilename
    from
        %s\n""" % tablename

            print(q)
            data=conn.quick(q, array=True)
            
        print("    writing to:",fname)
        fitsio.write(fname, data, clobber=True)
    else:
        print("    reading:",fname)
        data=fitsio.read(fname)

    return data


def get_tilename_cache_dir():
    """
    directory to hold cache
    """
    d=get_des_lensdir()
    return os.path.join(d, 'tilename_cache')

def get_tilename_cache_file(tablename):
    """
    file to hold tilenames
    """
    dir=get_tilename_cache_dir()
    fname='%s-tilenames.fits' % tablename
    return os.path.join(dir, fname)


