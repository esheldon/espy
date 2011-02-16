import os
from sys import stdout
import lensing
import numpy

import esutil as eu
from esutil.ostools import path_join

finfo={}
finfo['config']={'ext':'.dat'}
finfo['lcat']={'ext':'.bin'}
finfo['scat']={'ext':'.bin'}
finfo['lensout']={'ext':'.bin'}
finfo['pbslens']={'ext':'.pbs'}


def lensdir():
    if 'LENSDIR' not in os.environ:
        raise ValueError("LENSDIR is not set")
    return os.environ['LENSDIR']

def catalog_dir():
    """
    The root of the catalog directory. ${LENSDIR}/catalogs
    """
    catdir = path_join(lensdir(), 'catalogs')
    return catdir


#
# generic for lcat,scat,config,pbs
#
def sample_dir(type,sample):
    if type not in finfo:
        raise ValueError("Unknown file type: '%s'" % type)
    d = lensdir()
    d = path_join(d,type,sample)
    return d

def sample_file(type,sample,split=None):
    d = sample_dir(type,sample)
    fname=path_join(d, '%s-%s' % (type,sample))
    if split is not None:
        fname += '-%03d' % split
    fname += finfo[type]['ext']
    return fname


#
# fits. Go in same dir as bins
#

def lensfit_write(data, run, name, extra=None, verbose=True):
    d = lensfit_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    f = lensfit_file(run,name,extra=extra)
    eu.io.write(f,data,verbose=verbose)

def lensfit_print(data):
    f = '%0.7f %0.7f'
def lensfit_read(run, name, extra=None,verbose=True):
    f = lensfit_file(run,name,extra=extra)
    return eu.io.read(f,verbose=verbose)

def lensfit_dir(run,name):
    return lensbin_dir(run,name)
def lensfit_file(run,name,extra=None):
    d=lensfit_dir(run,name)
    if extra is not None:
        ex='-'+extra
    else:
        ex=''
    return path_join(d, 'lensfit-%s-%s%s.rec' % (run,name,ex))


#
# inversions
#

def lensinv_write(data, run, name, verbose=True):
    d = lensinv_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    f = lensinv_file(run,name)
    eu.io.write(f,data,verbose=verbose)

def lensinv_read(run, name, verbose=True):
    f = lensinv_file(run,name)
    return eu.io.read(f,verbose=verbose)

def lensinv_dir(run,name):
    return lensbin_dir(run,name)
def lensinv_file(run,name):
    d=lensinv_dir(run,name)
    return path_join(d, 'lensinv-%s-%s.rec' % (run,name))



# binned outputs
#
def lensbin_write(data, run, name, verbose=True):
    d = lensbin_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    f = lensbin_file(run,name)
    eu.io.write(f,data,verbose=verbose)

def lensbin_read(run, name, verbose=True):
    f = lensbin_file(run,name)
    return eu.io.read(f,verbose=verbose)

def lensbin_dir(run,name):
    """
    e.g. lensbin_dir('02','m12z3')
    """
    dir = sample_dir('lensout',run)
    return path_join(dir,'lensbin-'+name)
def lensbin_file(run,name):
    d=lensbin_dir(run,name)
    return path_join(d, 'lensbin-%s-%s.rec' % (run,name))

def lensbin_plot_dir(run,name):
    d = lensbin_dir(run,name)
    d=path_join(d,'plots')
    return d

#
# read objshear outputs
#

def lensout_collate(run):
    """
    Collate the lensing output with the original catalog
    """

    conf = json_read(run)
    cat = lensing.lcat.read_catalog(conf['lens_catalog'],conf['lens_version'])
    lout = lensout_read(run=run)

    if cat.size != lout.size:
        raise ValueError("catalog and lensout are not the same size")


    # create a new struct with out outputs at the end
    stdout.write('collating...')
    data = eu.numpy_util.add_fields(cat, lout.dtype)
    eu.numpy_util.copy_fields(lout, data)
    stdout.write('done\n')

    return data


def lensout_read(file=None, run=None, split=None, silent=False):
    if file is None and run is None:
        raise ValueError("usage: lensout_read(file=, run=)")
    if file is None:
        return lensout_read_byrun(run)

    if not silent:
        stdout.write('Opening lensout file: %s\n' % file)
    fobj = open(file,'r')

    nlens = numpy.fromfile(fobj,dtype='i4',count=1)
    if not silent:
        stdout.write('  nlens: %s\n' % nlens[0])
    nbin  = numpy.fromfile(fobj,dtype='i4',count=1)
    if not silent:
        stdout.write('  nbin: %s\n' % nbin[0])

    dt = lensout_dtype(nbin[0])
    if not silent:
        stdout.write('Reading with dtype: %s\n' % str(dt))
    data = numpy.fromfile(fobj,dtype=dt,count=nlens[0])

    fobj.close()

    return data

def lensout_read_byrun(run):
    conf = json_read(run)
    nsplit = conf['nsplit']
    if nsplit == 0:
        file = sample_file('lensout', run)
        return lensout_read(file)

    # count number in each first
    ntot = 0
    for i in xrange(nsplit):
        file = sample_file('lensout', run, split=i)
        fobj = open(file,'r')
        nlens = numpy.fromfile(fobj,dtype='i4',count=1)
        nbin  = numpy.fromfile(fobj,dtype='i4',count=1)
        fobj.close()
        ntot += nlens[0]

    stdout.write("Reading %s from %s files\n" % (ntot,nsplit))
    stdout.write("First file is: %s\n" % sample_file('lensout',run,split=0))
    dt = lensout_dtype(nbin[0])
    data = numpy.zeros(ntot, dtype=dt)
    beg=0
    for i in xrange(nsplit):
        file = sample_file('lensout', run, split=i)
        t = lensout_read(file,silent=True)

        end = beg+t.size

        data[beg:end] = t

        beg = end

    return data


def lensout_dtype(nbin):

    nbin = int(nbin)
    dt=[('zindex','i4'),
        ('weight','f8'),
        ('npair','i8',nbin),
        ('rsum','f8',nbin),
        ('wsum','f8',nbin),
        ('dsum','f8',nbin),
        ('osum','f8',nbin)]
    return numpy.dtype(dt)

#
# read/write objshear lens input catalogs
#


def lcat_write(data, file=None, sample=None, split=None):
    if file is None and sample is None:
        raise ValueError("usage: lcat_write(data, file=, sample=)")
    if file is None:
        file = sample_file('lcat',sample, split=split)

    stdout.write("Writing %d to lens cat: '%s'\n" % (data.size, file))

    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)

    fobj = open(file,'w')

    narr =  numpy.array([data.size],dtype='i4')
    narr.tofile(fobj)

    data.tofile(fobj)

    fobj.close()

def lcat_read(file=None, sample=None, split=None):
    """
    import lensing
    d = lcat_read(file='somefile')
    d = lcat_read(sample='03')
    d = lcat_read(sample='03', split=27)
    """

    if file is None and sample is None:
        raise ValueError("usage: lcat_write(data, file=, sample=)")

    if file is None:
        file = sample_file('lcat',sample, split=split)

    stdout.write('Reading lens cat: %s\n' % file)
    fobj = open(file,'r')

    narr = numpy.fromfile(fobj, dtype='i4', count=1)

    stdout.write('Reading %d lenses\n' % narr[0])
    dt = lcat_dtype()

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    return data

def lcat_dtype():
    dt=[('ra','f8'),
        ('dec','f8'),
        ('z','f4'),
        ('dc','f4'),
        ('zindex','i4'),
        ('padding','i4')]
    return dt


#
# source catalogs
#


def scat_write(data, file=None, sample=None):
    if file is None and sample is None:
        raise ValueError("usage: scat_write(data, file=, sample=)")
    if file is None:
        file = sample_file('scat',sample, split=split)


    stdout.write("Writing binary source file:'%s'\n" % file)

    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)

    fobj = open(file,'w')

    narr =  numpy.array([data.size],dtype='i4')
    narr.tofile(fobj)

    data.tofile(fobj)

    fobj.close()

def scat_read(file=None, sample=None, interp_scinv=False):
    if file is None and sample is None:
        raise ValueError("usage: scat_write(data, file=, sample=)")
    if file is None:
        file = sample_file('scat',sample)

    stdout.write('Reading source_ztrue to: %s\n' % file)
    fobj = open(file,'r')

    narr = numpy.fromfile(fobj, dtype='i4', count=1)

    stdout.write('Reading %d sources\n' % narr[0])
    dt = scat_dtype(interp_scinv=interp_scinv)

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    return data

def scat_dtype(interp_scinv=False, nz=20):
    if interp_scinv:
        dt=[('ra','f8'),
            ('dec','f8'),
            ('g1','f4'),
            ('g2','f4'),
            ('err','f4'),
            ('hpixid','i4'),
            ('mean_scinv','f4',nz)]

    else:
        dt=[('ra','f8'),
            ('dec','f8'),
            ('g1','f4'),
            ('g2','f4'),
            ('err','f4'),
            ('hpixid','i4'),
            ('z','f4'),
            ('dc','f4')]
    return dt




#
# hand written json run config files
#

def json_read(run):
    """
    You can also do:
        rc = lensing.config.RunConfig(run)
    and rc is inherited from a dict
    """
    return eu.io.read(json_file(run))

def json_file(run):
    dir = json_dir()
    fname = 'lconfig-%s.json' % run
    fname = path_join(dir,fname)
    return fname

def json_dir():
    if 'ESPY_DIR' not in os.environ:
        raise ValueError("ESPY_DIR is not set")
    dir = os.environ['ESPY_DIR']
    dir = path_join(dir,'lensing','config')
    return dir



