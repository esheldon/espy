from __future__ import print_function
import os
from sys import stdout
import lensing
import numpy

import esutil as eu
from esutil.ostools import path_join, expand_path
from esutil.numpy_util import where1

finfo={}
finfo['lcat']     = {'subdir':'lcat',    'front':'lcat',    'ext':'.bin'}
finfo['scat']     = {'subdir':'scat',    'front':'scat',    'ext':'.bin'}
finfo['lensout']  = {'subdir':'lensout', 'front':'lensout', 'ext':'.rec'}
finfo['lensred']  = {'subdir':'lensout', 'front':'lensred', 'ext':'.rec'}
finfo['lensred-collate'] = {'subdir':'lensout', 'front':'lensred-collate', 'ext':'.rec'}
finfo['config']   = {'subdir':'proc',    'front':'run',     'ext':'.config'}
finfo['condor']   = {'subdir':'proc',    'front':'run',     'ext':'.condor'}
finfo['script']   = {'subdir':'proc',    'front':'run',     'ext':'.sh'}

finfo['pbslens']  = {'subdir':'pbslens','ext':'.pbs'}


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
    dsub = finfo[type]['subdir']
    d = path_join(d,dsub,sample)
    return d

def sample_file(type, sample, split=None, extra=None):
    d = sample_dir(type,sample)
    front = finfo[type]['front']
    fname=path_join(d, '%s-%s' % (front,sample))
    if split is not None:
        fname += '-%03d' % split
    if extra is not None:
        fname += '-%s' % extra
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

def reduced_collated_lensout_read(run):
    f = sample_file('lensred-collate',run)
    return eu.io.read(f)

def reduced_lensout_read(run):
    file = sample_file('lensred', run)
    return eu.io.read(file)
    
def lensout_read(file=None, run=None, split=None, silent=False, old=False):
    if old:
        return lensout_read_old(file=file,run=run,split=split,silent=silent)

    if file is None and run is None:
        raise ValueError("usage: lensout_read(file=, run=, split=, silent=False)")
    if file is None:
        if split is not None:
            file=sample_file('lensout', run, split=split)
            return lensout_read(file=file)

        return lensout_read_byrun(run)

    if not silent:
        stdout.write('Reading lensout file: %s\n' % file)
    file=expand_path(file)

    return eu.io.read(file)

def lensout_read_old(file=None, run=None, split=None, silent=False, old=False):
    '''
    Note old means something different here
    '''
    if file is None and run is None:
        raise ValueError("usage: lensout_read(file=, run=)")
    if file is None:
        if split is not None:
            file=sample_file('lensout', run, split=split)
            return lensout_read(file=file)

        return lensout_read_byrun(run)

    if not silent:
        stdout.write('Opening lensout file: %s\n' % file)
    file=expand_path(file)
    fobj = open(file,'r')

    if old:
        idt='i4'
    else:
        idt='i8'
    nlens = numpy.fromfile(fobj,dtype=idt,count=1)
    if not silent:
        stdout.write('  nlens: %s\n' % nlens[0])
    nbin  = numpy.fromfile(fobj,dtype=idt,count=1)
    if not silent:
        stdout.write('  nbin: %s\n' % nbin[0])

    dt = lensout_dtype(nbin[0], old=old)
    if not silent:
        stdout.write('Reading with dtype: %s\n' % str(dt))
    data = numpy.fromfile(fobj,dtype=dt,count=nlens[0])

    fobj.close()

    return data

def lensout_read_byrun(run, old=False):
    conf = cascade_config(run)
    nsplit = conf['src_config']['nsplit']

    stdout.write("Combining %s splits from run %s\n" % (nsplit,run))

    for i in xrange(nsplit):
        tdata = lensout_read(run=run, split=i, old=old)
        if i == 0:
            data = tdata
        else:
            # note zindex match occurs in add_lensums
            lensing.outputs.add_lensums(data, tdata)

    return data


# the old split of the lenses: we now split the sources
def lensout_read_byrun_lensplit(run):
    runconf = read_config('run',run)
    conf = read_config('lcat', runconf['lens_sample'])
    nsplit = conf['nsplit']
    if nsplit == 0:
        file = sample_file('lensout', run)
        return lensout_read(file)

    # count number in each first
    ntot = 0
    for i in xrange(nsplit):
        file = sample_file('lensout', run, split=i)
        fobj = open(file,'r')
        nlens = numpy.fromfile(fobj,dtype='i8',count=1)
        nbin  = numpy.fromfile(fobj,dtype='i8',count=1)
        fobj.close()
        print("    %s:  nlens: %i   nbin: %i" % (file,nlens,nbin) )
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


def lensout_dtype(nbin, old=False):

    if old:
        zidt='i4'
    else:
        zidt='i8'
    nbin = int(nbin)
    dt=[('zindex',zidt),
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


def lcat_write(data, file=None, sample=None):
    if file is None and sample is None:
        raise ValueError("usage: lcat_write(data, file=, sample=)")
    if file is None:
        file = sample_file('lcat',sample)

    stdout.write("Writing %d to lens cat: '%s'\n" % (data.size, file))

    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)

    fobj = open(file,'w')

    narr =  numpy.array([data.size],dtype='i8')
    narr.tofile(fobj)

    data.tofile(fobj)

    fobj.close()

def lcat_read(sample=None, file=None):
    """
    import lensing
    d = lcat_read(file='somefile')
    d = lcat_read(sample='03')
    """

    if file is None and sample is None:
        raise ValueError("usage: lcat_read(data, file=, sample=)")

    if file is None:
        file = sample_file('lcat', sample)

    stdout.write('Reading lens cat: %s\n' % file)
    fobj = open(file,'r')

    narr = numpy.fromfile(fobj, dtype='i8', count=1)

    stdout.write('Reading %d lenses\n' % narr[0])
    dt = lcat_dtype()

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    return data

def lcat_dtype():
    dt=[('zindex','i8'),
        ('ra','f8'),
        ('dec','f8'),
        ('z','f8')]
    return dt


#
# source catalogs
#


def scat_write(sample, data,split=None):
    from . import sigmacrit
    file = sample_file('scat',sample, split=split)
    conf = read_config('scat', sample)
    style=conf['sigmacrit_style']
    if style not in [1,2]:
        raise ValueError("sigmacrit_style should be in [1,2]")


    if style == 2:
        if 'zlmin' not in conf or 'zlmax' not in conf or 'dzl' not in conf:
            raise ValueError("You must have zlmin,zlmax,dzl in config")

        zlvals=sigmacrit.make_zlvals(conf['dzl'], conf['zlmin'], conf['zlmax'])
        nzl = zlvals.size

        data_nzl = data['scinv'].shape[1]
        if nzl != data_nzl:
            raise ValueError("Calculated nzl of %d but data has nzl of %d" % (nzl,data_nzl))

    print("Writing binary source file:",file)
    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)
    fobj = open(file,'w')


    print("Writing style:",style,"to file")
    stylewrite = numpy.array([style],dtype='i8')
    stylewrite.tofile(fobj)

    if style == 2:
        nzlwrite = numpy.array([nzl],dtype='i8')
        zlwrite = numpy.array(zlvals,dtype='f8',copy=False)

        print("Writing nzl:",nzl,"to file")
        nzlwrite.tofile(fobj)
        print("Writing zl to file")
        zlwrite.tofile(fobj)

    narr =  numpy.array([data.size],dtype='i8')
    print("Writing narr:",narr[0],"to file")
    narr.tofile(fobj)

    data.tofile(fobj)

    fobj.close()

def scat_read(sample=None, file=None, split=None):

    if file is None and sample is None:
        raise ValueError("usage: scat_write(data, file=, sample=)")
    
    if file is None:
        file = sample_file('scat',sample, split=split)


    stdout.write('Reading sources: %s\n' % file)
    fobj = open(file,'r')

    style = numpy.fromfile(fobj, dtype='i8', count=1)[0]
    print("Found sigmacrit_style",style)

    if style == 2:
        nzl = numpy.fromfile(fobj, dtype='i8', count=1)[0]
        print("found nzl =",nzl)
        print("reading zl values")
        zl = numpy.fromfile(fobj, dtype='f8', count=nzl)
    elif style == 1:
        nzl=None
    else:
        raise ValueError("sigmacrit_style should be in [1,2]")
    dt = scat_dtype(style, nzl=nzl)

    narr = numpy.fromfile(fobj, dtype='i8', count=1)

    stdout.write('Reading %d sources\n' % narr[0])

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    if style == 2:
        return data, zl
    else:
        return data

def scat_dtype(sigmacrit_style, nzl=None):
    dt=[('ra','f8'),
        ('dec','f8'),
        ('g1','f8'),
        ('g2','f8'),
        ('err','f8')]
        #('hpixid','i8')]

    if sigmacrit_style == 1:
        #dt += [('z','f8'), ('dc','f8')]
        dt += [('z','f8')]
    elif sigmacrit_style == 2:
        if nzl == None:
            raise ValueError('you must send nzl for sigmacrit_style of 2')
        dt += [('scinv','f8',nzl)]
    else:
        raise ValueError("sigmacrit_style should be in [1,2]")

    return dt




#
# hand written json run config files
#

def cascade_config(run):
    conf=read_config('run',run)

    ls = conf['lens_sample']
    conf['lens_config'] = read_config('lcat',ls)
    ss = conf['src_sample']
    conf['src_config'] = read_config('scat',ss)
    cs=conf['cosmo_sample']
    conf['cosmo_config'] = read_config('cosmo',cs)

    return conf



def read_config(type,id):
    return eu.io.read(config_file(type,id))

def config_file(type, id):
    dir = json_dir()
    fname = '%s-%s.json' % (type,id)
    fname = path_join(dir,fname)
    return fname

def json_dir():
    if 'ESPY_DIR' not in os.environ:
        raise ValueError("ESPY_DIR is not set")
    dir = os.environ['ESPY_DIR']
    dir = path_join(dir,'lensing','config')
    return dir



