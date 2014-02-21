"""

Work with lensing files.

Many of the files can be dealt with using the generic functions
    sample_dir
    sample_file
    sample_read
    sample_write

The exceptions are the 

    lcat/scat: 
        simple binary format easily read by objshear

    config:
        These are the yaml config files.  These reside under
        a different directory structure.  
        
        These *could* be put under the other system if we move from 'subdir' to
        'dir' and use environment variables, e.g.

            ${ESPY_DIR}/lensing/config
            ${LENSDIR}/lcat

        And run os.path.expandvars() on each path *before* expanding
        the local variables.

"""

from __future__ import print_function
import os
from sys import stdout, stderr
import lensing
import numpy

import esutil as eu
from esutil.ostools import path_join, expand_path
from esutil.numpy_util import where1
from . import sigmacrit

try:
    import recfile
except:
    pass

finfo={}
#finfo['lcat']     = {'subdir':'lcat/{sample}',  'name':'lcat-{sample}{extra}.{ext}', 'default_ext':'dat'}
finfo['scat']     = {'subdir':'scat/{sample}',  'name':'scat-{sample}.{ext}','default_ext':'dat'}

finfo['lcat-split']  = {'subdir':'lcat/{sample}',  
                        'name':'lcat-{sample}-{lens_split}.{ext}', 
                        'default_ext':'dat'}
finfo['scat-split']  = {'subdir':'scat/{sample}',
                        'name':'scat-{sample}-{src_split}.{ext}','default_ext':'dat'}
finfo['pofz-split']  = {'subdir':'scat/{sample}',
                        'name':'pofz-{sample}-{src_split}.{ext}','default_ext':'dat'}

# usually we only have the split versions of these
finfo['lensout']  = {'subdir':'lensout/{sample}',
                     'name':'lensout-{sample}.{ext}','default_ext':'dat'}
finfo['lensout-split']  = {'subdir':'lensout/{sample}',
                           'name':'lensout-{sample}-{lens_split}-{src_split}.{ext}',
                           'default_ext':'dat'}

# reduced across sources, but still split on lenses
finfo['src-reduced-split']  = {'subdir':'lensout/{sample}',
                               'name':'lensout-{sample}-src-reduce-{lens_split}.{ext}',
                               'default_ext':'dat'}


# this is the final reduced one over all lens and source splits
finfo['reduced']  = {'subdir':'lensout/{sample}',
                     'name':'reduced-{sample}.{ext}','default_ext':'dat'}

finfo['collated']  = {'subdir':'lensout/{sample}',
                      'name':'collated-{sample}.fits'}
finfo['collated-split']  = {'subdir':'lensout/{sample}',
                            'name':'collated-{sample}-{lens_split}.fits'}

# here {extra} is really only used for the matched random sums
finfo['binned']       = {'subdir':'lensout/{sample}/binned-{name}',
                         'name':'binned-{sample}-{name}{extra}.fits'}
                         #'name':'binned-{sample}-{name}{extra}.rec'}

finfo['weights'] = {'subdir':'lensout/{sample}/binned-{name}',
                    'name':'weights-{sample}-{name}{extra}.fits'}
                    #'name':'weights-{sample}-{name}{extra}.rec'}

# probably always want to send extra here
finfo['binned-plots']       = {'subdir':'lensout/{sample}/binned-{name}/plots',
                               'name':'binned-{sample}-{name}{extra}.{ext}'}

# currently these must also be binned...
finfo['corrected']  = {'subdir':'lensout/{sample}/binned-{name}',
                       'name':'corrected-{sample}-{name}.fits'}
                       #'name':'corrected-{sample}-{name}.rec'}
finfo['corrected-plots']       = {'subdir':'lensout/{sample}/binned-{name}/plots',
                                  'name':'corrected-{sample}-{name}{extra}.{ext}'}

finfo['corrected-ssh']  = {'subdir':'lensout/{sample}/binned-{name}',
                           'name':'corrected-ssh-{sample}-{name}.fits'}
finfo['corrected-ssh-plots']       = {'subdir':'lensout/{sample}/binned-{name}/plots',
                                      'name':'corrected-ssh-{sample}-{name}{extra}.{ext}'}




finfo['invert']       = {'subdir':'lensout/{sample}/binned-{name}',
                         'name':'invert-{sample}-{name}.fits'}
                         #'name':'invert-{sample}-{name}.rec'}

# note if extra is not None/'' in a call to sample_file, it gets a '-'
# prepended
finfo['fit']       = {'subdir':'lensout/{sample}/binned-{name}',
                      'name':'fit-{sample}-{name}{extra}.fits'}



# this config is objshear_config not the yaml files
finfo['config']   = {'subdir':'proc/{sample}', 'name':'run-{sample}.cfg'}

finfo['wq-genrand-split'] = {'subdir':'proc/{sample}/genrand',
                             'name':'run-{sample}-{lens_split}-genrand.yaml'}
finfo['wq-split']   = {'subdir':'proc/{sample}', 
                       'name':'run-{sample}-{lens_split}-{src_split}.yaml'}

# reduce *all* files, even if split on both lenses and sources
finfo['wq-reduce']   = {'subdir':'proc/{sample}', 'name':'run-{sample}-reduce-all.yaml'}

# reduce over lenses at fixed source split
finfo['wq-src-reduce-split']   = {'subdir':'proc/{sample}', 
                                  'name':'run-{sample}-src-reduce-{lens_split}.yaml'}
#finfo['wq-lens-reduce-split']   = {'subdir':'proc/{sample}', 
#                                   'name':'run-{sample}-lens-reduce-{src_split}.yaml'}

# concat the src-reduce outputs over the lens splits
finfo['wq-lens-concat']   = {'subdir':'proc/{sample}', 'name':'run-{sample}-lens-concat.yaml'}

finfo['wq-collate-split']   = {'subdir':'proc/{sample}', 
                               'name':'run-{sample}-collate-{lens_split}.yaml'}
#
# base directories
#

def lensdir():
    """
    This is the root dir under which all data are stored
    """
    if 'LENSDIR' not in os.environ:
        raise ValueError("LENSDIR is not set")
    return os.environ['LENSDIR']

def hdfs_dir():
    """
    my area in hdfs
    """
    return 'hdfs://astro0034.rcf.bnl.gov/user/esheldon/lensing'

def local_dir():
    """
    This is a temporary "local" area on each node
    """
    return '/data/objshear/lensing'



def sample_dir(**keys):
    """

    Generic routine to get the directory for a sample of a given type, e.g.
    for lcat,scat,proc files etc

    See finfo for a list of types

    """

    type=keys.get('type',None)
    sample=keys.get('sample',None)
    name=keys.get('name',None)
    fs=keys.get('fs','nfs')

    if type not in finfo:
        if type+'-split' in finfo:
            type=type+'-split'
        else:
            raise ValueError("Unknown file type: '%s'" % type)
    if fs == 'nfs':
        d = lensdir()
    elif fs == 'hdfs':
        d = hdfs_dir()
    elif fs == 'local':
        d = local_dir()
    else:
        raise ValueError("file system not recognized: %s" % fs)

    dsub = finfo[type]['subdir'].format(sample=sample, name=name)
    d = path_join(d,dsub)
    return d


def sample_file(**keys):
    """
    Generic routine to get the file for a sample of a given type, e.g.  for
    lcat,scat files etc

    All parameters are named keywords.

    See finfo for a list of types

    parameters
    ----------
    type:
        type name, e.g. lcat,scat
    sample:
        sample name
    lens_split:
        lens split number
    src_split:
        src split number
    name:
        name for binned
    extra:
        extra name
    ext:
        extension
    """

    type=keys.get('type',None)
    sample=keys.get('sample',None)
    if type is None or sample is None:
        raise ValueError("send type= and sample=")

    lens_split=keys.get('lens_split',None)
    src_split=keys.get('src_split',None)

    name=keys.get('name',None)
    extra=keys.get('extra',None)
    ext=keys.get('ext',None)

    d = sample_dir(**keys)

    if src_split is not None or lens_split is not None:
        if 'split' not in type:
            type = type+'-split'

    if src_split is not None:
        src_split = '%03d' % src_split

    if lens_split is not None:
        lens_split = '%03d' % lens_split

    if extra is None:
        extra=''
    if extra != '':
        extra = '-'+extra

    if ext is None:
        ext=finfo[type].get('default_ext',None)

    f = finfo[type]['name'].format(sample=sample, 
                                   lens_split=lens_split, 
                                   src_split=src_split, 
                                   extra=extra, 
                                   ext=ext, 
                                   name=name)
    f = os.path.join(d,f)
    return f


def sample_read(**keys):
    """

    Generic reader.

    This works for most things, but won't work for reading multiple splits and
    combining which is done by /bin/reduce-lensout.py and doesn't work for the
    config, lcat and scat special file types.

    parameters
    -----------
    type:
        file type, e.g. lcat
    sample:
        sample name
    lens_split:
        lens split number
    src_split:
        src split number
    name:
        extra name
    extra:
        extra string
    fs:
        file system
    """
    type=keys['sample']
    if type == 'lcat':
        return lcat_read(**keys)
    elif type == 'scat':
        return scat_read_ascii(**keys)

    f=sample_file(**keys)
    return eu.io.read(f, verbose=True)

def sample_write(**keys):
    """

    Generic writer, just gets filename and runs eu.io.write

    clobber=True is default for this function.

    parameters
    ----------
    data:
        Keyword holding the data to write
    keys:
        other keywords for the sample_file(**keys) function.
    """

    data=keys.get('data',None)

    f=sample_file(**keys)
    if 'hdfs://' not in f:
        make_dir_from_path(f)

    print('writing file:',f)
    return eu.io.write(f, data, clobber=True)


#
# these are the basic catalogs we use as inputs.  they get converted to 'lcat'
# or 'scat' files for input to objshear, and collated with the reduced
# catalogs.
#
# Note these call into lensing.lcat and lensing.scat where specialization
# occurs
#

def catalog_dir():
    """
    The root of the catalog directory. ${LENSDIR}/catalogs
    """
    catdir = path_join(lensdir(), 'catalogs')
    return catdir


#
# read/write objshear lens and source input catalogs.  These are special
# formats understood by the objshear code.
#


def lcat_file(**keys):
    """
    import lensing
    lcat_file(sample='rm03')
    """
    keys['ext'] = 'dat'
    keys['type'] = 'lcat'
    return sample_file(**keys)

def lcat_write(**keys):

    sample=keys.get('sample',None)
    data=keys.get('data',None)

    if sample is None or data is None:
        raise ValueError("usage: lcat_write(sample=, data=)")

    keys['ext'] = 'dat'
    keys['type'] = 'lcat'
    fs=keys.get('fs','nfs')

    fname = lcat_file(**keys)

    stdout.write("Writing %d to %s: '%s'\n" % (data.size, 'lcat', fname))

    if fs=='nfs':
        if os.path.exists(fname):
            os.remove(fname)
        else:
            d=os.path.dirname(fname)
            if not os.path.exists(d):
                os.makedirs(d)
        with recfile.Open(fname,'w',delim=' ') as robj:
            robj.fobj.write('%d\n' % data.size)
            robj.write(data)
    else:
        with eu.hdfs.HDFSFile(fname,verbose=True) as fobj:
            with recfile.Open(fobj.localfile,'w',delim=' ') as rec:
                rec.fobj.write('%d\n' % data.size)
                rec.write(data)

            fobj.put(clobber=True)


def lcat_read(**keys):
    """
    import lensing
    d = lcat_read(sample='rm03')

    """
    file = lcat_file(**keys)

    stdout.write('Reading lens cat: %s\n' % file)
    with eu.hdfs.HDFSFile(file,verbose=True) as hf:
        hf.stage()

        with open(hf.localfile) as fobj:
            import recfile
            nlens=int(fobj.readline())
            print('Reading',nlens,'lenses',file=stderr)
            dt = lcat_dtype()
            with recfile.Open(fobj, mode='r', dtype=dt, delim=' ',nrows=nlens) as robj:
                data = robj[:]

    return data

def collated_file(**keys):
    keys['type'] = 'collated'
    return sample_file(**keys)

def collated_read(**keys):
    fname=collated_file(**keys)
    verbose=keys.get('verbose',True)
    print('Reading collated:',fname,file=stderr)
    return eu.io.read(fname,verbose=verbose)

def lensout_file(**keys):
    """
    import lensing
    lensout_file(sample='rm03')
    """
    keys['ext'] = 'dat'
    keys['type'] = 'lensout'
    return sample_file(**keys)

def lensout_read(**keys):

    sample=keys.get('sample',None)
    if sample is None:
        raise ValueError("send sample=")
    conf = cascade_config(sample)
    nbin = conf['lens_config']['nbin']

    dtype=lensout_dtype(nbin)

    fname=lensout_file(**keys)
    print('Reading lensout:',fname,file=stderr)
    return eu.io.read(fname, dtype=dtype, delim=' ', type='rec')

def reduced_file(**keys):
    keys['ext'] = 'dat'
    keys['type'] = 'reduced'
    return sample_file(**keys)

def reduced_read(**keys):
    sample=keys.get('sample',None)
    if sample is None:
        raise ValueError("send sample=")

    conf = cascade_config(sample)
    nbin = conf['lens_config']['nbin']

    dtype=lensout_dtype(nbin)

    fname = reduced_file(**keys)
    print('Reading reduced:',fname,file=stderr)
    return eu.io.read(fname, dtype=dtype, delim=' ', type='rec')

def lensout_dtype(nbin):
    dtype=[('index','i8'),
           ('zindex','i8'),
           ('weight','f8'),
           ('totpairs','i8'),
           ('sshsum','f8'),
           ('npair','i8',nbin),
           ('rsum','f8',nbin),
           ('wsum','f8',nbin),
           ('dsum','f8',nbin),
           ('osum','f8',nbin)]
    return numpy.dtype(dtype)

def lensout_im3shape_dtype(nbin):
    dtype=[('index','i8'),
           ('zindex','i8'),
           ('weight','f8'),
           ('totpairs','i8'),
           ('npair','i8',nbin),
           ('rsum','f8',nbin),
           ('wsum','f8',nbin),
           ('dsum','f8',nbin),
           ('osum','f8',nbin)]
    return numpy.dtype(dtype)


def lcat_read_old(sample=None, extra=None, file=None, old=False):
    """
    import lensing
    d = lcat_read(file='somefile')
    d = lcat_read(sample='03')

    """

    if file is None and sample is None:
        raise ValueError("usage: lcat_read(data, file=, sample=)")

    if file is None:
        file = sample_file(type='lcat', sample=sample, extra=extra)

    stdout.write('Reading lens cat: %s\n' % file)
    fobj = open(file,'r')

    narr = numpy.fromfile(fobj, dtype='i8', count=1)

    stdout.write('Reading %d lenses\n' % narr[0])
    dt = lcat_dtype(old=old)

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    return data

def lcat_dtype(**keys):
    dt=[('zindex','i4'),
        ('ra','f8'),
        ('dec','f8'),
        ('z','f8'),
        ('maskflags','i8')]

    return dt


def scat_file(**keys):
    """
    scat_file(sample=)
    """
    keys['type'] = 'scat'
    return sample_file(**keys)

def scat_write_ascii(**keys):
    from sys import stdout

    sample=keys.get('sample',None)
    data=keys.get('data',None)
    fs=keys.get('fs','nfs')

    if sample is None or data is None:
        raise ValueError("usage: scat_write_ascii(sample=, data= [, src_split=]")

    fname=scat_file(**keys)

    conf = read_config('scat', sample)

    print("Writing ascii source file:",fname)

    if fs=='nfs':
        if os.path.exists(fname):
            os.remove(fname)
        else:
            d=os.path.dirname(fname)
            if not os.path.exists(d):
                os.makedirs(d)
        with recfile.Open(fname,'w',delim=' ') as robj:
            robj.write(data)
    else:
        if eu.hdfs.exists(file):
                eu.hdfs.rm(file)
        with eu.hdfs.HDFSFile(file) as hdfs_file:
            with recfile.Open(hdfs_file.localfile,'w',delim=' ') as robj:
                robj.write(data)
            hdfs_file.put()

def scat_read_ascii(**keys):

    sample=keys.get('sample',None)

    if sample is None:
        raise ValueError("usage: data=scat_read_ascii(sample=, [, src_split=]")

    conf = read_config('scat', sample)
    style=conf['sigmacrit_style']
    if style not in [1,2]:
        raise ValueError("sigmacrit_style should be in [1,2]")

    if style == 2:
        sconf = read_config('scinv', conf['scinv_sample'])
        zlvals=sigmacrit.make_zlvals(sconf['dzl'],sconf['zlmin'],sconf['zlmax'])
        nzl = zlvals.size
    else:
        nzl=None
    if 'gmix' in sample:
        dt = scat_gmix_dtype(style, nzl=nzl)
    elif 'im3' in sample:
        dt = scat_im3shape_dtype(style, nzl=nzl)
    else:
        dt = scat_dtype(style, nzl=nzl)

    file=scat_file(**keys)
    print("reading scat file:",file,file=stderr)

    data = eu.io.read(file, delim=' ', dtype=dt, type='rec')

    #with eu.hdfs.HDFSFile(file) as hdfs_file:
    #    hdfs_file.stage()
    #    with recfile.Open(hdfs_file.localfile,'r',delim=' ',dtype=dt) as robj:
    #        data = robj[:]

    return data

'''
def scat_read_binary(sample=None, file=None, lens_split=None):
    """
    Not used
    """

    if file is None and sample is None:
        raise ValueError("usage: scat_write(data, file=, sample=)")
    
    if file is None:
        file = sample_file(type='scat',sample=sample, lens_split=lens_split)


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

    return data
'''


def scat_dtype(sigmacrit_style, nzl=None):
    dt=[('ra','f8'),
        ('dec','f8'),
        ('g1','f8'),
        ('g2','f8'),
        ('err','f8'),
        ('mag','f8'),
        ('R','f8')]
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


def scat_gmix_dtype(sigmacrit_style, nzl=None):
    dt=[('ra','f8'),
        ('dec','f8'),
        ('g1','f8'),
        ('g2','f8'),
        ('gcov11','f8'),
        ('gcov12','f8'),
        ('gcov22','f8'),
        ('gsens1','f8'),
        ('gsens2','f8')]

    if sigmacrit_style == 1:
        dt += [('z','f8')]
    elif sigmacrit_style == 2:
        if nzl == None:
            raise ValueError('you must send nzl for sigmacrit_style of 2')
        dt += [('scinv','f8',nzl)]
    else:
        raise ValueError("sigmacrit_style should be in [1,2]")

    return dt


def scat_im3shape_dtype(sigmacrit_style, nzl=None):
    dt=[('ra','f8'),
        ('dec','f8'),
        ('g1','f8'),
        ('g2','f8'),
        ('weight','f8')]

    if sigmacrit_style == 1:
        dt += [('z','f8')]
    elif sigmacrit_style == 2:
        if nzl == None:
            raise ValueError('you must send nzl for sigmacrit_style of 2')
        dt += [('scinv','f8',nzl)]
    else:
        raise ValueError("sigmacrit_style should be in [1,2]")

    return dt




#
# hand written yaml config files
#

def cascade_config(run):
    conf=read_config('run',run)

    ls = conf['lens_sample']
    conf['lens_config'] = read_config('lcat',ls)

    ss = conf['src_sample']
    conf['src_config'] = read_config('scat',ss)

    if 'scinv_sample' in conf['src_config']: 
        conf['scinv_config'] = read_config('scinv',
                                           conf['src_config']['scinv_sample'])

        scosmo = conf['scinv_config']['cosmo_sample']
    else:
        scosmo=None

    lcosmo = conf['lens_config']['cosmo_sample']
    if scosmo is not None and scosmo != lcosmo:
        mess="""
        mismatch between source scinv cosmo and lens cosmo
            scinv cosmo: %s
            lens cosmo:  %s
        """ % (scosmo, lcosmo)
        raise ValueError(mess)

    conf['cosmo_config'] = read_config('cosmo', lcosmo)
    return conf



def read_config(type,id):
    f = config_file(type, id, ext='yaml')
    if not os.path.exists(f):
        f = config_file(type, id, ext='json')
        if not os.path.exists(f):
            raise ValueError("No config file found for type %s id %s\n" % (type,id))

    return eu.io.read(f)

def config_file(type, id, ext='yaml'):
    """
    Want to move over to yaml files
    """
    dir = config_dir()
    fname = '%s-%s.%s' % (type,id,ext)
    fname = path_join(dir,fname)
    return fname

def config_dir():
    if 'ESPY_DIR' not in os.environ:
        raise ValueError("ESPY_DIR is not set")
    dir = os.environ['ESPY_DIR']
    dir = path_join(dir,'lensing','config')
    return dir





#
# old read objshear outputs
#


def lensout_read_old(file=None, run=None, lens_split=None, silent=False, old=False):
    '''
    Note old means something different here
    '''
    if file is None and run is None:
        raise ValueError("usage: lensout_read(file=, run=)")
    if file is None:
        if lens_split is not None:
            file=sample_file(type='lensout', sample=run, lens_split=lens_split)
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

    dt = lensout_dtype_old(nbin[0], old=old)
    if not silent:
        stdout.write('Reading with dtype: %s\n' % str(dt))
    data = numpy.fromfile(fobj,dtype=dt,count=nlens[0])

    fobj.close()

    return data


def lensout_dtype_old(nbin, old=False):

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

def make_dir_from_path(fname):
    d = os.path.dirname(fname)
    if not os.path.exists(d):
        print("Making dir:",d)
        os.makedirs(d)


