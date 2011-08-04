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
from sys import stdout
import lensing
import numpy

import esutil as eu
from esutil.ostools import path_join, expand_path
from esutil.numpy_util import where1


finfo={}
finfo['lcat']     = {'subdir':'lcat/{sample}',  'name':'lcat-{sample}.bin'}
finfo['scat']     = {'subdir':'scat/{sample}',  'name':'scat-{sample}.bin'}

finfo['scat-split']  = {'subdir':'scat/{sample}',
                            'name':'scat-{sample}-{split}.bin'}

finfo['lensout-split']  = {'subdir':'lensout/{sample}',
                               'name':'lensout-{sample}-{split}.rec'}

finfo['reduced']  = {'subdir':'lensout/{sample}',
                         'name':'reduced-{sample}.fits'}
finfo['collated']  = {'subdir':'lensout/{sample}',
                                 'name':'collated-{sample}.fits'}

# here {extra} is really only used for the matched random sums
# these need to be .rec because we have matrix columns
finfo['binned']       = {'subdir':'lensout/{sample}/binned-{name}',
                              'name':'binned-{sample}-{name}{extra}.rec'}

finfo['weights'] = {'subdir':'lensout/{sample}/binned-{name}',
                    'name':'weights-{sample}-{name}{extra}.rec'}

# probably always want to send extra here
finfo['binned-plots']       = {'subdir':'lensout/{sample}/binned-{name}/plots',
                                'name':'binned-{sample}-{name}{extra}.{ext}'}

# currently these must also be binned...
finfo['corrected']  = {'subdir':'lensout/{sample}/binned-{name}',
                       'name':'corrected-{sample}-{name}.rec'}

finfo['invert']       = {'subdir':'lensout/{sample}/binned-{name}',
                         'name':'invert-{sample}-{name}.rec'}

# note if extra is not None/'' in a call to sample_file, it gets a '-'
# prepended
finfo['fit']       = {'subdir':'lensout/{sample}/binned-{name}',
                      'name':'fit-{sample}-{name}{extra}.fits'}



# this config is objshear_config not the yaml files
finfo['config-split']   = {'subdir':'proc/{sample}', 'name':'run-{sample}-{split}.config'}

finfo['script-split']   = {'subdir':'proc/{sample}', 'name':'run-{sample}-{split}.sh'}
finfo['condor-split']   = {'subdir':'proc/{sample}', 'name':'run-{sample}-{split}.condor'}
finfo['condor']   = {'subdir':'proc/{sample}', 'name':'run-{sample}.condor'}


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
    return 'hdfs:///user/esheldon/lensing'

def local_dir():
    """
    This is a temporary "local" area on each node
    """
    return '/data/objshear/lensing'



def sample_dir(type, sample, name=None, fs='nfs'):
    """

    Generic routine to get the directory for a sample of a given type, e.g.
    for lcat,scat,proc files etc

    See finfo for a list of types

    """
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

def sample_file(type, sample, split=None, name=None, extra=None, ext=None, fs='nfs'):
    """

    Generic routine to get the file for a sample of a given type, e.g.  for
    lcat,scat files etc

    See finfo for a list of types

    """

    d = sample_dir(type, sample, name=name, fs=fs)
    if split is not None:
        split = '%03d' % split
        type = type+'-split'

    if extra is None:
        extra=''
    if extra != '':
        extra = '-'+extra

    f = finfo[type]['name'].format(sample=sample, split=split, extra=extra, ext=ext, name=name)
    f = os.path.join(d,f)
    return f


def sample_read(type, sample, split=None, name=None, extra=None, fs='nfs'):
    """

    Generic reader.

    This works for most things, but won't work for reading multiple splits and
    combining which is done by /bin/reduce-lensout.py and doesn't work for the
    config, lcat and scat special file types.

    """
    f=sample_file(type, sample, split=split, name=name, extra=extra, fs=fs)
    return eu.io.read(f, verbose=True)

def sample_write(data, type, sample, split=None, name=None, extra=None, fs='nfs', **keys):
    """

    Generic writer, just gets filename and runs eu.io.write

    """
    f=sample_file(type, sample, split=split, name=name, extra=extra, fs=fs)
    print("writing",type,"file:",f)
    #keys['verbose'] = True
    return eu.io.write(f, data, **keys)


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

def original_catalog_file(type, sample):
    """
    This is the original catalog from which we generate the objshear inputs.
    Note for dr8 we read from columns, so there is no "original file"
    """
    if type == 'lens':
        return lensing.lcat.original_file(sample)
    elif type == 'source':
        return lensing.scat.original_file(sample)

def read_original_catalog(type, sample):
    """
    This is the original catalog from which we generate the objshear inputs.
    Note for dr8 we read from columns, so there is no "original file"
    """
    if type == 'lens':
        return lensing.lcat.read_original(sample)
    elif type == 'source':
        return lensing.scat.read_original(sample)
 

#
# read/write objshear lens and source input catalogs.  These are special
# formats understood by the objshear code.
#


def lcat_write(sample, data):
    file = sample_file('lcat',sample)

    stdout.write("Writing %d to %s: '%s'\n" % (data.size, 'lcat', file))

    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)

    fobj = open(file,'w')

    narr =  numpy.array([data.size],dtype='i8')
    narr.tofile(fobj)

    data.tofile(fobj)

    fobj.close()

def lcat_read(sample=None, file=None, old=False):
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
    dt = lcat_dtype(old=old)

    data = numpy.fromfile(fobj, dtype=dt, count=narr[0])
    fobj.close()

    return data

def lcat_dtype(old=False):
    if not old:
        dt=[('zindex','i8'),
            ('ra','f8'),
            ('dec','f8'),
            ('z','f8'),
            ('maskflags','i8')]
    else:
        dt=[('zindex','i8'),
            ('ra','f8'),
            ('dec','f8'),
            ('z','f8')]

    return dt


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

def scat_write_ascii(sample, data,split=None):
    from esutil import recfile
    from sys import stdout
    from . import sigmacrit

    file = sample_file('scat',sample, split=split, ext='.dat')

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

    print("Writing ascii source file:",file)
    d = os.path.dirname(file)
    if not os.path.exists(d):
        stdout.write("Making dir: '%s'\n" % d)
        os.makedirs(d)

    with recfile.Open(file,'w',delim=' ') as robj:

        #robj = recfile.Open(file,'w',delim=' ')

        print("Writing header\n")
        #print("  Writing size of sources:",data.size)
        robj.fobj.write('%-15s = %d\n' % ("size",data.size))
        stdout.write('%-15s = %d\n' % ("size",data.size))

        #print("  Writing style:",style,"to file")
        robj.fobj.write('%-15s = %d\n' % ("sigmacrit_style",style))
        stdout.write('%-15s = %d\n' % ("sigmacrit_style",style))

        if style == 2:
            #print("  Writing nzl:",nzl,"to file")
            robj.fobj.write('%-15s = %d\n' % ("nzlens",nzl))
            stdout.write('%-15s = %d\n' % ("nzlens",nzl))

            #print("  Writing zl to file")
            #robj.fobj.write("%-15s = " % "zlens")
            stdout.write("%-15s = " % "zlens")
            zlwrite = numpy.array(zlvals,dtype='f8',copy=False)
            zlwrite.tofile(robj.fobj, sep=' ')
            zlwrite.tofile(stdout, sep=' ')
            robj.fobj.write('\n')
            stdout.write('\n')

        robj.fobj.write("END\n\n")
        stdout.write("END\n\n")

        print("Writing data")
        robj.write(data)


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
# hand written yaml config files
#

def cascade_config(run):
    conf=read_config('run',run)

    ls = conf['lens_sample']
    conf['lens_config'] = read_config('lcat',ls)
    ss = conf['src_sample']
    conf['src_config'] = read_config('scat',ss)
    cs=conf['lens_config']['cosmo_sample']
    conf['cosmo_config'] = read_config('cosmo',cs)

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

