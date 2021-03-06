"""
Procedure to do weighting

The basic idea is the following:

    * define a photometric sample.

    * Define a matched training sample that references the photometric sample
    and some other parameters. Match the training sample to that sample.  This
    simplifies everything because we only have one set of photometric cuts
    except in special cases (need to make sure e.g. SDSS mag range is preserved
    in the match, which isn't the case due to some deblending issues).   

    * Generate inputs for weighting for these samples.

Details:

First you need a sweep reduce columns database.  e.g. from run
prim04.

Create a photometric sample, both the one big file for doing the overall
weights and the smaller ones for the individual p(z) calculations, which can be
done in parallel.

Create the config file in 
    ${ESPY_DIR}/zphot/config/zinput-{photo_sample}.yaml
e.g.

photo_sample:   "06"
source:         dr8_final
primary:        true
cmodel_rmin:    15.0
cmodel_rmax:    21.8
modelmag_min:   15.0
modelmag_max:   29.0
weighting:      true
masktype:       [good,star]
filetype:       rec

Then run these scripts

    bin/create-input.py photo_sample
    bin/create-input.py -n 50 photo_sample

Then match the training samples.  You'll need to create a config file under 

    ${ESPY_DIR}/zphot/config/train-{train_sample}.yaml

e.g.

train_sample:   "13"
photo_sample:   "06"
sdss_rmax:      19.6
zmin:           -0.001
zmax:           2.0
comment:        photo sample 06, r < 21.8

NOTE: photo_sample is specified here, so train_sample is the only identifier
you'll need after this.

First do a match and then make the versions formatted for input to weighting:

    bin/match-training.py train_sample
    bin/make-weighting-train.py train_sample
    


Then run to get the overall weights.  Create the wq file.  

    /bin/create-wq.py weights wrun

For that you'll need to create a weights-{wrun}.yaml file:

wrun: "13"
train_sample: "13"
n_near1: 5
n_near2: 100


Then make some plots comparing the distributions of the variables:
    /bin/weighting-qaplots.py train_sample

To just do the z hist since varhist is slow:
    /bin/weighting-qaplots.py -z train_sample


To make p(z) for each object
    /bin/create-wq.py -n 50 pofz pzrun

Which will require a pofz run defined in pofz-{pzrun}.yaml in the
config directory, e.g.

pzrun:  15
wrun:   13
nz:     35
zmin:   0.0
zmax:   1.1


Then you can sum up the p(z) to get an overall histogram
    /bin/sum-pofz.py {pzrun}


Plot with the p(z) from summed individual overplotted
    /bin/weighting-qaplots.py -z --pzrun pzrun wrun

You can correct the individual p(z) using the ratio
of N(z) to summed p(z)

    /bin/correct-pofz.py pzrun

You can make a columns database after the corrected p(z) are created

    /bin/make-pofz-columns.py pzrun

Finally, after the columns db is made, you can make an official release split
up by run,camcol

    /bin/make-release.py pzrun
"""
from __future__ import print_function

import os
from sys import stdout,stderr
import glob
import pprint

import numpy
from numpy import where

import esutil as eu
from esutil.recfile import Recfile
from esutil.ostools import path_join, expand_path
from esutil.numpy_util import ahelp
from esutil.numpy_util import where1
from esutil.stat import histogram
from esutil.plotting import make_hist_curve

try:
    import es_sdsspy
    import sdsspy
except:
    pass

import biggles
from biggles import FramedPlot, FramedArray, PlotKey, Histogram, \
        SymmetricErrorBarsY, PlotLabel, Curve, Points
import converter

import copy

import pbs

import zphot

try:
    import columns
except:
    pass

_sample='dr8'
_basedir = '~esheldon/photoz/weighting/%s' % _sample
_basedir = os.path.expanduser(_basedir)


def pofz_columns_dir(pzrun):
    dir = pofz_dir(pzrun)
    dir += '.cols'
    #return path_join(dir,'pofz-%s.cols' % pzrun)
    return dir

def open_pofz_columns(pzrun):
    d=pofz_columns_dir(pzrun)
    print("Opening poz columns:",d)
    return columns.Columns(d)





def read_z(filename):
    stdout.write("Reading z file: '%s'\n" % filename)
    dt=[('zmin','f4'),('zmax','f4')]
    r = eu.recfile.Open(filename,'r',dtype=dt,delim=' ')
    data = r.read()
    r.close()
    return data

def photo_dtype(photo_sample):
    dtype = [('photoid','i8'),
             ('cmodelmag_dered_r','f4'),
             ('model_umg','f4'),
             ('model_gmr','f4'),
             ('model_rmi','f4'),
             ('model_imz','f4')]
    conf = zphot.read_config('zinput',photo_sample)
    if 'extra_columns' in conf:
        for cdict in conf['extra_columns']:
            dtype += [(cdict['name'], cdict['dtype'])]

    return dtype

def read_photo(photo_sample,rows=None):
    dt=photo_dtype(photo_sample)
    zs = zphot.select.ColumnSelector(photo_sample)
    filename = zs.filename()
    rec_filename=filename.replace('.dat','.rec')
    if not os.path.exists(rec_filename):
        print("Reading weighting photo file:",filename)
        r=eu.recfile.Open(filename,dtype=dt,delim=' ')
        data = r.read(rows=rows)
        r.close()

        print("caching binary version:",rec_filename)
        eu.io.write(rec_filename, data)
    else:
        print("Reading weighting photo rec file:",rec_filename)
        data = eu.io.read(rec_filename) 
    return data



def training_dtype(train_sample):
    out_dtype = [('z','f4'),
                 ('catid','i4'),
                 ('weight','f4'),
                 ('cmodelmag_dered_r','f4'),
                 ('model_umg','f4'),
                 ('model_gmr','f4'),
                 ('model_rmi','f4'),
                 ('model_imz','f4')]

    conf = zphot.read_config('train',train_sample)
    pconf = zphot.read_config('zinput',conf['photo_sample'])
    if 'extra_columns' in pconf:
        for cdict in pconf['extra_columns']:
            out_dtype += [(cdict['name'], cdict['dtype'])]

    return out_dtype
 
def read_training(train_sample, type):
    
    train = zphot.training.Training(train_sample)
    filename = train.fname_matched(type)
    print("Reading training file:", filename)
    dt = training_dtype(train_sample)
    return _read_training(filename, dt)
    '''
    r = eu.recfile.Open(filename, dtype=dt, delim=' ')
    data = r.read()
    r.close()
    return data
    '''

def _read_training(filename, dtype):
    r = eu.recfile.Open(filename, dtype=dtype, delim=' ')
    data = r.read()
    r.close()
    return data



def pofz_dir(pzrun):
    dir = weights_basedir()
    return path_join(dir, 'pofz-%s' % pzrun)

def pofz_file(pzrun, chunk=None, with_rmag=False):
    dir = pofz_dir(pzrun)
    f = 'pofz-%s' % pzrun
    if chunk is not None:
        if isinstance(chunk,basestring):
            if chunk == '*':
                f = f+'-chunk*'
            elif chunk == 'rand':
                f = f + '-'+chunk
            else:
                raise ValueError("Expected string chunk of '*' or 'rand', "
                                 "got '%s'" % chunk)
        else:
            f = f + '-chunk%03i' % chunk
    
    if with_rmag:
        f += '-withrmag'
    f = f+'.dat'
    return path_join(dir, f)

def pofz_release_dir(pzrun):
    dir = weights_basedir()
    return path_join(dir, 'pofz-%s-release' % pzrun)

def zbins_release_file(pzrun, ext='fits'):
    dir = pofz_release_dir(pzrun)
    f = 'zbins-%s.%s' % (pzrun,ext)
    return path_join(dir, f)

def pofz_release_file(pzrun, run, ext='fits'):
    dir = pofz_release_dir(pzrun)
    f = 'pofz-%s-%06d.%s' % (pzrun,run,ext)
    return path_join(dir, f)


def corrected_pofz_file(pzrun, chunk):
    """
    Files corrected so mean is same as overall
    N(z)
    """
    f=pofz_file(pzrun, chunk)
    f = f.replace('.dat','-corr.dat')
    return f

def pofz_correction_file(pzrun):
    """
    This holds the correction factor for p(z)s.
    """
    dir = pofz_dir(pzrun)
    f = 'pofz-correction-%s.rec' % pzrun
    return path_join(dir, f)

def pofz_hist_file(pzrun):
    """
    Summed zhist
    """
    dir = pofz_dir(pzrun)
    f = 'zhist-%s.rec' % pzrun
    return path_join(dir, f)


def z_file(pzrun, chunk=None):
    dir = pofz_dir(pzrun)
    f = 'z-%s' % pzrun
    if chunk is not None:
        f = f + '-chunk%03i' % chunk
    f = f+'.dat'
    return path_join(dir, f)


def read_weights(wrun, iteration, nozero=False):
    """
    Read weights output of calcweights code
    """
    f = weights_file(wrun, iteration, nozero=nozero)

    conf = zphot.read_config('weights', wrun)
    dt = training_dtype(conf['train_sample'])

    print("reading weights output:",f)
    r = eu.recfile.Open(f, dtype=dt, delim=' ')
    data = r.read()
    r.close()
    return data


def weights_basedir():
    pzdir = zphot.photoz_dir()
    dir = path_join(pzdir, 'weighting', _sample)
    return dir

def weights_dir(wrun):
    """
    loation of the output file
    """
    dir = weights_basedir()
    return path_join(dir, 'weights-%s' % wrun)

def weights_file(wrun, iteration, nozero=False):
    """
    Output of the calcweights code
    """
    conf = zphot.read_config('weights',wrun)
    dir = weights_dir(wrun)
    
    f = 'weights'
    if iteration == 1:
        if nozero:
            f += '-nozero'
        f += '-%s-%s.dat' % (wrun,conf['n_near1'])
    elif iteration == 2:
        f += '-%s-%s.dat' % (wrun,conf['n_near2'])
    else:
        raise ValueError("only support iteration 1 or 2")

    return path_join(dir, f)

def read_num(wrun, iteration, rows=None):
    """
    Read Output of the calcweights code.  You can
    call _read_num to read a known file name
    """
    f = num_file(wrun, iteration)
    return _read_num(f, rows=rows)

def num_dtype():
    return [('photoid','i8'),('num','i4')]

def _read_num(filename, rows=None):
    dt = num_dtype()
    stdout.write("Reading photo num file: '%s'\n" % filename)
    r=eu.recfile.Open(filename,dtype=dt,delim=' ')
    data = r.read(rows=rows)
    r.close()
    return data



def num_file(wrun, iteration):
    """
    num output file of calcweights code
    """
    dir = weights_dir(wrun)
    
    conf = zphot.read_config('weights',wrun)
    if iteration == 1:
        n_near = conf['n_near1']
    elif iteration == 2:
        n_near = conf['n_near2']
    else:
        raise ValueError("only support iteration 1 or 2")

    f = 'num-%s-%s.dat' % (wrun,n_near)
    return path_join(dir, f)



def pofz_dtype(nz, with_rmag=False):
    if with_rmag:
        dtype=[('photoid','i8'),('rmag','f4'),('pofz','f4',nz)]
    else:
        dtype=[('photoid','i8'),('pofz','f4',nz)]
    return dtype


def read_pofz(pzrun, chunk, with_rmag=False):
    conf = zphot.read_config('pofz',pzrun)
    nz = conf['nz']
    filename = pofz_file(pzrun, chunk, with_rmag=with_rmag)
    return read_pofz_file(filename, nz, with_rmag=with_rmag)

def read_corrected_pofz(pzrun, chunk, with_rmag=False):
    conf = zphot.read_config('pofz',pzrun)
    nz = conf['nz']
    filename = corrected_pofz_file(pzrun, chunk, with_rmag=with_rmag)
    return read_pofz_file(filename, nz, with_rmag=with_rmag)



def read_pofz_file(filename, nz, with_rmag=False):
    dt = pofz_dtype(nz, with_rmag=with_rmag)

    stdout.write("Reading pofz file: '%s' with dtype:\n" % filename)
    pprint.pprint(dt)
    r = eu.recfile.Open(filename, dtype=dt, delim=' ')
    data = r.read()
    r.close()
    stdout.write("    read %s\n" % data.size)
    return data

def wq_dir(type, runid):
    dir = path_join(zphot.photoz_dir(), 'wq')
    dir = path_join(dir, '%s-%s' % (type,runid))
    return dir

def wq_file(type, runid, chunk=None):
    dir = wq_dir(type, runid)

    if type == 'weights':
        fname = 'weights-%s.yaml' % runid
        fname = path_join(dir,fname)
    elif type == 'pofz':
        if chunk is None:
            raise ValueError("for pofz send chunk")
        fname = 'pofz-%s-%03i.yaml' % (runid,chunk)
        fname = path_join(dir, fname)
    else:
        raise ValueError("'weights' or 'pofz'")

    return fname

def create_weights_wq(wrun):

    conf = zphot.read_config('weights',wrun)
    pprint.pprint(conf)

    train_sample = conf['train_sample']

    # get the input training file
    wt = WeightedTraining(train_sample)
    train_file = wt.filename('all')
    photo_sample = wt.conf['photo_sample']

    # get the photo input file
    zs = zphot.select.ColumnSelector(photo_sample)

    ndim = 5
    if 'extra_columns' in zs.conf:
        ndim += len(zs.conf['extra_columns'])

    photo_file = zs.filename()


    setups = 'setup weighting -r ~esheldon/exports/weighting-work'



    script = """
mode: bynode
group: [new,new2]
job_name: {job_name}
command: |
    source ~/.bashrc

    logf="{logf}"
    photo_file="{photo_file}"
    train_file="{train_file}"

    n_near1={n_near1}
    n_near2={n_near2}

    weights_file1="{weights_file1}"
    weights_file_nozero1="{weights_file_nozero1}"
    weights_file2="{weights_file2}"

    num_file1="{num_file1}"
    num_file2="{num_file2}"

    # first run

    echo "
    Running calcweights

    First run with n_near=$n_near1
    " > "$logf"

    calcweights{ndim}         \\
        "$train_file"    \\
        "$photo_file"    \\
        "$n_near1"       \\
        "$weights_file1" \\
        "$num_file1"  >> "$logf" 2>&1

    if [ "$?" != "0" ]; then
        echo Halting >> "$logf"
        exit 45
    fi


    # remove training set objects with zero weight at first
    # n_near
    echo "
    removing zero weight training objects to
        "$weights_file_nozero1"
    " >> "$logf"
    awk '$3>0' < "$weights_file1" 1> "$weights_file_nozero1" 2>> "$logf"

    if [ "$?" != "0" ]; then
        echo Error running awk.  Halting >> "$logf"
        exit 45
    fi

    # second run
    echo "
    Second run with n_near=$n_near2
    " >> "$logf"

    calcweights{ndim}                \\
        "$weights_file_nozero1" \\
        "$photo_file"           \\
        "$n_near2"              \\
        "$weights_file2"        \\
        "$num_file2" >>"$logf" 2>&1

    if [ "$?" != "0" ]; then
        echo Halting >> "$logf"
        exit 45
    fi
    """.format(job_name             = 'weights-%s' % wrun,
               ndim                 = ndim,
               logf                 = wq_file('weights',wrun).replace('yaml','log'),
               train_file           = train_file,
               photo_file           = photo_file,
               n_near1              = conf['n_near1'],
               n_near2              = conf['n_near2'],
               weights_file1        = weights_file(wrun,1),
               num_file1            = num_file(wrun,1),
               weights_file_nozero1 = weights_file(wrun,1,nozero=True),
               weights_file2        = weights_file(wrun,2),
               num_file2            = num_file(wrun,2))

    # make output dir for calcweights
    wdir = weights_dir(wrun)
    if not os.path.exists(wdir):
        os.makedirs(wdir)

    fname=wq_file('weights',wrun)
    eu.ostools.makedirs_fromfile(fname)
    print("Writing wq config:",fname)
    with open(fname,'w') as fobj:
        fobj.write(script)
    

def create_pofz_wq(pzrun, nchunk):
    """
    For calculating the individual p(z)
    """

    pzdir = pofz_dir(pzrun)
    if not os.path.exists(pzdir):
        os.makedirs(pzdir)

    conf = zphot.read_config('pofz',pzrun)
    pprint.pprint(conf)

    # get output weights file from last round.
    wrun = conf['wrun']
    iteration = 2
    wfile = weights_file(wrun,iteration)

    # to get the photo input file
    wconf = zphot.read_config('weights',wrun)
    train_sample = wconf['train_sample']
    tconf = zphot.read_config('train',train_sample)
    photo_sample = tconf['photo_sample']
    zs = zphot.select.ColumnSelector(photo_sample, nchunk=nchunk)

    ndim = 5
    if 'extra_columns' in zs.conf:
        ndim += len(zs.conf['extra_columns'])

 
    template="""
group: [new,new2]
job_name: {job_name}

command: |
    source ~/.bashrc

    logf="{logf}"
    weights_file="{weights_file}"
    photo_file="{photo_file}"

    n_near={n_near}
    nz={nz}
    zmin={zmin}
    zmax={zmax}
    pofz_file="{pofz_file}"
    z_file="{z_file}"

    echo "
    Running calcpofz
    " > "$logf"

    calcpofz{ndim}           \\
        "$weights_file" \\
        "$photo_file"   \\
        "$n_near"       \\
        "$nz"           \\
        "$zmin"         \\
        "$zmax"         \\
        "$pofz_file"    \\
        "$z_file"  >>"$logf" 2>&1

    if [ "$?" != "0" ]; then
        echo Halting >> "$logf"
        exit 45
    fi
    """
    setups = 'setup weighting -r ~esheldon/exports/weighting-work'

    for chunk in xrange(nchunk):


        script=template.format(job_name = 'pofz-%s-%03i' % (pzrun,chunk),
                               ndim=ndim,
                               logf = wq_file('pofz',pzrun,chunk=chunk).replace('yaml','log'),
                               weights_file=wfile,
                               photo_file=zs.filename(chunk),
                               n_near=wconf['n_near2'],
                               nz=conf['nz'],
                               zmin=conf['zmin'],
                               zmax=conf['zmax'],
                               pofz_file=pofz_file(pzrun,chunk),
                               z_file=z_file(pzrun,chunk))



        fname = wq_file('pofz',pzrun,chunk=chunk)
        eu.ostools.makedirs_fromfile(fname)
        print("Writing wq file:",fname)
        with open(fname,'w') as fobj:
            fobj.write(script)


  
def pbs_dir(type, runid, create=False):
    dir = path_join(zphot.photoz_dir(), 'pbs')
    dir = path_join(dir, '%s-%s' % (type,runid))

    if create:
        if not os.path.exists(dir):
            os.makedirs(dir)
    return dir

def submit_file(type, runid):
    dir = pbs_dir(type, runid)
    f='submit-%s-%s.sh' % (type,runid)
    return path_join(dir,f)

def pbs_file(type, runid, chunk=None, create=False):
    dir = pbs_dir(type, runid, create=create)

    if type == 'weights':
        fname = 'weights-%s.pbs' % runid
        fname = path_join(dir,fname)
    elif type == 'pofz':
        if chunk is None:
            raise ValueError("for pofz send chunk")
        fname = 'pofz-%s-%03i.pbs' % (runid,chunk)
        fname = path_join(dir, fname)
    else:
        raise ValueError("'weights' or 'pofz'")

    return fname

def create_weights_pbs(wrun):

    conf = zphot.read_config('weights',wrun)
    pprint.pprint(conf)

    train_sample = conf['train_sample']

    # get the input training file
    wt = WeightedTraining(train_sample)
    train_file = wt.filename('all')
    photo_sample = wt.conf['photo_sample']

    # get the photo input file
    zs = zphot.select.ColumnSelector(photo_sample)

    ndim = 5
    if 'extra_columns' in zs.conf:
        ndim += len(zs.conf['extra_columns'])

    photo_file = zs.filename()


    setups = 'setup weighting -r ~esheldon/exports/weighting-work'

    wo = WeightedOutputs()

    script = """
photo_file="{photo_file}"
train_file="{train_file}"

n_near1={n_near1}
n_near2={n_near2}

weights_file1="{weights_file1}"
weights_file_nozero1="{weights_file_nozero1}"
weights_file2="{weights_file2}"

num_file1="{num_file1}"
num_file2="{num_file2}"

# first run

echo "
Running calcweights

First run with n_near=$n_near1
" > "$logf"

calcweights{ndim}         \\
    "$train_file"    \\
    "$photo_file"    \\
    "$n_near1"       \\
    "$weights_file1" \\
    "$num_file1"  >> "$logf" 2>&1

if [ "$?" != "0" ]; then
    echo Halting >> "$logf"
    exit 45
fi


# remove training set objects with zero weight at first
# n_near
echo "
removing zero weight training objects to
    "$weights_file_nozero1"
" >> "$logf"
awk '$3>0' < "$weights_file1" 1> "$weights_file_nozero1" 2>> "$logf"

if [ "$?" != "0" ]; then
    echo Error running awk.  Halting >> "$logf"
    exit 45
fi

# second run
echo "
Second run with n_near=$n_near2
" >> "$logf"

calcweights{ndim}                \\
    "$weights_file_nozero1" \\
    "$photo_file"           \\
    "$n_near2"              \\
    "$weights_file2"        \\
    "$num_file2" >>"$logf" 2>&1

if [ "$?" != "0" ]; then
    echo Halting >> "$logf"
    exit 45
fi
    """.format(ndim=ndim,
               train_file          = train_file,
               photo_file          = photo_file,
               n_near1             = conf['n_near1'],
               n_near2             = conf['n_near2'],
               weights_file1        = weights_file(wrun,1),
               num_file1           = wo.num_file(wrun,1),
               weights_file_nozero1 = weights_file(wrun,1,nozero=True),
               weights_file2        = weights_file(wrun,2),
               num_file2           = wo.num_file(wrun,2))

    wdir = weights_dir(wrun)
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    

    # we will use tons of memory
    ppn=8
    job_name='w-%s' % train_sample
    pbsfile = pbs_file('weights',wrun,create=True)
    p =pbs.PBS(pbsfile,
               script,
               job_name=job_name,
               ppn=ppn,
               setups=setups,
               no_redirect=True)
    stdout.write("Writing pbs file: %s\n" % pbsfile)

    p.write()




def create_pofz_pbs(pzrun, nchunk):
    """
    For calculating the individual p(z)
    """

    wo = WeightedOutputs()
    pofz_dir = pofz_dir(pzrun)
    if not os.path.exists(pofz_dir):
        os.makedirs(pofz_dir)

    # make sure pbsdir exists
    pbsdir = pbs_dir('pofz',pzrun,create=True)
    conf = zphot.read_config('pofz',pzrun)
    pprint.pprint(conf)

    # get output weights file from last round.
    wrun = conf['wrun']
    iteration = 2
    weights_file = weights_file(wrun,iteration)

    # to get the photo input file
    wconf = zphot.read_config('weights',wrun)
    train_sample = wconf['train_sample']
    tconf = zphot.read_config('train',train_sample)
    photo_sample = tconf['photo_sample']
    zs = zphot.select.ColumnSelector(photo_sample, nchunk=nchunk)

    subfile = open(submit_file('pofz',pzrun), 'w')
    stdout.write("Writing to submit file: '%s'\n" % subfile.name)

    ndim = 5
    if 'extra_columns' in zs.conf:
        ndim += len(zs.conf['extra_columns'])

 
    script_template="""
weights_file="{weights_file}"
photo_file="{photo_file}"

n_near={n_near}
nz={nz}
zmin={zmin}
zmax={zmax}
pofz_file="{pofz_file}"
z_file="{z_file}"

echo "
Running calcpofz
" > "$logf"

calcpofz{ndim}           \\
    "$weights_file" \\
    "$photo_file"   \\
    "$n_near"       \\
    "$nz"           \\
    "$zmin"         \\
    "$zmax"         \\
    "$pofz_file"    \\
    "$z_file"  >>"$logf" 2>&1

if [ "$?" != "0" ]; then
    echo Halting >> "$logf"
    exit 45
fi
    """
    setups = 'setup weighting -r ~esheldon/exports/weighting-work'

    for chunk in xrange(nchunk):


        script=script_template.format(ndim=ndim,
                                      weights_file=weights_file,
                                      photo_file=zs.filename(chunk),
                                      n_near=wconf['n_near2'],
                                      nz=conf['nz'],
                                      zmin=conf['zmin'],
                                      zmax=conf['zmax'],
                                      pofz_file=pofz_file(pzrun,chunk),
                                      z_file=z_file(pzrun,chunk))



        job_name='pofz-%s-%03i' % (pzrun,chunk)

        pbsfile = pbs_file('pofz',pzrun,chunk=chunk,create=True)
        p =pbs.PBS(pbsfile,
                   script,
                   job_name=job_name,
                   setups=setups,
                   no_redirect=True)
        stdout.write("Writing pbs file: %s\n" % pbsfile)

        p.write()

        b_pbsfile=os.path.basename(pbsfile)
        subfile.write('echo -n "%s "\n' % b_pbsfile)
        subfile.write('qsub "%s"\n' % b_pbsfile)

    stdout.write("Closing submit file: '%s'\n" % subfile.name)
    subfile.close()




class WeightedTraining:
    """
    Create training samples in the right format for the weighting code
    No additional cuts are made here, just formatting.
    """
    def __init__(self, train_sample):
        self.train_sample = train_sample
        self.zt = zphot.training.Training(train_sample)
        self.types = self.zt.types
        self.conf = self.zt.conf
        self.photo_conf = zphot.read_config('zinput', self.conf['photo_sample'])

    def convert_training(self):
        """
        Convert the training sets into an input usable
        by the weighting code
        """

        
        dt = self.get_dtype()

        # this will hold the sum of all
        allfile = self.filename('all')
        stdout.write("Combined file: %s\n" % allfile)
        rall=eu.recfile.Open(allfile, 'w', delim=' ')

        catid = 0
        for type in self.types:

            t=self.zt.read_matched(type)

            newt = numpy.zeros(t.size, dtype=dt)

            newt['z'] = t['z']
            newt['catid'] = catid
            newt['weight'] = 1.0

            newt['cmodelmag_dered_r'] = t['cmodelmag_dered_r']
            newt['model_umg'] = t['modelmag_dered_u']-t['modelmag_dered_g']
            newt['model_gmr'] = t['modelmag_dered_g']-t['modelmag_dered_r']
            newt['model_rmi'] = t['modelmag_dered_r']-t['modelmag_dered_i']
            newt['model_imz'] = t['modelmag_dered_i']-t['modelmag_dered_z']

            if 'extra_columns' in self.photo_conf:
                for cdict in self.photo_conf['extra_columns']:
                    n=cdict['name']
                    newt[n] = t[n]


            newf = self.filename(type)
            print('Writing weighting input:',newf)
            r=eu.recfile.Open(newf, 'w', delim=' ')
            r.write(newt)
            r.close()

            rall.write(newt)

            stdout.write('-'*70 + '\n')
            catid += 1

        stdout.write("\nClosing combined file: %s\n" % allfile)
        rall.close()

    def read(self, type):
        filename = self.filename(type)
        stdout.write("Reading weighting train file: '%s'\n" % filename)
        return self.readfile(filename)

    def readfile(self, filename):
        return read_training(filename)

    def filename(self, type):
        fname=self.zt.fname_matched(type)
        fname = fname.replace('.rec','-wtrain.dat')
        return fname

    def get_dtype(self):
        return training_dtype(self.train_sample)


def simulated_dir(pzrun):
    d=pofz_dir(pzrun)
    d=path_join(d, 'simulated')
    return d

def simulated_nofz_file(pzrun):
    d=simulated_dir(pzrun)
    f='nzsigf-edges-%s.rec' % pzrun
    f = path_join(d,f)
    return f

def read_simulated_nofz(pzrun):
    f=simulated_nofz_file(pzrun)
    data = eu.io.read(f, verbose=True)
    return data


def combine_nofz_and_simulated(pzrun):
    """

    Combined plot of the simulated sample variance errors with our estimate of
    N(z).  Only currently have this for pzrun 12

    Note we use pzrun, but that is just a marker.  We actually plot N(z)
    """
    # first read the weights data
    conf = zphot.read_config('pofz',pzrun)
    iteration=2
    wtrain = read_weights(conf['wrun'],iteration)

    # now get the simulated data
    simdata = read_simulated_nofz(pzrun)
    simnorm = simdata['pofz'].sum()
    simhist = simdata['pofz']/simnorm
    sample_variance = simdata['sample_variance']/simnorm

    # normalize the result

    # perform the weighted histogram on using the bins
    # from the simulated data

    zmin = simdata['zmin'][0]
    zmax = simdata['zmax'][-1]
    binsize = simdata['zmax'][0] - simdata['zmin'][0]
    wdict = eu.stat.histogram(wtrain['z'], 
                              weights=wtrain['weight'],
                              min=zmin, 
                              max=zmax,
                              binsize=binsize)
    """
    eu.misc.colprint(wdict['low'], wdict['high'], wdict['whist'])
    if wdict['whist'].size != simdata['zmin'].size:
        raise ValueError("Expected %d bins but "
                         "got %d" % (simdata['zmin'].size,wdict['whist'].size))
    """

    # we always get an extra bin
    wnorm = wdict['whist'][0:-1].sum()
    whist = wdict['whist'][0:-1]/wnorm
    
    out = numpy.zeros(simdata['zmin'].size,
                      dtype=[('zmin','f8'),('zmax','f8'),('zmid','f8'),
                             ('nofz','f8'),('sample_variance','f8')])
    out['zmin'] = simdata['zmin']
    out['zmax'] = simdata['zmax']
    out['zmid'] = (out['zmax']+out['zmin'])/2.
    out['nofz'] = whist
    out['sample_variance'] = sample_variance
    simf = simulated_nofz_file(pzrun)
    outf = simf.replace('.rec','-with-nofz.rec')
    print("Writing output file:",outf)
    eu.io.write(outf, out)

    # plot with N(z) and sample variance errors
    width=3
    aspect_ratio=0.7

    epsfile=simf.replace('.rec', '-nofz-errors.eps')
    nzeplt=FramedPlot()
    nzeplt.xlabel = 'z'
    nzeplt.ylabel = 'N(z)'
    nzeplt.aspect_ratio=aspect_ratio

    wH = Histogram(whist, x0=out['zmin'][0], binsize=out['zmax'][0]-out['zmin'][0],
                   width=width)
    wH.label = 'N(z)'
    wHerr = SymmetricErrorBarsY(out['zmid'],
                                out['nofz'], 
                                out['sample_variance'], width=width)
    nzerrlab = PlotLabel(0.9,0.9,'Sample Variance Errors', halign='right')


    nzeplt.add(wH,wHerr,nzerrlab)
    nzeplt.show()
    print("Writing plot file:",epsfile)
    nzeplt.write_eps(epsfile)


    # now add simulation
    epsfile=simf.replace('.rec', '-nofz-sim-errors.eps')
    plt=FramedPlot()
    plt.xlabel = 'z'
    plt.aspect_ratio=aspect_ratio

    simH = Points(out['zmid'], simhist, color='blue',type='filled circle',size=2)
    """
    simH = Curve(out['zmid'], simhist, color='blue',type='shortdashed',
                 width=width)
    """
    """
    simH = Histogram(simhist, x0=out['zmin'][0], 
                     binsize=out['zmax'][0]-out['zmin'][0],
                     color='blue', width=width, type='dotted')
    """
    simHerr = SymmetricErrorBarsY((out['zmin']+out['zmax'])/2., 
                                  simhist,
                                  out['sample_variance'], 
                                  color='blue', width=width)
    simH.label = 'simulation'
    key = PlotKey(0.9,0.9,[wH,simH], halign='right')
    plt.add(wH, simH, simHerr, key)
    plt.show()
    print("Writing plot file:",epsfile)
    plt.write_eps(epsfile)


    # two-panel plot with both N(z) and simulated
    epsfile=simf.replace('.rec', '-nofz-sim-errors-2panel.eps')
    biggles.configure('fontsize_min', 2.5)
    arr = FramedArray(2,1)
    arr.aspect_ratio=1
    arr.uniform_limits=1
    arr.xlabel = 'z'
    arr[0,0].add(simH, simHerr, wH)
    arr[0,0].add(key)
    #arr[0,0].add(PlotLabel(0.9,0.9,"simulated",halign='right'))
    arr[1,0].add(wH,wHerr)
    arr[1,0].add(PlotLabel(0.9,0.9,"N(z)+simulated error",halign='right'))

    arr.show()
    print("Writing plot file:",epsfile)
    arr.write_eps(epsfile)

    biggles.configure('fontsize_min', 0.1)



def plot_simulated(pzrun, errcolor='black'):
    """

    Plot carlos' simulated cosmic variance errors for the given sample.

    """

    data = read_simulated_nofz(pzrun)

    plt=FramedPlot()
    plt.aspect_ratio=0.7
    plt.xlabel = 'z'
    plt.ylabel = 'Simulated N(z)'

    epsfile = simulated_nofz_file(pzrun).replace('.rec','.eps')
    if errcolor != 'black':
        epsfile = epsfile.replace('.eps', '-err'+errcolor+'.eps')

    width=3
    h = Histogram(data['pofz'], x0=data['zmin'][0],
                  binsize=data['zmax'][0]-data['zmin'][0],
                  width=width,
                  line='solid')

    zmean = (data['zmax']+data['zmin'])/2.
    err = SymmetricErrorBarsY(zmean, data['pofz'], data['sample_variance'], 
                              width=width,color=errcolor)
    plt.add(h)
    plt.add(err)
    plt.show()

    print("Writing to file:",epsfile)
    plt.write_eps(epsfile)

def plot_6rand_pofz(pzstruct, zmin, binsize, seed=25):
    """

    plot 6 random p(z) from the input struct.  Note typically
    zmin=0.0 and binsize=0.031429

    after each plot is made, raw_input is read.  If the input is 'q' return, if
    the input is 'p' a plot is written to disk in ~/tmp/id-pofz.eps where i
    is the index.

    The input struct must have photoid, p(z) and rmag.  I've been using the
    random sample in pofz-12, which I created by running

        pv pofz-12-chunk0*.dat 0.034 > pofz-12-rand.dat

    then matched to sweep columns to get rmag using

        es_sdsspy.sweeps_collate.match_columns(photoid, 'cmodelmag_dered_r')
        esutil.numpy_util.add_fields(struct, ...)
        copy_fields
        recfile write file

    """
    
    biggles.configure('fontsize_min', 1.5)

    width=5
    w1=where1(pzstruct['rmag'] < 18.)
    w2=where1((pzstruct['rmag'] > 18.) & (pzstruct['rmag'] < 19) )
    w3=where1((pzstruct['rmag'] > 19.) & (pzstruct['rmag'] < 20) )
    w4=where1((pzstruct['rmag'] > 20.) & (pzstruct['rmag'] < 21) )
    w5=where1((pzstruct['rmag'] > 21.) & (pzstruct['rmag'] < 21.5) )
    w6=where1((pzstruct['rmag'] > 21.5) )

    labf = 'r: %0.2f'
    iloop=0
    while True:
        iloop +=1
        # grab 6 random
        if iloop == 1:
            numpy.random.seed(seed)
        i1=w1[numpy.random.randint(0,w1.size)]
        i2=w2[numpy.random.randint(0,w2.size)]
        i3=w3[numpy.random.randint(0,w3.size)]
        i4=w4[numpy.random.randint(0,w4.size)]
        i5=w5[numpy.random.randint(0,w5.size)]
        i6=w6[numpy.random.randint(0,w6.size)]
        
        tab=biggles.FramedArray(2,3)
        tab.uniform_limits=1
        tab.aspect_ratio=0.65
        tab.xlabel='z'
        tab.ylabel='P(z)'
        tab.label_offset=2

        hist = pzstruct['pofz'][i1]/pzstruct['pofz'][i1].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i1], 
                            halign='right')
        tab[0,0].add(h)
        tab[0,0].add(k)

        hist = pzstruct['pofz'][i2]/pzstruct['pofz'][i2].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i2], 
                            halign='right')
        tab[0,1].add(h)
        tab[0,1].add(k)


        hist = pzstruct['pofz'][i3]/pzstruct['pofz'][i3].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i3], 
                            halign='right')
        tab[0,2].add(h)
        tab[0,2].add(k)

        hist = pzstruct['pofz'][i4]/pzstruct['pofz'][i4].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i4], 
                            halign='right')
        tab[1,0].add(h)
        tab[1,0].add(k)

        hist = pzstruct['pofz'][i5]/pzstruct['pofz'][i5].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i5], 
                            halign='right')
        tab[1,1].add(h)
        tab[1,1].add(k)

        hist = pzstruct['pofz'][i6]/pzstruct['pofz'][i6].max()
        h=Histogram(hist, x0=zmin, binsize=binsize, width=width)
        k=biggles.PlotLabel(0.9,0.9,
                            labf % pzstruct['rmag'][i6], 
                            halign='right')
        tab[1,2].add(h)
        tab[1,2].add(k)



        tab.show()


        key=raw_input('hit a key (q quit, p plot): ')
        if key == 'q':
            biggles.configure('fontsize_min',0.1)
            return
        if key == 'p':
            outdir=expand_path('~/tmp/pofz-plots')
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            epsfile=path_join(outdir, 'seed%d-%d-6pofz.eps' % (seed,iloop))

            print("plotting to:",epsfile)
            tab.write_eps(epsfile)

            key=raw_input('hit a key again (q quit): ')
            if key == 'q':
                biggles.configure('fontsize_min',0.1)
                return

    biggles.configure('fontsize_min',0.1)

class WeightedOutputs:
    """
    Just a simple class to work with the outputs
    """
    def __init__(self):
        self.overall_sample = 'dr8'

        pzdir = zphot.photoz_dir()
        self.basedir = path_join(pzdir, 'weighting', self.overall_sample)

    def make_columns(self, pzrun):
        """
        Collate the outputs into a columns database
        Note we use the corrected p(z) here
        """
        conf = zphot.cascade_config(pzrun)

        pprint.pprint(conf)

        coldir = pofz_columns_dir(pzrun)
        if os.path.exists(coldir):
            raise ValueError("coldir exists, start fresh: " + coldir)
        cols = columns.Columns(coldir)
        cols.create()

        # get the config
        nz = conf['pofz']['nz']

        #
        # metadata
        #

        cols.write_column('conf', conf, type='json')

        #
        # the z values of the histogram
        #

        zf = z_file(pzrun,chunk=0)
        zdata = read_z(zf)

        print("Writing zbins column")
        cols.write_column('zbins', zdata)

        pofz_dt = [('photoid','i8'),('pofz','f4',nz)]

        # the num file, to determine of an object is "recoverable"
        iteration=2
        numdata = read_num( conf['pofz']['wrun'],iteration )
        

        pzdir=pofz_dir(pzrun)
        pattern=path_join(pzdir,'pofz-%s-chunk*-corr.dat' % pzrun)
        flist = glob.glob(pattern)
        ntot = len(flist)
        if ntot == 0:
            raise ValueError("No files matched pattern: '%s'" % pattern)

        # ii keeps track of location in the overall num file
        ii = 0
        zhist=None
        for f in sorted(flist):
            print(f)
            r=eu.recfile.Open(f,'r',delim=' ',dtype=pofz_dt)
            data = r.read()
            r.close()
            n = data.size

            # make sure they line up
            wbad,=where(numdata['photoid'][ii:ii+n] != data['photoid'])
            if wbad.size != 0:
                raise ValueError("photoid mismatch for file %s" % f)

            tnum = numdata['num'][ii:ii+n]

            cols.write_column('photoid', data['photoid'])
            cols.write_column('num', tnum)
            cols.write_column('pofz', data['pofz'])


            ii += n

        #cols['photoid'].create_index(verbose=True)
        #cols['num'].create_index(verbose=True)

    def add_columns(self, pzrun):
        """
        Add columns to the collated photoz table
        """

        print("Adding columns\n")
        cols = open_pofz_columns(pzrun)

        # we will also do colors and objid below
        newcols = ['run','rerun','camcol','field','id',
                   'cmodelmag_dered_r','ra','dec','inbadfield']

        print("reading pofz photoid")
        photoid = cols['photoid'][:]

        print("reading sweeps photoid")
        scols = es_sdsspy.sweeps_collate.open_columns('primgal')
        sphotoid = scols['photoid'][:]

        print("matching")
        mcols, mscols = eu.numpy_util.match(photoid, sphotoid)
        if mcols.size != photoid.size:
            raise ValueError("Some did not match: %d/%d" % (mcols.size,photoid.size))

        # first write the index into the primgal sweeps
        matchcol='primgal_id'
        print("    writing",matchcol)
        cols.write_column(matchcol,mscols)

        for col in newcols:
            print("Reading data for: ",col)
            data = scols[col][mscols]
            print("    Writing data")
            if col == 'cmodelmag_dered_r':
                col = 'cmodelmag_r'
            cols.write_column(col, data, create=True)


        print("creating umg")
        u=scols['modelmag_dered_u'][mscols]
        g=scols['modelmag_dered_g'][mscols]
        umg = u-g
        cols.write_column('modelmag_umg',umg,create=True)
        del u

        print("creating gmr")
        r=scols['modelmag_dered_r'][mscols]
        gmr = g-r
        cols.write_column('modelmag_gmr',gmr,create=True)
        del g

        print("creating rmi")
        i=scols['modelmag_dered_i'][mscols]
        rmi = r-i
        cols.write_column('modelmag_rmi',rmi,create=True)
        del r

        print("creating imz")
        z=scols['modelmag_dered_z'][mscols]
        imz = i-z
        cols.write_column('modelmag_imz',imz,create=True)
        del i
        del z

        # don't forget objid
        print("creating objid")
        run=cols['run'][:]
        rerun=cols['rerun'][:]
        camcol=cols['camcol'][:]
        field=cols['field'][:]
        id=cols['id'][:]
        objid = sdsspy.objid(run,rerun,camcol,field,id)
        cols.write_column('objid', objid, create=True)


    def make_release(self, pzrun):
        d = pofz_release_dir(pzrun)
        if not os.path.exists(d):
            print("making release dir:",d)
            os.makedirs(d)

        cols = open_pofz_columns(pzrun)
        print("getting zbins")
        zbindata = cols['zbins'][:]

        outf=zbins_release_file(pzrun)
        print("writing zbins file:",outf)
        eu.io.write(outf, zbindata, clobber=True)

        outf=zbins_release_file(pzrun, ext='dat')
        with Recfile(outf,'w',delim=' ') as fobj:
            print("writing zbins file:",outf)
            fobj.write(zbindata)

        print("reading runs")
        runs = cols['run'][:]
        urun = numpy.unique(runs)

        columns=['objid',
                 'run','rerun','camcol','field','id',
                 'ra','dec',
                 'inbadfield',
                 'cmodelmag_r',
                 'modelmag_umg',
                 'modelmag_gmr',
                 'modelmag_rmi',
                 'modelmag_imz',
                 'pofz']
        i=1
        for run in urun:
            print("run: %d  (%d/%d)" % (run,i,urun.size))
            w,=where(runs == run)

            data = cols.read_columns(columns, rows=w)
            
            fitsfile=pofz_release_file(pzrun, run)
            print("writing fits file:",fitsfile)
            eu.io.write(fitsfile, data, clobber=True)

            datfile=pofz_release_file(pzrun, run, ext='dat')
            with Recfile(datfile,'w',delim=' ') as fobj:
                print("writing dat file:",datfile)
                fobj.write(data)
            del data

            i+=1


    def make_pofz_hist(self, pzrun, num=None, numdata=None, keepflags=None):
        output = self.calculate_pofz_hist(pzrun,num,numdata,keepflags)
        outfile = pofz_hist_file(pzrun)
        stdout.write("Writing summed p(z) hist: '%s'\n" % outfile)
        eu.io.write(outfile, output, delim=' ')

    def calculate_pofz_hist(self, pzrun, num=None, numdata=None, keepflags=None):
        """

        make a combined p(z) histogram from the chunk
        outputs

        send num to only read through files [0,num)

        keepflags is a flag the same length as the sample with 1 for
        keep 0 for not keep

        """
        
        conf = zphot.read_config('pofz',pzrun)
        pprint.pprint(conf)

        wrun = conf['wrun']

        outfile = pofz_hist_file(pzrun)
        out_dt = [('zmin','f8'),('zmax','f8'),('pofz','f8')]

        # get the config
        nz = conf['nz']
        pofz_dt = [('photoid','i8'),('pofz','f8',nz)]

        # first the num file
        iteration=2
        if numdata is None and keepflags is None:
            numdata = read_num(wrun, iteration)

        # the z values for the histogram
        zf = z_file(pzrun,chunk=0)
        zdata = read_z(zf)

        output = numpy.zeros(zdata.size, dtype=out_dt)
        output['zmin'] = zdata['zmin']
        output['zmax'] = zdata['zmax']
        print('nz:',nz)
        print('n(zmin):',output['zmin'].size)
        print('n(pofz):',output['pofz'].size)

        dir=pofz_dir(pzrun)

        pattern=path_join(dir,'pofz-%s-chunk*.dat' % pzrun)
        flist = glob.glob(pattern)
        ntot = len(flist)
        if ntot == 0:
            raise ValueError("No files matched pattern: '%s'" % pattern)
        
        if num is not None:
            if num > ntot:
                raise ValueError("num must be <= %s" % ntot)
        else:
            num=ntot

        flist.sort()

        # ii keeps track of location in the overall num file
        ii = 0
        zhist=None
        for i in xrange(num):
            f=flist[i]
            stdout.write("%s\n" % f)
            r=eu.recfile.Open(f,'r',delim=' ',dtype=pofz_dt)
            data = r.read()
            r.close()
            n = data.size

            if keepflags is not None:
                tflags = keepflags[ii:ii+n]
                w, = where(tflags > 0)
                stdout.write("    keeping: %s/%s\n" % (w.size,data.size))
            else:
                tnum = numdata['num'][ii:ii+n]

                w,=where(tnum != 0)
                stdout.write("    recoverable: %s/%s\n" % (w.size,data.size))

            if w.size > 0:
                if zhist is None:
                    zhist=data['pofz'][0] 
                    if zhist.size != output['pofz'].size:
                        raise ValueError("zhist and pofz different size: "
                                         "%s/%s" % (zhist.size,output['pofz'].size))

                for wi in xrange(w.size):
                    ind = w[wi]
                    zhist[:] += data['pofz'][ind][:]

            ii += n

        zhist = zhist/zhist.sum()

        output['pofz'] = zhist
        return output

    def plot_pofz_hist(self, pzrun):

        outfile = pofz_hist_file(pzrun)
        if not os.path.exists(outfile):
            raise ValueError("Run make_pofz_hist first")
        stdout.write("Reading zhist file: '%s'\n" % outfile)
        output = eu.io.read(outfile)

        epsfile = outfile.replace('.rec','.eps')

        plt=FramedPlot()
        binsize = output['zmax'][1]-output['zmin'][1]
        zhp = biggles.Histogram(output['pofz'],
                                x0=output['zmin'][0],
                                binsize=binsize,
                                width=2)
        plt.add(zhp)
        plt.xlabel='z'
        plt.ylabel='p(z)'
        plt.aspect_ratio=1
        plt.write_eps(epsfile)
        converter.convert(epsfile,dpi=90,verbose=True)


class CompareWeighting:
    """
    Compare the distribution of mag/color for the training
    set, the weighted training set, and the photometric set.
    """

    def __init__(self, wrun, pzrun=None):
        self.wrun=wrun
        self.pzrun=pzrun

        self.conf = zphot.read_config('weights',wrun)
        self.train_conf = zphot.read_config('train',self.conf['train_sample'])
        self.photo_conf = zphot.read_config('zinput', self.train_conf['photo_sample'])

        # ranges for plotting
        self.magmin = 15.0
        self.magmax = 22.2
        self.magbin = 0.1

        self.cmin = -1.0
        self.cmax = 3.0
        self.cbin = 0.05

        # this will probably change with sample
        self.wmin = 0.0
        self.wmax = 1.e-7
        self.wbin = self.wmax/50.0

        self.zmin = 0.0
        self.zmax = 1.1
        self.znbin = 35

        self.weights_data_loaded = False
        self.photo_data_loaded = False


        self.line_width=3

    def zhist(self, dopng=True):

        if not self.weights_data_loaded:
            self.load_weights_data()

        # rows,columns
        tab = biggles.Table( 2, 1 )
            

        wh_plt = FramedPlot()

        biggles.configure('fontsize_min',1.5)

        whist = eu.stat.histogram(self.wtrain['weight'], 
                                  min=self.wmin, 
                                  max=self.wmax, 
                                  binsize=self.wbin)

        pwhist = biggles.Histogram(whist, 
                                   x0=self.wmin, 
                                   binsize=self.wbin, 
                                   width=self.line_width)
        wh_plt.add( pwhist )
        wh_plt.xlabel = 'weights'
        wh_plt.aspect_ratio=1

        psfile = self.psfile('whist')
        print("writing epsfile:",psfile)
        wh_plt.write_eps(psfile)
        if dopng:
            converter.convert(psfile,dpi=90,verbose=True)


        z_plt = FramedPlot()

        # raw zspec of training set
        zbin = (self.zmax-self.zmin)/self.znbin
        htrain = eu.stat.histogram(self.wtrain['z'], 
                                   min=self.zmin, 
                                   max=self.zmax, 
                                   binsize=zbin)
        htrain = htrain/float(htrain.sum())
        p_htrain = biggles.Histogram(htrain, 
                                     x0=self.zmin, 
                                     binsize=zbin,
                                     width=self.line_width, 
                                     color='blue')
        p_htrain.label = 'train'

        # weighted histogram of z from training set
        wdict = eu.stat.histogram(self.wtrain['z'], 
                                  weights=self.wtrain['weight'],
                                  min=self.zmin, 
                                  max=self.zmax, 
                                  binsize=zbin)
        whtrain = wdict['whist']
        whtrain = whtrain/float(whtrain.sum())
        p_whtrain = biggles.Histogram(whtrain, 
                                      x0=self.zmin, 
                                      binsize=zbin, 
                                      color='black',
                                      width=self.line_width)
        p_whtrain.label = 'weighted train'


        # fine weighted histogram of z from training set

        fine_zmin = -0.001
        fine_zmax = 0.015
        fine_zbin = 0.0005
        wdict_fine = eu.stat.histogram(self.wtrain['z'], 
                                       weights=self.wtrain['weight'],
                                       min=fine_zmin,
                                       max=fine_zmax, 
                                       binsize=fine_zbin)
        whtrain_fine = wdict_fine['whist']
        whtrain_fine = whtrain_fine/float(whtrain.sum())
        p_whtrain_fine = biggles.Histogram(whtrain_fine, 
                                           x0=fine_zmin, 
                                           binsize=fine_zbin, 
                                           color='red',
                                           width=self.line_width*2)

        z_plt_inset = FramedPlot()
        z_plt_inset.xlabel = 'z'
        z_plt_inset.add( p_whtrain_fine )



        z_plt.add( p_whtrain )
        inset_lowleft = (0.6,0.6)
        inset_upright = (0.95,0.95)
        z_plt.add( biggles.Inset(inset_lowleft, inset_upright, z_plt_inset) )

        z_plt.xlabel = 'z'
        z_plt.aspect_ratio = 1


        psfile = self.psfile('zhist')
        print("writing epsfile:",psfile)
        z_plt.write_eps(psfile)
        if dopng:
            converter.convert(psfile,dpi=90,verbose=True)


        # also make one with the original included
        keyx = 0.65
        keyy = 0.4
        pk = biggles.PlotKey(keyx,keyy,[p_htrain,p_whtrain])
        z_plt.add( p_htrain )
        z_plt.add( pk )
        psfile = self.psfile('zhist-withorig')
        print("writing epsfile:",psfile)
        z_plt.write_eps(psfile)
        if dopng:
            converter.convert(psfile,dpi=90,verbose=True)

    
        # this version includes the summed p(z) from
        # individual objects
        if self.pzrun is not None:

            z_plt = FramedPlot()
            z_plt.xlabel = 'z'
            z_plt.aspect_ratio = 1

            pzfile = pofz_hist_file(self.pzrun)
            stdout.write("Reading summed zhist file: '%s'\n" % pzfile)
            pzdata = eu.io.read(pzfile)
            binsize = pzdata['zmax'][1]-pzdata['zmin'][1]
            sum_zvals = (pzdata['zmax']+pzdata['zmin'])/2

            pzhist = pzdata['pofz']/pzdata['pofz'].sum()
            p_hsum = Points(sum_zvals, pzhist, type='filled circle',
                            color='red', size=2)
            """
            p_hsum = biggles.Histogram(pzhist, 
                                       x0=0.0, 
                                       binsize=binsize, 
                                       color='grey',
                                       width=self.line_width,
                                       type='solid')
            """
            p_hsum.label = 'summed p(z)'
            pk = biggles.PlotKey(keyx,keyy,[p_htrain,p_whtrain,p_hsum])

            z_plt.add( biggles.Inset(inset_lowleft, inset_upright, z_plt_inset) )
            #z_plt.add(p_htrain, p_whtrain, p_hsum, pk)
            z_plt.add(p_htrain, p_hsum, p_whtrain, pk)
            psfile = self.psfile('zhist-withorig-withsum-%s' % self.pzrun)
            print("writing epsfile:",psfile)
            z_plt.write_eps(psfile)
            if dopng:
                converter.convert(psfile,dpi=90,verbose=True)



        #tab[0,0] = wh_plt
        #tab[1,0] = z_plt

        #stdout.write("Writing zhist eps file: %s\n" % psfile)
        #tab.write_eps(psfile)

        biggles.configure('fontsize_min',0.1)

    def varhist(self, subphoto=None):

        biggles.configure('fontsize_min',1.5)
        if not self.weights_data_loaded:
            self.load_weights_data()

        if not self.photo_data_loaded:
            self.load_photo_data(subphoto=subphoto)

        extra_colname=None
        if 'extra_columns' in self.photo_conf:
            extra_columns = self.photo_conf['extra_columns']
            if len(extra_columns) > 1:
                raise ValueError('can only support 1 extra column for now')
            extra_colname = extra_columns[0]['name']

        a=biggles.Table(3,2)
        a[0,0] = self.varhist1('cmodelmag_dered_r',
                               'r',self.magmin,self.magmax,self.magbin,
                               dokey=True,keyx=0.5,)
        a[0,1] = self.varhist1('model_umg','u-g',-1.0,3.5,0.05)

        a[1,0] = self.varhist1('model_gmr','g-r',-0.2,2.4,0.04)
        a[1,1] = self.varhist1('model_rmi','r-i',-0.4,1.4,0.025)

        a[2,0] = self.varhist1('model_imz','i-z',-0.8,1.1,0.025)#,
                               #dokey=True)

        if extra_colname is None:
            print("whist")
            a[2,1] = self.whist()
        else:
            if extra_colname == 'psf_fwhm_r':
                xmin=0.7
                xmax=2.5
                binsize=0.02
                a[2,1] = self.varhist1(extra_colname,r'$seeing_{r}$',xmin,xmax,
                                       binsize)
            else:
                raise ValueError("only support psf_fwhm_r for extra column for now")

        psfile = self.psfile('varhist')
        stdout.write("Writing var compare eps file: %s\n" % psfile)
        a.write_eps(psfile)
        converter.convert(psfile,dpi=120,verbose=True)
        biggles.configure('fontsize_min',0.1)

    def whist(self):

        biggles.configure('_HalfAxis','ticks_size',2)
        biggles.configure('_HalfAxis','subticks_size',1)

        whdict = eu.stat.histogram(self.wtrain['weight'], 
                                   min=self.wmin, 
                                   max=self.wmax, 
                                   binsize=self.wbin,
                                   more=True)
        #whist = whist/float( whist.sum() )
        whist = whdict['hist']/float( whdict['hist'].sum() )

        #x0 = self.wmin/self.wmax
        #binsize = self.wbin/self.wmax
        #x0 = self.wmin/1.e-7
        #binsize = self.wbin/1.e-7

        pwhist = Curve(whdict['center']/1.e-7, whist, 
                       width=self.line_width*1.5)
        """
        pwhist = biggles.Histogram(whist, 
                                   x0=x0, 
                                   binsize=binsize, 
                                   width=self.line_width*1.5)
        """

        plt = FramedPlot()
        plt.add(pwhist)
        plt.xlabel = r'weight/$10^{-7}$'
        return plt

    def varhist1(self, tagname, xlabel, xmin, xmax, binsize, dokey=False, keyx=0.1):

        print("varhist:",tagname)
        biggles.configure('_HalfAxis','ticks_size',2)
        biggles.configure('_HalfAxis','subticks_size',1)
        thdict = eu.stat.histogram(self.wtrain[tagname],
                                   min=xmin,
                                   max=xmax,
                                   binsize=binsize,
                                   more=True)

        wdict = eu.stat.histogram(self.wtrain[tagname],
                                  min=xmin,
                                  max=xmax,
                                  binsize=binsize,
                                  weights=self.wtrain['weight'])
        phdict = eu.stat.histogram(self.photo[tagname],
                                   min=xmin,
                                   max=xmax,
                                   binsize=binsize,
                                   more=True)

        th=thdict['hist']/float(thdict['hist'].sum())
        wh=wdict['whist']/float(wdict['whist'].sum())
        ph=phdict['hist']/float(phdict['hist'].sum())

        width=self.line_width*1.5
        """
        p_ph = biggles.Histogram(ph, x0=xmin, binsize=binsize,
                                 width=width)
        p_th = biggles.Histogram(th, x0=xmin, color='blue',type='longdashed', 
                                 binsize=binsize,width=width)
        p_wh = biggles.Histogram(wh, x0=xmin, color='red', 
                                 binsize=binsize,width=width )
        """
        """
        p_ph = make_hist_curve(phdict['low'],phdict['high'],ph,
                               width=width)
        p_th = make_hist_curve(thdict['low'],thdict['high'],th,
                               width=width,color='blue',type='longdashed')
        p_wh = make_hist_curve(wdict['low'],wdict['high'],wh,
                               width=width,color='red')
        """
        p_ph = Curve(phdict['center'],ph,width=width,type='solid')
        p_th = Curve(thdict['center'],th,width=width,type='longdashed',color='blue')
        p_wh = Curve(wdict['center'],wh,width=width,type='solid',color='red')
        plt = FramedPlot()
        plt.add(p_ph,p_th,p_wh)
        plt.xlabel = xlabel

        if dokey:
            p_ph.label = 'photo'
            p_th.label = 'train'
            p_wh.label = 'weighted train'
            key=PlotKey(keyx,0.9, [p_ph,p_th,p_wh])
            plt.add(key)
        return plt

    def load_weights_data(self):

        stdout.write('\n')

        iteration=2
        self.wtrain = read_weights(self.wrun,iteration)
        self.weights_data_loaded=True
 

    def load_photo_data(self, subphoto=None):

        iteration=1
        if subphoto is not None:
            stdout.write("    Reading subset of size: %s\n" % len(subphoto))
        num = read_num(self.wrun,iteration, rows=subphoto)

        stdout.write("Selecting photo with num > 0\n")
        w,=where(num['num'] > 0)
        stdout.write("    keeping %s/%s\n" % (w.size,num.size))
        self.num = num[w]

        train_sample = self.conf['train_sample']
        tconf = zphot.read_config('train',train_sample)
        photo_sample = tconf['photo_sample']
        photo = read_photo(photo_sample,rows=subphoto) 
        print("trimming")
        self.photo = photo[w]

        self.photo_data_loaded=True
    
    def psfile(self, type):
        psfile = 'zweight-%s-%s.eps' % (self.wrun, type)
        psfile = path_join(weights_dir(self.wrun),psfile)
        return psfile


def compare_two_zhist(wrun1, label1, wrun2, label2):
    wo=WeightedOutputs()
    iteration=2
    wtrain1 = read_weights(wrun1,iteration)
    wtrain2 = read_weights(wrun2,iteration)
    zmin = 0.0
    zmax = 1.1
    nz = 35

    zbin = (zmax-zmin)/nz

    line_width=2

    plt = FramedPlot()

    color_train='magenta'
    color1='blue'
    color2='red'

    #
    # sample 1
    #

    # unweighted hist
    hist1 = eu.stat.histogram(wtrain1['z'], 
                              min=zmin, 
                              max=zmax,
                              binsize=zbin)
    hist1 = hist1/float(hist1.sum())
    phist1 = biggles.Histogram(hist1, 
                               x0=zmin, 
                               binsize=zbin, 
                               color=color_train,
                               width=line_width)
    phist1.label = 'unweighted'

    # unweighted hist
    wdict1 = eu.stat.histogram(wtrain1['z'], 
                               weights=wtrain1['weight'],
                               min=zmin, 
                               max=zmax,
                               binsize=zbin)
    whist1 = wdict1['whist']
    whist1 = whist1/float(whist1.sum())
    pwhist1 = biggles.Histogram(whist1, 
                                x0=zmin, 
                                binsize=zbin, 
                                color=color1,
                                width=line_width)
    pwhist1.label = 'weighted '+label1

    #
    # sample 2
    #

    # unweighted hist
    junk="""
    hist2 = eu.stat.histogram(wtrain2['z'], 
                              min=zmin, 
                              max=zmax,
                              binsize=zbin)
    hist2 = hist2/float(hist2.sum())
    phist2 = biggles.Histogram(hist2, 
                               x0=zmin, 
                               binsize=zbin, 
                               type='dotted',
                               color=color2,
                               width=line_width)
    phist2.label = 'unweighted '+label2
    """

    # unweighted hist
    wdict2 = eu.stat.histogram(wtrain2['z'], 
                               weights=wtrain2['weight'],
                               min=zmin, 
                               max=zmax,
                               binsize=zbin)
    whist2 = wdict2['whist']
    whist2 = whist2/float(whist2.sum())
    pwhist2 = biggles.Histogram(whist2, 
                                x0=zmin, 
                                binsize=zbin, 
                                color=color2,
                                width=line_width)
    pwhist2.label = 'weighted '+label2

    pk = PlotKey(0.5,0.8,[phist1,pwhist1,pwhist2])
    plt.add( phist1, pwhist1, pwhist2, pk )

    plt.xlabel = 'z'
    plt.ylabel = 'P(z)'

    plt.aspect_ratio = 1


    dir = wo.weights_dir('compare-'+wrun1+'-'+wrun2)
    if not os.path.exists(dir):
        os.makedirs(dir)

    epsfile = 'zweight-compare-%s-%s.eps' % (wrun1,wrun2)
    epsfile = path_join(dir,epsfile)


    plt.write_eps(epsfile)
    converter.convert(epsfile,dpi=90,verbose=True)

def reconstruct_from_superset():
    """
    Take the r < 22.0 individual p(z), limit to r < 21.5 and
    compare the resulting summed p(z)
    """

    wo = WeightedOutputs()

    #
    # the sample limited to r < 21.5
    #

    maglim215 = 21.5
    pzrun215 = '03'
    zhist215_file = pofz_hist_file(pzrun215)
    zhist215 = eu.io.read(zhist215_file)

    #
    # the sample limited to 22.5
    #

    pzrun220 = '04'
    pzconf220 = zphot.read_config('pofz',pzrun220)
    wrun220 = pzconf220['wrun']

    wconf220 = zphot.read_config('weights',wrun220)
    train_sample220 = wconf220['train_sample']
    tconf220 = zphot.read_config('train',train_sample220)
    photo_sample220 = tconf220['photo_sample']
    

    # read in the num 5 data and the original input data
    # and limit to those with n > 0 and r < 21.5

    iteration=1
    num220 = read_num(wrun220,iteration)
    ahelp(num220)
    
    zs = zphot.select.ColumnSelector(photo_sample220)
    f = zs.filename()
    photo220 = read_photo(photo_sample220)
    ahelp(photo220)


    stdout.write("Selecting r < 21.5\n")
    w,=where( (num220['num'] > 0) & (photo220['cmodelmag_dered_r'] < maglim215))
    stdout.write("keeping %s/%s\n" % (w.size, photo220.size) )

    nobjects = photo220.size

    del num220
    del photo220

    keepflags = numpy.zeros(nobjects, dtype='i4')
    keepflags[w] = 1
    zhist220_recon = wo.calculate_pofz_hist(pzrun220, keepflags=keepflags)

    line_width=2
    binsize=zhist215['zmax'][1] - zhist215['zmin'][1]
    phist215 = biggles.Histogram(zhist215['pofz'], 
                                 x0=zhist215['zmin'][0], 
                                 color='blue',
                                 binsize=binsize,
                                 width=line_width)
    phist215.label = 'r < 21.5'

    binsize=zhist220_recon['zmax'][1] - zhist220_recon['zmin'][1]
    phist220_recon = biggles.Histogram(zhist220_recon['pofz'], 
                                       x0=zhist220_recon['zmin'][0], 
                                       binsize=binsize,
                                       color='red',
                                       width=line_width)
    phist220_recon.label = 'reconstructed from r < 22'

    pk = PlotKey(0.5,0.8,[phist215,phist220_recon])
    plt = FramedPlot()
    plt.add( phist215, phist220_recon, pk )

    dir = wo.weights_dir('recon-pofz-'+pzrun215+'-'+pzrun220)
    if not os.path.exists(dir):
        os.makedirs(dir)

    epsfile = 'zweight-reconstruct-pofz-%s-%s.eps' % (pzrun215,pzrun220)
    epsfile = path_join(dir,epsfile)

    plt.write_eps(epsfile)
    converter.convert(epsfile,dpi=90,verbose=True)

def get_frac_lowz(wrun):
    
    wo = WeightedOutputs()
    iter=2
    data = read_weights(wrun,iter)

    zmin = -0.001
    zmax = 1.1
    zbin = 0.0005

    d = eu.stat.histogram(data['z'],
                          weights=data['weight'],
                          min=zmin,
                          max=zmax,
                          binsize=zbin)

    totsum = d['whist'].sum()
    zlowvals = [0.002,0.004,0.008,0.01]
    for zlow in zlowvals:

        w,=where(d['high'] < zlow)

        lowsum = d['whist'][w].sum()

        fraclow = lowsum/totsum
        stdout.write("Fraction with z < %s: %s\n" % (zlow, fraclow))


def make_pofz_correction(pzrun, plot_only=False):
    """
    Get the "correction", the ratio of recovered weighted
    N(z) to the summed p(z)

               N(z)
            -----------
            <sum(p(z))>

    """
    import converter
    from biggles import FramedPlot, Histogram, PlotKey, Table, Points, Curve

    biggles.configure('fontsize_min', 2.0)

    pzconf = zphot.cascade_config(pzrun)

    wo = zphot.weighting.WeightedOutputs()

    zhfile = pofz_hist_file(pzrun)
    if not os.path.exists(zhfile):
        print("generating overall pofz summed hist")
        wo.make_pofz_hist(pzrun)
    sumpofz_struct = eu.io.read(zhfile)

    minz = sumpofz_struct['zmin'][0]
    maxz = sumpofz_struct['zmax'][-1]
    nbin = sumpofz_struct['zmin'].size

    iteration=2
    wtrain = read_weights(pzconf['weights']['wrun'], iteration)

    bs = eu.stat.Binner(wtrain['z'], weights=wtrain['weight'])

    bs.dohist(nbin=nbin, min=minz, max=maxz)
    bs.calc_stats()

    if bs['whist'].size != nbin:
        raise ValueError("whist not same size as summed pofz: %s/%s" % (bs['whist'].size,nbin))


    #eu.numpy_util.aprint(sumpofz_struct, format='%15.8f')


    binsize = bs['high'][0] - bs['low'][0]
    zmin = bs['low'][0]
    sumpofz = sumpofz_struct['pofz']
    wnofz = bs['whist']
    sum_zvals = (sumpofz_struct['zmin']+sumpofz_struct['zmax'])/2
    width=3
    #sumpofz_h = Histogram(sumpofz/sumpofz.sum(), x0=zmin, binsize=binsize, width=width, color='red')
    sumpofz_h = Points(sum_zvals, sumpofz/sumpofz.sum(),
                       type='filled circle',size=2, color='red')
    sumpofz_h.label = 'summed p(z)'
    #wnofz_h = Histogram(wnofz/wnofz.sum(), x0=zmin, binsize=binsize, width=width, color='blue', type='shortdashed')
    wnofz_h = Histogram(wnofz/wnofz.sum(), x0=zmin, binsize=binsize, width=width)
    wnofz_h.label = 'Weighted N(z)'


    sumpofz_norm = sumpofz/sumpofz.sum()
    wnofz_norm = wnofz/wnofz.sum()

    corrall = wnofz_norm/sumpofz_norm

    w=where((bs['center'] > 0.85))
    corr = corrall.copy()
    corr[w] = corrall[w].mean()



    key = PlotKey(0.6,0.9, [sumpofz_h, wnofz_h])
    hplt=FramedPlot()
    hplt.xlabel = 'z'
    hplt.ylabel = 'N(z)'
    hplt.add(sumpofz_h, wnofz_h, key)

    size=3
    cplt = FramedPlot()
    c_corrall = Curve(bs['center'], corrall, type='shortdashed', width=width)
    p_corrall = Points(bs['center'], corrall, type='circle',size=size)
    c_corr = Curve(bs['center'], corr, color='blue', width=width)
    p_corr = Points(bs['center'], corr, color='blue', type='filled diamond',
                    size=size)
    c_corrall.label = 'ratio'
    c_corr.label = 'applied correction'
    key = PlotKey(0.1,0.9,[c_corrall,c_corr])
    cplt.add(c_corrall, p_corrall, c_corr, p_corr, key)
    cplt.add( )
    cplt.xlabel = 'z'
    cplt.ylabel = r'$N(z)/\Sigma p(z)$'


    tab = Table(2,1)
    tab[0,0] = hplt
    tab[1,0] = cplt
    tab.show()

    wo=WeightedOutputs()
    dir=pofz_dir(pzrun)
    epsfile=path_join(dir,'pofz-correct-%s.eps' % pzrun)
    print("Writing eps file:",epsfile)
    tab.write_eps(epsfile)
    converter.convert(epsfile, dpi=120, verbose=True)

    out_dt = [('zmin','f8'),('zmax','f8'),('corr','f8')]
    output = numpy.zeros(corr.size, dtype=out_dt)
    output['zmin'] = sumpofz_struct['zmin']
    output['zmax'] = sumpofz_struct['zmax']
    output['corr'] = corr

    if not plot_only:
        outfile = pofz_correction_file(pzrun)
        print("writing correction to:",outfile)
        eu.io.write(outfile, output, delim=' ')

