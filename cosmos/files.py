from sys import stderr
import os

# number of splits to do in processing
NPERSPLIT=30

VERSION='23.5'
def get_dir():
    d=os.environ['COSMOS_DIR']
    vdir='COSMOS_%s_training_sample_sub' % VERSION
    return os.path.join(d, vdir)

def get_cat_path():
    """
    Get the path to the main catalog (not the one with fit parameters)
    """
    d=get_dir()
    path=os.path.join(d, 'real_galaxy_catalog_%s_sub.fits' % VERSION)
    return path

def get_fits_cat_path():
    """
    Get the path to the fit parameters
    """
    d=get_dir()
    path=os.path.join(d,'real_galaxy_catalog_%s_sub_fits.fits' % VERSION)
    return path

def read_cat():
    import fitsio
    path=get_cat_path()
    print 'reading:',path
    return fitsio.read(path,lower=True)

def read_fits_cat():
    import fitsio
    path=get_fits_cat_path()
    print 'reading:',path
    return fitsio.read(path,lower=True)

def get_galsim_catalog():
    import galsim
    path=get_cat_path()
    return galsim.RealGalaxyCatalog(path)

def get_gmix_dir(version):
    """
    The base dir
    """
    d=os.environ['COSMOS_DIR']
    d=os.path.join(d,'gmix-fits', version)
    return d

def get_gmix_output_dir(version):
    """
    Outputs from the code
    """
    d=get_gmix_dir(version)
    d=os.path.join(d, 'outputs')
    return d

def get_gmix_output_path(version, obj_range=None):
    """
    The output file from processing with gmix fits
    """

    d=get_gmix_output_dir(version)
    fname='gmix-cosmos-%s' % version
    if obj_range is not None:
        fname = '%s-%05d-%05d' % (fname, obj_range[0], obj_range[1])

    fname = '%s.fits' % fname
    path=os.path.join(d, fname)
    return path

def get_gmix_condor_path(version):
    """
    The condor submit file path
    """

    d=get_gmix_output_dir(version)
    fname='gmix-cosmos-%s.condor' % version

    path=os.path.join(d, fname)
    return path

def get_splits():
    """
    get the beg:end inclusive splits
    """
    cat=read_cat()
    nobj=cat.size

    nsplit = nobj//NPERSPLIT
    nleft  = nobj % NPERSPLIT

    if nleft > 0:
        nsplit += 1

    beglist=[]
    endlist=[]

    for i in xrange(nsplit):
        beg=i*NPERSPLIT

        if i == nsplit-1:
            end=nobj-1
        else:
            end=(i+1)*NPERSPLIT-1
        
        beglist.append(beg)
        endlist.append(end)
    return beglist, endlist

def get_all_gmix_paths(version):
    """
    Get a list of all split output files
    """
    pathlist=[]

    beglist,endlist = get_splits()

    for rng in zip(beglist,endlist):
        path=get_gmix_output_path(version, obj_range=rng)

        pathlist.append(path)

    return pathlist

def combine_outputs(version):
    """
    Combine into a single output file
    """
    import fitsio
    outfile=get_gmix_output_path(version)
    paths=get_all_gmix_paths(version)

    print >>stderr,'writing:',outfile
    with fitsio.FITS(outfile,'rw',clobber=True) as fobj:
        for i,path in enumerate(paths):
            data=fitsio.read(path)
            if i==0:
                fobj.write(data)
            else:
                fobj[-1].append(data)


def make_condor_text(version):
    """
    Make the condor script for processing
    """

    gmix_dir=get_gmix_dir(version)

    head="""
Universe        = vanilla

Notification    = Never 

# Run this exe with these args
Executable      = {gmix_dir}/master.sh

Image_Size      = 700000

GetEnv = True

kill_sig        = SIGINT

+Experiment     = "astro"

Output          = ./gmix-cosmos-{version}.$(cluster).out
Error           = ./gmix-cosmos-{version}.$(cluster).err
Log             = /data/esheldon/tmp/gmix-cosmos-{version}.$(cluster).log
    """.format(gmix_dir=gmix_dir,
               version=version)

    beglist,endlist = get_splits()
    text_list=[head]

    arg_template="""
+job_name = "gc-%(beg)05d-%(end)05d"
Arguments = %(beg)d %(end)d %(output_file)s
Queue\n"""
    for rng in zip(beglist,endlist):
        output_file=get_gmix_output_path(version, obj_range=rng)

        arg = arg_template % {'beg':rng[0],
                              'end':rng[1],
                              'output_file':output_file}
        text_list.append(arg)

    text = ''.join(text_list) 
    return text

def make_condor_script(version):
    """
    Write the condor script
    """
    gmix_dir=get_gmix_dir(version)
    fname='gmix-cosmos-{version}.condor'.format(version=version)
    path=os.path.join(gmix_dir, fname)

    if not os.path.exists(gmix_dir):
        print >>stderr,'making dir:',gmix_dir
        os.makedirs(gmix_dir)

    text=make_condor_text(version)
    print 'writing:',path
    with open(path,'w') as fobj:
        fobj.write(text)

def make_master(version):
    """
    Make the master executable script
    """
    gmix_dir=get_gmix_dir(version)
    fname='master.sh'
    path=os.path.join(gmix_dir, fname)

    text="""#!/bin/bash
source ~/.bashrc

beg=$1
end=$2
output_file=$3
log_file=${output_file}.log

python ~/python/cosmos/bin/gmix-cosmos.py --obj-range $beg,$end $output_file &> $log_file
    \n"""

    if not os.path.exists(gmix_dir):
        print >>stderr,'making dir:',gmix_dir
        os.makedirs(gmix_dir)

    print 'writing:',path
    with open(path,'w') as fobj:
        fobj.write(text)
    
    os.system('chmod 755 %s' % path)
