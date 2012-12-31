import os
from numpy import zeros, where, sqrt

default_version='2012-10-16'
psfnums=[1,2,3,4,5,6]
shnums=[1,2,3,4,5,6,7,8]
ccds=range(1,62+1)

def get_basedir(**keys):
    fs=keys.get('fs','nfs')
    if fs=='nfs':
        return os.environ['CLUSTERSTEP']
    else:
        return os.environ['CLUSTERSTEP_HDFS']

def get_version_dir(**keys):
    version=keys.get('version',None)
    if not version:
        version=default_version

    bdir=get_basedir()
    return os.path.join(bdir, version)

def get_input_path(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        The file type, 'image', 'cat', 'seg'

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ccd=keys['ccd']
    ftype=keys['ftype']

    if ftype=='cat':
        subdir='cats'
        ext='.cat'
    elif ftype=='image':
        subdir='fits'
        ext='.fits'
    elif ftype=='seg':
        subdir='seg.fits'
        ext='.seg.fits'
    else:
        raise ValueError("bad ftype: '%s'" % ftype)

    name='psf{psfnum}/{subdir}/im_p{psfnum}_s{shnum}_{ccd}{ext}'
    name=name.format(psfnum=psfnum,subdir=subdir,shnum=shnum,
                     ccd=ccd,ext=ext)

    vdir=get_version_dir(**keys)

    return os.path.join(vdir, name)


def read_image(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number

    ftype: optional
        send ftype='seg' to read the segmentation map

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio

    if 'ftype' not in keys:
        keys['ftype']='image'
    path=get_input_path(**keys)
    data,hdr=fitsio.read(path,header=True)
    return data, hdr

def read_cat(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear.

    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import recfile
    keys['ftype']='cat'
    path=get_input_path(**keys)

    # 1 ID 
    # 2 XCENTROID_IMAGE 
    # 3 YCENTROID_IMAGE 
    # 4 RA_IMAGE 
    # 5 DEC_IMAGE 
    # 6 R-MAG-AUTO 
    # 9 FLAGS 
    # 7 CLASS 
    # 8 SIMID  
    #2        1653.320        27.031        334.8027607        -41.4594665        21.7761        0        1.0        7254574.0  

    dt0=[('id','i4'),('col','f8'),('row','f8'),('ra','f8'),('dec','f8'),
        ('mag_auto_r','f8'),('flags','i4'),('class','f4'),('simid','f4')]
    skiplines=9
    with recfile.Open(path,delim=' ',dtype=dt0,skiplines=skiplines) as fobj:
        #data=fobj[:]
        data0=fobj.read()

    # now fix the types
    dt=[('id','i4'),('col','f8'),('row','f8'),('ra','f8'),('dec','f8'),
        ('mag_auto_r','f8'),('flags','i4'),('class','i4'),('simid','i4')]
    data=zeros(data0.size, dtype=dt)
    for n in data0.dtype.names:
        data[n][:] = data0[n][:].astype(data[n].dtype)

    data['row'] -= 1
    data['col'] -= 1
    return data

def get_config_dir():
    dir=os.environ['ESPY_DIR']
    return os.path.join(dir, 'cluster_step', 'config')

def get_config_path(run):
    dir=get_config_dir()
    return os.path.join(dir, '%s.yaml' % run)

def read_config(run):
    import yaml
    path=get_config_path(run)
    conf=yaml.load(open(path))
    if conf['run'] != run:
        mess="run does not match itself in config: %s instead of  %s"
        mess=mess % (conf['run'],run)
        raise ValueError(mess)
    return conf

def write_fits_output(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    data:
        The data to write
    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...

    header: optional
        optional header to write

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio
    data=keys['data']
    header=keys.get('header',None)

    path=get_output_path(**keys)
    dir=os.path.dirname(path)
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except:
            pass
    print 'writing:',path
    with fitsio.FITS(path,mode='rw',clobber=True) as fobj:
        fobj.write(data, header=header)


def read_output_set(run, psfnums, shnums, 
                    objtype=None, 
                    s2n_min=None,
                    s2_max=None,
                    gsens_min=None,
                    gerr_max=None,
                    columns=None,
                    subtract_mean=False,
                    progress=False):
    """
    Read some data based on the input.
    
    Multiple files may be read. If files are missing they will be skipped

    Note only a single shear number is expected but many psfnums can
    be sent.  
    
    Only those with flags==0 are kept.

    parameters
    ----------
    run: string
        run id
    psfnums: integers
        the psf numbers to read
    shnums: integers
        The shear numbers to read.
    objtype: string, optional
        optionally select only objects with this best-fit model
    columns: optional
        only return these columns
    subtract_mean: bool, optional
        Calculate the mean g and subtract it
    """
    from esutil.numpy_util import strmatch, combine_arrlist
    psfnums=get_psfnums(psfnums)
    shnums=get_shnums(shnums)

    ntot=len(shnums)*len(psfnums)*62

    itot=1
    if progress:
        from progressbar import ProgressBar
        prog=ProgressBar(width=70, color='green')

    datalist=[]
    for shnum in shnums:
        shlist=[]
        for psfnum in psfnums:
            for ccd in xrange(1,62+1):
                if progress:
                    #prog.update(frac=float(itot)/ntot,
                    #            message='%s/%s' % (itot,ntot))
                    prog.update(frac=float(itot)/ntot)
                    itot += 1

                fname=get_output_path(run=run, psfnum=psfnum, shnum=shnum, 
                                      ccd=ccd, ftype='shear')
                if os.path.exists(fname):
                    data0=read_fits_output(run=run, psfnum=psfnum, 
                                           shnum=shnum, ccd=ccd, 
                                           ftype='shear',
                                           columns=columns, 
                                           verbose=False)

                    logic=data0['flags']==0
                    if objtype:
                        logic=logic & strmatch(data0['model'],objtype)
                    if s2n_min is not None:
                        logic=logic & (data0['s2n_w'] > s2n_min)
                    if s2_max is not None:
                        logic=logic & (data0['s2'] < s2_max)
                    if gsens_min is not None:
                        logic=logic \
                            & (data0['gsens'][:,0] > gsens_min) \
                            & (data0['gsens'][:,1] > gsens_min)
                    if gerr_max is not None:
                        g1err=sqrt(data0['gcov'][:,0,0])
                        g2err=sqrt(data0['gcov'][:,1,1])
                        logic=logic \
                            & (g1err < gerr_max) & (g2err < gerr_max)



                    wkeep,=where(logic)
                    data0=data0[wkeep]
                    shlist.append(data0)
        shdata=combine_arrlist(shlist)

        if subtract_mean:
            g1mean = shdata['g'][:,0].mean()
            g2mean = shdata['g'][:,1].mean()
            shdata['g'][:,0] -= g1mean
            shdata['g'][:,1] -= g2mean
        datalist.append(shdata)

    if len(datalist)==0:
        raise RuntimeError("no outputs were found")
    data=combine_arrlist(datalist)
    return data
 

def read_fits_output(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio

    path=get_output_path(**keys)
    verbose=keys.get('verbose',True)
    if verbose:
        print 'reading:',path
    return fitsio.read(path, **keys)



def get_output_path(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...
    ext: optional
        The extension,will default to fits
        for the appropriate ftypes

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ccd='%02d' % int(keys['ccd'])
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'psf%s' % psfnum)

    if ftype in ['admom','psf','shear']:
        ext='fits'
    else:
        ext=keys.get('ext',None)
        if ext is None:
            raise ValueError("send ext= for non-standard file types")

    name='{run}-p{psfnum}-s{shnum}-{ccd}-{ftype}.{ext}'
    name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                     ccd=ccd,ftype=ftype,ext=ext)

    return os.path.join(dir,name)

def get_summary_plot_dir(**keys):
    """
    This is opposed to the plots generated to go with
    the exposure outputs
    """
    run=keys['run']
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'plots', ftype)
    return dir

def get_summary_plot_path(**keys):
    """
    Specifically plots that are from summary data, averaged in some set
    """
    dir=get_summary_plot_dir(**keys)

    run=keys['run']
    ftype=keys['ftype']
    extra=keys.get('extra',None)
    ext=keys.get('ext','eps')

    psfnum=keys.get('psfnum',None)
    shnum=keys.get('shnum',None)

    name='{run}'
    if psfnum is not None:
        name += '-p{psfnum}'
    if shnum is not None:
        name += '-s{shnum}'

    name += '-{ftype}'

    if extra is not None:
        name += '-{extra}'

    name += '.%s' % ext

    name=name.format(run=run,
                     psfnum=psfnum,
                     shnum=shnum,
                     ftype=ftype,
                     extra=extra)

    path=os.path.join(dir, name)
    return path

def get_wq_dir(**keys):
    run=keys['run']
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'wq',ftype)
    return dir

def get_wq_path(**keys):
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ftype=keys['ftype']

    dir=get_wq_dir(**keys)

    if 'ccd' in keys:
        ccd='%02d' % int(keys['ccd'])

        name='{run}-p{psfnum}-s{shnum}-{ccd}-{ftype}.yaml'
        name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                         ccd=ccd,ftype=ftype)
    else:
        name='{run}-p{psfnum}-s{shnum}-{ftype}.yaml'
        name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                         ftype=ftype)


    return os.path.join(dir,name)

def get_psfnums(psfnum=None):
    return get_nums(psfnum, 1, 6)

def get_shnums(shnum=None):
    return get_nums(shnum, 1, 8)

def get_nums(nums, nmin, nmax):
    if nums is None:
        nums=range(nmin, nmax+1)
    elif isinstance(nums,basestring):
        nums=nums.split(',')

    if not isinstance(nums,list):
        nums=[nums]

    nums=[int(s) for s in nums]
    for n in nums:
        if n < nmin or n > nmax:
            raise ValueError("number %d out of range: [%d,%d]" % (n,nmin,nmax))
    return nums


