import os
from numpy import zeros

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

def read_fits_output(**keys):
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

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio

    path=get_output_path(**keys)
    print 'reading:',path
    return fitsio.read(path)



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


