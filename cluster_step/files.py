import os

default_version='2012-10-16'

def get_basedir():
    return os.path.expanduser("~/oh/cluster-step")

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
    repeat:
        The repeat number
    ftype:
        The file type, 'image', 'cat', 'seg'

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    repeat=keys['repeat']
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

    name='psf{psfnum}/{subdir}/im_p{psfnum}_s{shnum}_{repeat}{ext}'
    name=name.format(psfnum=psfnum,subdir=subdir,shnum=shnum,
                     repeat=repeat,ext=ext)

    vdir=get_version_dir(**keys)

    return os.path.join(vdir, name)

def get_output_path(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    repeat:
        The repeat number

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    repeat=keys['repeat']

    vdir=get_version_dir(**keys)

    dir=os.path.join(vdir, 'shear', 'psf%s' % psfnum)

    name='mixmc_im_p{psfnum}_s{shnum}_{repeat}.fits'
    name=name.format(psfnum=psfnum,shnum=shnum, repeat=repeat)

    return os.path.join(dir,name)


def read_image(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    repeat:
        The repeat number

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
    repeat:
        The repeat number

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

    dt=[('id','i4'),('col','f8'),('row','f8'),('ra','f8'),('dec','f8'),
        ('mag_auto_r','f8'),('flags','i4'),('class','f4'),('simid','f4')]
    skiplines=9
    with recfile.Open(path,delim=' ',dtype=dt,skiplines=skiplines) as fobj:
        #data=fobj[:]
        data=fobj.read()

    data['row'] -= 1
    data['col'] -= 1
    return data

def get_config_dir():
    dir=os.environ['ESPY_DIR']
    return os.path.join(dir, 'cluster_step', 'config')

def get_config_path(run):
    dir=get_config_dir()
    return os.path.join(dir, 'run-%s.yaml' % run)

def read_config(run):
    import yaml
    path=get_config_path(run)
    return yaml.load(open(path))
