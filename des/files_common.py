from __future__ import print_function
import os

def config_dir():
    """
    dir holding config files
    """
    d=os.environ['ESPY_DIR']
    d=os.path.join(d, 'des', 'config')
    return d

def config_file(name):
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
    fname=config_file(name)
    print("reading:",fname)
    with open(fname) as fobj:
        data=yaml.load(fobj)
    return data

def cascade_config(run):
    """
    read the run config as well as the lens, source, and cosmo
    configs
    """
    conf=read_config(run)

    conf['lens_conf']=read_config(conf['lcat_vers'])
    conf['source_conf']=read_config(conf['scat_vers'])
    conf['cosmo_conf']=read_config(conf['lens_conf']['cosmo_vers'])

    lc=conf['lens_conf']['cosmo_vers']
    sc=conf['source_conf']['cosmo_vers']
    if lc != sc:
        raise ValueError("cosmo mismatch: '%s' '%s'" % (lc,sc))

    return conf

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
