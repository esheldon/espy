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

def get_lensdir():
    """
    This is the root dir under which all data are stored
    """
    if 'LENSDIR' not in os.environ:
        raise ValueError("LENSDIR is not set")
    return os.environ['LENSDIR']

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
