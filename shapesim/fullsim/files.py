import os
import esutil as eu

def get_config_dir():
    return os.path.join(os.environ['ESPY_DIR'],
                        'shapesim',
                        'fullsim',
                        'config')

def get_config_url(simname):
    d=get_config_dir()
    return os.path.join(d, 'sim-%s.yaml' % simname)

def read_config(simname):
    import yaml
    url=get_config_url(simname)
    return yaml.load(open(url))

def get_original_dir(vers):
    d='/astro/u/esheldon/lensing/catalogs/%s_truth'
    return d % vers

def get_original_url(vers, fnum):
    d=get_original_dir(vers)
    f='Aardvark_v0.5d_truth_des_masked.%s.fit' % fnum
    return os.path.join(d, f)

def read_original(vers, fnum, **keys):
    url=get_original_url(fnum)
    return eu.io.read(url, lower=True, **keys)

def read_original_list(vers, fnum_list=None, **keys):

    if fnum_list is not None:
        flist=[]
        for fnum in fnum_list:
            flist.append( get_original_url(vers, fnum) )
    else:
        import glob
        d=get_original_dir(vers)
        patter=os.path.join(d,'*.fit')
        flist=glob.glob(pattern)

    return eu.io.read(flist, verbose=True, combine=True, lower=True, **keys)

def get_basedir():
    return os.path.join(os.environ['LENSDIR'],
                        'shapesim',
                        'fullsim')

def get_simdir(simname):
    d=get_basedir()
    return os.path.join(d, simname)

def get_pointings_dir(simname):
    d=get_simdir(simname)

    return os.path.join(d, 'pointings')

def get_pointings_url(simname):
    d=get_pointings_dir(simname)
    return os.path.join(d, 'pointings.pickle')
