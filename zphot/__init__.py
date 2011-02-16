import os
import sys

import select
import training
import weighting
import validation
import rachel

try:
    import cas
except:
    sys.stdout.write("Could not import cas module\n")
    pass

import esutil as eu

def config_dir():
    if 'ESPY_DIR' not in os.environ:
        raise ValueError("espy is not set up")
    espy_dir = os.environ['ESPY_DIR']
    espy_dir = os.path.join(espy_dir,'zphot')
    espy_dir = os.path.join(espy_dir,'config')
    return espy_dir

def photoz_dir():
    if 'PHOTOZ_DIR' not in os.environ:
        raise ValueError("PHOTOZ_DIR not set")
    return os.environ['PHOTOZ_DIR']

def read_config(type, id):
    f = config_file(type,id)
    return eu.io.read(f)

def config_file(type, id):
    dir = config_dir()
    f = '%s-%s.json' % (type,id)
    return os.path.join(dir,f)
 
def cascade_config(pzrun):
    conf={}
    conf['pofz'] = read_config('pofz',pzrun)
    conf['weights'] = read_config('weights',conf['pofz']['wrun'])
    conf['train'] = read_config('train',conf['weights']['train_sample'])
    conf['photo'] = read_config('zinput',conf['train']['photo_sample'])

    return conf

