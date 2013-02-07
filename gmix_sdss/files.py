import os

def get_config_dir():
    dir=os.environ['ESPY_DIR']
    return os.path.join(dir, 'gmix_sdss', 'config')

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


def get_basedir():
    if 'GMIX_SDSS' not in os.environ:
        raise ValueError("set GMIX_SDSS environment variable")
    return os.environ['GMIX_SDSS']

def get_output_dir(**keys):
    if ('gmix_run' not in keys
            or 'run' not in keys
            or 'camcol' not in keys):
        raise ValueError("send gmix_run=,run=, camcol=")

    bdir=get_basedir()
    d=os.path.join(bdir, 
                   keys['gmix_run'],
                   str(keys['run']),
                   str(keys['camcol']))
    return d

def get_output_url(**keys):
    if ('gmix_run' not in keys
            or 'run' not in keys
            or 'camcol' not in keys
            or 'field' not in keys
            or 'filter' not in keys):
        raise ValueError("send gmix_run=,run=, camcol=, field=, filter=")

    d=get_output_dir(**keys)

    fname='%(gmix_run)s-%(run)06d-%(camcol)d-%(field)04d-%(filter)s.fits' % keys

    url=os.path.join(d,fname)
    return url

def read_output(**keys):
    import fitsio

    url=get_output_url(**keys)

    if not os.path.exists(url):
        raise IOError("file not found: %s" % url)

    fobj=fitsio.FITS(url)
    if len(fobj)==1:
        # No objects were selected for this field.  This is OK
        return numpy.array([])

    return fobj[1].read()

def write_output(**keys):
    import esutil as eu
    if 'data' not in keys:
        raise ValueError("send data=")

    d=get_output_dir(**keys)
    url=get_output_url(**keys)

    if not os.path.exists(d):
        try:
            os.makedirs(d)
        except:
            pass

    print 'writing:',url
    eu.io.write(url, keys['data'], clobber=True)
