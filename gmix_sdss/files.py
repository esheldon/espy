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

def get_wq_dir(**keys):
    if ('gmix_run' not in keys
            or 'run' not in keys):
        raise ValueError("send gmix_run=,run=")

    bdir=get_basedir()
    d=os.path.join(bdir, 
                   keys['gmix_run'],
                   'wq',
                   str(keys['run']))
    return d

def get_wq_url(**keys):
    if ('gmix_run' not in keys
            or 'run' not in keys
            or 'camcol' not in keys
            or 'field' not in keys):
        raise ValueError("send gmix_run=,run=, camcol=, field=")

    d=get_wq_dir(**keys)

    fname='%(gmix_run)s-%(run)06d-%(camcol)d-%(field)04d.yaml' % keys

    url=os.path.join(d,fname)
    return url


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
            or 'field' not in keys):
        raise ValueError("send gmix_run=,run=, camcol=, field=")

    d=get_output_dir(**keys)

    fname='%(gmix_run)s-%(run)06d-%(camcol)d-%(field)04d.fits' % keys

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


def get_primary_boss_fields(minscore, doplot=False):
    """
    Get primary fields and trim to the boss area
    """
    import numpy
    import sdsspy
    import es_sdsspy
    win=sdsspy.window.Window()
    mask=es_sdsspy.mangle_masks.load('boss','basic')

    print 'getting primary fields'
    flist = win.get_primary_fields(minscore=minscore)

    print 'getting field limits'
    ra1,dec1 = sdsspy.astrom.gc2eq(flist['mu_start'], flist['nu_start'],
                                   flist['node'],flist['incl'])
    ra2,dec2 = sdsspy.astrom.gc2eq(flist['mu_start'], flist['nu_end'],
                                   flist['node'],flist['incl'])
    ra3,dec3 = sdsspy.astrom.gc2eq(flist['mu_end'], flist['nu_start'],
                                   flist['node'],flist['incl'])
    ra4,dec4 = sdsspy.astrom.gc2eq(flist['mu_end'], flist['nu_end'],
                                   flist['node'],flist['incl'])

    print 'checking mask'
    c1=mask.contains(ra1,dec1)
    c2=mask.contains(ra2,dec2)
    c3=mask.contains(ra3,dec3)
    c4=mask.contains(ra4,dec4)
    
    w,=numpy.where( (c1==1) | (c2==1) | (c3==1) | (c4==1) )

    flist=flist[w]

    if doplot:
        import esutil as eu
        import os
        png_name=os.path.expanduser('~/tmp/boss-primary-flist.png')
        eps_name=os.path.expanduser('~/tmp/boss-primary-flist.eps')

        print 'making plot'
        plt=eu.plotting.bscatter(flist['ra'], flist['dec'], type='dot',
                                 xlabel='RA',ylabel='DEC',
                                 show=False)

        print eps_name
        plt.write_eps(eps_name)
        print png_name
        plt.write_img(2048,1024,png_name)

    return flist
