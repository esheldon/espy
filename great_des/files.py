import os

def get_shear_name_dict(model=None):
    # add as many as you need here
    names=['nuse','shear','shear_cov','shear_err',
           'P','Q','R','flags']

    ndict={}
    if model is not None:
        for n in names:
            name='%s_%s' % (model,n)
            ndict[n] = name
    else:
        for n in names:
            ndict[n] = n
    return ndict

def get_config_file(**keys):
    run=keys['run']

    fname='%s.yaml' % run

    d=os.environ['ESPY_DIR']
    fname=os.path.join(d, 'great_des','config',fname)
    return fname

def read_config(**keys):
    import yaml
    fname=get_config_file(**keys)
    data= yaml.load( open(fname) )

    if data['run'] != keys['run']:
        raise ValueError("run '%s' doesn't "
                         "match '%s'" % (data['run'],keys['run']))
    return data

def get_input_dir(**keys):
    """
    parameters
    ----------
    gdtype: keyword
        The type of run, e.g. nbc2
    gdrun: keyword
        The great des run e.g. 140307_first_test
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des',keys['gdtype'], keys['gdrun'])
    return d

def get_input_file(**keys):
    """
    parameters
    ----------
    gdtype: keyword
        The type of gdrun, e.g. nbc2
    gdrun: keyword
        The gdrun e.g. 140307_first_test
    ftype: keyword
        The file type, e.g. 'meds' 'truth'
    fnum: keyword
        The file number within given g set
    gnum: keyword
        The g (shear) number set
    """
    d=get_input_dir(**keys)

    fname='%(gdtype)s.%(ftype)s.%(fnum)03i.g%(gnum)02i.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)
    return fname

def get_psf_file(**keys):
    """
    parameters
    ----------
    gdtype: keyword
        The type of gdrun, e.g. nbc2
    gdrun: keyword
        The gdrun e.g. 140307_first_test
    res: keyword
        The res, default 'lores'
    """

    d=get_input_dir(**keys)

    if 'res' not in keys:
        keys['res'] = 'lores'

    fname='%(gdtype)s.psf.%(res)s.fits'
    fname=fname % keys

    return os.path.join(d,fname)

def get_output_dir(**keys):
    """
    parameters
    ----------
    run: keyword
        The processing run
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des',keys['run'], 'output')
    return d

def get_wq_dir(**keys):
    """
    parameters
    ----------
    run: keyword
        The processing run
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des',keys['run'],'wq')
    return d

def get_wq_file(**keys):
    d=get_wq_dir(**keys)

    start=keys['start']
    end=keys['end']

    fname='%(run)s-%(fnum)03i-g%(gnum)02i-%(start)05d-%(end)05d.yaml'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def get_output_file(**keys):
    d=get_output_dir(**keys)

    if 'start' in keys:
        fname='%(run)s-%(fnum)03i-g%(gnum)02i-%(start)05d-%(end)05d.fits'
    else:
        fname='%(run)s-g%(gnum)02i.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def read_output(**keys):
    import fitsio
    fname=get_output_file(**keys)
    print 'reading:',fname
    data = fitsio.read(fname)
    return data


def get_averaged_file(**keys):
    d=get_output_dir(**keys)
    
    if 'gnum' in keys:
        fname='%(run)s-g%(gnum)02i-avg.fits'
    else:
        fname='%(run)s-avg.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def read_averaged(**keys):
    import fitsio
    fname=get_averaged_file(**keys)
    print 'reading:',fname
    data = fitsio.read(fname)
    return data


def get_chunk_ranges(**keys):
    import meds
    nper=int(keys['nper'])
    
    keys['ftype']='meds'
    meds_file=get_input_file(**keys)

    m=meds.MEDS(meds_file)

    ntotal=m.size

    nchunks = ntotal/nper
    nleft = ntotal % nper

    if nleft != 0:
        nchunks += 1

    low=[]
    high=[]

    for i in xrange(nchunks):

        low_i = i*nper

        # minus one becuase it is inclusive
        if i == (nchunks-1) and nleft != 0:
            high_i = low_i + nleft -1
        else:
            high_i = low_i + nper  - 1

        low.append( low_i )
        high.append( high_i )

    return low,high
