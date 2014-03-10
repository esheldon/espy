import os

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

    start=keys['start']
    end=keys['end']

    fname='%(run)s-%(fnum)03i-g%(gnum)02i-%(start)05d-%(end)05d.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname



def get_chunk_ranges(**keys):
    import meds
    nper=int(keys['nper'])
    
    keys['ftype']='meds'
    meds_file=get_input_file(**keys)

    m=meds.MEDS(meds_file)

    ntotal=m.size

    nchunks = ntotal/nper
    if (nchunks % nper) != 0:
        nchunks += 1

    low=[]
    high=[]

    for i in xrange(nchunks):

        low.append( i*nper )

        # minus one becuase it is inclusive
        high.append( (i+1)*nper -1 )

    return low,high
