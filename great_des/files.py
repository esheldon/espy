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
    gdrun: keyword
        The great des run e.g. nbc-sva1-001
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des', keys['gdrun'], 'data')
    return d

def get_input_file(**keys):
    """
    parameters
    ----------
    gdrun: keyword
        The gdrun e.g. nbc-sva1-001
    ftype: keyword
        The file type, e.g. 'meds' 'truth'
    fnum: keyword
        The file number within given g set
    gnum: keyword
        The g (shear) number set

    noisefree: bool
        If true, return path to noisefree data; meds only.
    """
    d=get_input_dir(**keys)

    noisefree=keys.get("noisefree",False)
    ftype=keys['ftype']

    if noisefree and ftype=='meds':
        fname='nbc2.%(ftype)s.%(fnum)03i.g%(gnum)02i.noisefree.fits'
    else:
        fname='nbc2.%(ftype)s.%(fnum)03i.g%(gnum)02i.fits'

    fname=fname % keys

    fname=os.path.join(d, fname)
    return fname

def get_psf_file(**keys):
    """
    parameters
    ----------
    gdrun: keyword
        The gdrun e.g. nbc-sva1-001
    res: keyword
        The res, default 'lores'
    """

    d=get_input_dir(**keys)

    if 'res' not in keys:
        keys['res'] = 'lores'

    fname='nbc2.psf.%(res)s.fits'
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

def get_collated_dir(**keys):
    """
    parameters
    ----------
    run: keyword
        The processing run
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des',keys['run'], 'collated')
    return d


def get_condor_dir(**keys):
    """
    parameters
    ----------
    run: keyword
        The processing run
    """
    h=os.environ['HOME']
    d=os.path.join(h, 'lensing','great-des',keys['run'],'condor')
    return d

def get_condor_master(**keys):
    """
    parameters
    ----------
    run
    """
    d=get_condor_dir(**keys)

    fname='master.sh'
    fname=os.path.join(d, fname)

    return fname


def get_condor_file(**keys):
    """
    parameters
    ----------
    run
    fnum
    gnum
    start
    end
    """
    d=get_condor_dir(**keys)

    missing=keys.get('missing',False)
    if missing:
        fname='%(run)s-%(filenum)05d-missing.condor'
    else:
        fname='%(run)s-%(filenum)05d.condor'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname



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
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    fnum: int
        File number for given gnum
    start: int
        Integer representing starting object number
    end: int
        Integer representing ending object number
    """
    d=get_output_dir(**keys)

    fname='%(run)s-%(fnum)03i-g%(gnum)02i-%(start)05d-%(end)05d.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def get_collated_file(**keys):
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    """
    d=get_collated_dir(**keys)

    fname='%(run)s-g%(gnum)02i.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def read_collated(**keys):
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    """
    import fitsio

    fname=get_collated_file(**keys)
    print "reading:",fname
    return fitsio.read(fname)

def read_output(**keys):
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    fnum: int, keyword
        File number for given gnum
    start: int, optional keyword
        Integer representing starting object number
    end: int, optional keyword
        Integer representing ending object number
    """

    import fitsio
    fname=get_output_file(**keys)
    print 'reading:',fname
    data = fitsio.read(fname)
    return data


def get_averaged_file(**keys):
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    """

    d=get_collated_dir(**keys)
    
    if 'gnum' in keys:
        fname='%(run)s-g%(gnum)02i-avg.fits'
    else:
        fname='%(run)s-avg.fits'
    fname=fname % keys

    fname=os.path.join(d, fname)

    return fname

def read_averaged(**keys):
    """
    parameters
    ----------
    run: string, keyword
        String representing the run, e.g. nfit-noisefree-04
    gnum: int, keyword
        Integer representing the shear number.
    """

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


