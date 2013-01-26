# -*- coding: utf-8 -*-
import os
from numpy import zeros, where, sqrt, array, arange

DEFAULT_VERSION='2012-10-16'
PSFNUMS=array([1,2,3,4,5,6])
SHNUMS=array([1,2,3,4,5,6,7,8])
CCDS=arange(1,62+1)



def get_basedir(**keys):
    fs=keys.get('fs','nfs')
    if fs=='nfs':
        return os.environ['CLUSTERSTEP']
    else:
        return os.environ['CLUSTERSTEP_HDFS']

def get_version_dir(**keys):
    version=keys.get('version',None)
    if not version:
        version=DEFAULT_VERSION

    bdir=get_basedir()
    return os.path.join(bdir, version)

def get_prior_dir(**keys):
    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'pofe')
    return dir

def get_gprior_path(**keys):
    dir=get_prior_dir(**keys)
    nosplit=keys.get('nosplit',False)
    if nosplit:
        name='pofe-fits-nosplit.fits' % objtype
    else:
        objtype=keys['type']
        name='pofe-fits-%s.fits' % objtype
    return os.path.join(dir, name)

def read_gprior(**keys):
    import fitsio
    path=get_gprior_path(**keys)
    print 'reading:',path
    return fitsio.read(path)

def get_sprior_path(**keys):
    dir=get_prior_dir(**keys)
    objtype=keys['type']
    name='pofs-fits-%s.fits' % objtype
    return os.path.join(dir, name)

def read_sprior(**keys):
    import fitsio
    path=get_sprior_path(**keys)
    print 'reading:',path
    return fitsio.read(path)



def read_prior_original(**keys):
    from esutil import recfile

    dir=get_prior_dir()
    old=keys.get('old',False)

    if old:
        name='pe_dist.dat'
        # julia writes integers as floating format
        dt0=[('chip','f8'),
             ('simid','f8'),
             ('mag','f8'),
             ('ellipb','f8'),
             ('phib','f8'),
             ('g','f8',2)]
        skiplines=7
    else:
        # julia writes integers as floating format
        dt0=[('chip','f8'),
             ('simid','f8'),
             ('mag','f8'),
             ('ellipb','f8'),
             ('phib','f8'),
             ('g','f8',2),
             ('n','f8'),
             ('scale','f8'),
             ('snr','f8')]

        name='pe_dist_snr.dat'
        skiplines=10

    path=os.path.join(dir,name)

    with recfile.Recfile(path, mode='r', delim=' ', dtype=dt0, skiplines=skiplines) as fobj:
        data=fobj[:]

    return data


def get_input_path(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        The file type, 'image', 'cat', 'seg'

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ccd=keys['ccd']
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

    name='psf{psfnum}/{subdir}/im_p{psfnum}_s{shnum}_{ccd}{ext}'
    name=name.format(psfnum=psfnum,subdir=subdir,shnum=shnum,
                     ccd=ccd,ext=ext)

    vdir=get_version_dir(**keys)

    return os.path.join(vdir, name)


def read_image(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number

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
    ccd:
        The ccd number

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    from esutil import recfile
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

    dt0=[('id','i4'),('col','f8'),('row','f8'),('ra','f8'),('dec','f8'),
        ('mag_auto_r','f8'),('flags','i4'),('class','f4'),('simid','f4')]
    skiplines=9

    verbose=keys.get('verbose',True)
    if verbose:
        print 'reading:',path
    with recfile.Open(path,delim=' ',dtype=dt0,skiplines=skiplines) as fobj:
        #data=fobj[:]
        data0=fobj.read()

    # now fix the types
    dt=[('id','i4'),('col','f8'),('row','f8'),('ra','f8'),('dec','f8'),
        ('mag_auto_r','f8'),('flags','i4'),('class','i4'),('simid','i4')]
    data=zeros(data0.size, dtype=dt)
    for n in data0.dtype.names:
        data[n][:] = data0[n][:].astype(data[n].dtype)

    data['row'] -= 1
    data['col'] -= 1
    return data

def get_config_dir():
    dir=os.environ['ESPY_DIR']
    return os.path.join(dir, 'cluster_step', 'config')

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

def write_fits_output(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    data:
        The data to write
    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...

    header: optional
        optional header to write

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio
    data=keys['data']
    header=keys.get('header',None)

    path=get_output_path(**keys)
    dir=os.path.dirname(path)
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except:
            pass
    print 'writing:',path
    with fitsio.FITS(path,mode='rw',clobber=True) as fobj:
        fobj.write(data, header=header)


class Reader(dict):
    def __init__(self, **keys):
        """
        Default to reading everything

        parameters
        ----------
        run:
            run identifier

        psfnums: optional
            restrict to these psf numbers.  Can be scalar or sequence
        shnums: optional
            restrict to these shear numbers.  Can be scalar or sequence
        ccds: optional
            restrict to these ccd numbers.  Can be scalar or sequence
        setname: string, optional
            The set name, e.g. 'use3'.  Default 'use1'
            Set names represent a set of cuts.
        sratio_range: sequence, optional
            Min and max sratio values.
        s2n_range: sequence, optional
            Min and max s2n values.
        Ts2n_range: sequence, optional
            Min and max Ts2n values.
        Tmean: sequence, optional
            Min and max Tmean range
        mag: sequence, optional
            Min and max mag range
        prob: sequence, optional
            Min and max fit_prob range

        verbose: bool
            Set True to be more verbose
        progress: bool
            Set True to show progress bar
        ignore_missing: bool
            Set to True to ignore missing files.

        """
        
        self._check_keys(**keys)

    def _check_keys(self, **keys):
        if 'run' not in keys:
            raise ValueError("send run=")
        self['ignore_missing']=keys.get('ignore_missing',False)

        for k in keys:
            self[k] = keys[k]

        self['setname'] = keys.get('setname','use1')
        self['verbose']=keys.get('verbose',False)
        self['progress']=keys.get('progress',False)

        psfnums=keys.get('psfnums',None)
        shnums=keys.get('shnums',None)
        ccds=keys.get('ccds',None)

        if psfnums is None:
            self['psfnums']=PSFNUMS.copy()
        else:
            self['psfnums']=array(self['psfnums'],ndmin=1,dtype='i2')

        if shnums is None:
            self['shnums']=SHNUMS.copy()
        else:
            self['shnums']=array(self['shnums'],ndmin=1,dtype='i2')

        if ccds is None:
            self['ccds']=CCDS.copy()
        else:
            self['ccds']=array(self['ccds'],ndmin=1,dtype='i2')

        self['prob_range']=self.get('prob_range',[0,1])

    def _load_data(self):
        from esutil.numpy_util import combine_arrlist
        ntot=len(self['shnums'])*len(self['psfnums'])*len(self['ccds'])

        itot=1
        if self['progress']:
            from progressbar import ProgressBar
            prog=ProgressBar(width=70, color='green')

        datalist=[]
        for shnum in self['shnums']:
            shlist=[]
            for psfnum in self['psfnums']:
                for ccd in self['ccds']:
                    if self['progress']:
                        prog.update(frac=float(itot)/ntot)
                        itot += 1

                    data0=self.read_one(shnum,psfnum,ccd)
                    if data0 is not None:
                        data=self.select(data0)
                        if data is not None:
                            datalist.append(data)
                    del data0

        if len(datalist) == 0:
            if self['verbose']:
                print 'no data read or passed cuts'
            self._data=None
        else:
            self._data=combine_arrlist(datalist)

    def get_data(self):
        if not hasattr(self, '_data'):
            self._load_data()
        return self._data

    def select(self, data0):
        from esutil.numpy_util import strmatch
        from .select import Selector

        selector=Selector()

        logic=selector.get_logic(data0, self['setname'])


        # these can be applied in addition to the set logic
        if 'objtype' in self:
            if self['objtype'] is not None:
                logic=logic & strmatch(data0['model'],self['objtype'])

        for name in ['sratio','s2n','Ts2n','s2','Tmean','mag','prob']:
            rname='%s_range' % name
            if name=='mag':
                name='mag_auto_r'
            if name=='prob':
                name='fit_prob'
            if name=='s2n':
                name='s2n_w'
            if rname in self:
                r=self[rname]
                logic=logic & ((data0[name] > r[0]) & (data0[name] < r[1]))

        w,=where(logic)
        if w.size == 0:
            data=None
        else:
            data=data0[w]
        return data

    def read_one(self, shnum, psfnum, ccd):
        fname=get_output_path(run=self['run'], 
                              psfnum=psfnum, 
                              shnum=shnum, 
                              ccd=ccd, ftype='shear')
        columns=self.get('columns',None)
        if os.path.exists(fname):
            data0=read_fits_output(run=self['run'], 
                                   psfnum=psfnum, 
                                   shnum=shnum, 
                                   ccd=ccd, 
                                   ftype='shear',
                                   columns=columns, 
                                   verbose=self['verbose'])

            if 'ccd' not in data0.dtype.names:
                data=self.collate_one(data0, shnum, psfnum, ccd)
                del data0
            else:
                data=data0
        else:
            if not self['ignore_missing']:
                raise RuntimeError("missing file: %s" % fname)
            data=None

        return data           

    def collate_one(self, data0, shnum, psfnum, ccd):
        from esutil.numpy_util import match
        if self['verbose']:
            print '    collating with catalog'
        cat=read_cat(run=self['run'],
                     shnum=shnum,
                     psfnum=psfnum,
                     ccd=ccd,
                     verbose=self['verbose'])

        extra=[('mag_auto_r','f8'),
               ('shnum','i2'),
               ('psfnum','i2'),
               ('ccd','i2'),
               ('sratio','f8')]
        dt=data0.dtype.descr + extra

        data=zeros(data0.size, dtype=dt)
        for n in data0.dtype.names:
            data[n] = data0[n]

        mcat,mshear=match(cat['id'], data0['id'])
        if mshear.size != data0.size:
            raise ValueError("not all matched!")

        data['mag_auto_r'][mshear] = cat['mag_auto_r'][mcat]
        data['shnum'] = shnum
        data['psfnum'] = psfnum
        data['ccd'] = ccd
        data['sratio'] = sqrt(1./data0['s2'])
        return data

def read_output_set(run, psfnums, shnums, 
                    objtype=None, 
                    s2n_field='s2n_w',
                    s2n_min=None,
                    s2n_max=None,
                    s2_max=None,
                    gsens_min=None,
                    gerr_max=None,
                    columns=None,
                    subtract_mean=False,
                    progress=False):
    """
    Read some data based on the input.
    
    Multiple files may be read. If files are missing they will be skipped

    Note only a single shear number is expected but many psfnums can
    be sent.  
    
    Only those with flags==0 are kept.

    parameters
    ----------
    run: string
        run id
    psfnums: integers
        the psf numbers to read
    shnums: integers
        The shear numbers to read.
    objtype: string, optional
        optionally select only objects with this best-fit model
    columns: optional
        only return these columns
    subtract_mean: bool, optional
        Calculate the mean g and subtract it
    """
    from esutil.numpy_util import strmatch, combine_arrlist
    psfnums=get_psfnums(psfnums)
    shnums=get_shnums(shnums)

    ntot=len(shnums)*len(psfnums)*62

    itot=1
    if progress:
        from progressbar import ProgressBar
        prog=ProgressBar(width=70, color='green')
        #prog=ProgressBar(width=70, color='green', block='▣', empty='□')
        #prog=ProgressBar(width=70, color='green', block='◧', empty='◫')
        #prog=ProgressBar(width=70, color='green', block='■', empty='□')
        #prog=ProgressBar(width=70, color='green', block='=', empty='-')

    datalist=[]
    for shnum in shnums:
        shlist=[]
        for psfnum in psfnums:
            for ccd in xrange(1,62+1):
                if progress:
                    #prog.update(frac=float(itot)/ntot,
                    #            message='%s/%s' % (itot,ntot))
                    prog.update(frac=float(itot)/ntot)
                    itot += 1

                fname=get_output_path(run=run, psfnum=psfnum, shnum=shnum, 
                                      ccd=ccd, ftype='shear')
                if os.path.exists(fname):
                    data0=read_fits_output(run=run, psfnum=psfnum, 
                                           shnum=shnum, ccd=ccd, 
                                           ftype='shear',
                                           columns=columns, 
                                           verbose=False)

                    logic=(data0['flags']==0) | (data0['flags']==65536)
                    if objtype:
                        logic=logic & strmatch(data0['model'],objtype)

                    if s2n_min is not None:
                        logic=logic & (data0[s2n_field] > s2n_min)
                    if s2n_max is not None:
                        logic=logic & (data0[s2n_field] < s2n_max)

                    if s2_max is not None:
                        logic=logic & (data0['s2'] < s2_max)
                    if gsens_min is not None:
                        logic=logic \
                            & (data0['gsens'][:,0] > gsens_min) \
                            & (data0['gsens'][:,1] > gsens_min)
                    if gerr_max is not None:
                        g1err=sqrt(data0['gcov'][:,0,0])
                        g2err=sqrt(data0['gcov'][:,1,1])
                        logic=logic \
                            & (g1err < gerr_max) & (g2err < gerr_max)



                    wkeep,=where(logic)
                    if wkeep.size==0:
                        print 'No objects passed cuts'
                    else:
                        data0=data0[wkeep]
                        shlist.append(data0)

        shdata=combine_arrlist(shlist)

        if subtract_mean:
            g1mean = shdata['g'][:,0].mean()
            g2mean = shdata['g'][:,1].mean()
            shdata['g'][:,0] -= g1mean
            shdata['g'][:,1] -= g2mean
        datalist.append(shdata)

    if len(datalist)==0:
        raise RuntimeError("no outputs were found")
    data=combine_arrlist(datalist)
    return data
 

def read_fits_output(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    import fitsio

    path=get_output_path(**keys)
    verbose=keys.get('verbose',True)
    if verbose:
        print 'reading:',path
    return fitsio.read(path, **keys)



def get_julia_collate_path(**keys):
    """
    All ccds for psf/shear are combined
    """

    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'collate')

    name='%(run)s-p%(psfnum)d-s%(shnum)d.fits'
    name = name % {'run':run,
                   'psfnum':psfnum,
                   'shnum':shnum}

    return os.path.join(dir, name)

def get_output_path(**keys):
    """
    parameters
    ----------

    All keywords to keep things clear

    run:
        run identifier
    psfnum:
        psf number
    shnum:
        The shear number
    ccd:
        The ccd number
    ftype:
        e.g. shear, admom, psf, sizemag, ...
    ext: optional
        The extension,will default to fits
        for the appropriate ftypes

    version: optional
        The version of cluster step, defaults
        to global variable default_version
    """
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ccd='%02d' % int(keys['ccd'])
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'psf%s' % psfnum)

    if ftype in ['admom','psf','shear']:
        ext='fits'
    else:
        ext=keys.get('ext',None)
        if ext is None:
            raise ValueError("send ext= for non-standard file types")

    name='{run}-p{psfnum}-s{shnum}-{ccd}-{ftype}.{ext}'
    name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                     ccd=ccd,ftype=ftype,ext=ext)

    return os.path.join(dir,name)

def get_summary_plot_dir(**keys):
    """
    This is opposed to the plots generated to go with
    the exposure outputs
    """
    run=keys['run']
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'plots', ftype)
    return dir

def get_summary_plot_path(**keys):
    """
    Specifically plots that are from summary data, averaged in some set
    """
    dir=get_summary_plot_dir(**keys)

    run=keys['run']
    ftype=keys['ftype']
    extra=keys.get('extra',None)
    ext=keys.get('ext','eps')

    psfnum=keys.get('psfnum',None)
    shnum=keys.get('shnum',None)

    name='{run}'
    if psfnum is not None:
        name += '-p{psfnum}'
    if shnum is not None:
        name += '-s{shnum}'

    name += '-{ftype}'

    if extra is not None:
        name += '-{extra}'

    name += '.%s' % ext

    name=name.format(run=run,
                     psfnum=psfnum,
                     shnum=shnum,
                     ftype=ftype,
                     extra=extra)

    path=os.path.join(dir, name)
    return path

def get_script_dir(**keys):
    run=keys['run']
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'script', ftype)
    return dir

def get_script_path(**keys):
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ftype=keys['ftype']
    ccd=keys['ccd']

    dir=get_script_dir(**keys)

    ccd='%02d' % int(ccd)

    name='{run}-p{psfnum}-s{shnum}-{ccd}-{ftype}.sh'
    name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                     ccd=ccd,ftype=ftype)

    return os.path.join(dir,name)

def get_pbs_dir(**keys):
    vdir=get_version_dir(**keys)
    run=keys['run']
    ftype=keys['ftype']
    dir=os.path.join(vdir, 'shear', run, 'script')
    return dir

def get_pbs_path(**keys):
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ftype=keys['ftype']


    dir=get_pbs_dir(**keys)

    name='{run}-p{psfnum}-s{shnum}-{ftype}.pbs'
    name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                     ftype=ftype)

    return os.path.join(dir,name)

def get_pbs_all_path(**keys):
    run=keys['run']
    ftype=keys['ftype']


    dir=get_pbs_dir(**keys)

    name='{run}-{ftype}.pbs'
    name=name.format(run=run, ftype=ftype)
    return os.path.join(dir,name)




def get_wq_dir(**keys):
    run=keys['run']
    ftype=keys['ftype']

    vdir=get_version_dir(**keys)
    dir=os.path.join(vdir, 'shear', run, 'wq',ftype)
    return dir


def get_wq_path(**keys):
    run=keys['run']
    psfnum=keys['psfnum']
    shnum=keys['shnum']
    ftype=keys['ftype']

    dir=get_wq_dir(**keys)

    ccd=keys.get('ccd',None)
    if ccd is not None:
        ccd='%02d' % int(ccd)

        name='{run}-p{psfnum}-s{shnum}-{ccd}-{ftype}.yaml'
        name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                         ccd=ccd,ftype=ftype)
    else:
        name='{run}-p{psfnum}-s{shnum}-{ftype}.yaml'
        name=name.format(run=run,psfnum=psfnum,shnum=shnum,
                         ftype=ftype)


    return os.path.join(dir,name)

def get_psfnums(psfnum=None):
    return get_nums(psfnum, PSFNUMS[0], PSFNUMS[-1])

def get_shnums(shnum=None):
    return get_nums(shnum, SHNUMS[0], SHNUMS[-1])

def get_nums(nums, nmin, nmax):
    if nums is None:
        nums=range(nmin, nmax+1)
    elif isinstance(nums,basestring):
        nums=nums.split(',')

    if not isinstance(nums,list):
        nums=[nums]

    nums=[int(s) for s in nums]
    for n in nums:
        if n < nmin or n > nmax:
            raise ValueError("number %d out of range: [%d,%d]" % (n,nmin,nmax))
    return nums


