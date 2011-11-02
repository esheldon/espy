import os
from sys import stderr
import desdb

def get_url(type, root=None, **keys):
    df=DESFiles(root=root)
    return df.url(type, **keys)

get_name=get_url
get_path=get_url

def des_rootdir():
    if 'DESDATA' not in os.environ:
        raise ValueError("The DESDATA environment variable is not set")
    return os.environ['DESDATA']

def des_net_rootdir():
    if 'DESREMOTE' not in os.environ:
        raise ValueError("The DESREMOTE environment variable is not set")
    return os.environ['DESREMOTE']


class Coadd(dict):
    def __init__(self, 
                 id=None, 
                 run=None, band=None, 
                 dataset=None, tilename=None,
                 verbose=False, 
                 user=None, password=None,
                 expandroot=True,
                 conn=None):
        """
        Construct either with
            c=Coadd(id=)
        or
            c=Coadd(run=, band=)
        or
            c=Coadd(dataset=, tilename=, band=)

        The tilename can be inferred (at least for now) from the run

        Sending a connection can speed things up greatly.
        """
        if id is not None:
            self.method='id'
        elif run is not None and band is not None:
            self.method='runband'
        elif (dataset is not None 
              and tilename is not None 
              and band is not None):
            self.method='dataset'
        else:
            raise ValueError("Send id= or (run=,band=) or (dataset=,tilename=,band=")

        self['image_id'] = id
        self['cat_id']   = None
        self['run']      = run
        self['band']     = band
        self['dataset']  = dataset
        self['tilename'] = tilename
        self.expandroot=expandroot

        self.verbose=verbose


        if conn is None:
            self.conn=desdb.Connection(user=user,password=password)
        else:
            self.conn=conn

    def load(self, srclist=False):

        if self.method == 'id':
            self._get_info_by_id()
        elif self.method == 'runband':
            self._get_info_by_runband()
        else:
            self._get_info_by_dataset()

        df=DESFiles()
        self['image_url'] = df.url('coadd_image', 
                                   run=self['run'], 
                                   tilename=self['tilename'], 
                                   band=self['band'],
                                   expandroot=self.expandroot)
        self['cat_url'] = df.url('coadd_cat', 
                                 run=self['run'], 
                                 tilename=self['tilename'], 
                                 band=self['band'],
                                 expandroot=self.expandroot)

        if srclist:
            self._load_srclist()
        

    def _get_info_by_runband(self):
        query="""
        select
            im.id as image_id,
            cat.id as cat_id,
            im.tilename
        from
            coadd im,
            catalog cat 
        where
            cat.catalogtype='coadd_cat'
            and cat.parentid = im.id
            and im.run = '%(run)s'
            and im.band = '%(band)s'\n""" % {'run':self['run'],
                                             'band':self['band']}

        res=self.conn.quick(query,show=self.verbose)

        for key in res[0]:
            self[key] = res[0][key]

    def _get_info_by_id(self):
        query="""
        select
            cat.id as cat_id,
            im.run,
            im.band,
            im.tilename
        from
            coadd im,
            catalog cat
        where
            cat.catalogtype='coadd_cat'
            and cat.parentid = im.id
            and im.id = %(id)s\n""" % {'id':self['image_id']}

        res=self.conn.quick(query,show=self.verbose)

        if len(res) > 1:
            raise ValueError("Expected a single result, found %d")

        for key in res[0]:
            self[key] = res[0][key]

    def _get_info_by_dataset(self):
        query="""
        select
            im.id as image_id,
            cat.id as cat_id,
            im.run
        from
            %(release)s_files cat,
            %(release)s_files im
        where
            cat.filetype='coadd_cat'
            and cat.catalog_parentid = im.id
            and cat.tilename = '%(tile)s'
            and cat.band='%(band)s'\n""" % {'tile':self['tilename'],
                                            'band':self['band'],
                                            'release':self['dataset']}

        res=self.conn.quick(query,show=self.verbose)
        if len(res) != 1:
            raise ValueError("Expected a single result, found %d")

        for key in res[0]:
            self[key] = res[0][key]


    def _load_srclist(self):
        query="""
        SELECT
            image.parentid
        FROM
            image,coadd_src
        WHERE
            coadd_src.coadd_imageid = %d
            AND coadd_src.src_imageid = image.id\n""" % self['image_id']

        res = self.conn.quick(query, show=self.verbose)

        idlist = [str(d['parentid']) for d in res]

        ftype=None
        itmax=5

        i=0 
        while ftype != 'red' and i < itmax:
            idcsv = ', '.join(idlist)

            query="""
            SELECT
                id,
                imagetype,
                parentid
            FROM
                image
            WHERE
                id in (%s)\n""" % idcsv

            res = self.conn.quick(query)
            idlist = [str(d['parentid']) for d in res]
            ftype = res[0]['imagetype']
            
            if self.verbose: stderr.write('ftype: %s\n' % ftype)
            i+=1

        if ftype != 'red':
            raise ValueError("Reach itmax=%s before finding 'red' images. last is %s" % (itmax, ftype))

        if self.verbose: stderr.write("Found %d red images after %d iterations\n" % (len(idlist),i))

        query="""
        select 
            id,run,exposurename,ccd
        from 
            location 
        where 
            id in (%(idcsv)s) 
        order by id\n""" % {'idcsv':idcsv}

        res = self.conn.quick(query)

        df=DESFiles(expandroot=self.expandroot)
        srclist=[]
        for r in res:
            url=df.url('red_image',
                       run=r['run'],
                       expname=r['exposurename'],
                       ccd=r['ccd'])
            r['url'] = url
            srclist.append(r)

        self.srclist=srclist


class DESFiles:
    """
    Generate file urls/paths from filetype, run, etc.

    The returned name is a local path or web url.  The generic
    name "url" is used for both.

    parameters
    ----------
    root: string, optional
        The root for filenames.  Defaults to the DESDATA environment
        variable.

        If you send root='net', the following is used:
            ftp://desar.cosmology.illinois.edu/DESFiles/desardata/DES

    Notes
    -----
    - Currently / is used for all path separators.  Good on unix and the web.
    
    """
    def __init__(self, root=None, expandroot=True):

        if root == 'net':
            self._root = des_net_rootdir()
        elif root is None:
            self._root = des_rootdir()
        else:
            self._root=root

        self.expandroot=expandroot

        if isinstance(root,(str,unicode)):
            self.root_type = 'remote'
        else:
            self.root_type = 'local'

    def root(self):
        return self._root
    
    def dir(self, type=None, **keys):
        if type is None:
            return self.root()

        if type not in _fs:
            raise ValueError("Unsupported path type '%s'" % type)
        
        url = _fs[type]['dir']

        url = self._expand_desvars(url, **keys)
        return url

    def url(self, type=None, **keys):
        if type is None:
            return self.root()

        url = self.dir(type, **keys)
        if 'name' in _fs[type]:
            url = os.path.join(url, _fs[type]['name'])
        url = self._expand_desvars(url, **keys)
        return url
    name=url

    def _expand_desvars(self, url, **keys):
        if 'root' not in keys and 'desdata' not in keys:
            keys['root'] = self._root
        if 'expandroot' not in keys:
            keys['expandroot'] = self.expandroot
        return expand_desvars(url, **keys)


# notes 
#   - .fz might not always hold
#   - EXPOSURENAME can also be built from POINTING-BAND-VISIT
_fs={}
_fs['red_run']   = {'remote_dir':'$DESREMOTE/red/$RUN/red',
                    'dir':         '$DESDATA/red/$RUN/red'}

_fs['red_exp']   = {'remote_dir':'$DESREMOTE/red/$RUN/red/$EXPOSURENAME',
                    'dir':         '$DESDATA/red/$RUN/red/$EXPOSURENAME'}

_fs['red_image'] = {'remote_dir':_fs['red_exp']['remote_dir'],
                    'dir':       _fs['red_exp']['dir'], 
                    'name':'$EXPOSURENAME_$CCD.fits.fz'}

_fs['red_cat']   = {'remote_dir':_fs['red_exp']['remote_dir'],
                    'dir':       _fs['red_exp']['dir'], 
                    'name':'$EXPOSURENAME_$CCD_cat.fits'}

_fs['coadd_run']   = {'remote_dir': '$DESREMOTE/coadd/$RUN/coadd',
                      'dir':        '$DESDATA/coadd/$RUN/coadd'}
_fs['coadd_image'] = {'remote_dir': _fs['coadd_run']['remote_dir'],
                      'dir':        _fs['coadd_run']['dir'], 
                      'name':       '$TILENAME_$BAND.fits.fz'}
_fs['coadd_cat']   = {'remote_dir': _fs['coadd_run']['remote_dir'],
                      'dir':_fs['coadd_run']['dir'], 
                      'name':'$TILENAME_$BAND_cat.fits'}

def expand_desvars(string_in, **keys):

    string=string_in

    expandroot=keys.get('expandroot',True)
    # allow keyword root=, desdata=, or fall back to environment variable
    if string.find('$DESDATA') != -1 and expandroot:
        root=keys.get('root', None)
        if root is None:
            root=keys.get('desdata',None)
        if root is None:
            if 'DESDATA' not in os.environ:
                raise ValueError("send desdata keyword or have DESDATA set for '%s'" % string_in)
            root=os.environ['DESDATA']
        string = string.replace('$DESDATA', str(root))

    if string.find('$DESREMOTE') != -1:
        root=keys.get('root', None)
        if root is None:
            root=keys.get('desremote',None)
        if root is None:
            if 'DESREMOTE' not in os.environ:
                raise ValueError("send desdata keyword or have DESREMOTE set for '%s'" % string_in)
            root=os.environ['DESREMOTE']
        string = string.replace('$DESREMOTE', str(root))



    if string.find('$RUN') != -1:
        run=keys.get('run', None)
        if run is None:
            raise ValueError("run keyword must be sent: '%s'" % string_in)
        string = string.replace('$RUN', str(run))

    if string.find('$EXPOSURENAME') != -1:
        expname=keys.get('expname', None)
        if expname is None:
            if 'pointing' in keys and 'band' in keys and 'visit' in keys:
                expname='%s-%s-%s' % (keys['pointing'],keys['band'],keys['visit'])
        if expname is None:
            raise ValueError("expname keyword or pointing,band,visit keywords "
                             "must be sent: '%s'" % string_in)

        string = string.replace('$EXPOSURENAME', str(expname))

    if string.find('$CCD') != -1:
        ccd=keys.get('ccd', None)
        if ccd is None:
            raise ValueError("ccd keyword must be sent: '%s'" % string_in)
        ccd=int(ccd)

        string = string.replace('$CCD', '%02i' % ccd)

    if string.find('$BAND') != -1:
        band=keys.get('band', None)
        if band is None:
            raise ValueError("band keyword must be sent: '%s'" % string_in)

        string = string.replace('$BAND', str(band))


    if string.find('$TILENAME') != -1:
        run=keys.get('tilename', None)
        if run is None:
            raise ValueError("run keyword must be sent: '%s'" % string_in)
        string = string.replace('$TILENAME', str(run))



    # see if there are any leftover un-expanded variables.  If so
    # raise an exception
    if string.find('$') != -1 and expandroot:
        raise ValueError("There were unexpanded variables: '%s'" % string_in)

    return string

