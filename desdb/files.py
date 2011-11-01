import os
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
    def __init__(self, id=None, run=None, band=None, verbose=False, 
                 user=None, password=None):
        """
        Construct either with
            c=Coadd(id=id)
        or with both
            c=Coadd(run=run, band=band)

        The tilename can be inferred (at least for now) from the run
        """
        if (id is None and (run is None or band is None)):
            raise ValueError("Send id= or both run= and band=")
        self['image_id'] = id
        self['cat_id']   = None
        self['run']      = run
        self['band']     = band

        self['verbose']=verbose

        self.user=user
        self.password=password

        self._get_info()
    
    def _get_info(self):
        if self['image_id'] is not None:
            self._get_info_by_id()
        else:
            self._get_info_by_runband()

        df=DESFiles()
        self['image_url'] = df.url('coadd_image', 
                                   run=self['run'], 
                                   tilename=self['tilename'], 
                                   band=self['band'])
        self['cat_url'] = df.url('coadd_cat', 
                                 run=self['run'], 
                                 tilename=self['tilename'], 
                                 band=self['band'])

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

        if self['verbose']:
            show=True
        else:
            show=False
        conn=desdb.Connection(user=self.user,password=self.password)
        res=conn.quick(query,show=show)

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

        if self['verbose']:
            show=True
        else:
            show=False
        conn=desdb.Connection(user=self.user,password=self.password)
        res=conn.quick(query,show=show)

        if len(res) > 1:
            raise ValueError("Expected a single result, found %d")

        for key in res[0]:
            self[key] = res[0][key]

class CoaddList:
    """
    Similar to Coadd but only takes ids and can take a list

    Retrieves the relevant info for the input list.  Get the
    result as a dictionary using get_info(), keyed by id
    """
    def __init__(self, ids, user=None, password=None, verbose=False):

        if not isinstance(ids,(type,list)):
            ids=[ids]

        self.image_ids = ids
        self.idcsv = ','.join([str(id) for id in ids])
        self.verbose=verbose

        self.user=user
        self.password=password

        self._get_info()
    
    def get_info(self):
        return self._info

    def _get_info(self):
        import desdb

        query="""
        select
            im.id,
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
            and im.id in (%(idcsv)s)\n""" % {'idcsv':self.idcsv}

        if self.verbose:
            show=True
        else:
            show=False
        conn=desdb.Connection(user=self.user,password=self.password)
        res=conn.quick(query,show=show)

        df=DESFiles()

        info={}
        for r in res:
            r['image_url'] = df.url('coadd_image', 
                                    run=r['run'], 
                                    tilename=r['tilename'], 
                                    band=r['band'])
            r['cat_url'] = df.url('coadd_cat', 
                                  run=r['run'], 
                                  tilename=r['tilename'], 
                                  band=r['band'])
            info[r['id']] = r

        self._info = info



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
    def __init__(self, root=None):

        if root == 'net':
            self._root = des_net_rootdir()
        elif root is None:
            self._root = des_rootdir()
        else:
            self._root=root

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
        return expand_desvars(url, **keys)


# notes 
#   - .fz might not always hold
#   - EXPOSURENAME can also be built from POINTING-BAND-VISIT
_fs={}
_fs['red_run']   = {'dir':'$DESDATA/red/$RUN/red'}
_fs['red_exp']   = {'dir':'$DESDATA/red/$RUN/red/$EXPOSURENAME'}
_fs['red_image'] = {'dir':_fs['red_exp']['dir'], 'name':'$EXPOSURENAME_$CCD.fits.fz'}
_fs['red_cat']   = {'dir':_fs['red_exp']['dir'], 'name':'$EXPOSURENAME_$CCD_cat.fits'}
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

    # allow keyword root=, desdata=, or fall back to environment variable
    if string.find('$DESDATA') != -1:
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

        string = string.replace('$CCD', str(ccd))

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
    if string.find('$') != -1:
        raise ValueError("There were unexpanded variables: '%s'" % string_in)

    return string

