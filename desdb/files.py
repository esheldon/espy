import os

def get_url(type, root=None, **keys):
    df=DESFiles(root=root)
    return df.url(type, **keys)

get_name=get_url
get_path=get_url

def des_rootdir():
    if 'DESDATA' not in os.environ:
        raise ValueError("The DESDATA environment variable is not set")
    return os.environ['DESDATA']

def des_web_rootdir():
    return 'https://desar.cosmology.illinois.edu/DESFiles/desardata/DES'

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

        If you send root='web', the following is used:
            https://desar.cosmology.illinois.edu/DESFiles/desardata/DES


    Notes
    -----

    - Currently / is used for all path separators.  Good on unix and the web.
    - Because jython is stuck in the dark ages, we can't use nice template
    strings with curly bracked named specifiers.
    
    """
    def __init__(self, root=None):

        if root == 'web':
            self._root = des_web_rootdir()
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
_fs['red_run']   = {'dir':'$DESDATA/red/$RUN/red/'}
_fs['red_exp']   = {'dir':'$DESDATA/red/$RUN/red/$EXPOSURENAME'}
_fs['red_image'] = {'dir':_fs['red_exp']['dir'], 'name':'$EXPOSURENAME_$CCD.fits.fz'}
_fs['red_cat']   = {'dir':_fs['red_exp']['dir'], 'name':'$EXPOSURENAME_$CCD_cat.fits'}
_fs['coadd_run']   = {'dir':'$DESDATA/coadd/$RUN/coadd/'}
_fs['coadd_image'] = {'dir':_fs['coadd_run']['dir'], 'name':'$TILENAME-$BAND.fits.fz'}
_fs['coadd_cat']   = {'dir':_fs['coadd_run']['dir'], 'name':'$TILENAME-$BAND_cat.fits'}

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

    # see if there are any leftover un-expanded variables.  If so
    # raise an exception
    if string.find('$') != -1:
        raise ValueError("There were unexpanded variables: '%s'" % string_in)

    return string

