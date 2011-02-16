"""
    Module:
        sdsspy.files

    Description:
        A set of functions for dealing with SDSS files, including reading files,
        finding file lists, getting run metadata, converting SDSS id numbers 
        to strings, etc.

    Classes:
        FileList:
            A class for genering file lists.  This reads actual
            files from disk.  If you want to do more complex
            selections, such as using the resolve and calibration
            information, see the sdsspy.window module.


    Documented Functions:
            read(type, 
                 run=, 
                 rerun=301, 
                 camcol='*', 
                 field='*', 
                 band='*', 
                 type=,
                 nocombine=)

                Read from SDSS files.

            filename(type, 
                     run=, 
                     rerun=301, 
                     camcol='*', 
                     field='*', 
                     band='*', 
                     add_dir=)

                Generate SDSS file names from id info.

            filedir(subdir, 
                    run=, 
                    rerun=301, 
                    camcol='*', 
                    field='*',
                    band='*')

                    
                Generate an SDSS directory name from id info.

        Routines for converting SDSS id numbers to strings:
            These functions return scalars for scalar input, otherwise lists. 

            int2string(numbers, min, max)
                Convert numbers to strings with zfill between min and max
            stripe2string(stripes)
                Convert input stripe(s) to string(s)
            run2string(runs)
                Convert input run(s) to string(s)
            rerun2string(reruns)
                Convert input rerun(s) to string(s)
            camcol2string(camcols, band=None)
                Convert input camcol(s) to string(s). If band is sent it
                is incorporated into the string as in SDSS files.
            field2string(fields)
                Convert input field(s) to string(s)
            id2string(ids)
                Convert input id(s) to string(s)
            band2string(bands)
                Return alphabetic representation of band; works on both string
                and numerical input.

    Dependencies:
        numpy
        esutil


    Modification History:
        Documented: 2007-06-01, Erin Sheldon
"""
import os
import sys
from sys import stdout

import re
import glob

import numpy

import esutil
from esutil.ostools import path_join


import sdsspy

_DEFAULT_RERUN=301

__gflcache = {}
debug=False

# list the required elements for each file type
filespec ={}


filespec['calibobj'] = \
    {'dir':'$PHOTO_SWEEP/$RERUN',
     'check_compress': True,
     'byfield': False,
     'hdu': 1,
     'name':'calibObj-$RUN-$CAMCOL-$TYPE.fits'}





       
def read(filetype, 
         run=None, 
         rerun=_DEFAULT_RERUN, 
         camcol=None,
         field='*', 
         **keywords):
    """
    res = sdsspy.files.read(filetype, 
                            run=None, 
                            rerun=301,
                            camcol=None, 
                            field='*',
                            columns=None,
                            combine=True, 
                            verbose=False, 
                            header=False, 
                            ensure_native=False,
                            **file_keywords)

    Inputs:
        filetype: File type, e.g. 'calibobj'
    Keywords:
        run: 
            The SDSS run.  Most files require this.
        rerun: 
            The SDSS rerun.  Default is 301
        camcol: 
            Camera column.  Required for many SDSS files.  Can be '*'
        field: 
            The field to read.  Can be '*' or a sequence or array.  Default is
            '*' meaning "read all fields"
        band: 
            bandpass. Required for some files.  Can also be '*'
        type: 
            The type, e.g. 'gal' for calibObj.  Can also be '*'
        combine: 
            Results are initially kept in a list of arrays, with each list
            element corresponding to a single file read.  By default these are
            combined into one big array.  If this keyword is False this
            combination step is not performed, saving some memory but making
            use difficult.

        rows: 
            A subset of rows to read from a file.  Ignored if reading multiple
            files.
        columns: 
            A subset of columns to read.
        ensure_native:
            For some file formats, e.g. .fits and .rec, you can force the 
            byte ordering to native by setting this keyword.
        verbose:  
            Print the names of files being read.
    Revision History:
        Documented: 2007-06-01, Erin Sheldon, NYU
    """

    # get some keywords
    rows=keywords.get('rows',None)
    columns=keywords.get('columns',None)
    verbose=keywords.get('verbose',False)
    combine=keywords.get('combine',True)
    ensure_native=keywords.get('ensure_native',False)
    
    # Note, we always get the full field, if this filetype is split by fields.
    # then we can extract a subset later

    ftype=filetype.lower()
    fl = FileList(ftype, run, rerun=rerun, camcol=camcol, **keywords)
    flist = fl.flist(subfields=field)
    
    data = esutil.io.read(flist, 
                          ext=filespec[ftype]['hdu'],
                          rows=rows,
                          columns=columns, 
                          view=numpy.ndarray, 
                          verbose=verbose, 
                          lower=True,
                          ensure_native=ensure_native,
                          combine=combine)

    # for some types we can split by field after the fact,e.g. calibobj
    # only try this if the data on disk are not split by field already.
    post_fieldsplit_types = ['calibobj']
    if field != '*' \
            and not filespec[ftype]['byfield'] \
            and ftype in post_fieldsplit_types:

        data = extract_subfields(data, field, verbose=verbose)

    return data


def extract_subfields(data, fields, verbose=False):
    if 'field' in data.dtype.names:
        if verbose:
            stdout.write("Extracting a subset of fields\n")
        h,rev=esutil.stat.histogram(data['field'],min=0,rev=True)

        f2get = ensure_sequence(fields)
        keep = numpy.zeros(data.size)
        for f in f2get:
            if rev[f] != rev[f+1]:
                w=rev[ rev[f]:rev[f+1] ]
                keep[w] = 1
        wkeep, = numpy.where(keep != 0)
        if wkeep.size == 0:
            raise ValueError("None of the requested fields matched")

        data_keep = data[wkeep]
    else:
        data_keep = data

    return data_keep



def filedir(filetype, run=None, rerun=_DEFAULT_RERUN, camcol='*', **keywords):
    """
    dir = sdsspy.files.filedir(filetype, run=None, **keywords)

    E.g. For calibobj, all that is required is the filetype
    and rerun.

       filedir = filedir('calibobj', rerun=301)

    Note also, the rerun can be omitted if using the default
    which is 301

       filedir = filedir('calibobj')

    Keywords:
        run
        rerun
        camcol
        field
        id
        band
        type

    """

    ftype=filetype.lower()
    if ftype not in filespec:
        raise ValueError("Unsupported file type: '%s'" % ftype)

    dir = expand_sdssvars(filespec[ftype]['dir'], 
                          run=run, 
                          rerun=rerun, 
                          camcol=camcol,
                          **keywords)
    return dir




def filename(filetype_input, 
             run=None,
             rerun=_DEFAULT_RERUN,
             camcol='*', 
             field='*',
             add_dir=True, 
             **keywords):
    """
    name = sdsspy.files.filename(filetype, run, rerun=301, add_dir=True, **keywords)

    Keywords:
        run (can either be argument or keyword)
        rerun (default 301)
        camcol
        field
        id
        band
        type

    Get an SDSS file name, such as calibobj, etc, given the relevant info.
    run,rerun,camcol,band must be scalars.  
    
    Certain file names requre more info than others. asTrans only requires the
    run, but calibObj needs run,camcol (assuming default rerun). The directory
    is also added unless add_dir = False


    # example: calibobj for run 756, camcol 1 type 'gal' and default
    #     rerun of 301
    >>> filename('calibobj', 756, camcol=1, type='gal')
    '/global/data/boss/sweeps/2010-01-11/301/calibObj-000756-1-gal.fits.gz'
    """

    ftype=filetype_input.lower()

    if ftype not in filespec:
        raise ValueError("Unsupported file type: '%s'" % ftype)

    name = filespec[ftype]['name']
    name = expand_sdssvars(name, 
                           run=run, 
                           rerun=rerun, 
                           camcol=camcol,
                           field=field,
                           **keywords)

    if add_dir:
        dir = filedir(ftype, 
                      run=run, 
                      rerun=rerun, 
                      camcol=camcol,
                      field=field, 
                      **keywords)
        name = os.path.join(dir, name)
    
    return name



                 

__files2fields_pattern_compiled = re.compile('-[0-9][0-9][0-9][0-9]$')
__files2fields_errmess = \
      'Incorrectly formatted name %s.  Must end in -iiii.extension\n'
def files2fields(filenames):
    """
    field = sdsspy.files.files2fields(filenames)
    Extract the field number from field-by-field type of files such
    as a tsObj file that end with front-iiii.extension
    """
    if numpy.isscalar(filenames):
        filenames = [filenames]

    fields = []
    search = __files2fields_pattern_compiled.search

    for name in filenames:
        if name.find('.') != -1:
            sp = name.split('.')[0]
            mo = search(sp)
            if mo:
                field = int( sp[mo.start()+1:] )
                fields.append(field)
            else:
                sys.stdout.write(__files2fields_errmess % s)
                return []
        else:
            sys.stdout.write(__files2fields_errmess % s)
            return []

    return fields

class FileList():
    """
    Class:
        FileList
    Purpose:
        A class to find sdss file lists on disk.  Lists for previous calls
        are cached internally for speed.

    Construction and loading of lists:
        By default only the file type is needed, in which case a glob pattern
        with '*' for all file elements is created.  Subsets are determined
        through the keyword.

    Required Arguments for Construction or load():
        filetype:  An SDSS filetype.  Supported types:
            calibObj

    Keywords for Construction and load():
        run
        rerun (default 301)
        camcol
        field
        id
        band
        type (e.g. 'gal' for calibObj)
        
    Examples:
        import sdsspy
        ft=sdsspy.files.FileList('calibobj', run=756, type='gal')
        flist = ft.flist()
        for f in flist:
            print f

        ft.load('fpAtlas', run=94, camcol=3)
        for f in ft.flist():
            print f
    """
    def __init__(self, 
                 filetype, 
                 run=None, 
                 rerun=_DEFAULT_RERUN, 
                 camcol='*', 
                 field='*', 
                 id='*', 
                 band='*',
                 type='*',
                 **keywords):

        if not hasattr(FileList, '_flist'):
            # this class variable will hold all previously requested file 
            # lists
            FileList._flist={} 

        self.load(filetype, run, rerun, camcol, field, id, band, type)

    def load(self, 
             filetype, 
             run=None, 
             rerun=301, 
             camcol='*', 
             field='*', 
             id='*', 
             band='*', 
             type='*'):

        self._filetype=filetype.lower()
        self._run=run
        self._rerun=rerun
        self._camcol=camcol
        self._field=field
        self._id=id
        self._band=band
        self._type=type
        self._key = self.makekey(filetype,run,rerun,camcol,field,id,band,type)

        if self._key not in FileList._flist:
            
            # now load the file list
            pattern = filename(filetype,
                               run=run,rerun=rerun,camcol=camcol,field=field,id=id,
                               band=band,type=type)
            if filespec[self._filetype]['check_compress']:
                pattern += '*'

            flist = glob.glob(pattern)
            if len(flist) == 0:
                raise RuntimeError("no files matched pattern: %s" % pattern)
            FileList._flist[self._key] = flist

            self._used_cache = False
        else:
            self._used_cache = True

    def used_cache(self):
        return self._used_cache

    def reload(self):
        if self._key in FileList._flist:
            del FileList._flist[self._key]
        self.load(self._filetype,
                  self._run,
                  self._rerun,
                  self._camcol,
                  self._field,
                  self._id,
                  self._band,
                  self._type)

    def flist(self, subfields='*'):
        if self._key not in FileList._flist:
            raise ValueError("File has not been loaded")

        if subfields == '*' or not filespec[self._filetype]['byfield']:
            # just return everything if '*' or if the data is not split
            # by field anyway
            return FileList._flist[self._key]
        else:
            # if this file type is split by fields, extract the desired subset
            if isinstance(subfields,(str,unicode)):
                if subfields != '*':
                    raise ValueError("String not supported for field keyword "
                                     "except for '*'")
            else:
                if isinstance(subfields,(list,tuple,numpy.ndarray)):
                    f2get = subfields
                else:
                    f2get = [subfields]

                keep = []
                this_flist=FileList._flist[self._key]
                for f in f2get:
                    tname = filename(self._filetype,self._run,
                                     rerun=self._rerun,
                                     camcol=self._camcol,
                                     field=f,
                                     id=self._id,
                                     band=self._band,
                                     type=self._type)

                    if tname in this_flist:
                        keep.append(tname)

                if len(keep) == 0:
                    raise ValueError("None of the requested fields were found")

                return keep 



    def clear(self):
        FileList._flist={}

    def makekey(self,filetype,run,rerun,camcol,field,id,band,type):
        key='%s-%s-%s-%s-%s-%s-%s-%s' % (filetype,run,rerun,camcol,field,id,band,type)
        return key



# convert sdss numbers to strings in file names and such
def stripe2string(stripes):
    """
    ss = stripe2String(stripes)
    Return the string version of the stripe.  9->'09'
    Range checking is applied.
    """
    return tostring(stripes, 0, 99)

def run2string(runs):
    """
    rs = run2string(runs)
    Return the string version of the run.  756->'000756'
    Range checking is applied.
    """
    return tostring(runs,0,999999)


def rerun2string(reruns):
    """
    rrs = rerun2string(reruns)
    Return the string version of the rerun.  No zfill is used.
    """
    return tostring(reruns)

def camcol2string(camcols, band=None):
    """
    cs = camcol2string(camcols)
    Return the string version of the camcol.  1 -> '1'
    Range checking is applied.
    """
    cstr=tostring(camcols,1,6)
    if band is not None:
        # make sure it is a scalar
        if numpy.isscalar(band):
            buse=band
        else:
            buse=band[0]
        bstr=band2string(buse)
        if numpy.isscalar(cstr):
            cstr = bstr+cstr
        else:
            cstr=[bstr+c for c in cstr]
    return cstr
    
def field2string(fields):
    """
    fs = field2string(field)
    Return the string version of the field.  25->'0025'
    Range checking is applied.
    """
    return tostring(fields,0,9999)

def id2string(ids):
    """
    istr = id2string(ids)
    Return the string version of the id.  25->'00025'
    Range checking is applied.
    """
    return tostring(ids,0,99999)



band_dict = {0: 'u',
             1: 'g',
             2: 'r',
             3: 'i',
             4: 'z',
             'u':'u',
             'g':'g',
             'r':'r',
             'i':'i',
             'z':'z'}

def band2string(band):
    """
    bstr = band2string(bands)
    Return alphabetic representation of band
      bpstr = band2string(2)      # returns 'r'
      bpstr = band2string('2')    # returns 'r'
      bpstr = band2string('r')    # returns 'r'
      bpstr = band2string('R')    # returns 'r'
    """

    if not numpy.isscalar(band):
        # Scalar pars cannot be modified
        return [band2string(b) for b in band]


    if band not in band_dict:
        raise ValueError("bad bandpass indicator: %s" % bint)

    return band_dict[band]

def tostring(val, nmin=None, nmax=None):
    
    if not numpy.isscalar(val):
        return [tostring(v,nmin,nmax) for v in val]

    if isinstance(val, (str,unicode)):
        return val

    if nmin is not None:
        if val < nmin:
            raise ValueError("Number ranges below min value of %s\n" % nmin)
    if nmax is not None:
        if val > nmax:
            raise ValueError("Number ranges higher than max value of %s\n" % nmax)


    if nmax is not None:
        nlen = len(str(nmax))
        vstr = str(val).zfill(nlen)
    else:
        vstr = str(val)

    return vstr






def expand_sdssvars(string_in, **keywords):

    string = string_in

    # this will expand all environment variables, e.g. $PHOTO_SWEEP
    # if they don't exist, the result will be incomplete

    string = os.path.expandvars(string)

    if string.find('$RUN') != -1:
        if 'run' not in keywords:
            raise ValueError("run keyword must be sent")
        run=keywords['run']
        if run is None:
            raise ValueError("run keyword must be sent")
        
        string = string.replace('$RUN', run2string(run))

    if string.find('$RERUN') != -1:
        if 'rerun' not in keywords:
            raise ValueError("rerun keyword must be sent")
        rerun=keywords['rerun']
        if rerun is None:
            raise ValueError("rerun keyword must be sent")
        string = string.replace('$RERUN', rerun2string(rerun))


    if string.find('$CAMCOL') != -1:
        if 'camcol' not in keywords:
            raise ValueError("camcol keyword must be sent")
        camcol=keywords['camcol']
        if camcol is None:
            raise ValueError("camcol keyword must be sent")

        string = string.replace('$CAMCOL', camcol2string(camcol))

    if string.find('$FIELD') != -1:
        if 'field' not in keywords:
            raise ValueError("field keyword must be sent")
        field=keywords['field']
        if field is None:
            raise ValueError("field keyword must be sent")
        string = string.replace('$FIELD', field2string(field))

    if string.find('$ID') != -1:
        if 'id' not in keywords:
            raise ValueError("id keyword must be sent")
        id=keywords['id']
        if id is None:
            raise ValueError("id keyword must be sent")
        string = string.replace('$ID', id2string(id))

    if string.find('$BAND') != -1:
        if 'band' not in keywords:
            raise ValueError("band keyword must be sent")
        band=keywords['band']
        if band is None:
            raise ValueError("band keyword must be sent")
        string = string.replace('$BAND', band2string(band))

    if string.find('$TYPE') != -1:
        if 'type' not in keywords:
            raise ValueError("type keyword must be sent")
        type=keywords['type']
        if type is None:
            raise ValueError("type keyword must be sent")
        string = string.replace('$TYPE', type)



    return string

def ensure_sequence(var):
    if isinstance(var, (list,tuple,numpy.ndarray)):
        return var
    else:
        return [var]

