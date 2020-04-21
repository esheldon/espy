import numpy
import sys
import os
try:
    import recfile
    recfile_found=True
except:
    recfile_found=False
    sys.stdout.write('recfile could not be imported.  Only simple access '+\
                     'supported\n')

def write(array, outfile, delim=None, append=False, header={}):
    """
    sfile.write()

    Write a numpy array into a simple self-describing file format with an 
    ascii header.  See the read() function for information
    about reading this file format.

    Calling Sequence:
        sfile.write(array, outfile, delim=None, append=False, header=)

    Inputs:
        array: A numpy array.
        outfile: A string or file pointer for the output file.
    Optional Inputs:
        delim: Delimiter between fields.  Default is None for binary.
        append=False: Append to the file. Default is False. If set to True,
            then what happens is situation dependent:
                1) if the input is a file object then it is assumed there is 
                    an existing header.  The header is updated to reflect the
                    new appended rows after writing.
                2) if the input is a string and the file exists, then the file
                    is opened with mode "r+", it is assumed the header
                    exists, and the header is updated to reflext the new
                    appended rows after writing.
                3) if the input is a string and the file does *not* exist,
                    then the file is opened with mode "w" and the request
                    to append is ignored.
        header=: A dictionary containing keyword-value pairs to be added to
            the header.

    Examples:
        import sfile
        hdr={'date': '2007-05-12','age': 33} 
        sfile.write(arr, 'test.pya', header=hdr)

        sfile.write(arr2, 'test.pya', append=True)

    The file format:
      First line:
          NROWS = ----number---------
              where the number is formatted as %20d
      The header must also contain:
          DTYPE = array data type description in list of tuples form. For
              example DTYPE = [('field1', 'f8'), ('f2','2i4')]
      Case doesn't matter in keyword names. There can be other header 
      keyword-value pairs, e.g.
          x = 3
          y = 'hello world'
      Last two lines of header:
          END
          blank line
      Then the entire array is written as a single chunk of data.


    Modification History:
        Created: 2007-05-25, Erin Sheldon, NYU
        Ignore append=True when file does not yet exist.
        Allow continuation characters "\\" to continue header keywords
        onto the next line.  2009-05-05

    """

    # Get the file object.  If append is requested and it is found 
    # appropriate, doappend will be True, else False
    fobj, f_isstring, doappend = _GetFobj(outfile, 'write', append=append)

    # Write the header
    nrows = array.size
    descr = array.dtype.descr
    _WriteHeader(fobj, nrows, descr, delim=delim, append=doappend, 
                 header=header)

    if (delim is not None) and (recfile_found is not True):
        sys.stdout.write("delim='"+delim+"'\n")
        sys.stdout.write("recfile_found=%s\n" % recfile_found)
        raise ValueError('recfile module not found for ascii writing')
    elif (delim is not None):
        # let recfile deal with delimiter writing
        r = recfile.Open(fobj, mode='u', delim=delim, dtype=array.dtype)
        #r.Write(array)
        r.Write(array.view(numpy.ndarray))
    else:
        # Write data out as a binary chunk
        array.tofile(fobj)

    # Need to close if string was input
    if f_isstring:
        fobj.close()

def _remove_byteorder(descr):
    import copy
    new_descr = []
    for d in descr:
        # d is a tuple, make it a list
        newd = list( copy.copy(d) )

        tdef = newd[1]
        tdef = tdef[1:]
        newd[1] = tdef
        newd = tuple(newd)
        new_descr.append(newd)

    return new_descr

def _get_rows2read(rows):
    if rows is None:
        r=None
        l=0
    else:
        r = numpy.array(rows,ndmin=1)
        l = r.size
    return r,l

def _get_fields2read(fields):
    if fields is None:
        f=None
        l=0
    elif type(fields) == list:
        f=fields
        l=len(f)
    elif type(fields) == numpy.ndarray:
        f=fields
        l=len(f)
    elif type(fields) == str:
        f=[fields]
        l=len(f)
    else:
        raise ValueError('fields must be list,string or array')
    return f,l

def read(infile, rows=None, fields=None, header=False):
    """
    sfile.read()

    Read a numpy array from a simple self-describing file format with an 
    ascii header.  See the write() function for information
    about reading this file format.

    Calling Sequence:
        arr = sfile.read(infile, rows=None, fields=None, header=False)

    Inputs:
        infile: A string or file pointer for the input file.
    Optional Inputs:
        rows=: A scalar, array, list or tuple of row numbers to read.  
            Default is all.
        fields=: A scalar, list or tuple of strings naming fields to read 
            from the file.  Default is all fields.
        header=False:  If True, return both the array and the header dict in
            a tuple.

    Examples:
        import sfile
        arr=sfile.read('test.pya')

        arr2, hdr = sfile.read('test.pya', 
                               rows=[3,4], fields=['ra', 'dec'], header=True)

    The file format:
      First line:
          NROWS = ----number---------
              where the number is formatted as %20d
      The header must also contain:
          DTYPE = array data type description in list of tuples form. For
              example DTYPE = [('field1', 'f8'), ('f2','2i4')]
      Case doesn't matter in keyword names. There can be other header 
      keyword-value pairs, e.g.
          x = 3
          y = 'hello world'
      Last two lines of header:
          END
          blank line
      Then the entire array is written as a single chunk of data.

    Modification History:
        Created: 2007-05-25, Erin Sheldon, NYU
        Allow continuation characters "\\" to continue header keywords
        onto the next line.  2009-05-05
    """


    # Get the file object
    fobj,f_isstring,junk = _GetFobj(infile,'read')

    # Get the header
    hdr=read_header(fobj)

    delim = _MatchKey(hdr,'delim')

    # read the header
    # number of rows
    nrows = _GetNrows(hdr)

    # The dtype
    dt = _GetDtype(hdr)
    
    rows2read, nr =_get_rows2read(rows)
    fields2read, nf =_get_fields2read(fields)
    if (nf == 0) & (nr == 0) & (delim is None):
        # Its binary and all, just use fromfile
        result = numpy.fromfile(fobj,dtype=dt)
    else:
        if recfile_found:
            dtype=numpy.dtype(dt)
            if delim is None:
                send_delim = ""
            else:
                send_delim = delim
            robj = recfile.Open(fobj, nrows=nrows, mode='r', dtype=dtype,
                                delim=send_delim)
            result = robj.Read(rows=rows2read, fields=fields2read)
        else:
            if delim is not None:
                raise ValueError("Cannot extract rows/fields for text files without the recfile package\n")
            # Do it the brute force way
            result = numpy.fromfile(fobj,dtype=dt)
            result = result[rows2read]
            result = extract_fields(result, fields2read)


    # close file if was opened locally
    if f_isstring:
        fobj.close()

    if header:
        return result, hdr
    else:
        return result


def copy_fields(arr1, arr2):
    """
    Copy the fields that match
    """
    if arr1.size != arr2.size:
        raise ValueError('arr1 and arr2 must be the same size')

    names1=arr1.dtype.names
    names2=arr2.dtype.names
    for name in names1:
        if name in names2:
            arr2[name] = arr1[name]

def extract_fields(arr, keepnames):
    if type(keepnames) != list and type(keepnames) != numpy.ndarray:
        keepnames=[keepnames]
    arrnames = list( arr.dtype.names )
    new_descr = []
    for d in arr.dtype.descr:
        name=d[0]
        if name in keepnames:
            new_descr.append(d)
    if len(new_descr) == 0:
        raise ValueError('No field names matched')

    shape = arr.shape
    new_arr = numpy.zeros(shape,dtype=new_descr)
    copy_fields(arr, new_arr)
    return new_arr




def read_header(infile):
    """
    sfile..read_header()

    Read the header from a simple self-describing file format with an 
    ascii header.  See the write() function for information
    about reading this file format, and read() for reading.

    Calling Sequence:
        sfile.read_header(infile)

    Inputs:
        infile: A string or file pointer for the input file.

    Examples:
        import sfile
        hdr=sfile.read_header('test.pya')

    The file format:
      First line:
          NROWS = ----number---------
              where the number is formatted as %20d
      The header must also contain:
          DTYPE = array data type description in list of tuples form. For
              example DTYPE = [('field1', 'f8'), ('f2','2i4')]
      Case doesn't matter in keyword names. There can be other header 
      keyword-value pairs, e.g.
          x = 3
          y = 'hello world'
      For lines other than the first and last, continuation marks \\ 
      can be used.  e.g
          dtype = [('x','f4'), \\
                   ('y','i8')]
      Last two lines of header:
          END
          blank line
      Then the entire array is written as a single chunk of data
      or rows of ascii.

    Modification History:
        Created: 2007-05-25, Erin Sheldon, NYU
        Allow continuation characters "\\" to continue header keywords
        onto the next line.  2009-05-05
    """


    # Get the file object
    fobj,f_isstring,junk = _GetFobj(infile,'read')

    # Read through the header until we hit "END"
    lineold=None
    hdr={}
    line=fobj.readline().strip()
    while line.upper() != 'END':

        if lineold is not None:
            line=lineold+line
            lineold=None

        # see if this line is a continuation line
        if line[-1] == '\\':
            lineold=line[0:len(line)-1]
        else:
            ls = line.split('=')
            if len(ls) < 2:
                sys.stdout.write('Skipping malformed header line: %s\n' % line)
            else:
                key=ls[0].strip().lower()
                command="hdr['"+key+"'] = "+"=".join(ls[1:])
                exec(command)
        line=fobj.readline().strip()

    # read one more line, which should be blank
    line = fobj.readline()
    if line.strip() != '':
        raise RuntimeError('Header should end END and then a blank line')

    # close file if was opened locally
    if f_isstring:
        fobj.close()

    return hdr


def _WriteHeader(fobj, nrows, descr, delim=None, header={}, append=False):
    if append:
        # Just update the nrows and move to the end
        _UpdateNrows(fobj, nrows)
        fobj.seek(0,2) # Seek to end of file
    else:
        # Write a full new header

        if delim is not None:
            # Text file
            # In this case byte order in the dtype could cause errors
            descr = _remove_byteorder(descr)

        fobj.seek(0)

        _WriteNrows(fobj, nrows)
        if delim is not None:
            _WriteDelim(fobj, delim)
        _WriteDescr(fobj, descr)

        if isinstance(header, dict):
            for k,v in header.items():
                if not isinstance(k,str):
                    mess='All header keyword names must be strings'
                    raise RuntimeError(mess)
                if isinstance(v,str):
                    v = "'"+v+"'"
                # Force to be lower
                k=k.lower()
                if k not in ['dtype','nrows','delim']:
                    fobj.write(k+' = '+str(v)+'\n')

        fobj.write('END\n')
        fobj.write('\n')


def _WriteNrows(fobj, nrows):
    # Go to the beginning of the file
    fobj.seek(0)
    # Specially formatted for updating later
    out = 'NROWS = %20d\n' % (nrows,)
    fobj.write(out)

def _WriteDelim(fobj, delim):
    if delim == '\t':
        fobj.write("DELIM = '\\t'\n")
    else:
        fobj.write("DELIM = '%s'\n" % delim)

def _WriteDescr(fobj, descr):
    descr_strings = ['         '+str(d) for d in descr]
    descr_head = ', \\\n'.join(descr_strings)
    descr_head = '['+descr_head.strip()+']'
    fobj.write('DTYPE = '+descr_head+'\n')

def _WriteDescrOld(fobj, descr):
    fobj.write('DTYPE = '+str(descr)+'\n')

def _UpdateNrows(fobj, nrows_add):
    fobj.seek(0)
    hdr=read_header(fobj)
    nrows_current = _GetNrows(hdr)
    if nrows_current is None:
        raise RuntimeError('Attempting to update nrows but not found in header')
    nrows_new = nrows_current + nrows_add
    _WriteNrows(fobj, nrows_new)

def _MatchKey(d,key):
    """
    Match the key in a case-insensitive way and return the value.
    """
    if not isinstance(d,dict):
        raise RuntimeError('Input object must be a dict')
    keys = list( d.keys() )
    keyslow = [k.lower() for k in keys]
    keylow = key.lower()
    if keyslow.count(keylow) != 0:
        ind = keyslow.index(keylow)
        return d[keys[ind]]
    else:
        return None
        

def _GetNrows(header):
    if not isinstance(header,dict):
        raise RuntimeError('header must be a dict')

    nrows = _MatchKey(header,'nrows')
    if nrows is not None:
        return int(nrows)
    else:
        raise RuntimeError('Header must contain the nrows keyword')

def _GetDtype(header):
    if not isinstance(header,dict):
        raise RuntimeError('header must be a dict')

    dtype = _MatchKey(header,'dtype')
    if dtype is not None:
        return dtype
    else:
        raise RuntimeError('Header must contain the DTYPE keyword')


def _GetFobj(outfile, major_mode, append=False):
    doappend=False

    if isinstance(outfile,str):
        if major_mode == 'read':
            mode='r'
        else:
            if append:
                # make sure it actually exists first before trying to 
                # append
                if os.path.exists(outfile):
                    doappend=True
                    mode="r+"
                else:
                    mode="w"
            else:
                mode="w"
        f_isstring = True
        fname_full = os.path.expanduser(outfile)
        fobj = open(fname_full,mode)
    elif isinstance(outfile,file):
        f_isstring = False
        fobj=outfile
        if append:
            doappend=True
    else:
        raise RuntimeError('File must be a string or a file object')

    return fobj, f_isstring, doappend


def test():
    """
    very simple test
    """
    import tempfile
    tmpfile = tempfile.TemporaryFile(suffix='.sf')
    x=numpy.array([(1.0, 3),(4.5,2)], dtype=[('fcol','f4'),('icol','i4')])

    write(x, tmpfile)
    tmpfile.seek(0)
    read(tmpfile)

    tmpfile.close()
