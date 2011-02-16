"""
NAME:
  es_util

PURPOSE:
  A set of utility modules for array manipulation.

MODULES:
    arrscl
    make_xy_grid
    wmom
    sigma_clip
    interplin
    gauleg
    qgauss
    combine_arrlist

    histogram
    objhist

    copy_fields
    extract_fields
    remove_fields
    copy_field_by_name
    add_fields
    compare_arrays
    expand_filename
    fits2array
    array2fits (now obsolete since I've added this
        functionality to my version of pyfits)
    _extract_format
    pad_string
    unique
    unique_seq
    match
  
REVISION HISTORY:
  Began converting routines from IDL: 2006-10-23. Erin Sheldon, NYU

"""

license="""
  Copyright (C) 2009  Erin Sheldon

    This program is free software; you can redistribute it and/or modify it
    under the terms of version 2 of the GNU General Public License as
    published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""

import sys
from sys import stdout,stderr
import os

try:
    import numpy as numpy
except:
    stderr.write('Failed to import numpy\n')
try:
    import scipy
    from scipy.interpolate import interp1d
    from scipy import weave
except:
    stderr.write('Failed to import scipy\n')


def arrscl(arr, minval, maxval, arrmin=None, arrmax=None):
    """
    NAME:
      arrscl()

    PURPOSE:
      Rescale the range of an array.

    CALLING SEQUENCE:
      newarr = arrscl(arr, minval, maxval, arrmin=None, arrmax=None)
    
    INPUTS:
      arr: An array
      minval: The minimum value for the output array
      maxval: The maximum value for the output array

    OPTIONAL OUTPUTS:
      arrmin=None: An number to use for the min range of the input array. By
        default it is taken from the input array.
      arrmax=None: An number to use for the max range of the input array. By
        default it is taken from the input array.
      * arrmin,arrmax are useful if you know the array is a sample of a
        particular range, for example of they are random numbers drawn
        from [0,1] you would send arrmin=0., arrmax=1.

    OUTPUTS:
      The new array.

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU
      
    """

    # makes a copy either way (asarray would not if it was an array already)
    output = numpy.array(arr)
    
    if arrmin == None: arrmin = output.min()
    if arrmax == None: arrmax = output.max()
    
    if output.size == 1:
        return output
    
    if (arrmin == arrmax):
        raise ValueError('arrmin must not equal arrmax')

    #try:
    a = (maxval - minval)/(arrmax - arrmin)
    b = (arrmax*minval - arrmin*maxval)/(arrmax - arrmin)
    #except:
    #print "Error calculating a,b: ", \
    #      sys.exc_info()[0], sys.exc_info()[1]
    #return None

    # in place
    numpy.multiply(output, a, output)
    numpy.add(output, b, output)
    
    return output

def make_xy_grid(n, xrang, yrang):
    # Create a grid on input ranges
    rng = numpy.arange(n, dtype='f8')
    ones = numpy.ones(n, dtype='f8')

    x = arrscl(rng, xrang[0], xrang[1])
    y = arrscl(rng, yrang[0], yrang[1])

    x= numpy.outer(x, ones)
    y= numpy.outer(ones, y)
    x = x.flatten(1)
    y = y.flatten(1)

    return x,y



def wmom(arrin, weightsin, inputmean=None, calcerr=False, sdev=False):
    """
    NAME:
      wmom()
      
    PURPOSE:
      Calculate the weighted mean, error, and optionally standard deviation
      of an input array.

    CALLING SEQUENCE:
     wmean,werr = wmom(arr, weights, inputmean=None, calcerr=False, sdev=False)
    
    INPUTS:
      arr: A numpy array or a sequence that can be converted.
      weights: A set of weights for each elements in array.
    OPTIONAL INPUTS:
      inputmean: An input mean value, around which them mean is calculated.
      calcerr=False: Calculate the weighted error.  By default the error
        is calculated as 1/sqrt( weights.sum() ).  In this case it is
        calculated as sqrt( (w**2 * (arr-mean)**2).sum() )/weights.sum()
      sdev=False: If True, also return the weighted standard deviation 
        as a third element in the tuple.

    OUTPUTS:
      wmean, werr: A tuple of the weighted mean and error. If sdev=True the
         tuple will also contain sdev: wmean,werr,wsdev

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU

   """
    from numpy import float64
    
    # no copy made if they are already arrays
    arr = numpy.array(arrin, ndmin=1, copy=False)
    weights = numpy.array(weightsin, ndmin=1, copy=False)
    
    # Weights is forced to be type double. All resulting calculations
    # will also be double
    if weights.dtype != float64:
        weights = numpy.array(weights, dtype=float64)
  
    wtot = weights.sum()
        
    # user has input a mean value
    if inputmean is None:
        wmean = ( weights*arr ).sum()/wtot
    else:
        wmean=float(inputmean)

    # how should error be calculated?
    if calcerr:
        werr2 = ( weights**2 * (arr-wmean)**2 ).sum()
        werr = numpy.sqrt( werr2 )/wtot
    else:
        werr = 1.0/numpy.sqrt(wtot)

    # should output include the weighted standard deviation?
    if sdev:
        wvar = ( weights*(arr-wmean)**2 ).sum()/wtot
        wsdev = numpy.sqrt(wvar)
        return wmean,werr,wsdev
    else:
        return wmean,werr




def sigma_clip(arrin, niter=4, nsig=4, extra={}, verbose=False):
    """
    NAME:
      sigma_clip()
      
    PURPOSE:
      Calculate the mean/stdev of an array with sigma clipping. Iterate
      niter times, removing elements that are outside nsig, and recalculating
      mean/stdev.

    CALLING SEQUENCE:
      mean,stdev = sigma_clip(arr, niter=4, nsig=4, extra={})
    
    INPUTS:
      arr: A numpy array or a sequence that can be converted.

    OPTIONAL INPUTS:
      niter: number of iterations, defaults to 4
      nsig: number of sigma, defaults to 4

    OUTPUTS:
      mean,stdev: A tuple containing mean and standard deviation.
    OPTIONAL OUTPUTS
      extra={}: Dictionary containing the array of used indices in
         extra['index']

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU

   """
    arr = numpy.array(arrin, ndmin=1, copy=False)

    index = numpy.arange( arr.size )

    for i in numpy.arange(niter):
        m = arr[index].mean()
        s = arr[index].std()

        if verbose:
            sys.stdout.write('iter %s\tnuse: %s\tmean %s\tstdev %s\n' % \
                (i+1, index.size,m,s))

        clip = nsig*s

        w, = numpy.where( (numpy.abs(arr[index]) - m) < clip )

        if w.size == 0:
            sys.stderr.write("nsig too small. Everything clipped on iteration %d" % i+1)
            return m,s


        index = index[w]

    # Calculate final stats
    amean = arr[index].mean()
    asig = arr[index].std()

    extra['index'] = index
    return amean, asig
        


def interplin(vin, xin, uin):
    """
    NAME:
      interplin()
      
    PURPOSE:
      Perform 1-d linear interpolation.  Values outside the bounds are
      permitted unlike the scipy.interpolate.interp1d module. The are
      extrapolated from the line between the 0,1 or n-2,n-1 entries.
      This program is not as powerful as interp1d but it does provide
      this which makes it compatible with the IDL interpol() function.

    CALLING SEQUENCE:
      yint = interplin(y, x, u)

    INPUTS:
      y, x:  The y and x values of the data.
      u: The x-values to which will be interpolated.

    REVISION HISTORY:
      Created: 2006-10-24, Erin Sheldon, NYU
    """
    # Make sure inputs are arrays.  Copy only made if they are not.
    v=numpy.array(vin, ndmin=1, copy=False)
    x=numpy.array(xin, ndmin=1, copy=False)
    u=numpy.array(uin, ndmin=1, copy=False)

    # Find closest indices
    xm = x.searchsorted(u) - 1
    
    # searchsorted returns size(array) when the input is larger than xmax
    # Also, we need the index to be less than the last since we interpolate
    # *between* points.
    w, = numpy.where(xm >= (x.size-1))
    if w.size > 0:
        xm[w] = x.size-2

    w, = numpy.where(xm < 0)
    if w.size > 0:
        xm[w] = 0
        
    xmp1 = xm+1
    return (u-x[xm])*(v[xmp1] - v[xm])/(x[xmp1] - x[xm]) + v[xm]



def gauleg(x1, x2, npts):
    """
    NAME:
      gauleg()
      
    PURPOSE:
      Calculate the weights and abscissa for Gauss-Legendre integration.
    
    CALLING SEQUENCE:
      x,w = gauleg(x1,x2,npts)

    INPUTS:
      x1,x2: The range for the integration.
      npts: Number of points to use in the integration.

    REVISION HISTORY:
      Created: 2006-10-24. Adapted from Numerial recipes in C. Uses
        scipy.weave.inline for the C loops.  2006-10-24 Erin Sheldon NYU
    """

    # outputs
    x = numpy.zeros(npts, dtype=numpy.float64)
    w = numpy.zeros(npts, dtype=numpy.float64)

    # Code string for weave
    code = \
         """
         int i, j, m;
         double xm, xl, z1, z, p1, p2, p3, pp, pi, EPS, abszdiff;
         
         EPS = 3.e-11;
         pi=3.1415927;

         m = (npts + 1)/2;

         xm = (x1 + x2)/2.0;
         xl = (x2 - x1)/2.0;
         z1 = 0.0;

         for (i=1; i<= m; ++i) 
         {
      
           z=cos( pi*(i-0.25)/(npts+.5) );

           abszdiff = fabs(z-z1);

           while (abszdiff > EPS) 
           {
             p1 = 1.0;
             p2 = 0.0;
             for (j=1; j <= npts;++j)
             {
                p3 = p2;
                p2 = p1;
                p1 = ( (2.0*j - 1.0)*z*p2 - (j-1.0)*p3 )/j;
             }
             pp = npts*(z*p1 - p2)/(z*z -1.);
             z1=z;
             z=z1 - p1/pp;

             abszdiff = fabs(z-z1);

           }
      
           x(i-1) = xm - xl*z;
           x(npts+1-i-1) = xm + xl*z;
           w(i-1) = 2.0*xl/( (1.-z*z)*pp*pp );
           w(npts+1-i-1) = w(i-1);


         }

         return_val = 0;
         

         """
    
    weave.inline(code, ['x1', 'x2', 'npts', 'x', 'w'],
                       type_converters = weave.converters.blitz)

    return x,w


    
# Global variables.  We can save time if we don't need to generate
# absissa and weights.
XXi=numpy.array([])
WWi=numpy.array([])

def _run_gauleg(npts):
    if XXi.size != npts:
        globals()['XXi'], globals()['WWi'] = gauleg(-1.0,1.0,npts)
    
def qgauss(yin, xin, npts):
    """
    NAME:
      qgauss()
      
    PURPOSE:
      Calculate the integral of data y,x using the Gauss-Legendre method.
    
    CALLING SEQUENCE:
      int = qgauss(y, x, npts)

    INPUTS:
      y, x: The y and x values of the data.
      npts: Number of points to use in the integration.

    REVISION HISTORY:
      Created: 2006-10-24. Converted from my IDL code Erin Sheldon, NYU
    """
    y = numpy.array(yin, ndmin=1, copy=False)
    x = numpy.array(xin, ndmin=1, copy=False)

    _run_gauleg(npts)

    # Note, x must be ordered
    x1 = x[0]
    x2 = x[-1]

    # these needed for coordinate transformation
    f1 = (x2-x1)/2.
    f2 = (x2+x1)/2.

    # interpolate the data to XXi
    xvals = XXi*f1 + f2
    yvals = interplin(y, x, xvals)

    return f1 * ((yvals*WWi).sum())




def combine_arrlist(arrlist, keep=False):
    """
    Combined the list of arrays into one big array.  The arrays must all
    be the same data type.

    By default the elements are deleted as they are added to the big array.
    Turn this off with keep=True
    """
    import numpy
    if not isinstance(arrlist,list):
        raise RuntimeError('Input must be a list of arrays')

    isarray = isinstance(arrlist[0], numpy.ndarray)
    isrec = isinstance(arrlist[0], numpy.recarray)
        
    if not isarray:
        mess = 'Input must be a list of arrays or recarrays. Found %s' % \
                type(arrlist[0])
        raise RuntimeError(mess)

    # loop and get total number of entries
    counts=0
    for data in arrlist:
        counts = counts+data.size

    output = numpy.zeros(counts, dtype=arrlist[0].dtype)
    if isrec:
        output = output.view(numpy.recarray)

    beg=0
    if keep:
        for data in arrlist:
            num = data.size
            output[beg:beg+num] = data
            beg=beg+num
    else:
        while len(arrlist) > 0:
            data = arrlist.pop(0)
            num = data.size
            output[beg:beg+num] = data
            beg=beg+num

    return output


def copy_fields(arr1, arr2):
    """
    Copy common fields from one numpy array or recarray to another.
    """
    if arr1.size != arr2.size:
        raise ValueError('arr1 and arr2 must be the same size')

    names1=arr1.dtype.names
    names2=arr2.dtype.names
    for name in names1:
        if name in names2:
            arr2[name] = arr1[name]

def extract_fields(arr, keepnames):
    """
    Extract a set of fields from a numpy array or recarray.
    """
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



def remove_fields(arr, rmnames):
    """
    Remove a set of fields from a numpy array or recarray
    """
    if type(rmnames) != list:
        rmnames=[rmnames]
    descr = arr.dtype.descr
    new_descr = []
    for d in descr:
        name=d[0]
        if name not in rmnames:
            new_descr.append(d)

    if len(new_descr) == 0:
        raise ValueError('Error: All fields would be removed')

    shape = arr.shape
    new_arr = numpy.zeros(shape, dtype=new_descr)
    copy_fields(arr, new_arr)
    return new_arr

def copy_fields_by_name(arr, names, vals):
    """
    Copy values into an array with fields, or recarray, by name.
    """
    if type(names) != list and type(names) != numpy.ndarray:
        names=[names]
    if type(vals) != list and type(vals) != numpy.ndarray:
        vals=[vals]
    if len(names) != len(vals):
        raise ValueError('Length of names and values must be the same')

    arrnames = list(arr.dtype.names)
    for name,val in zip(names,vals):
        if name in arrnames:
            arr[name] = val

def add_fields(arr, add_dtype_or_descr, defaults=[]):
    """
    Add new fields to a numpy array or recarray, with an optional
    set of default values
    """
    # the descr is a list of tuples
    import copy
    old_descr = arr.dtype.descr
    add_dtype = numpy.dtype(add_dtype_or_descr)
    add_descr = add_dtype.descr

    new_descr = copy.deepcopy(old_descr)

    old_names = list(arr.dtype.names)
    new_names = list(add_dtype.names)
    for d in add_descr:
        name=d[0]
        if old_names.count(name) ==0:
            new_descr.append(d)
        else:
            raise ValueError( 'field '+str(name)+' already exists')

    shape = arr.shape
    new_arr = numpy.zeros(shape, dtype=new_descr)
    
    copy_fields(arr, new_arr)
    
    # See if the user has indicated default values for the new fields
    if type(defaults) != list:
        defaults=[defaults]
    ldef=len(defaults)
    if ldef > 0:
        if ldef != len(add_descr):
            raise ValueError('defaults must be same length as new dtype')
        copy_fields_by_name(new_arr, list(add_dtype.names), defaults)

    return new_arr


def compare_arrays(arr1, arr2, verbose=False):
    """
    Compare the values field-by-field in two sets of numpy arrays or
    recarrays.
    """

    nfail = 0
    for n2 in arr2.dtype.names:
        n1 = n2
        if n1 not in arr1.dtype.names:
            n1 = n1.lower()
            if n1 not in arr1.dtype.names:
                n1 = n1.upper()
                if n1 not in arr1.dtype.names:
                    raise ValueError('field name %s not found in array 1' % n2)
            
        if verbose:
            sys.stdout.write("    testing field: '%s'\n" % n2)
            sys.stdout.write('        shape...........')
        if arr2[n2].shape != arr1[n1].shape:
            nfail += 1
            if verbose:
                sys.stdout.write('shapes differ\n')
        else:
            if verbose:
                sys.stdout.write('OK\n')
                sys.stdout.write('        elements........')
            w,=numpy.where(arr1[n1].ravel() != arr2[n2].ravel())
            if w.size > 0:
                nfail += 1
                if verbose:
                    sys.stdout.write('\n        '+\
                        '%s elements in field %s differ\n' % (w.size,n2))
            else:
                if verbose:
                    sys.stdout.write('OK\n')

    if nfail == 0:
        if verbose:
            sys.stdout.write('All tests passed\n')
        return True
    else:
        if verbose:
            sys.stdout.write('%d differences found\n' % nfail)
        return False







def expand_filename(filename):
    """
    expand all user info such as ~userid and environment
    variables such as $SOMEVAR.
    """
    fname = os.path.expanduser(filename)
    fname = os.path.expandvars(fname)
    return fname

    
def fits2array(fobj, ext=0, fields=None, view=numpy.ndarray, 
               verbose=False, combine=False):
    """
    A wrapper for the pyfits.getdata() method.  This extracts the header
    as well and if the input is a list of files then they are all read

    By default the first extension (ext=0) is read.  

    If the result is table data the fields keyword may be sent to 
    extract a subset of the fields.

    If input is a list of files and combine=True then an attempt is made
    to combine them into one big array.  Only makes sense if they are the
    same data type
    """
    import pyfits

    if type(fobj) == type([]):
        # the input is a list of files
        res, hdrlist = _fits2array_list(fobj, ext=ext, fields=fields, 
                                        view=view,
                                        verbose=verbose)
        if combine:
            res = combine_arrlist(res)
        return res, hdrlist

    if verbose and (type(fobj) == type('')):
        stdout.write('Reading: %s\n' % fobj)

    hdr = pyfits.getheader(fobj, ext=ext)
    tdata,hdr = pyfits.getdata(fobj, ext=ext, header=True)
    if view is None:
        data=tdata
    else:
        data = tdata.view(view)
    if fields is not None:
        data = extract_fields(data, fields)
    return data,hdr

def _fits2array_list(flist, ext=0, fields=None, view=numpy.ndarray):
    datalist = []
    hdrlist=[]
    for f in flist:
        data,hdr = fits2array(f, ext=ext, fields=fields, view=view)
        datalist.append(data)
        hdrlist.append(hdr)
    return datalist, hdrlist



def bigendian(array):
    machine_little_endian = numpy.little_endian

    byteorder = array.dtype.base.byteorder
    return (byteorder == '>') \
            or (machine_little_endian and byteorder== '=')




def to_bigendian(array, inplace=False, keep_dtype=False):
    """
    Convert an array to big endian byte order, updating the dtype
    to reflect this.  Send keep_dtype=True to prevent the dtype
    from being updated.
    """

    doswap=False
    if array.dtype.names is None:
        if not bigendian(array):
            doswap=True
    else:
        # assume all are same byte order: we only need to find one with
        # little endian
        for fname in array.dtype.names:
            if not bigendian(array[fname]):
                doswap=True
                break

    outdata = array
    if doswap:
        outdata = byteswap(outdata, inplace, keep_dtype=keep_dtype)

    return outdata

def byteswap(array, inplace=False, keep_dtype=False):
    """
    byteswap an array, updating the dtype to reflect this

    If you *don't* want the dtype changed, simply use the
    built-int array method.  E.g.  array.byteswap()
    """

    outdata = array.byteswap(inplace)
    if not keep_dtype:
        outdata.dtype = outdata.dtype.newbyteorder()

    return outdata
      

def _extract_format(dtype):
    """
    the pyfits format extractor chokes on multidimensional arrays.

    This replacement returns the format string for any type and this can be
    sent to the pyfits format extractor.  If multi-d it converts to the
    equivalent 1-d but also returns a dimensions string that can be placed in
    a TDIM keyword field.
    
    the input should be the dtype of an individual field, e.g. 
        arr.dtype[3]
    """
    shape = dtype.shape
    tname = dtype.base.kind + str(dtype.base.itemsize)
    byteorder = dtype.base.byteorder

    ndims = len(shape)
    if ndims > 0:
        nel = numpy.array(shape, dtype='i8').prod()
        if nel > 1:
            tname = str(nel)+tname
    else:
        nel=1

    # will only use this for ndims > 1
    dimstring = [str(d) for d in shape] 
    dimstring = '('+','.join(dimstring)+')'

    return tname, byteorder, dimstring



def array2fits(data, filename, header=None, clobber=True, makecopy=False, 
               noswap=False):
    """
    This is now obsolete since I've addes this functionality to pyfits

    Dump an array to a fits file.  If the data has fields it is put into
    a binary table in the first extension.

    Image data:
        In this case the data is simply written to the file using the
        pyfits.writeto() function.

    Table data:

        Input data must be a numpy ndarray, recarray, pyfits.rec.recarray, or
        pyfits.FITS_rec 

        Currently requires a patched pyfits because of some errors involving
        using buffers

        In order to save memory, by default no copies are made.  Before
        writing the data are byteswapped if the user's machine is little
        endian. By default, the data are then byteswapped back again. Because
        this can be slow the noswap=True keyword can be sent.  Then the data
        returns with big-endian ordering and the dtype reflects this.  Note
        for noswap=True to work properly all fields in the array must have the
        same byte ordering; this is because of a limitation of numpy arrays.

        WARNING:  By default, if the code crashes the byte ordering may not be
        correct on return from the program.

        If you are worried the program might crash and the data will end up in
        an odd state, send the makecopy=True keyword.  Be warned that this
        will use twice the memory.

        On finish the hdu that was created is returned.

    Todo:
        Allow appending a new hdu

        Want to allow header inputs of other than type Header.  This would
        involve copying cardlists or dicts to a Header type.

    """

    import numpy
    import pyfits

    if data.dtype.fields is None:
        # This is simple array data
        pyfits.writeto(filename, data, header=header, clobber=clobber)

        return


    # the data has fields

    if makecopy:
        tmp = data.copy()
        fitsrec = tmp.view(pyfits.FITS_rec)
    else:
        # no copy is made here
        fitsrec = data.view(pyfits.FITS_rec)

    byteorder_list=[]
    if fitsrec._coldefs == None:
        #
        # The data does not have a _coldefs attribute so 
        # create one from the underlying recarray.
        #
        columns = []
        formats = []

        for i in range(len(data.dtype.names)):
            cname = data.dtype.names[i]

            if data.dtype.fields[cname][0].type == numpy.string_:
                # Still need to deal with strings
                format = \
                        'A'+str(data.dtype.fields[cname][0].itemsize)
                # STILL NEED TO put this into extract_format
            else:
                fmt, byteorder, dimstring = _extract_format(data.dtype[i])
                format = \
                        pyfits.NP_pyfits._convert_format(fmt, reverse=True)
                #stdout.write("npy: '%s' fits: '%s'\n" % (fmt, format))

            #stdout.write('%s format = %s\n' % (cname, format))
            formats.append(format)
            
            # This was making a copy of the data!
            #c = pyfits.Column(name=cname,format=format,array=data[cname])
            c = pyfits.Column(name=cname,format=format)
            columns.append(c)

            byteorder_list.append(byteorder)

        try:
            tbtype = 'BinTableHDU'

            if fitsrec._xtn == 'TABLE':
                tbtype = 'TableHDU'
        except AttributeError:
            pass

        fitsrec.formats = formats
        fitsrec._coldefs = pyfits.ColDefs(columns,tbtype=tbtype)


    hdu = pyfits.BinTableHDU(fitsrec, header=header)
    hdu.writeto(filename, clobber=clobber)

    if not makecopy:
        is_little_endian = numpy.little_endian
        if noswap:
            # in this case just determine if the thing was byteswapped
            # and just view it differently on return.  We have to assume
            # the whole array has the same byte ordering since we can only
            # to the switch in place for the whole array
            change_view=False
            for i in range(len(data.dtype)):
                byetorder=byteorder_list[i]
                name=data.dtype.names[i]
                if ((byteorder == '<') 
                    or (byteorder == '=' and is_little_endian) ):
                    change_view=True
            if change_view:
                data.dtype = data.dtype.newbyteorder()
        else:
            inplace=True
            for i in range(len(data.dtype)):
                # if the old byteorder was '<' or if it was '=' and this
                # arch is little endian, then we need to byteswap back
                if ((byteorder_list[i] == '<') 
                    or (byteorder_list[i] == '=' and is_little_endian) ):
                    name=data.dtype.names[i]
                    #stdout.write('Swapping %s\n' % name)
                    #if noswap:
                        #    data[name] = data[name].newbyteorder()
                        #else:
                    data[name].byteswap(inplace)

    return hdu

def dict_select(input_dict, keep=[], remove=[]):
    if len(keep) == 0 and len(remove) == 0:
        stderr.write('Enter one or both of keep and remove keywords\n')

    outdict={}

    if len(keep) == 0:
        keep = input_dict.keys()

    for key in keep:
        if key in input_dict and key not in remove:
            outdict[key] = input_dict[key]

    return outdict

def dict2fitsheader(d, prompt_comments=False):
    import pyfits
    h=pyfits.Header()

    for key in d:
        if prompt_comments:
            comment=raw_input('Enter a comment for field %s\n' % key)
        else:
            comment=''
        h.update(key, d[key], comment)
    return h

def dict2array(d, sort=False, keys=None):
    """

    Convert a dictionary to an array with fields (recarray, structured array).
    This works for simple types e.g.  strings, integers, floating points.

    You can send the keys= keyword to provide a sequence of keys to copy.
    This can be used to order the fields (standard dictionary keys are
    unordered) or copy only a subset of keys.

    You can also send sort=True to sort the keys.

    In python >= 3.1 you can used ordered dictionaries.

    """
    desc=[]

    if keys is None:
        if sort:
            keys=sorted(d)
        else:
            keys=d.keys()

    for key in keys:
        # check key existence in case a set of keys was sent
        if key not in d:
            raise KeyError("Requested key %s not in dictionary" % key)

        if isinstance(d[key], int):
            dt=int
        elif isinstance(d[key], float):
            dt=float
        elif isinstance(d[key], str):
            dt='S%s' % len(d[key])
        else:
            raise ValueError("Only support int, float, string currently")

        desc.append( (key, dt) )

    a=numpy.zeros(1, dtype=desc)

    for key in keys:
        a[key] = d[key]

    return a



fhc_ignores = ['comment','end']
def FieldHeaderCard2Dict(fhc):
    import pyfits
    card = pyfits.Card()

    d={}
    clen=80
    ncards = len(fhc)/clen
    for i in range(ncards):
        card.fromstring(fhc[i*80:(i+1)*80])
        d[card.key.lower()] = card.value
    return d

def FieldHeaderCard2Array(fhc):
    import pyfits
    card = pyfits.Card()

    descr = []
    cdict = {}

    clen=80
    ncards = len(fhc)/clen
    for i in range(ncards):
        card.fromstring(fhc[i*80:(i+1)*80])
        key = card.key.lower()
        value = card.value
        vtype = type(value)
        asarray = numpy.array(value)

        if key not in fhc_ignores:
            if key not in cdict:
                if asarray.dtype == type(True):
                    asarray = numpy.array('False')
                    value = str(value)

                d = (key,asarray.dtype)
                descr.append(d)
                cdict[key] = value

    arr = numpy.zeros(1, dtype=descr)
    for key in arr.dtype.names:
        arr[key] = cdict[key]

    return arr




def pad_string(instring, length, val=' ', left=False):
    """
    Pad the input string to the specified length, by default with spaces
    on the right.  You can specify the value and left padding with the
    keywords
    """
    # the sense is reversed here.  left padding is the same as right 
    # justifying
    if not left:
        return instring.ljust(length,val)
    else:
        return instring.rjust(length,val)


def unique_seq(seq, idfun=None, values=False): 
    """
    Get an order-preserving list of indices of unique values in the 
    sequence.  Send values=True to get the values instead of the 
    indicies.  Send idfun= keyword in order to specify a method for
    identifying, or marking, the object, default is the value itself.
    """
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []

    i=0
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        if values:
            result.append(item)
        else:
            result.append(i)

        i+= 1
    return result

def unique(arr, values=False):
    """
    un = unique(arr, values=False)

    return indices of unique elements of a numpy array, or the 
    unique values if values=True.  This is not order preserving.

    """
    n = arr.size
    keep = numpy.zeros(n, dtype='i8')

    s = arr.argsort()

    val = arr[0]
    i=1
    nkeep = 0
    while i < n:
        ind = s[i]
        if arr[ind] != val:
            val = arr[ind]
            nkeep += 1
            keep[nkeep] = ind
        i += 1

    keep = keep[0:nkeep+1]
    if values:
        return arr[keep]
    else:
        return keep


def match(arr1, arr2):
    """
    ind1,ind2 = match(arr1, arr2)
    match two numpy arrays.  Return the indices of the matches or [-1] if no
    matches are found.  This means
        arr1[ind1] == arr2[ind2]
    is true for all corresponding pairs

    Arrays must contain only unique elements
    """
    dtype = 'i8'
    n1 = len(arr1)
    n2 = len(arr2)

    if (n1 == 1) or (n2 == 1):
        # one of the arrays is length one
        if n2 > 1:
            sub2, = numpy.where(arr2 == arr1[0])
            if sub2.size > 0:
                sub1 = numpy.array([0], dtype=dtype)
            else:
                sub1 = numpy.array([-1], dtype=dtype)
        else:
            sub1, = numpy.where(arr1 == arr2[0])
            if sub1.size > 0:
                sub2 = numpy.array([0], dtype=dtype)
            else:
                sub2 = numpy.array([-1], dtype=dtype)

        return sub1, sub2


    # make a combined set
    tmp = numpy.zeros(n1+n2, dtype=arr1.dtype)
    tmp[0:n1] = arr1[:]
    tmp[n1:] = arr2[:]

    ind = numpy.zeros(n1+n2, dtype=dtype)
    ind[0:n1] = numpy.arange(n1)
    ind[n1:] = numpy.arange(n2)

    vec = numpy.zeros(n1+n2, dtype='b1')
    vec[n1:] = 1

    # sort combined list
    sortind = tmp.argsort()
    tmp = tmp[sortind]
    ind = ind[sortind]
    vec = vec[sortind]

    # this finds adjacent dups but only if they are not from the
    # same array.  Since we demand unique arrays I'm not sure why
    # the second check is needed
    firstdup, = numpy.where((tmp == numpy.roll(tmp,-1)) &
                            (vec != numpy.roll(vec,-1)) )
    if firstdup.size == 0:
        sub1 = numpy.array([-1], dtype=dtype)
        sub2 = numpy.array([-1], dtype=dtype)
        return sub1, sub2

    # both duplicate values...?
    dup = numpy.zeros(firstdup.size*2, dtype=dtype)

    even = numpy.arange( firstdup.size, dtype=dtype)*2
    dup[even] = firstdup
    dup[even+1] = firstdup+1

    # indices of duplicates
    ind = ind[dup]
    # vector id of duplicates
    vec = vec[dup]

    # now subscripts
    sub1 = ind[ numpy.where( vec == 0 ) ]
    sub2 = ind[ numpy.where( vec != 0 ) ]
    return sub1, sub2

def _match_test_write(arr1, arr2):
    stdout.write('Error matching arrays: %s : %s\n' % 
                 (str(arr1),str(arr2)))
def _match_test():
    nerror = 0

    arr1 = numpy.array([3])
    arr2 = numpy.array([1,3,4])

    m1,m2 = match(arr1, arr2)
    if (m1[0] != 0) or (m2[0] != 1):
        nerror += 1
        _match_test_write(arr1,arr2)

    
    arr1 = numpy.array([3,8,-15, 6])
    arr2 = numpy.array([1,7,4, 16, 8, 25, -15])

    m1,m2 = match(arr1, arr2)
    if m1.size != 2:
        nerror += 1
        _match_test_write(arr1,arr2)
    m1.sort()
    m2.sort()
    if (m1[0] != 1) or (m1[1] != 2) or (m2[0] != 4) or (m2[1] != 6):
        nerror += 1
        _match_test_write(arr1,arr2)


    arr1 = numpy.array([3,-100,188, 6])
    arr2 = numpy.array([1,7,4, 16, 8, 25, -15])

    m1,m2 = match(arr1, arr2)
    if (m1.size != 1) or (m2.size != 1) or (m1[0] != -1) or (m2[0] != -1):
        nerror += 1
        _match_test_write(arr1,arr2)

    if nerror > 0:
        stdout.write('Errors: %s\n' % nerror)
    else:
        stdout.write('All tests passed\n')

def list2file(inlist, fname):
    fobj = open(fname, 'w')
    for l in inlist:
        fobj.write(str(l))
        fobj.write('\n')
    fobj.close()


def intersect_sequences(seq1, seq2):
    """
    find the intersection between the sequences
    not finished
    """
    narg = len(arrays)
    stdout.write("narg = %s\n" % narg)

    arr0 = arrays[0]

#def _rev_hist(data, s, binsize, hist, rev):
def _weave_rev_hist(data, s, binsize, hist, rev):
    """
    Weave version of histogram with reverse_indices
    """
    code = """

    int64_t nbin = hist.size();
    int64_t binnum_old = -1;

    // index of minimum value
    int64_t imin = s(0);
    for (int64_t i=0; i<s.size(); i++) {

        int64_t offset = i+nbin+1;
        int64_t data_index = s(i);


        rev(offset) = data_index;

        int64_t binnum = (int64_t) ( (data(data_index)-data(imin))/binsize);

        if (binnum >= 0 && binnum < nbin) {
            if (binnum > binnum_old) {
                int64_t tbin = binnum_old + 1;
                while (tbin <= binnum) {
                    rev(tbin) = offset;
                    tbin++;
                }
            }
            hist(binnum) = hist(binnum) + 1;
            binnum_old = binnum;
        }
    }

    int64_t tbin = binnum_old + 1;
    while (tbin <= nbin) {
        rev(tbin) = rev.size();
        tbin++;
    }

    """

    weave.inline(code, ['data','s','binsize','hist','rev'],
                 type_converters = weave.converters.blitz)
    return


def histogram(data, binsize=1., min=None, max=None, rev=False, use_weave=False):
    """
    simplified interface for now, just binsize
    """
     
    dowhere=False
    s = data.argsort()
    if min is not None:
        dmin = min
        dowhere=True
    else:
        dmin = data[s[0]]

    if max is not None:
        dmax = max
        dowhere=True
    else:
        dmax = data[s[-1]]

    bsize = float(binsize)

    if dowhere:
        # where will preserve order, so subscript with s
        w,=numpy.where( (data[s] >= dmin) & (data[s] <= dmax) )
        if w.size == 0:
            raise ValueError("No data in specified min/max range\n")
        s = s[w]

    nbin = long( (dmax-dmin)/bsize ) + 1
    revsize = s.size + nbin+1

    if rev:
        revind = numpy.zeros(revsize, dtype='i8')
    hist = numpy.zeros(nbin, dtype='i8')

    # populate the array from nbin+1:nbin+1+s.size
    # with the sort indices.  Simultaneosly record bin
    # edges at the beginning of reverse indices

    if use_weave:
        _weave_rev_hist(data, s, bsize, hist, revind)
        return hist, revind

    offset = nbin+1
    i=0
    binnum_old = -1
    while i < s.size:
        data_index = s[i]
        if rev:
            revind[offset] = data_index

        val = data[data_index]

        binnum = long( (val-dmin)/bsize )
        #print 'binnum:',binnum,' binnum old:',binnum_old, 'val:',val
        if binnum >= 0 and binnum < nbin:
        #if binnum >= 0:
            if binnum > binnum_old:
                tbin = binnum_old + 1
                while tbin <= binnum:
                    if rev:
                        revind[tbin] = offset
                        #print '\t\trevind[%d] = %d' % (tbin,offset)
                    tbin += 1

            hist[binnum] += 1
            binnum_old = binnum
        #print 'rev:',revind[binnum]

        i += 1
        offset += 1

    if rev:
        #pass
        # Fill in the last ones
        tbin = binnum_old + 1
        while tbin <= nbin:
            revind[tbin] = revsize
            tbin += 1
        #revind[nbin] = revsize

    if rev:
        return hist, revind
    else:
        return hist

def objhist(objects):
    """
    Count the occurences of objects in an interable input.  The result
    is a dictionary keyed by the objects with values the count
    """
    objdict={}
    for s in objects:
        if s in objdict:
            objdict[s] += 1
        else:
            objdict[s] = 1
    return objdict

def test_histogram(nrand=None, fixed=False, use_weave=False, display=True):

    min=None
    max=None
    if fixed:
        data= \
            numpy.array([11, 62,  4, 69, 56, 
                         32, 80, 31, 19,  4, 
                         19, 81, 37, 89, 14, 
                         55,  2, 36, 81, 76])
        binsize=1
        min = -10
        max=85
    elif nrand is not None:
        data = numpy.random.normal(size=nrand)
        binsize=0.1
    else:
        data = \
            numpy.array([0.01233105,0.1967937,0.56444893,0.57024105,0.95079757,
                         0.80643415,0.36878545,0.5514516,0.87020452,0.27447491,
                         0.30478068,0.10748157,0.12471763,0.45594969,0.80562055,
                         0.62194449,0.38645207,0.54947444,0.56945579,0.44712302])*0.7
        binsize=0.05

    crap="""
    data = \
        numpy.array([ 0.39518813,  0.73823606,  0.55035441,  0.89738103,  0.81981696,
           0.5737272 ,  0.62388723,  0.60248498,  0.58692789, 0.25261468,
           0.0967087 ,  0.96778357,  0.9853306 ,  0.47983496, 0.63522547,
           0.51335563,  0.42617909,  0.53522043,  0.22824262, 0.16393504])
           """


    hist,rev = histogram(data, binsize=binsize, min=min, max=max,
                         rev=True,use_weave=use_weave)

    #print hist
    #print rev

    if display:
        print hist
        print rev[0:hist.size+1]

        #for r in rev: print r

        for i in range(hist.size):
            if rev[i] != rev[i+1]:
                w = rev[rev[i]:rev[i+1]]
                w.sort()
                #print 'rev[%d]: %d rev[%d]: %d' % (i,rev[i],i+1,rev[i+1])
                print '\tbin:',i,' w:',w, 'vals:',data[w]
        print 'nbin =',hist.size



def print_human_readable_filesize(filename, sz_in):
    sz=long(sz_in) 
    if sz > 1024*1024*1024*1024:
        psize="%sT" % (sz/(1024*1024*1024*1024))
    elif sz > 1024*1024*1024:
        psize="%sG" % (sz/(1024*1024*1024))
    elif sz > 1024*1024:
        psize="%sM" % (sz/(1024*1024))
    elif sz > 1024:
        psize="%sK" % (sz/1024)
    else:
        psize="%s" % sz

    stdout.write( "%-10s %s\n" % (psize, filename) )
    return psize


def du(filelist, verbose=True):
    import os
    sztot = 0
    for f in filelist:
        # size in bytes
        sz=os.path.getsize(f)
        sztot += sz

        if verbose:
            print_human_readable_filesize(f, sz)


    if verbose:
        print_human_readable_filesize('total', sztot)

    return sztot

def dict_compare_field(fieldname):
    """
    For sorting a list of dictionaries by field

    list_of_dicts.sort( dict_compare_fields('somename'))
    """
    def compare_two_dicts_byfield(a, b):
        return cmp(a[fieldname], b[fieldname])
    return compare_two_dicts_byfield

def dict_compare_fields(fieldnames):
    """
    For sorting a list of dictionaries by more than one field

    list_of_dicts.sort( dict_compare_fields('somename'))
    """
    def compare_two_dicts_byfields(a, b):
        sa = []
        sb = []
        for f in fieldnames:
            sa.append(a[f])
            sb.append(b[f])
        sa='-'.join(sa)
        sb='-'.join(sb)
        return cmp(sa, sb)
    return compare_two_dicts_byfields



def dict_print(d, keys=None):
    if keys is None:
        print d
        return
    newd={}
    for k in keys:
        newd[k] = d[k]

    dict_print(newd)



