import numpy
import pyfits
from sys import stdout


def compare_arrays(arr1in, arr2in, verbose=False):
    """
    Compare the values field-by-field in two sets of numpy arrays or
    recarrays.
    """

    arr1 = arr1in.view(numpy.ndarray)
    arr2 = arr2in.view(numpy.ndarray)

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
            stdout.write("    testing field: '%s'\n" % n2)
            stdout.write('        shape...........')
        if arr2[n2].shape != arr1[n1].shape:
            nfail += 1
            if verbose:
                stdout.write('shapes differ\n')
        else:
            if verbose:
                stdout.write('OK\n')
                stdout.write('        elements........')
            w,=numpy.where(arr1[n1].ravel() != arr2[n2].ravel())
            if w.size > 0:
                nfail += 1
                if verbose:
                    stdout.write('\n        '+\
                        '%s elements in field %s differ\n' % (w.size,n2))
            else:
                if verbose:
                    stdout.write('OK\n')

    if nfail == 0:
        if verbose:
            stdout.write('All tests passed\n')
        return True
    else:
        if verbose:
            stdout.write('%d differences found\n' % nfail)
        return False



def get_test_data(verbose=False):
    st = numpy.zeros(3, [('f1','f4'),('f2','S6'),('f3','2f8')])

    numpy.random.seed(35)
    st['f1'] = [1,3,5]
    st['f2'] = ['hello','world','byebye']
    st['f3'] = numpy.random.random(st['f3'].shape)
    print st.dtype.descr
    print st

    return st


def test(verbose=False, copy=True, noswap=False):

    fname='qsodata.fits'

    print 'Reading from ',fname
    data1,h1 = pyfits.getdata(fname, ext=1, header=True)
    data2,h2 = pyfits.getdata(fname, ext=2, header=True) 

    st = get_test_data()

    outfile = 'test.fits'
    print 'Writing to file data1:',outfile
    pyfits.writeto(outfile, data1, clobber=True, noswap=noswap, copy=copy)
    print 'Appending to file: data2',outfile
    pyfits.writeto(outfile, data2, append=True, noswap=noswap, copy=copy)

    print 'Appending to file: st',outfile
    pyfits.writeto(outfile, st, append=True, noswap=noswap, copy=copy)
    print st.dtype.descr
    print st


    print 'Reading data back'
    data1check, h1check = pyfits.getdata(outfile, ext=1, header=True)
    data2check, h2check = pyfits.getdata(outfile, ext=2, header=True)
    stcheck, sthcheck = pyfits.getdata(outfile, ext=3, header=True)

    if verbose:
        print acheck
        print stcheck

    if not compare_arrays(data1, data1check, verbose=True):
        raise ValueError,'Fail'
    if not compare_arrays(data2, data2check, verbose=True):
        raise ValueError,'Fail'
    if not compare_arrays(st, stcheck, verbose=True):
        raise ValueError,'Fail'
