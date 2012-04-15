from sys import stderr
import numpy
from numpy import where, sqrt, zeros
import esutil as eu
from esutil.numpy_util import where1
import fitsio
import admom

def process_image_cat(imfile, catfile, get_meta=False):
    """
    Process the input image and cat files.

    If get_meta=True, a tuple of (array, meta) is returned.  meta is a
    dictionary with 'shiftmax' and 'nsub'
    """
    im, ivar_im, cat = read_image_cat(imfile, catfile)

    col=cat['x_image']-1
    row=cat['y_image']-1
    sky=cat['background']

    colint=col.astype('i4')
    rowint=row.astype('i4')
    ivar = ivar_im[rowint,colint] 
    sigsky = zeros(cat.size) 
    w=where1(ivar > 0)
    wbad=where1(ivar <= 0)

    if w.size > 0:
        sigsky[w] = 1.0/sqrt(ivar[w])

        print >>stderr,'running admom'
        tmp=admom.admom(im, row[w], col[w], sky=sky[w], sigsky=sigsky[w])
        print >>stderr,'copying to array output'
        a,h=as_array(cat, w, wbad, tmp)
        if get_meta:
            out=(a,h)
        else:
            out=a
    else:
        out=None

    return out

def read_image_cat(imfile, catfile):
    im_ext = 1
    ivar_ext = 3
    cat_ext = 2
    pos_offset=1

    with eu.hdfs.HDFSFile(imfile,verbose=True) as fobj:
        fobj.stage()
        f = fitsio.FITS(fobj.localfile)
        im = f[im_ext][:,:]
        ivar_im = f[ivar_ext][:,:]

    cat = eu.io.read(catfile, type='fits', verbose=True, 
                     lower=True, ext=cat_ext)

    return im, ivar_im, cat

def as_array(cat, w, wbad, odict):
    dt=[]
    for n in ['row','col','Irr','Irc','Icc','e1','e2',
              'rho4','a4','s2','uncer','s2n','numiter',
              'wrow','wcol','whyflag','whystr',
              'sky','sigsky','guess']:
        if n in odict:
            dt += [(n,odict[n].dtype.str)]
    for n in cat.dtype.names:
        # only copy scalar fields
        if len(cat[n].shape) == 1:
            dt += [(n,cat[n].dtype.str)]

    arr=zeros(cat.size, dtype=dt) 

    eu.numpy_util.copy_fields(cat, arr)

    for n in odict:
        if n in arr.dtype.names:
            arr[n][w] = odict[n]

    if wbad.size > 0:
        arr['whyflag'][wbad] = 2**31
        arr['whystr'][wbad] = 'badivar'
        arr['Irr'][wbad] = -9999.
        arr['Irc'][wbad] = -9999.
        arr['Icc'][wbad] = -9999.
        arr['e1'][wbad] = -9999.
        arr['e2'][wbad] = -9999.
        arr['rho4'][wbad] = -9999.
        arr['uncer'][wbad] = 9999.
        arr['wrow'][wbad] = -9999.
        arr['wcol'][wbad] = -9999.
        arr['a4'][wbad] = -9999.
        arr['s2n'][wbad] = -9999.

    h={'shiftmax':odict['shiftmax'],'nsub':odict['nsub']}
    return arr, h
