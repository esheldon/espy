from __future__ import print_function
import os
import numpy
from sys import stdout
import copy

def multiview(image, **keys):
    """
    View the image and also some cross-sections through it.  Good for
    postage stamp type images
    """

    import biggles

    cen = keys.get('cen', None)
    if cen is None:
        # use the middle as center
        cen = [(image.shape[0]-1)/2., (image.shape[1]-1)/2.]

    keys2 = copy.copy(keys)
    keys2['show'] = False
    imp = view(image, **keys2)

    # cross-section across rows
    imrows = image[:, cen[1]]
    imcols = image[cen[0], :]
    

    crossplt = biggles.FramedPlot()
    hrows = biggles.Histogram(imrows, color='blue')
    hrows.label = 'Center rows'
    hcols = biggles.Histogram(imcols, color='red')
    hcols.label = 'Center columns'

    key = biggles.PlotKey(0.1, 0.9, [hrows, hcols])

    crossplt.add(hrows, hcols, key)
    crossplt.aspect_ratio=1
    yr = crossplt._limits1().yrange()
    yrange = (yr[0], yr[1]*1.2)
    crossplt.yrange = yrange


    tab = biggles.Table( 1, 2 )

    tab[0,0] = imp
    tab[0,1] = crossplt

    title=keys.get('title',None)
    if title:
        tab.title=title

    show = keys.get('show', True)
    if show:
        tab.show()
    return tab

def view(image, **keys): 
    import biggles

    # image_scale does not reverse by default, but for display purposes we
    # will want to
    if 'reverse' not in keys:
        keys['reverse'] = True
    im = scale(image, **keys)

    # we need to transpose for biggles to display properly
    trans = keys.get('transpose',True)
    if trans:
        im = im.transpose()

    plt = biggles.FramedPlot()

    s = im.shape
    # this is a difference from Contours which can be alarming
    #ranges = ((0, 0), (s[0], s[1]))
    ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
    d = biggles.Density(im, ranges)
    plt.add(d)

    levels=keys.get('levels',None)
    if levels is not None:
        ccolor = keys.get('ccolor','grey')
        c = biggles.Contours(im, color=ccolor)
        c.levels = levels
        plt.add(c)

    # make sure the pixels look square
    plt.aspect_ratio = im.shape[1]/float(im.shape[0])

    if 'epsfile' in keys:
        epsfile=os.path.expandvars(keys['epsfile'])
        epsfile=os.path.expanduser(epsfile)
        plt.write_eps(epsfile)
    show=keys.get('show',True)
    if show:
        plt.show()

    return plt

def scale(image, **keys):
    reverse=keys.get('reverse', False)
    scaling = keys.get('scaling', 'asinh')

    if 'min' in keys or 'max' in keys:
        doclip=True
    else:
        doclip=False

    minval = keys.get('min', None)
    maxval = keys.get('max', None)

    thismax = image.max()

    if minval is None:
        minval = image.min()

    if maxval is None:
        setmax=False
        maxval = thismax
    else:
        setmax=True

    im = numpy.array(image, copy=True, dtype='f4')
    if doclip:
        im.clip(minval, maxval, im)

    # subtract so min is zero
    im -= minval
    # and make sure in [0,1]
    im /= (maxval-minval)

    if scaling == 'asinh':
        #alpha        = keys.get('alpha',0.02)
        #nonlinearity = keys.get('nonlinearity',8.0)

        # image scaled [0-1], make The crossover point to log scaling 0.001,
        # which is alpha=1000 this won't work for everything of course.

        alpha        = keys.get('alpha',1000.0)
        nonlinearity = keys.get('nonlinearity',1.0)

        # do asinh scaling
        # arcsinh( alpha*nonlinearity*image )/nonlinearity
        #im = asinh_scale(im, alpha=alpha, nonlinearity=nonlinearity)

        im *= alpha*nonlinearity
        numpy.arcsinh(im, im)
        im /= nonlinearity

        if setmax:
            tmv = asinh_scale(thismax, alpha=alpha, nonlinearity=nonlinearity)
            mv = asinh_scale(maxval, alpha=alpha, nonlinearity=nonlinearity)


    elif scaling == 'linear':
        pass
    else:
        raise ValueError('Uknown scaling: %s' % scaling)

    im /= im.max()
    if setmax:
        # e.g. if we set the max to a value *higher* than our max input value, 
        # the max now will be less than one
        
        im *= (tmv/mv)
        # in case max value was less than input max
        im.clip(0.,1.,im)

    if reverse:
        # reverse the color. For biggles this is what we want since background
        # can only be white
        # im = 1.0-im
        im *= -1
        im += 1

    return im

def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)

    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = numpy.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]


def rebin_old(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)

    To preserve the counts, multiple by the factors as needed.  For example, 
    if reducing from (4,4) to (2,2) you need to multiply by 4
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = numpy.asarray(shape)/numpy.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in xrange(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in xrange(lenShape)] + \
             ['/factor[%d]'%i for i in xrange(lenShape)]
    print(''.join(evList))
    return eval(''.join(evList))


def asinh_scale(image, alpha=0.02, nonlinearity=8.0, dtype='f4'):
    #image_out = numpy.array( image, dtype='f4')
    image_out = \
        numpy.arcsinh( alpha*nonlinearity*image )/nonlinearity

    return image_out


def imprint(im):
    if len(im.shape) != 2:
        raise ValueError("image must be 2-dimensional")

    if im.dtype.char == 'f':
        f = '%+e'
    else:
        maxint = im.max()
        minint = im.min()
        l = max( len(str(maxint)), len(str(minint)) )
        f = '%'+str(l)+'d'

    nrow = im.shape[0]
    ncol = im.shape[1]

    for row in xrange(nrow):
        for col in xrange(ncol):
            stdout.write(f % im[row,col] )
            if col < (ncol-1):
                stdout.write(" ")
        stdout.write("\n")


def expand(image, new_dims, padval=0, verbose=False):
    """

    Expand an image to the specified size.  Pad with the specified value
    (default 0).  If the new dims are all less than the existing image, the
    original image is returned.

    """

    sz = image.shape
    if new_dims[0] > sz[0] or new_dims[1] > sz[1]:
        srow = max(sz[0], new_dims[0])
        scol = max(sz[1], new_dims[1])

        if verbose:
            print("  expanding image from",sz,"to:",[srow,scol])

        new_image = numpy.empty( (srow, scol), dtype=image.dtype )
        new_image[:] = padval

        new_image[0:sz[0], 0:sz[1]] = image[0:sz[0], 0:sz[1]]
        
        return new_image
    else:
        return image

def compare_images(im1, im2, cen=None, minval=None, 
                   label1='im1',label2='im2',
                   show=True,
                   skysig=None, **keys):
    import biggles

    labelres='%s-%s' % (label2,label1)

    biggles.configure( 'default', 'fontsize_min', 1.)

    if im1.shape != im2.shape:
        raise ValueError("images must be the same shape")

    if cen is None:
        cen = [(im1.shape[0]-1)/2., (im1.shape[1]-1)/2.]

    resid = im2-im1

    maxval = max( im1.max(), im2.max() )
    if minval is None:
        minval = min( im1.min(), im2.min() )

    levels=keys.get('levels',7)
    tab=biggles.Table(2,3)
    if 'title' in keys:
        tab.title=keys['title']
    #tab=biggles.Table(3,2)

    im1plt=view(im1, levels=levels, show=False, min=minval, max=maxval)

    im2plt=view(im2, levels=levels, show=False, min=minval, max=maxval)
    residplt=view(resid, show=False, min=minval, max=maxval)

    if skysig is not None:
        chi2perpix = (resid**2).sum()/skysig**2/im1.size
        lab = biggles.PlotLabel(0.1,0.1,
                    r'$\chi^2/N$: %0.2e' % chi2perpix,halign='left')
    else:
        chi2perpix = (resid**2).sum()/im1.size
        lab = biggles.PlotLabel(0.1,0.1,
                    r'noerr $\chi^2/N$: %0.2e' % chi2perpix,halign='left')
    residplt.add(lab)

    im1plt.title=label1
    im2plt.title=label2
    residplt.title=labelres


    # cross-sections
    im1rows = im1[:,cen[1]]
    im1cols = im1[cen[0],:]
    im2rows = im2[:,cen[1]]
    im2cols = im2[cen[0],:]
    resrows = resid[:,cen[1]]
    rescols = resid[cen[0],:]

    him1rows = biggles.Histogram(im1rows, color='blue')
    him1cols = biggles.Histogram(im1cols, color='blue')
    him2rows = biggles.Histogram(im2rows, color='orange')
    him2cols = biggles.Histogram(im2cols, color='orange')
    hresrows = biggles.Histogram(resrows, color='red')
    hrescols = biggles.Histogram(rescols, color='red')

    him1rows.label = label1
    him2rows.label = label2
    hresrows.label = labelres
    key = biggles.PlotKey(0.1,0.9,[him1rows,him2rows,hresrows]) 

    rplt=biggles.FramedPlot()
    rplt.add( him1rows, him2rows, hresrows,key )

    cplt=biggles.FramedPlot()
    cplt.add( him1cols, him2cols, hrescols )

    rplt.aspect_ratio=1
    cplt.aspect_ratio=1


    tab[0,0] = im1plt
    tab[0,1] = im2plt
    tab[0,2] = residplt
    tab[1,0] = rplt
    tab[1,1] = cplt

    if show:
        tab.show()

    biggles.configure( 'default', 'fontsize_min', 1.25)

    return tab
