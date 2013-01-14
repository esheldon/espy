from __future__ import print_function
import os
import numpy
from numpy import array, arcsinh, zeros,where
from sys import stdout, stderr
import copy

def view(image, **keys): 
    """
    View the image and return the biggles plot object

    If show=False just return the plot object.

    Values below zero are clipped, so pre-subtract any background as
    necessary.

    parameters
    ----------
    image or r,g,b: ndarray
        The image(s) as a 2-d array or 
    type: string
        Type of plot, 'dens', 'cont', 'dens-cont'.  (cont is short
        for contour).  Default is 'dens'.
    nonlinear:
        Non-linear scale for an asinh scaling.  If not sent a linear scale is
        used. See the asinh_scale() function.  For asinh scaling you must scale
        your image so that it does not saturate.
    autoscale:
        For linear scaling, re-scale the image so that the maximum value is
        1.0.  This guarantees all pixels are shown without saturation.  

        For asinh scaling you are responsible for performing such scalings
        Subtract this value from the image and clip values below zero.
    xlabel: string
        Label for the X axis
    ylabel: string
        Label for the Y axis
    show: bool
        Set False to not show the image in an X window.
    levels: int
        Number of levels for contours. Default 8
    ccolor:
        color for contour, default 'white' unless type is dens-cont
        then grey.
    transpose: bool
        Transpose the image.  Default True because of biggles' conventions.
    file: string
        Write the image to the intput file name.  .png will be written
        as a png file, else an eps file.
    dims:
        [width,height] for png output, default 800x800
    """
    import biggles

    # we need to transpose for biggles to display properly
    trans = keys.get('transpose',True)

    im=scale_image(image, **keys)

    if trans:
        im = im.transpose()

    type=keys.get('type','dens')

    plt = biggles.FramedPlot()
    if 'title' in keys:
        plt.title=keys['title']

    s = im.shape
    if 'xdr' in keys and 'ydr' in keys:
        xdr=keys['xdr']
        ydr=keys['ydr']
        ranges = ((xdr[0], ydr[0]), (xdr[1], ydr[1]))

        x=numpy.linspace(xdr[0],xdr[1],s[0])
        y=numpy.linspace(ydr[0],ydr[1],s[1])
    else:
        # this is a difference from Contours which can be alarming
        ranges = ((-0.5, -0.5), (s[0]-0.5, s[1]-0.5))
        x=None
        y=None

    def_contour_color='white'
    if 'dens' in type:
        d = biggles.Density(im, ranges)
        plt.add(d)
        def_contour_color='grey'

    if 'cont' in type:
        levels=keys.get('levels',8)
        if levels is not None:
            ccolor = keys.get('ccolor',def_contour_color)
            c = biggles.Contours(im,x=x,y=y,color=ccolor)
            c.levels = levels
            plt.add(c)

    # make sure the pixels look square
    plt.aspect_ratio = im.shape[1]/float(im.shape[0])

    if 'xlabel' in keys:
        plt.xlabel=keys['xlabel']
    if 'ylabel' in keys:
        plt.ylabel=keys['ylabel']

    if 'file' in keys and keys['file'] is not None:
        file=os.path.expandvars(keys['file'])
        file=os.path.expanduser(file)
        if '.png' in file:
            png_dims=keys.get('dims',[800,800])
            plt.write_img(dims[0],dims[1],pngfile)
        else:
            plt.write_eps(file)

    show=keys.get('show',True)
    if show:
        plt.show()

    return plt

def write_jpg(filename, image, **keys):
    """
    Write the scaled image to a jpeg file.
    """
    from PIL import Image
    #from numpy import flipud

    im=scale_image(image, **keys)
    #im=flipud(im)

    #im = (im*256).astype('u1')
    im = (im*255).astype('u1')
    print('max val:',im.max())

    pim=Image.fromarray(im)

    fname=os.path.expandvars(filename)
    fname=os.path.expanduser(fname)
    pim.save(fname, **keys)


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

def compare_images(im1, im2, **keys):
    import biggles

    show=keys.get('show',True)
    skysig=keys.get('skysig',None)
    dof=keys.get('dof',None)

    label1=keys.get('label1','im1')
    label2=keys.get('label2','im2')

    cen=keys.get('cen',None)
    if cen is None:
        cen = [(im1.shape[0]-1)/2., (im1.shape[1]-1)/2.]

    labelres='%s-%s' % (label2,label1)

    biggles.configure( 'default', 'fontsize_min', 1.)

    if im1.shape != im2.shape:
        raise ValueError("images must be the same shape")


    resid = im2-im1

    # will only be used if type is contour
    tab=biggles.Table(2,3)
    if 'title' in keys:
        tab.title=keys['title']
    #tab=biggles.Table(3,2)

    tkeys=copy.deepcopy(keys)
    tkeys['show']=False
    im1plt=view(im1, **tkeys)
    im2plt=view(im2, **tkeys)

    tkeys['nonlinear']=None
    tkeys['min'] = resid.min()
    tkeys['max'] = resid.max()
    residplt=view(resid, **tkeys)

    if skysig is not None:
        if dof is None:
            dof=im1.size
        chi2per = (resid**2).sum()/skysig**2/dof
        lab = biggles.PlotLabel(0.1,0.1,
                    r'$\chi^2/dof$: %0.2f' % chi2per,halign='left')
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
    rplt.xlabel = 'Center Rows'

    cplt=biggles.FramedPlot()
    cplt.add( him1cols, him2cols, hrescols )
    cplt.xlabel = 'Center Columns'

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

def image_read_text(fname):
    """
    Read the simple text image format:
        nrows ncols
        im00 im01 im02...
        im10 im11 im12...
    """

    with open(fname) as fobj:
        ls=fobj.readline().split()
        nrows=int(ls[0])
        ncols=int(ls[1])


        image=numpy.fromfile(fobj, 
                             sep=' ', 
                             count=(nrows*ncols), 
                             dtype='f8').reshape(nrows,ncols)
    return image


def _get_max_image(im1, im2, im3):
    maximage=im1.copy()

    w=where(im2 > maximage)
    if w[0].size > 1:
        maximage[w] = im2[w]

    w=where(im3 > maximage)
    if w[0].size > 1:
        maximage[w] = im3[w]

    return maximage

def _fix_hard_satur(r, g, b, satval):
    """
    Clip to satval but preserve the color
    """

    # make sure you send scales such that this occurs at 
    # a reasonable place for your images

    maximage=_get_max_image(r,g,b)

    w=where(maximage > satval)
    if w[0].size > 1:
        # this preserves color
        fac=satval/maximage[w]
        r[w] *= fac
        g[w] *= fac
        b[w] *= fac
        maximage[w]=satval

    return maximage

def _fix_rgb_satur(r,g,b,fac):
    """
    Fix the factor so we don't saturate the
    RGB image (> 1)

    maximage is the
    """
    maximage=_get_max_image(r,g,b)

    w=where( (r*fac > 1) | (g*fac > 1) | (b*fac > 1) )
    if w[0].size > 1:
        # this preserves color
        fac[w]=1.0/maximage[w]


def get_color_image(imr, img, imb, **keys):
    """
    Create a color image.

    The idea here is that, after applying the asinh scaling, the color image
    should basically be between [0,1] for all filters.  Any place where a value
    is > 1 the intensity will be scaled back in all but the brightest filter
    but color preserved.

    In other words, you develaop a set of pre-scalings the images so that after
    multiplying by

        asinh(I/nonlinear)/(I/nonlinear)

    the numbers will be mostly between [0,1].  You can send scales using the
    scale= keyword

    It can actually be good to have some color saturation so don't be too
    agressive.  You'll have to play with the numbers for each image.

    Note also the image is clipped at zero.

    TODO:
        Implement a "saturation" level in the raw image values as in
        djs_rgb_make.  Even better, implement an outside function to do this.
    """

    nonlinear=keys.get('nonlinear',1.0)
    scales=keys.get('scales',None)
    satval=keys.get('satval',None)

    r = imr.astype('f4')
    g = img.astype('f4')
    b = imb.astype('f4')

    r.clip(0.,r.max(),r)
    g.clip(0.,g.max(),g)
    b.clip(0.,b.max(),b)

    if scales is not None:
        r *= scales[0]
        g *= scales[1]
        b *= scales[2]

    if satval is not None:
        # note using rescaled images so the satval
        # means the same thing (e.g. in terms of real flux)
        maximage=_fix_hard_satur(r,g,b,satval)


    # average images and divide by the nonlinear factor
    fac=1./nonlinear/3.
    I = fac*(r + g + b)

    # make sure we don't divide by zero
    # due to clipping, average value is zero only if all are zero
    w=where(I <= 0)
    if w[0].size > 0:
        I[w] = 1./3. # value doesn't matter images are zero

    f = arcsinh(I)/I

    # limit to values < 1
    # make sure you send scales such that this occurs at 
    # a reasonable place for your images
    _fix_rgb_satur(r,g,b,f)

    R = r*f
    G = g*f
    B = b*f

    st=R.shape
    colorim=zeros( (st[0], st[1], 3) )

    colorim[:,:,0] = R[:,:]
    colorim[:,:,1] = G[:,:]
    colorim[:,:,2] = B[:,:]

    return colorim




def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)

    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = numpy.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]


def scale_image(im, **keys):
    nonlinear=keys.get('nonlinear',None)
    if nonlinear is not None:
        return asinh_scale(im, nonlinear)
    else:
        return linear_scale(im, **keys)

def linear_scale(im, **keys):
    autoscale=keys.get('autoscale',True)

    I=im.astype('f4')

    if autoscale:
        maxval=I.max()
        I  *= (1.0/I.max())

    I.clip(0.0, 1.0, I)

    return I

def asinh_scale(im, nonlinear):
    """
    Scale the image using and asinh stretch

        I = image*f
        f=asinh(image/nonlinear)/(image/nonlinear)

    Values greater than 1.0 after this scaling will be shown as white, so you
    are responsible for pre-scaling your image.  Values < 0 are clipped.

    parameters
    ----------
    image:
        The image.
    nonlinear: keyword
        The non-linear scale.
    """
    from numpy import arcsinh
    I=im.astype('f4')

    I *= (1./nonlinear)

    # make sure we don't divide by zero
    w=where(I <= 0)
    if w[0].size > 0:
        I[w] = 1. # value doesn't matter since images is zero

    f = arcsinh(I)/I

    imout = im*f

    imout.clip(0.0, 1.0, imout)

    return imout



def imprint(im, stream=stdout, fmt=None):
    if len(im.shape) != 2:
        raise ValueError("image must be 2-dimensional")

    if fmt is None:
        if im.dtype.char in ['f','d']:
            fmt = '%+e'
        else:
            maxint = im.max()
            minint = im.min()
            l = max( len(str(maxint)), len(str(minint)) )
            fmt = '%'+str(l)+'d'

    nrow = im.shape[0]
    ncol = im.shape[1]

    for row in xrange(nrow):
        for col in xrange(ncol):
            stream.write(fmt % im[row,col] )
            if col < (ncol-1):
                stream.write(" ")
        stream.write("\n")


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


