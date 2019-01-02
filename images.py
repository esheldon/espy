from __future__ import print_function
try:
    xrange
except:
    xrange=range

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
        for contour, dens for density).  Default is 'dens'.
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
        Transpose the image.  Default False.
    file: string
        Write the image to the intput file name.  .png will be written
        as a png file, else an eps file.
    dims:
        [width,height] for png output, default 800x800
    """
    import biggles

    # we need to transpose for biggles to display properly
    trans = keys.pop('transpose',True)
    doscale = keys.pop('scale',False)

    if len(image.shape) == 2:
        im=scale_image(image, **keys)
        if trans:
            im = im.transpose()
    else:
        # for 3 d we need to trust the input to
        # be properly scaled and transposed
        im=image


    type=keys.get('type','dens')

    plt = biggles.FramedPlot()

    if 'title' in keys:
        plt.title=keys['title']

    x, y, ranges = _extract_data_ranges(im.shape[0:0+2], **keys)

    def_contour_color='black'
    if 'dens' in type:
        d = biggles.Density(im, ranges)
        plt.add(d)
        def_contour_color='grey'

    if 'cont' in type:
        levels=keys.get('levels',8)
        if levels is not None:
            if 'zrange' in keys:
                zrange=keys['zrange']
            else:
                zrange=numpy.percentile(im, [1.0,100.0])

            ccolor = keys.get('ccolor',def_contour_color)
            c = biggles.Contours(im,x=x,y=y,color=ccolor,zrange=zrange)
            c.levels = levels
            plt.add(c)

    if 'xrange' in keys:
        plt.xrange=keys['xrange']
    if 'yrange' in keys:
        plt.yrange=keys['yrange']

    # make sure the pixels look square
    aratio=im.shape[1]/float(im.shape[0])
    plt.aspect_ratio = aratio

    if 'xlabel' in keys:
        plt.xlabel=keys['xlabel']
    if 'ylabel' in keys:
        plt.ylabel=keys['ylabel']

    if doscale and 'dens' in type:
        # the following all works as long as the font size doesn't
        # get to small on the scale.  Making sure the biggles
        # fontsize_min is less than about 0.75 works even in
        # some bad cases
        immin, immax = image.min(), image.max()
        num=100

        tab_aratio = (im.shape[1] * 0.8)/float(im.shape[0])
        tab=biggles.Table(
            1,2,
            aspect_ratio=tab_aratio,
            col_fractions=[0.1,0.9],
        )
        tab.cellpadding=0.0
        tab.cellspacing=0.0

        scale_plt=biggles.FramedPlot(
            yrange=[immin,immax],
            xrange=[0.0,1],
        )

        scale_im = numpy.linspace(0, 1, num).reshape(1,num)
        drngs = ((0,immin),(1,immax))
        d = biggles.Density(scale_im, drngs)
        scale_plt.add(d)

        scale_plt.x1.draw_ticks=False
        scale_plt.x1.draw_ticklabels=False
        scale_plt.x2.draw_ticks=False
        scale_plt.x2.draw_ticklabels=False

        scale_plt.y1.draw_ticks=False
        scale_plt.y1.draw_ticklabels=False
        scale_plt.y2.draw_ticks=False
        scale_plt.y2.draw_ticklabels=True


        tab[0,0] = plt
        tab[0,1] = scale_plt
        plt=tab
        
    _writefile_maybe(plt, **keys)
    _show_maybe(plt, **keys)

    return plt


def get_profile(image, cen=None):
    if cen is None:
        cen=(numpy.array(image.shape)-1.0)/2.0
    else:
        assert len(cen)==2,"cen must have two elements"

    rows,cols=numpy.mgrid[
        0:image.shape[0],
        0:image.shape[1],
    ]

    rows = rows.astype('f8')-cen[0]
    cols = cols.astype('f8')-cen[1]

    r = numpy.sqrt(rows**2 + cols**2).ravel()
    s = r.argsort()
    r=r[s]
    pim = image.ravel()[s]

    return r, pim

def view_profile(image, **keys): 
    """
    View the image as a radial profile vs radius from the center.

    If show=False just return the plot object.

    Values below zero are clipped, so pre-subtract any background as
    necessary.

    parameters
    ----------
    image or r,g,b: ndarray
        The image(s) as a 2-d array or 
    cen: [c1,c2], optional
        Optional center
    show: bool
        Set False to not show the image in an X window.
    file: string
        Write the image to the intput file name.  .png will be written
        as a png file, else an eps file.
    width, height: integers
        Size for output
    **kw:
        keywords for the FramedPlot and for the output image dimensions
        width
    """
    import biggles

    plt=keys.pop('plt',None)
    show=keys.pop('show',True)
    cen=keys.pop('cen',None)
    keys['xlabel']=keys.get('xlabel','radius')

    r, pim = get_profile(image, cen=cen)

    keys['visible']=False
    plt = biggles.plot(r, pim, plt=plt, **keys)

    _writefile_maybe(plt, **keys)
    _show_maybe(plt, show=show, **keys)

    return plt



def _writefile_maybe(plt, **keys):
    if 'file' in keys and keys['file'] is not None:
        file=os.path.expandvars(keys['file'])
        file=os.path.expanduser(file)
        if '.png' in file:
            dims=keys.get('dims',[800,800])
            plt.write_img(dims[0],dims[1],file)
        else:
            plt.write_eps(file)

def _show_maybe(plt, **keys):
    show=keys.get('show',None)
    if show is None:
        if 'file' in keys and keys['file'] is not None:
            # don't show anything if show not explicitly
            # sent and we are writing a file
            show=False
        else:
            show=True
    else:
        show=keys['show']

    if show:
        wkeys={}
        dims=keys.get('dims',None)
        if dims is None:
            width=keys.get('width',None)
            height=keys.get('height',None)
            if width is not None:
                dims=[width,height]

        if dims is not None:
            wkeys['width']=dims[0]
            wkeys['height']=dims[1]
        plt.show(**wkeys)


def _extract_data_ranges(imshape, **keys):
    if 'xdr' in keys and 'ydr' in keys:
        xdr=keys['xdr']
        ydr=keys['ydr']
        ranges = ((xdr[0], ydr[0]), (xdr[1], ydr[1]))

        x=numpy.linspace(xdr[0],xdr[1],imshape[0])
        y=numpy.linspace(ydr[0],ydr[1],imshape[1])
    elif 'ranges' in keys:
        ranges = keys['ranges']
        xmin,ymin = ranges[0]
        xmax,ymax = ranges[1]
        x=numpy.linspace(xmin, xmax, imshape[0])
        y=numpy.linspace(ymin, ymax, imshape[1])
    else:
        # this is a difference from Contours which can be alarming
        ranges = ((-0.5, -0.5), (imshape[0]-0.5, imshape[1]-0.5))
        x=None
        y=None

    return x, y, ranges
 
def make_combined_mosaic(imlist):
    """
    only works if all images are same size

    Also should be "sky subtracted" for best
    effect when the grid is not fully packed
    """
    import plotting
    nimage=len(imlist)
    grid=plotting.Grid(nimage)
    shape=imlist[0].shape

    imtot=numpy.zeros( (grid.nrow*shape[0], grid.ncol*shape[1]) )

    for i in xrange(nimage):
        im=imlist[i]
        row,col=grid(i)

        rstart = row*shape[0]
        rend   = (row+1)*shape[0]

        cstart = col*shape[1]
        cend   = (col+1)*shape[1]

        imtot[rstart:rend, cstart:cend] = im

    return imtot


def view_mosaic(imlist, titles=None, combine=False, **keys):
    import biggles
    import plotting

    tabtitle=keys.pop('title',None)

    if combine:
        imtot=make_combined_mosaic(imlist)
        return view(imtot, **keys)

    nimage=len(imlist)
    grid=plotting.Grid(nimage)

    tab=biggles.Table(grid.nrow,grid.ncol)

    tkeys={}
    tkeys.update(keys)
    tkeys['show']=False

    for i in xrange(nimage):
        im=imlist[i]

        row,col = grid(i)

        if titles is not None:
            title=titles[i]
        else:
            title=None

        implt=view(im, title=title, **tkeys)
        tab[row,col] = implt

    # aspect is ysize/xsize
    aspect=keys.get('aspect',None)
    if aspect is None:
        aspect = float(grid.nrow)/grid.ncol
    tab.aspect_ratio=aspect

    tab.title=tabtitle

    _writefile_maybe(tab, **keys)
    _show_maybe(tab, **keys)

    return tab

def bytescale(im):
    """ 
    The input should be between [0,1]

    output is [0,255] in a unsigned byte array
    """
    imout = (im*255).astype('u1')
    return imout

def write_image(filename, image, **keys):
    """
    Write an image to an image file.

    the image must be bytescaled between [0,255] and by of type 'u1'.  See
    the scale_image() and bytescale() functions.

    The file type, compression scheme, etc are determined by the file name.
    Extra keywords such as quality for jpegs are passed along.
    """
    from PIL import Image

    pim=Image.fromarray(image)

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
        cen = [
            int(round( (image.shape[0]-1)/2. )),
            int(round( (image.shape[1]-1)/2. )),
        ]

    assert len(cen)==2

    keys2 = copy.copy(keys)
    keys2['show'] = False
    keys2['file'] = None

    do_profile=keys2.pop('profile',False)

    tab = biggles.Table( 1, 2 )

    imp = view(image, **keys2)
    tab[0,0] = imp

    if do_profile:
        tab[0,1] = view_profile(image, **keys2)
        tab[0,1].aspect_ratio=1
    else:
        # cross-section across rows
        imrows = image[:, cen[1]]
        imcols = image[cen[0], :]
        
        # in case xdr, ydr were sent
        x, y, ranges = _extract_data_ranges(image.shape, **keys)
        if x is not None:
            x0 = ranges[0][0]
            y0 = ranges[0][1]
            xbinsize=x[1]-x[0]
            ybinsize=y[1]-y[0]
        else:
            x0=0
            y0=0
            xbinsize=1
            ybinsize=1


        crossplt = biggles.FramedPlot()
        hrows = biggles.Histogram(imrows, x0=y0, binsize=ybinsize, color='blue')
        hrows.label = 'Center rows'
        hcols = biggles.Histogram(imcols, x0=x0, binsize=xbinsize, color='red')
        hcols.label = 'Center columns'

        key = biggles.PlotKey(0.1, 0.9, [hrows, hcols])

        crossplt.add(hrows, hcols, key)
        crossplt.aspect_ratio=1
        yr = crossplt._limits1().yrange()
        yrange = (yr[0], yr[1]*1.2)
        crossplt.yrange = yrange

        tab[0,1] = crossplt


    _writefile_maybe(tab, **keys)
    _show_maybe(tab, **keys)
    return tab

def compare_images(im1, im2, **keys):
    import biggles

    show=keys.get('show',True)
    skysig=keys.get('skysig',None)
    dof=keys.get('dof',None)
    cross_sections=keys.get('cross_sections',True)
    ymin=keys.get('min',None)
    ymax=keys.get('max',None)

    color1=keys.get('color1','blue')
    color2=keys.get('color2','orange')
    colordiff=keys.get('colordiff','red')

    nrow=2
    if cross_sections:
        ncol=3
    else:
        ncol=2

    label1=keys.get('label1','im1')
    label2=keys.get('label2','im2')

    cen=keys.get('cen',None)
    if cen is None:
        cen = [(im1.shape[0]-1)/2., (im1.shape[1]-1)/2.]

    labelres='%s-%s' % (label1,label2)

    biggles.configure( 'default', 'fontsize_min', 1.)

    if im1.shape != im2.shape:
        raise ValueError("images must be the same shape")


    #resid = im2-im1
    resid = im1-im2

    # will only be used if type is contour
    tab=biggles.Table(nrow,ncol)
    if 'title' in keys:
        tab.title=keys['title']

    tkeys=copy.deepcopy(keys)
    tkeys['show']=False
    tkeys['file']=None
    im1plt=view(im1, **tkeys)
    im2plt=view(im2, **tkeys)

    tkeys['nonlinear']=None
    # this has no effect
    tkeys['min'] = resid.min()
    tkeys['max'] = resid.max()
    residplt=view(resid, **tkeys)

    if skysig is not None:
        if dof is None:
            dof=im1.size
        chi2per = (resid**2).sum()/skysig**2/dof
        lab = biggles.PlotLabel(0.1,0.1,
                                r'$\chi^2/dof$: %0.2f' % chi2per,
                                color='red',
                                halign='left')
    else:
        if dof is None:
            dof=im1.size
        chi2per = (resid**2).sum()/dof
        lab = biggles.PlotLabel(0.1,0.1,
                                r'$\chi^2/npix$: %.3e' % chi2per,
                                color='red',
                                halign='left')
    residplt.add(lab)

    im1plt.title=label1
    im2plt.title=label2
    residplt.title=labelres


    # cross-sections
    if cross_sections:
        cen0=int(cen[0])
        cen1=int(cen[1])
        im1rows = im1[:,cen1]
        im1cols = im1[cen0,:]
        im2rows = im2[:,cen1]
        im2cols = im2[cen0,:]
        resrows = resid[:,cen1]
        rescols = resid[cen0,:]

        him1rows = biggles.Histogram(im1rows, color=color1)
        him1cols = biggles.Histogram(im1cols, color=color1)
        him2rows = biggles.Histogram(im2rows, color=color2)
        him2cols = biggles.Histogram(im2cols, color=color2)
        hresrows = biggles.Histogram(resrows, color=colordiff)
        hrescols = biggles.Histogram(rescols, color=colordiff)

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
    else:
        tab[0,0] = im1plt
        tab[0,1] = im2plt
        tab[1,0] = residplt

    _writefile_maybe(tab, **keys)
    _show_maybe(tab, **keys)

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
        if maxval != 0.0:
            I  *= (1.0/maxval)

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

def rebin(im, factor, dtype=None):
    """
    Rebin the image so there are fewer pixels.  The pixels are simply
    averaged.
    """
    factor=int(factor)
    s = im.shape
    if ( (s[0] % factor) != 0
            or (s[1] % factor) != 0):
        raise ValueError("shape in each dim (%d,%d) must be "
                   "divisible by factor (%d)" % (s[0],s[1],factor))

    newshape=array(s)/factor
    if dtype is None:
        a=im
    else:
        a=im.astype(dtype)

    return a.reshape(newshape[0],factor,newshape[1],factor,).sum(1).sum(2)/factor/factor

def boost( a, factor):
    """
    Resize an array to larger shape, simply duplicating values.
    """
    from numpy import mgrid
    
    factor=int(factor)
    if factor < 1:
        raise ValueError("boost factor must be >= 1")

    newshape=array(a.shape)*factor

    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]

def expand(image, new_dims, padval=0, verbose=False):
    """

    Expand an image to the specified size.  The extra pixels are set to the
    padval value.  Note this does not change the pixel scale, just adds new
    pixels.  See boost to change the pixel scale.
    
    If the new dims are all less than the existing image, the original image is
    returned.

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



def ds9(im):
    """
    view the image in ds9
    """
    import fitsio
    import os
    import tempfile
    tmpdir=os.environ.get('TMPDIR','/tmp')
    tfile = tempfile.mktemp(suffix='.fits')
    tfile=os.path.join(tmpdir, tfile)
    fitsio.write(tfile, im, clobber=True)
    os.system('ds9 %s' % tfile)
    if os.path.exists(tfile):
        os.remove(tfile)
