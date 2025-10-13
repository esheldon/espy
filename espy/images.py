import os
import numpy as np
from sys import stdout
from .plotting import _show_andor_save


def view(
    image,
    nonlinear=None,
    autoscale=False,
    colorbar=False,
    figax=None,
    cmap='gray',
    dpi=None,
    title=None,
    **kw
):
    """
    View the image and return the plot object

    If show=False just return the plot object.

    Values below zero are clipped, so pre-subtract any background as
    necessary.

    parameters
    ----------
    image or r,g,b: ndarray
        The image(s) as a 2-d array or
    nonlinear:
        Non-linear scale for an asinh scaling.  If not sent a linear scale is
        used. See the asinh_scale() function.  For asinh scaling you must scale
        your image so that it does not saturate.
    autoscale:
        For linear scaling, re-scale the image so that the maximum value is
        1.0.  This guarantees all pixels are shown without saturation.

        For asinh scaling you are responsible for performing such scalings
        Subtract this value from the image and clip values below zero.
    colorbar: bool
        If set to True, a color bar will be drawn next to the image
    show: bool
        Set False to not show the image in an X window.
    file: string
        Write the image to the intput file name.  .png will be written
        as a png file, else an eps file.
    dpi: int
        Dots per inch for image files or display
    plt_kws: dict
        dict of keywords for creating the figure
    plt: figure
        If not sent, a new one is created using the plt_kws
    """

    fig, ax, file, show = _prep_plot(
        figax=figax, title=title, kw=kw,
    )

    if len(image.shape) == 2:
        # for 3D color images we need to trust the input to be properly scaled
        image = scale_image(
            image=image,
            nonlinear=nonlinear,
            autoscale=autoscale,
        )

    pim = ax.imshow(image, cmap=cmap)

    # print(type(super(plt)))
    if colorbar:
        ax.colorbar(pim)

    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)

    return fig, ax


def get_profile(image, cen=None):
    if cen is None:
        cen = (np.array(image.shape) - 1.0) / 2.0
    else:
        assert len(cen) == 2, "cen must have two elements"

    rows, cols = np.mgrid[
        0: image.shape[0],
        0: image.shape[1],
    ]

    rows = rows.astype("f8") - cen[0]
    cols = cols.astype("f8") - cen[1]

    r = np.sqrt(rows ** 2 + cols ** 2).ravel()
    s = r.argsort()
    r = r[s]
    pim = image.ravel()[s]

    return r, pim


def view_profile(
    image,
    cen=None,
    figax=None,
    xlim=None,
    ylim=None,
    xlabel='radius [pixels]',
    ylabel=None,
    xlog=False,
    ylog=False,
    title=None,
    dpi=None,
    **kw
):
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
    **kw:
        keywords for the FramedPlot and for the output image dimensions
    """

    fig, ax, file, show = _prep_plot(
        figax=figax, title=title, kw=kw,
        xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog,
        xlabel=xlabel, ylabel=ylabel,
    )

    r, pim = get_profile(image, cen=cen)

    ax.scatter(r, pim, **kw)

    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)

    return fig, ax


def make_combined_mosaic(imlist):
    """
    only works if all images are same size

    Also should be "sky subtracted" for best
    effect when the grid is not fully packed
    """
    from . import plotting

    nimage = len(imlist)
    grid = plotting.Grid(nimage)
    shape = imlist[0].shape

    imtot = np.zeros((grid.nrow * shape[0], grid.ncol * shape[1]))

    for i in range(nimage):
        im = imlist[i]
        row, col = grid(i)

        rstart = row * shape[0]
        rend = (row + 1) * shape[0]

        cstart = col * shape[1]
        cend = (col + 1) * shape[1]

        imtot[rstart:rend, cstart:cend] = im

    return imtot


def view_mosaic(
    imlist,
    colorbar=False,
    titles=None,
    suptitle=None,
    combine=False,
    plt_kws={},
    figax=None,
    dpi=None,
    nonlinear=None,
    figsize=None,
    **kws,
):
    from . import plotting

    if combine:
        imtot = make_combined_mosaic(imlist)
        return view(
            imtot, figax=figax, colorbar=colorbar, dpi=dpi, plt_kws=plt_kws,
            nonlinear=nonlinear,
        )

    nimage = len(imlist)
    grid = plotting.Grid(nimage)

    aratio = grid.nrow / grid.ncol

    if figsize is None:
        if aratio > 1:
            figsize = (8 / aratio, 8)
        else:
            figsize = (8, 8 * aratio)

    fig, axs, file, show = _prep_plot(
        figax=figax, kw=kws,
        nrows=grid.nrow, ncols=grid.ncol,
        figsize=figsize,
    )

    if titles is None:
        titles = ['im%d' for i in range(nimage)]

    add_plt_kws = {}

    if "constrained_layout" not in plt_kws:
        add_plt_kws["constrained_layout"] = False

    if len(add_plt_kws) > 0:
        plt_kws = _get_updated_keywords(plt_kws, **add_plt_kws)

    for i in range(nimage):
        ax = axs.ravel()[i]
        tmp_plt_kws = _get_updated_keywords(
            plt_kws,
            show=False,
            file=None,
            colorbar=colorbar,
        )
        ax.set_title(titles[i])
        view(
            imlist[i], figax=(fig, ax), nonlinear=nonlinear,
            **tmp_plt_kws,
        )

    nax = axs.size
    if nax > nimage:
        for i in range(nimage, axs.size):
            ax = axs.ravel()[i]
            ax.axis('off')

    if suptitle is not None:
        fig.suptitle(suptitle)

    fig.tight_layout()

    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    return fig, axs


def bytescale(im):
    """
    The input should be between [0,1]

    output is [0,255] in a unsigned byte array
    """
    imout = (im * 255).astype("u1")
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

    pim = Image.fromarray(image)

    fname = os.path.expandvars(filename)
    fname = os.path.expanduser(fname)
    pim.save(fname, **keys)


def multiview(
    image,
    nonlinear=None,
    autoscale=False,
    colorbar=False,
    cen=None,
    profile=False,
    xlim=None,
    ylim=None,
    xlabel=None,
    ylabel=None,
    xlog=False,
    ylog=False,
    title=None,
    cmap='gray',
    dpi=None,
    figax=None,
    **kw,
):
    """
    View the image and also some cross-sections through it.  Good for
    postage stamp type images
    """

    fig, axs, file, show = _prep_2plot(
        figax=figax, kw=kw,
        xlim=xlim, ylim=ylim, xlog=xlog, ylog=ylog,
        xlabel=xlabel, ylabel=ylabel,
    )

    if cen is None:
        # use the middle as center
        cen = [
            int(round((image.shape[0] - 1) / 2.0)),
            int(round((image.shape[1] - 1) / 2.0)),
        ]

    assert len(cen) == 2

    view(
        image, figax=(fig, axs[0]),
        nonlinear=nonlinear, autoscale=autoscale, colorbar=colorbar,
        cmap=cmap,
        show=False, file=None,
    )

    if profile:
        view_profile(
            image,
            show=False, file=None,
            figax=(fig, axs[1]),
            **kw
        )
    else:
        # cross-section across rows
        imrows = image[:, cen[1]]
        imcols = image[cen[0], :]

        yvals = np.arange(image.shape[0])
        xvals = np.arange(image.shape[1])

        # axs[1].plot(yvals, imrows)
        # axs[1].plot(xvals, imcols)

        axs[1].step(
            yvals, imrows,
            label='rows', linestyle='-', marker=None,
            where='mid',
            **kw
        )
        axs[1].step(
            xvals, imcols,
            label='cols', linestyle='-', marker=None,
            where='mid',
            **kw
        )
        axs[1].legend()

    if title is not None:
        fig.format(suptitle=title)

    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    return fig, axs


# def compare_images(
#     im1,
#     im2,
#     cross_sections=True,
#     colorbar=False,
#     color1="blue",
#     color2="orange",
#     label1="im1",
#     label2="im2",
#     colordiff="red",
#     plt_kws={},
#     file_kws={},
#     **kws
# ):
#
#     import hickory
#
#     nrows = 2
#     if cross_sections:
#         ncols = 3
#     else:
#         ncols = 2
#
#     aratio = nrows / ncols
#
#     add_plt_kws = {}
#
#     if "constrained_layout" not in plt_kws:
#         add_plt_kws["constrained_layout"] = False
#
#     if "figsize" not in plt_kws:
#         if aratio > 1:
#             add_plt_kws["figsize"] = (8 / aratio, 8)
#         else:
#             add_plt_kws["figsize"] = (8, 8 * aratio)
#
#     if len(add_plt_kws) > 0:
#         plt_kws = _get_updated_keywords(plt_kws, **add_plt_kws)
#
#     tab = hickory.Table(
#         nrows=nrows,
#         ncols=ncols,
#         **plt_kws,
#     )
#     ax1 = tab[0, 0]
#     ax2 = tab[0, 1]
#     if cross_sections:
#         axresid = tab[0, 2]
#         tab[1, 2].axis('off')
#     else:
#         axresid = tab[1, 0]
#         tab[1, 1].axis('off')
#
#     labelres = "%s-%s" % (label1, label2)
#
#     if im1.shape != im2.shape:
#         raise ValueError("images must be the same shape")
#
#     resid = im1 - im2
#
#     tmp_plt_kws = _get_updated_keywords(
#         plt_kws,
#         show=False,
#         file=None,
#         title=label1,
#         colorbar=colorbar,
#     )
#     view(im1, plt=ax1, **tmp_plt_kws)
#
#     tmp_plt_kws = _get_updated_keywords(
#         plt_kws,
#         show=False,
#         file=None,
#         title=label2,
#         colorbar=colorbar,
#     )
#
#     view(im2, plt=ax2, **tmp_plt_kws)
#
#     tmp_plt_kws = _get_updated_keywords(
#         plt_kws,
#         show=False,
#         file=None,
#         title=labelres,
#         colorbar=colorbar,
#     )
#     view(resid, plt=axresid, **tmp_plt_kws)
#     # if colorbar:
#     #     for ax in plt.axes[:3]:
#     #         ax.colorbar(
#
#     if cross_sections:
#         cen = (np.array(im1.shape) - 1) / 2
#         cen0 = int(cen[0])
#         cen1 = int(cen[1])
#         im1rows = im1[:, cen1]
#         im1cols = im1[cen0, :]
#         im2rows = im2[:, cen1]
#         im2cols = im2[cen0, :]
#         resrows = resid[:, cen1]
#         rescols = resid[cen0, :]
#
#         yvals = np.arange(im1.shape[0])
#         xvals = np.arange(im1.shape[1])
#         tab[1, 0].step(
#             yvals, im1rows,
#             # color=color1,
#             label=label1, linestyle='-', marker=None,
#         )
#         tab[1, 0].step(
#             yvals, im2rows,
#             # color=color2,
#             label=label1, linestyle='-', marker=None,
#         )
#         tab[1, 0].step(
#             yvals, resrows,
#             # color=colordiff,
#             label=labelres, linestyle='-', marker=None,
#         )
#         tab[1, 0].legend()
#         tab[1, 0].set(xlabel="center rows")
#
#         tab[1, 1].step(
#             xvals, im1cols, color=color1, label=label1,
#             linestyle='-', marker=None,
#         )
#         tab[1, 1].step(
#             xvals, im2cols, color=color2, label=label1,
#             linestyle='-', marker=None,
#         )
#         tab[1, 1].step(
#             xvals, rescols, color=colordiff, label=labelres,
#             linestyle='-', marker=None,
#         )
#         tab[1, 1].legend()
#         tab[1, 1].set(xlabel="center cols")
#
#     tab.tight_layout()
#
#     _writefile_maybe(plt=tab, **kws)
#     _show_maybe(plt=tab, **kws)
#
#     return tab


def image_read_text(fname):
    """
    Read the simple text image format:
        nrows ncols
        im00 im01 im02...
        im10 im11 im12...
    """

    with open(fname) as fobj:
        ls = fobj.readline().split()
        nrows = int(ls[0])
        ncols = int(ls[1])

        image = np.fromfile(
            fobj, sep=" ", count=(nrows * ncols), dtype="f8"
        ).reshape(nrows, ncols)
    return image


def _get_max_image(im1, im2, im3):
    maximage = im1.copy()

    w = np.where(im2 > maximage)
    if w[0].size > 1:
        maximage[w] = im2[w]

    w = np.where(im3 > maximage)
    if w[0].size > 1:
        maximage[w] = im3[w]

    return maximage


def _fix_hard_satur(r, g, b, satval):
    """
    Clip to satval but preserve the color
    """

    # make sure you send scales such that this occurs at
    # a reasonable place for your images

    maximage = _get_max_image(r, g, b)

    w = np.where(maximage > satval)
    if w[0].size > 1:
        # this preserves color
        fac = satval / maximage[w]
        r[w] *= fac
        g[w] *= fac
        b[w] *= fac
        maximage[w] = satval

    return maximage


def _fix_rgb_satur(r, g, b, fac):
    """
    Fix the factor so we don't saturate the
    RGB image (> 1)

    maximage is the
    """
    maximage = _get_max_image(r, g, b)

    w = np.where((r * fac > 1) | (g * fac > 1) | (b * fac > 1))
    if w[0].size > 1:
        # this preserves color
        fac[w] = 1.0 / maximage[w]


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

    nonlinear = keys.get("nonlinear", 1.0)
    scales = keys.get("scales", None)
    satval = keys.get("satval", None)
    clip = keys.get("clip", None)

    r = imr.astype("f4")
    g = img.astype("f4")
    b = imb.astype("f4")

    if clip is not None:
        r.clip(clip, r.max(), r)
        g.clip(clip, g.max(), g)
        b.clip(clip, b.max(), b)

    if scales is not None:
        r *= scales[0]
        g *= scales[1]
        b *= scales[2]

    if satval is not None:
        # note using rescaled images so the satval
        # means the same thing (e.g. in terms of real flux)
        _ = _fix_hard_satur(r, g, b, satval)

    # average images and divide by the nonlinear factor
    fac = 1.0 / nonlinear / 3.0
    I = fac * (r + g + b)  # noqa

    # make sure we don't divide by zero
    # due to clipping, average value is zero only if all are zero
    w = np.where(I <= 0)
    if w[0].size > 0:
        I[w] = 1.0 / 3.0  # value doesn't matter images are zero

    f = np.arcsinh(I) / I

    # limit to values < 1
    # make sure you send scales such that this occurs at
    # a reasonable place for your images
    _fix_rgb_satur(r, g, b, f)

    R = r * f
    G = g * f
    B = b * f

    st = R.shape
    colorim = np.zeros((st[0], st[1], 3))

    colorim[:, :, 0] = R[:, :]
    colorim[:, :, 1] = G[:, :]
    colorim[:, :, 2] = B[:, :]

    return colorim


def scale_image(*, image, nonlinear=None, autoscale=False):

    if nonlinear is not None:
        return asinh_scale(image=image, nonlinear=nonlinear)
    elif autoscale:
        return linear_autoscale(image=image)
    else:
        return image


def linear_autoscale(*, image):

    I = image.astype("f4")  # noqa

    maxval = I.max()
    if maxval != 0.0:
        I *= 1.0 / maxval  # noqa

    # I.clip(0.0, 1.0, I)

    return I


def asinh_scale(*, image, nonlinear):
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

    I = image.astype("f4")  # noqa

    I *= 1.0 / nonlinear  # noqa

    # make sure we don't divide by zero
    w = np.where(I <= 0)
    if w[0].size > 0:
        I[w] = 1.0  # value doesn't matter since images is zero

    f = np.arcsinh(I) / I

    imout = image * f

    # imout.clip(0.0, 1.0, imout)

    return imout


def imprint(im, stream=stdout, fmt=None):
    if len(im.shape) != 2:
        raise ValueError("image must be 2-dimensional")

    if fmt is None:
        if im.dtype.char in ["f", "d"]:
            fmt = "%+e"
        else:
            maxint = im.max()
            minint = im.min()
            ln = max(len(str(maxint)), len(str(minint)))
            fmt = "%" + str(ln) + "d"

    nrow = im.shape[0]
    ncol = im.shape[1]

    for row in range(nrow):
        for col in range(ncol):
            stream.write(fmt % im[row, col])
            if col < (ncol - 1):
                stream.write(" ")
        stream.write("\n")


def rebin(im, factor, dtype=None):
    """
    Rebin the image so there are fewer pixels.  The pixels are simply
    averaged.
    """
    factor = int(factor)
    s = im.shape
    if (s[0] % factor) != 0 or (s[1] % factor) != 0:
        raise ValueError(
            "shape in each dim (%d,%d) must be "
            "divisible by factor (%d)" % (s[0], s[1], factor)
        )

    newshape = np.array(s) // factor
    if dtype is None:
        a = im
    else:
        a = im.astype(dtype)

    return (
        a.reshape(
            newshape[0],
            factor,
            newshape[1],
            factor,
        )
        .sum(1)
        .sum(2)
        / factor
        / factor
    )


def boost(a, factor):
    """
    Resize an array to larger shape, simply duplicating values.
    """
    from numpy import mgrid

    factor = int(factor)
    if factor < 1:
        raise ValueError("boost factor must be >= 1")

    newshape = np.array(a.shape) * factor

    slices = [slice(0, old, float(old) / new) for old, new in zip(a.shape, newshape)]  # noqa
    coordinates = mgrid[slices]

    # choose the biggest smaller integer index
    indices = coordinates.astype("i")
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
            print("  expanding image from", sz, "to:", [srow, scol])

        new_image = np.empty((srow, scol), dtype=image.dtype)
        new_image[:] = padval

        new_image[0: sz[0], 0: sz[1]] = image[0: sz[0], 0: sz[1]]

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
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as tmpdir:
        tfile = tempfile.mktemp(suffix=".fits")
        tfile = os.path.join(tmpdir, tfile)
        fitsio.write(tfile, im)
        os.system("ds9 %s" % tfile)


def _get_updated_keywords(input_kws, **kws):

    new_kws = input_kws.copy()
    new_kws.update(kws)
    return new_kws


def _prep_plot(
    figax, kw,
    xlim=None, ylim=None, xlog=False, ylog=False,
    xlabel=None, ylabel=None, title=None,
    **subplots_kws
):
    import matplotlib.pyplot as mplt

    file = kw.pop('file', None)

    if file is not None:
        show = kw.pop('show', False)
    else:
        show = kw.pop('show', True)

    if figax is None:
        figax = mplt.subplots(**subplots_kws)
        fig, axs = figax

        if isinstance(axs, np.ndarray):
            axlist = axs.ravel()
        else:
            axlist = [axs]

        for ax in axlist:
            axis_kw = {
                'xlabel': xlabel,
                'ylabel': ylabel,
            }
            if xlim is not None:
                axis_kw['xlim'] = xlim

            if ylim is not None:
                axis_kw['ylim'] = ylim

            ax.set(**axis_kw)

            if title is not None:
                ax.set_title(title)

            if xlog:
                ax.set_xscale('log')

            if ylog:
                ax.set_yscale('log')

    else:
        fig, axs = figax

    return fig, axs, file, show


def _prep_2plot(
    figax, kw,
    xlim=None, ylim=None, xlog=False, ylog=False,
    xlabel=None, ylabel=None, title=None,
):
    import matplotlib.pyplot as mplt

    file = kw.pop('file', None)

    if file is not None:
        show = kw.pop('show', False)
    else:
        show = kw.pop('show', True)

    if figax is None:
        fig, axs = mplt.subplots(ncols=2)

        if xlabel is None:
            xlabel = 'radius [pixels]'

        axis_kw = {
            'xlabel': xlabel,
            'ylabel': ylabel,
            'title': title,
        }
        if xlim is not None:
            axis_kw['xlim'] = xlim

        if ylim is not None:
            axis_kw['ylim'] = ylim

        axs[1].set(**axis_kw)

        if xlog:
            axs[1].set_xscale('log')

        if ylog:
            axs[1].set_yscale('log')

    else:
        fig, axs = figax

    return fig, axs, file, show


def _get_demo_image():
    nx, ny = 32, 32

    cen = (np.array([ny, nx]) - 1)/2
    sigma = 3.0

    x, y = np.mgrid[
        0:ny,
        0:ny,
    ]
    y = y - cen[0]
    x = x - cen[1]

    arg = -0.5 * (x**2 + y**2)/sigma**2
    image = np.exp(arg)

    return image


def demo():
    colorbar = True

    image0 = _get_demo_image()

    image = image0 + np.random.normal(size=image0.shape, scale=0.05)

    # view_profile(image, title='profile')
    #
    # view(image, colorbar=colorbar, title='image view')
    # multiview(image, colorbar=colorbar, title='multiview')
    multiview(image, colorbar=colorbar, profile=True,
              title='multiview profile')

    # view_mosaic([image]*5, colorbar=colorbar, title='mosaic')
    # view_mosaic([image]*5, colorbar=colorbar, title='mosaic combined',
    #             combine=True)
