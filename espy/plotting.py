GOLDEN_RATIO = 1.61803398875
GOLDEN_ARATIO = 1.0/GOLDEN_RATIO

def plot(
    x, y, xerr=None, yerr=None,
    xlabel=None,
    ylabel=None,
    title=None,
    xlim=None,
    ylim=None,
    xlog=False,
    ylog=False,
    aspect=1.618,  # golden ratio
    legend=None,
    figax=None,

    figsize=None,
    dpi=None,
    **kw
):
    """
    make a plot

    Parameters
    ----------
    x: array or sequences
        Array of x values
    y: array or sequences
        Array of y values
    xerr: array or sequence, optional
        Optional array of x errors
    yerr: array or sequence, optional
        Optional array of y errors

    xlabel: str, optional
        Label for x axis
    ylabel: str, optional
        Label for y axis
    title: str, optional
        Title string for plot
    xlim: 2-element sequence, optional
        Optional limits for the x axis
    ylim: 2-element sequence, optional
        Optional limits for the y axis
    xlog: bool, optional
        If True, use log x axis
    ylog: bool, optional
        If True, use log y axis
    aratio: float, optional
        Axis ratio of plot, ysize/xsize
    legend: bool or Legend instance
        If True, a legend is created. You can also send a Legend() instance.
        If None or False, no legend is created
    plt: Plot instance, optional
        If sent, a new Plot is not created, the input one
        is reused

    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False
    file: str, optional
        Filename to write.

    Keywords for the Plot/matplotlib Figure.  See docs for
        the matplotlib Figure class

        figsize dpi facecolor edgecolor linewidth
        frameon subplotpars tight_layout constrained_layout

    Keywords for plot or errorbar, depending if xerr/yerr
    are sent.  See docs for matplotlib axes.plot and errorbar
    commands

    Returns
    -------
    Plot instance
    """
    file = kw.pop('file', None)
    if file is not None:
        show = kw.pop('show', False)
    else:
        show = kw.pop('show', config['show'])

    if figax is None:
        axis_kw = {
            'xlabel': xlabel,
            'ylabel': ylabel,
            'title': title,
        }
        if xlim is not None:
            axis_kw['xlim'] = xlim
        if ylim is not None:
            axis_kw['ylim'] = ylim
        if xlog:
            axis_kw['xscale'] = 'log'
        if ylog:
            axis_kw['yscale'] = 'log'

        plt = Plot(
            aratio=aratio,
            legend=legend,
            figsize=figsize,
            dpi=dpi,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            frameon=frameon,
            subplotpars=subplotpars,
            tight_layout=tight_layout,
            constrained_layout=constrained_layout,  # default to rc
            **axis_kw
        )

    if xerr is not None or yerr is not None:
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)
    else:
        plt.plot(x, y, **kw)

    if file is not None:
        plt.savefig(file)

    if show:
        plt.show()

    return plt



def plot_residuals(
    *, x, y, model, yerr=None, frac=0.2, pad=0,
    data_kw={}, model_kw={},
    resid_axis_kw={},
    no_resid_xticklabels=False,
    no_resid_yticklabels=False,
    plot_kw={},
    show=False,
    plt=None,
):
    """
    plot data and model with a residuals plot below

    Parameters
    ----------
    x: array
        Array of x values
    y: array
        Array of y values
    model: array
        Array of model values
    yerr: array, optional
        Optional array of errors for y
    frac: float, optional
        Fraction of figure taken up by residuals plot. Default 0.2
    pad: float, optional
        Optional padding, default 0
    data_kw: dict, optional
        Optional dict of keywords for the data points, e.g. color
        and label, etc.
    model_kw: dict, optional
        Optional dict of keywords for the model points, e.g. color
        and label, etc.  Note the model is always plotted as a curve
    plot_kw: dict, optional
        Keywords to use when constructing the plot object.
    resid_axis_kw: dict, optional
        Optional dict of keywords for the residuals axis set() method
    no_resid_xticklabels: bool, optional
        If True, don't add xticklabels for residuals plot
    no_resid_yticklabels: bool, optional
        If True, don't add yticklabels for residuals plot
    show: bool, optional
        If True, show the plot on the screen
    plt: hickory Plot object, optional
        Used for plotting if sent rather than creating a new one
    """
    import numpy as np
    import hickory
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if plt is None:
        plt = hickory.Plot(**plot_kw)

    divider = make_axes_locatable(plt)
    size_string = '%d%%' % (frac*100)
    ax2 = divider.append_axes("bottom", size=size_string, pad=pad)
    plt.figure.add_axes(ax2, label='%d' % np.random.randint(0, 2**15))

    residuals = model - y
    residuals_err = yerr

    ax2.set(**resid_axis_kw)
    ax2.set_xscale('log')
    ax2.axhline(0, color='black')

    if yerr is not None:
        plt.errorbar(x, y, yerr, **data_kw)
        ax2.errorbar(x, residuals, residuals_err, **data_kw)
    else:
        plt.plot(x, y, **data_kw)
        ax2.plot(x, residuals, **data_kw)

    plt.curve(x, model, **model_kw)

    if no_resid_xticklabels:
        ax2.set_xticklabels([])

    if no_resid_yticklabels:
        ax2.set_yticklabels([])

    if show:
        plt.show()

    return plt


def plot_hist2d(x, y, **kw):
    """
    create a 2-d histogram of the data and view it


    parameters
    ----------
    x: array
        x data
    y: array
        y data

    log: bool
        If true, plot the log of the histogram, default
        True
    xlabel: string
        label for x axis
    ylabel: string
        label for y axis

    **kw:
        keywords for the histogram2d routine

    dependencies
    ------------
    esutil
    images
    """
    import numpy as np
    import images
    import esutil as eu

    dolog = kw.pop('log', True)

    kw['more'] = True
    hdict = eu.stat.histogram2d(x, y, **kw)

    hist = hdict['hist']
    if dolog:
        hist = np.log10(hist.clip(min=0.1))

    kw['transpose'] = False
    kw['ranges'] = hdict['ranges']

    # we want exact ranges here, to avoid confusing areas
    # of black where there is no histogram
    kw['xrange'] = [hdict['ranges'][0][0], hdict['ranges'][1][0]]
    kw['yrange'] = [hdict['ranges'][0][1], hdict['ranges'][1][1]]

    plt = images.view(hist, **kw)

    return plt


def multihist(data, binfac=0.1, **kw):
    """
    plot a histogram for each dimension of the data

    parameters
    ----------
    data: array
        array with shape [npoints, ndim]
    binfac: float
        The binsize for each dimension will be chosen as binfac*std(dimdata)
    labels: optional
        A sequence of labels for each dimension
    """
    import biggles
    import esutil as eu

    if len(data.shape) != 2:
        raise ValueError("data should have shape [npoints,ndim]")

    ndim = data.shape[1]

    labels = kw.pop('labels', None)
    if labels is not None:
        nl = len(labels)
        assert len(labels) == ndim, "len(labels) = %d != %d" % (nl, ndim)

    grid = Grid(ndim)

    tab = kw.pop('plt', None)
    if tab is not None:
        add_to_existing_plots = True
    else:
        add_to_existing_plots = False
        tab = biggles.Table(grid.nrow, grid.ncol)
        tab.aspect_ratio = kw.pop('aspect_ratio', None)

        if tab.cols != grid.ncol or tab.rows != grid.nrow:
            m = "input table has wrong dims.  Expected %s got %s"
            tup = ((grid.nrow, grid.ncol), (tab.rows, tab.cols))
            raise ValueError(m % tup)

    for dim in range(ndim):

        ddata = data[:, dim]

        mn, std = eu.stat.sigma_clip(ddata)

        binsize = binfac*std

        hc = biggles.make_histc(
            ddata,
            binsize=binsize,
            **kw
        )

        row, col = grid(dim)
        if add_to_existing_plots:
            plt = tab[row, col]
            plt.add(hc)
        else:
            plt = biggles.FramedPlot(aspect_ratio=1, **kw)
            plt.add(hc)
            tab[row, col] = plt

        if labels is not None:
            lab = labels[dim]
        else:
            lab = 'dim %d' % dim
        tab[row, col].xlabel = lab

    return tab


class Grid(object):
    """
    represent plots in a grid.  The grid is chosen
    based on the number of plots

    example
    -------
    grid=Grid(n)

    for i in range(n):
        row,col = grid(i)

        # equivalently grid.get_rowcol(i)

        plot_table[row,col] = plot(...)
    """
    def __init__(self, nplot):
        self.set_grid(nplot)

    def set_grid(self, nplot):
        """
        set the grid given the number of plots
        """
        from math import sqrt

        self.nplot = nplot

        # first check some special cases
        if nplot == 8:
            self.nrow, self.ncol = 2, 4
        else:

            sq = int(sqrt(nplot))
            if nplot == sq*sq:
                self.nrow, self.ncol = sq, sq
            elif nplot <= sq*(sq+1):
                self.nrow, self.ncol = sq, sq+1
            else:
                self.nrow, self.ncol = sq+1, sq+1

        self.nplot_tot = self.nrow*self.ncol

    def get_rowcol(self, index):
        """
        get the grid position given the number of plots

        move along columns first

        parameters
        ----------
        index: int
            Index in the grid

        example
        -------
        nplot=7
        grid=Grid(nplot)
        arr=biggles.FramedArray(grid.nrow, grid.ncol)

        for i in range(nplot):
            row,col=grid.get_rowcol(nplot, i)
            arr[row,col].add( ... )
        """

        imax = self.nplot_tot-1
        if index > imax:
            raise ValueError("index too large %d > %d" % (index, imax))

        row = index//self.ncol
        col = index % self.ncol

        return row, col

    def __call__(self, index):
        return self.get_rowcol(index)


def get_corner_data(*args):
    """
    get corner data array from inputs
    """
    import numpy as np

    ndim = len(args)
    if ndim == 0:
        raise ValueError('no args sent')

    npts = len(args[0])

    data = np.zeros((npts, ndim))

    for dim, arg in enumerate(args):
        assert len(arg) == npts, 'all arrays must be same size'

        data[:, dim] = arg

    return data


def color_color(
    gmr, rmi, imz,
    marker='dot',
    show=True, file=None,
    dpi=100,
    figsize=(11, 5),
    plt=None,
    **kw
):
    """
    make a color color plot

    Parameters
    ----------
    gmr: array
        Array of g-r
    rmi: array
        Array of r-i
    imz: array
        Array of i-z
    show: bool
        If True, show a plot on the screen, default True unless file= is
        sent
    file: string
        Filename to write for figure
    dpi: int
        dpi for non-pdf/eps type figures
    figsize: sequence
        Default (11, 5)
    **kw: keywords
        Extra keywords for the plot markers

    Returns
    -------
    the hickory plot object
    """
    import hickory

    if plt is None:
        plt = hickory.Table(
            nrows=1, ncols=2,
            figsize=figsize,
        )

    plt[0].set(
        xlim=(-1, 3),
        ylim=(-1, 3),
        xlabel=r'$g - r$',
        ylabel=r'$r - i$',
    )
    plt[0].plot(gmr, rmi, marker=marker, **kw)
    plt[1].set(
        xlim=(-1, 3),
        ylim=(-1, 2),
        xlabel=r'$r - i$',
        ylabel=r'$i - z$',
    )
    plt[1].plot(rmi, imz, marker=marker, **kw)

    if file is not None:
        plt.savefig(file, **kw)
    elif show:
        plt.show()

    return plt


def color_color_hexbin(
    gmr, rmi, imz,
    weights=None,
    C=None,
    nbin=100,
    show=True, file=None,
    dpi=100,
    title=None,
    figsize=(12, 5),
    plt=None,
    **kw
):
    """
    make a color color hexbin plot

    Parameters
    ----------
    gmr: array
        Array of g-r
    rmi: array
        Array of r-i
    imz: array
        Array of i-z
    weights or C:
        Weights to accumulate in each bin.
    show: bool
        If True, show a plot on the screen, default True unless file= is
        sent
    file: string
        Filename to write for figure
    dpi: int
        dpi for non-pdf/eps type figures
    figsize: sequence
        Default (11, 5)
    **kw: keywords
        Extra keywords for the plot markers

    Returns
    -------
    the hickory plot object
    """
    import numpy as np
    from esutil.numpy_util import between
    import hickory

    if weights is not None:
        kw['C'] = weights
    elif C is not None:
        kw['C'] = C
    else:
        kw['C'] = gmr*0 + 1

    kw['reduce_C_function'] = np.sum

    if plt is None:
        plt = hickory.Table(
            nrows=1, ncols=2,
            figsize=figsize,
        )

    if title is not None:
        plt.suptitle(title, fontsize=16)

    gmr_min, gmr_max = -1, 3
    rmi_min, rmi_max = -1, 3
    imz_min, imz_max = -1, 2

    w, = np.where(
        between(gmr, -1, 3) &
        between(rmi, -1, 3) &
        between(imz, -1, 2)
    )

    plt[0].set(
        xlim=(gmr_min, gmr_max),
        ylim=(rmi_min, rmi_max),
        xlabel=r'$g - r$',
        ylabel=r'$r - i$',
    )

    hb0 = plt[0].hexbin(gmr[w], rmi[w], gridsize=nbin, **kw)

    plt[1].set(
        xlim=(rmi_min, rmi_max),
        ylim=(imz_min, imz_max),
        xlabel=r'$r - i$',
        ylabel=r'$i - z$',
    )
    hb1 = plt[1].hexbin(rmi[w], imz[w], gridsize=nbin, **kw)

    plt.colorbar(hb0, ax=plt.axes[0])
    plt.colorbar(hb1, ax=plt.axes[1])

    if file is not None:
        plt.savefig(file, **kw)
    elif show:
        plt.show()

    return plt


def whiskers(
    *,
    x, y,
    u=None, v=None,
    g1=None, g2=None,
    color='black',
    scale=1.0,
    linewidth=0.5,
    linestyle='-',
    show=False,
    plt=None,
    **plotting_keywords
):
    """
    Name:
        mwhiskers
    Calling Sequence:
        whiskers(plt, x, y, u, v, scale=1, **plotting_keywords)
    Plotting Context:
        matplotlib.  Do make whiskers using biggles use the bwhiskers function

    Purpose:

        Using matplotlib, draw lines centered a the input x,y positions, with
        length
            sqrt(u**2 + v**2)
        and angle
            arctan(v,u)

    plt could be an axes instance
        ax = pyplot.subplot(1,2,1)
    or could it self be pyplot or pylab

    """
    import numpy as np
    import esutil as eu

    if plt is None:
        import hickory
        plt = hickory.Plot()

    x = np.array(x, copy=False, ndmin=1)
    y = np.array(y, copy=False, ndmin=1)

    if g1 is not None and g2 is not None:
        u, v = eu.plotting.polar2whisker(g1, g2)
    elif u is not None and v is not None:
        u = np.array(u, copy=False, ndmin=1)
        v = np.array(v, copy=False, ndmin=1)
    else:
        raise ValueError('send u, v or g1, g2')

    if x.size != y.size or x.size != u.size or x.size != v.size:
        raise ValueError(
            "Sizes don't match: %s %s %s %s\n" % (x.size, y.size, u.size, v.size)  # noqa
        )

    for i in range(x.size):
        # create the line to draw.
        xvals = x[i] + np.array([-u[i]/2.0, u[i]/2.0], dtype='f4')*scale
        yvals = y[i] + np.array([-v[i]/2.0, v[i]/2.0], dtype='f4')*scale
        # print(xvals, yvals)

        plt.curve(
            xvals, yvals,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            **plotting_keywords
        )

    if show:
        plt.show()

    return plt
