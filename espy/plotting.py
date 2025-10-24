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
    width=3.5,
    legend=None,
    figax=None,

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
        Axis ratio of plot, ysize/xsize, default golden ratio 1.618
    width: float, optional
        Optional reference width, default 3.5
    legend: bool or dicdt
        If True, a legend is created. If a dict, then the keywords
        are sent to the legend call
    figax: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False
    file: str, optional
        Filename to write.
    dpi: int, optional
        dots-per-inch for a bitmap output

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """

    fig, ax, file, show = _prep_plot(
        figax=figax,
        xlim=xlim, ylim=ylim,
        title=title, xlabel=xlabel, ylabel=ylabel,
        xlog=xlog, ylog=ylog,
        aspect=aspect, width=width,
        kw=kw,
    )

    if xerr is not None or yerr is not None:
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kw)
    elif 'ls' in kw or 'linestyle' in kw:
        ax.plot(x, y, **kw)
    else:
        ax.scatter(x, y, **kw)

    _do_legend_maybe(ax=ax, legend=legend)
    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    return fig, ax


def plot_hist(
    x,
    xlabel=None,
    ylabel=None,
    title=None,
    xlim=None,
    ylim=None,
    ylog=False,
    aspect=1.618,  # golden ratio
    width=3.5,
    legend=None,
    figax=None,

    dpi=None,
    **kw
):
    """
    make a plot

    Parameters
    ----------
    x: array or sequences
        Array of x values
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
        Axis ratio of plot, ysize/xsize, default golden ratio 1.618
    width: float, optional
        Optional reference width, default 3.5
    legend: bool or dicdt
        If True, a legend is created. If a dict, then the keywords
        are sent to the legend call
    figax: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False
    file: str, optional
        Filename to write.
    dpi: int, optional
        dots-per-inch for a bitmap output

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """
    import numpy as np

    fig, ax, file, show = _prep_plot(
        figax=figax,
        xlim=xlim, ylim=ylim,
        title=title, xlabel=xlabel, ylabel=ylabel,
        xlog=False, ylog=ylog,
        aspect=aspect, width=width,
        kw=kw,
    )

    binsize = kw.pop('binsize', None)
    if binsize is not None:
        if 'range' in kw:
            xmin, xmax = kw['range']
        else:
            xmin, xmax = [x.min(), x.max()]

        nbin = int((xmax - xmin) / binsize)
        kw['bins'] = np.linspace(xmin, xmax, nbin)

    ax.hist(x, **kw)

    _do_legend_maybe(ax=ax, legend=legend)
    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    return fig, ax


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


def multihist(
    data,
    binfac=0.1,
    nsig=5,
    labels=None,
    ylog=False,
    aspect=None,
    width=6,
    figax=None,
    alpha=1,
    density=False,
    dpi=150,
    **kw,
    # file=None,
    # show=True,
):
    """
    plot a histogram for each dimension of the data

    For many cases, the corner code is a better option

    parameters
    ----------
    data: array
        array with shape [npoints, ndim]
    binfac: float, optional
        The binsize for each dimension will be chosen as binfac*std(dimdata)
        Default 0.1
    nsig: float, optional
        If sent, plot range will be [-nsig*sigma, nsig*sigma].  Otherwise
        the full range is shown
    labels: optional
        A sequence of labels for each dimension
    figax: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    """
    import numpy as np
    import matplotlib.pyplot as mplt
    import esutil as eu

    file = kw.pop('file', None)

    if file is not None:
        show = kw.pop('show', False)
    else:
        show = kw.pop('show', True)

    if len(data.shape) != 2:
        raise ValueError("data should have shape [npoints,ndim]")

    ndim = data.shape[1]

    if labels is not None:
        nl = len(labels)
        assert len(labels) == ndim, "len(labels) = %d != %d" % (nl, ndim)

    ndim = data.shape[1]

    grid = Grid(ndim)

    if aspect is None:
        aspect = grid.ncol / grid.nrow

    if figax is not None:
        fig, axs = figax
    else:
        height = width / aspect
        fig, axs = mplt.subplots(
            nrows=grid.nrow,
            ncols=grid.ncol,
            figsize=(width, height),
            layout='constrained',
        )
        for ax in axs.ravel():
            ax.axis('off')

    for dim in range(ndim):
        ax = axs.ravel()[dim]
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

        ddata = data[:, dim]

        mn, std, ind = eu.stat.sigma_clip(ddata, get_indices=True)

        if nsig is not None:
            low, high = mn - nsig * std, mn + nsig * std
        else:
            low, high = ddata.min(), ddata.max()

        binsize = binfac * std
        nbin = int((high - low) / binsize)
        bins = np.linspace(low, high, nbin)

        ax.axis('on')

        if labels is not None:
            lab = labels[dim]
        else:
            lab = 'dim %d' % dim

        ax.set(xlabel=lab)

        if ylog:
            ax.set_yscale('log')

        ax.hist(ddata, bins=bins, alpha=alpha, density=density)

    # fig.tight_layout()

    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    # if show:
    #     mplt.show()
    # if file is not None:
    #     fig.savefig(file, dpi=dpi)

    return fig, axs


def test_multihist(num=100_000, show=False):
    import numpy as np

    # sigmas = np.array([1, 2, 3])
    # locs = np.array([1, 2, 3])
    # labels = ['a', 'b', 'c']
    #
    # rng = np.random.default_rng()
    #
    # ndim = sigmas.size
    # data = np.zeros((num, ndim))
    #
    # for i in range(ndim):
    #     data[:, i] = rng.normal(scale=sigmas[i], loc=locs[i], size=num)
    #

    means = [1, 2]
    stds = [1, 2]
    corr = 0.8
    covs = [[stds[0]**2, stds[0]*stds[1]*corr],
            [stds[0]*stds[1]*corr, stds[1]**2]]

    data = np.random.multivariate_normal(means, covs, num)
    print(data.shape)
    multihist(data, labels=['a', 'b'], show=show)

    return data


class Grid(object):
    """
    represent plots in a grid.  The grid is chosen
    based on the number of plots

    example
    -------
    grid=Grid(n)

    for row, col in grid:
        axs[row, col].plot(...)

    for i in range(n):
        row, col = grid(i)

        axs[row, col].plot(...)

    for row, col in grid:
        index = grid.get_index(row, col)
        axs[row, col].plot(images[index])
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
            self.nrows, self.ncols = 2, 4
        else:

            sq = int(sqrt(nplot))
            if nplot == sq*sq:
                self.nrows, self.ncols = sq, sq
            elif nplot <= sq*(sq+1):
                self.nrows, self.ncols = sq, sq+1
            else:
                self.nrows, self.ncols = sq+1, sq+1

        self.nrow = self.nrows
        self.ncol = self.ncols
        self.nplot_tot = self.nrows * self.ncols

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
        nplot = 7
        grid = Grid(nplot)

        for i in range(nplot):
            row, col = grid.get_rowcol(nplot, i)
            axs[row, col].plot(...)
        """

        imax = self.nplot_tot-1
        if index > imax:
            raise ValueError("index too large %d > %d" % (index, imax))

        row = index//self.ncol
        col = index % self.ncol

        return row, col

    def __len__(self):
        return self.nplot

    def __iter__(self):
        for i in range(self.nplot):
            yield self(i)
    #     self._current = -1
    #     return self
    #
    # def __next__(self):
    #     self._current += 1
    #     if self._current >= self.nplot:
    #         raise StopIteration
    #     return self(self._current)

    def index(self, row, col):
        return row * self.ncols + col

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
    C=None,  # noqa
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


def plot_ranges(
    list_of_ranges,
    labels,
    xlabel=None,
    ylabel=None,
    title=None,
    xlim=None,
    ylim=None,
    xlog=False,
    ylog=False,
    aspect=1.618,  # golden ratio
    width=3.5,
    figax=None,
    textoff=0.1,
    buff=0.5,

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
        Axis ratio of plot, ysize/xsize, default golden ratio 1.618
    width: float, optional
        Optional reference width, default 3.5
    legend: bool or dicdt
        If True, a legend is created. If a dict, then the keywords
        are sent to the legend call
    figax: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False
    file: str, optional
        Filename to write.
    dpi: int, optional
        dots-per-inch for a bitmap output

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """
    import numpy as np

    fig, ax, file, show = _prep_plot(
        figax=figax,
        xlim=xlim, ylim=ylim,
        title=title, xlabel=xlabel, ylabel=ylabel,
        xlog=xlog, ylog=ylog,
        aspect=aspect, width=width,
        kw=kw,
    )

    nr = len(list_of_ranges)
    nlab = len(labels)
    if nlab != nr:
        raise ValueError(f'got {nlab} labels and {nr} ranges')

    lowest = np.inf
    highest = -np.inf

    biggest = 0.0
    for range_ in list_of_ranges:
        low, high = range_
        size = (high - low)
        if size > biggest:
            biggest = size
        if low < lowest:
            lowest = low
        if high > highest:
            highest = high

    # text_x = low - biggest * 0.1

    ystep = 1
    buff = 0.5
    ylim = [-(nr - 1) * ystep - buff, buff]
    fac = 0.2
    xlim = [lowest - fac * biggest, highest + fac * biggest]
    ax.set(xlim=xlim, ylim=ylim, yticklabels=[])
    ax.tick_params(left=False, right=False)
    ax.tick_params(which='minor', left=False, right=False)

    ystart = 0

    for i, range_ in enumerate(list_of_ranges):
        low, high = range_
        y = ystart - i * ystep

        mid = (high + low) * 0.5

        ax.plot(
            [low, high],
            [y, y],
            lw=10,
        )

        ax.text(
            mid,
            y + textoff,
            labels[i],
            horizontalalignment='center',
            fontsize='x-large',
        )

    fig.tight_layout()
    _show_andor_save(fig=fig, file=file, show=show, dpi=dpi)
    return fig, ax


def _prep_plot(
    figax, xlim, ylim, xlabel, ylabel, title, xlog, ylog, aspect, width, kw,
):
    import matplotlib.pyplot as plt

    file = kw.pop('file', None)

    if file is not None:
        show = kw.pop('show', False)
    else:
        show = kw.pop('show', True)

    if figax is None:
        height = width / aspect
        figax = plt.subplots(figsize=(width, height))
        fig, ax = figax
        axis_kw = {
            'xlabel': xlabel,
            'ylabel': ylabel,
            'title': title,
        }
        if xlim is not None:
            axis_kw['xlim'] = xlim

        if ylim is not None:
            axis_kw['ylim'] = ylim

        ax.set(**axis_kw)

        if xlog:
            ax.set_xscale('log')

        if ylog:
            ax.set_yscale('log')

    else:
        fig, ax = figax

    return fig, ax, file, show


def _show_andor_save(fig, file, show, dpi):
    import matplotlib.pyplot as plt

    if file is not None:
        fig.savefig(file, dpi=dpi)

    if show:
        plt.show()


def _do_legend_maybe(ax, legend):
    if legend is not None:
        if legend is True:
            ax.legend()
        else:
            ax.legend(**legend)


def zenburn():
    from cycler import cycler
    colors = [
        '#8cd0d3', '#cc9393', '#7f9f7f', '#f0dfaf', '#dcdccc',
        '#4A7274', '#466F46', '#A55D5D', '#A35E2E', 'white',
    ]

    colorcyc = cycler('color', colors)

    return {
        'axes.axisbelow': True,
        'axes.edgecolor': 'CFCFCF',
        'axes.facecolor': '3F3F3F',

        'axes.grid': True,

        'axes.labelcolor': 'white',

        'axes.labelsize': 'large',

        'axes.prop_cycle': colorcyc,
        'axes.titlesize': 'x-large',

        'figure.edgecolor': 'FFFFEF',
        'figure.facecolor': '3F3F3F',

        'grid.color': '6F6F6F',
        'grid.linestyle': '-',

        'legend.facecolor': '9F9F9F',
        'legend.fancybox': True,

        'lines.color': '3F3F3F',

        'patch.antialiased': True,
        'patch.facecolor': '8cd0d3',

        'savefig.edgecolor': 'CFCFCF',
        'savefig.facecolor': '3F3F3F',

        'text.color': 'white',

        'xtick.color': 'CFCFCF',
        'ytick.color': 'CFCFCF',
    }
