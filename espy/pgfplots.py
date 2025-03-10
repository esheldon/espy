def plot(
    x,
    y,
    xerr=None,
    yerr=None,

    # axis keywords
    minor_tick_num=3,
    scaled_ticks=True,
    xlabel=None,
    ylabel=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xlog=False,
    ylog=False,
    cycle_list_name='exotic',
    width=None,
    height=None,
    legend=None,
    grid=None,
    grid_style=None,
    enlarge_limits=True,
    axis_equal_image=False,
    colormap=None,

    # for saving png or show
    dpi=150,

    plt=None,
    file=None,
    show=None,

    # plot keywords
    label=None,
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

    minor_tick_num: int
        Default 3
    scaled_ticks: bool
        Default True
    xlabel: str, optional
        Label for x axis
    ylabel: str, optional
        Label for y axis
    title: str, optional
        Title string for plot
    xmin: float
        Optional lower limit
    xmax: float
        Optional lower limit
    ymin: float
        Optional lower limit
    ymax: float
        Optional lower limit
    xlog: bool, optional
        If True, use log x axis
    ylog: bool, optional
        If True, use log y axis
    cyhcle_list_name: str
        Default 'exotic'
    width: float, optional
        Optional width
    height: float
        Optional height
    grid: bool
        If set to True, add grid
    grid_style: dict
        Style for grid
    enlarge_limits: bool
        Defaul True
    axis_equal_image: bool
        Default False
    colormap: str
        Name for colormap, default None
    dpi: int, optional
        dots-per-inch for a bitmap output, default 150
    plt: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    file: str, optional
        Filename to write.
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """
    import pypgf

    if plt is None:
        plt = pypgf.Plot(
            minor_tick_num=minor_tick_num,
            scaled_ticks=scaled_ticks,
            xlabel=xlabel,
            ylabel=xlabel,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            xlog=xlog,
            ylog=ylog,
            cycle_list_name=cycle_list_name,
            width=width,
            height=height,
            legend=legend,
            grid=legend,
            grid_style=legend,
            enlarge_limits=enlarge_limits,
            axis_equal_image=axis_equal_image,
            colormap=colormap,
        )

    plt.plot(
        x=x,
        y=y,
        xerr=xerr,
        yerr=yerr,
        label=label,
        **kw
    )

    if file is None and show is None:
        show = True

    _show_andor_save(plt=plt, file=file, show=show, dpi=dpi)
    return plt


def scatter(
    x,
    y,

    # axis keywords
    minor_tick_num=3,
    scaled_ticks=True,
    xlabel=None,
    ylabel=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xlog=False,
    ylog=False,
    cycle_list_name='exotic',
    width=None,
    height=None,
    legend=None,
    grid=None,
    grid_style=None,
    enlarge_limits=True,
    axis_equal_image=False,
    colormap=None,

    # for saving png or show
    dpi=150,

    plt=None,
    file=None,
    show=None,

    # plot keywords
    label=None,
    only_markers=True,
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

    minor_tick_num: int
        Default 3
    scaled_ticks: bool
        Default True
    xlabel: str, optional
        Label for x axis
    ylabel: str, optional
        Label for y axis
    title: str, optional
        Title string for plot
    xmin: float
        Optional lower limit
    xmax: float
        Optional lower limit
    ymin: float
        Optional lower limit
    ymax: float
        Optional lower limit
    xlog: bool, optional
        If True, use log x axis
    ylog: bool, optional
        If True, use log y axis
    cyhcle_list_name: str
        Default 'exotic'
    width: float, optional
        Optional width
    height: float
        Optional height
    grid: bool
        If set to True, add grid
    grid_style: dict
        Style for grid
    enlarge_limits: bool
        Defaul True
    axis_equal_image: bool
        Default False
    colormap: str
        Name for colormap, default None
    dpi: int, optional
        dots-per-inch for a bitmap output, default 150
    plt: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    file: str, optional
        Filename to write.
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """
    import pypgf

    if plt is None:
        plt = pypgf.Plot(
            minor_tick_num=minor_tick_num,
            scaled_ticks=scaled_ticks,
            xlabel=xlabel,
            ylabel=xlabel,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            xlog=xlog,
            ylog=ylog,
            cycle_list_name=cycle_list_name,
            width=width,
            height=height,
            legend=legend,
            grid=legend,
            grid_style=legend,
            enlarge_limits=enlarge_limits,
            axis_equal_image=axis_equal_image,
            colormap=colormap,
        )

    plt.scatter(
        x=x,
        y=y,
        label=label,
        only_markers=only_markers,
        **kw
    )

    if file is None and show is None:
        show = True

    _show_andor_save(plt=plt, file=file, show=show, dpi=dpi)
    return plt


def plot_hist(
    data,

    # axis keywords
    minor_tick_num=3,
    scaled_ticks=True,
    xlabel=None,
    ylabel=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xlog=False,
    ylog=False,
    cycle_list_name='exotic',
    width=None,
    height=None,
    legend=None,
    grid=None,
    grid_style=None,
    enlarge_limits=True,
    axis_equal_image=False,
    colormap=None,

    # for saving png or show
    dpi=150,

    plt=None,
    file=None,
    show=None,

    # plot keywords
    label=None,
    bins=10,
    range=None,
    density=None,
    weights=None,
    no_markers=True,
    fill=True,
    draw=True,

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

    minor_tick_num: int
        Default 3
    scaled_ticks: bool
        Default True
    xlabel: str, optional
        Label for x axis
    ylabel: str, optional
        Label for y axis
    title: str, optional
        Title string for plot
    xmin: float
        Optional lower limit
    xmax: float
        Optional lower limit
    ymin: float
        Optional lower limit
    ymax: float
        Optional lower limit
    xlog: bool, optional
        If True, use log x axis
    ylog: bool, optional
        If True, use log y axis
    cyhcle_list_name: str
        Default 'exotic'
    width: float, optional
        Optional width
    height: float
        Optional height
    grid: bool
        If set to True, add grid
    grid_style: dict
        Style for grid
    enlarge_limits: bool
        Defaul True
    axis_equal_image: bool
        Default False
    colormap: str
        Name for colormap, default None
    dpi: int, optional
        dots-per-inch for a bitmap output, default 150
    plt: (fig, ax)
        If sent, a new figure and axis is not created, the input one is
        reused
    file: str, optional
        Filename to write.
    show: bool, optional
        If True, show the plot on the screen.  If the file= is
        not sent, this defaults to True.  If file= is sent
        this defaults to False

    **kw extra keywords for plot call

    Returns
    -------
    fig, ax
    """
    import pypgf

    if plt is None:
        plt = pypgf.Plot(
            minor_tick_num=minor_tick_num,
            scaled_ticks=scaled_ticks,
            xlabel=xlabel,
            ylabel=xlabel,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
            xlog=xlog,
            ylog=ylog,
            cycle_list_name=cycle_list_name,
            width=width,
            height=height,
            legend=legend,
            grid=legend,
            grid_style=legend,
            enlarge_limits=enlarge_limits,
            axis_equal_image=axis_equal_image,
            colormap=colormap,
        )

    plt.hist(
        data,
        label=label,
        bins=bins,
        range=range,
        density=density,
        weights=weights,
        no_markers=no_markers,
        fill=fill,
        draw=draw,
        **kw
    )

    if file is None and show is None:
        show = True

    _show_andor_save(plt=plt, file=file, show=show, dpi=dpi)
    return plt


def _show_andor_save(plt, file, show, dpi):
    if file is not None:
        plt.write(file, dpi=dpi)

    if show:
        plt.show(dpi=dpi)
