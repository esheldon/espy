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


def _show_andor_save(plt, file, show, dpi):
    if file is not None:
        plt.write(file, dpi=dpi)

    if show:
        plt.show(dpi=dpi)
