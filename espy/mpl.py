"""
TODO:
    - non-symmetric error bars
"""
import numpy as np
import tempfile
import subprocess
import matplotlib.pyplot as plt


class Points(object):
    def __init__(
        self,
        x, y, xerr=None, yerr=None,
        marker='o',
        size=None,
        linestyle=None,
        linewidth=None,
        color=None,
        edgecolor=None,
        edgewidth=None,
        alpha=None,
        capsize=2,
    ):

        self._set_points(x, y, xerr, yerr)
        self.marker = marker
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.color = color
        self.edgecolor = edgecolor
        self.edgewidth = edgewidth
        self.size = size
        self.alpha = alpha

        self.capsize = capsize

    def _set_points(self, x, y, xerr, yerr):
        self.x = _make_array_maybe(x)
        self.y = _make_array_maybe(y)
        self.xerr = _make_array_maybe(xerr)
        self.yerr = _make_array_maybe(yerr)

        if self.x.size != self.y.size:
            raise ValueError(
                "x and y must be same "
                "size, got %d and %d" % (self.x.size, self.y.size)
            )

        if self.xerr is not None and self.x.size != self.xerr.size:
            raise ValueError(
                "x and xerr must be same "
                "size, got %d and %d" % (self.x.size, self.xerr.size)
            )
        if self.yerr is not None and self.y.size != self.yerr.size:
            raise ValueError(
                "y and yerr must be same "
                "size, got %d and %d" % (self.y.size, self.yerr.size)
            )

    def _add_to_axes(self, ax):
        linestyle = 'none' if self.linestyle is None else self.linestyle
        ax.errorbar(
            self.x, self.y,
            xerr=self.xerr,
            yerr=self.yerr,
            marker=self.marker,
            markersize=self.size,
            linestyle=linestyle,
            linewidth=self.linewidth,
            color=self.color,
            capsize=self.capsize,
            ecolor=self.color,
            elinewidth=self.linewidth,
            alpha=self.alpha,
            # zorder=3,
            # barsabove=True,
        )


class Curve(Points):
    """
    Same as Points but defaults to marker None
    """
    def __init__(
        self, x, y,
        marker=None,
        linestyle='-',
        linewidth=None,
        color=None,
        alpha=None,
    ):

        super().__init__(
            x, y,
            marker=marker,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            alpha=alpha,
        )


class HLine(object):
    def __init__(self, y=0,
                 xmin=0, xmax=1,
                 linestyle='-', linewidth=None,
                 color=None,
                 alpha=None):
        self.y = y
        self.xmin = xmin
        self.xmax = xmax
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.color = color
        self.alpha = alpha

    def _add_to_axes(self, ax):
        linestyle = 'none' if self.linestyle is None else self.linestyle
        ax.axhline(
            y=self.y,
            xmin=self.xmin,
            xmax=self.xmax,
            linestyle=linestyle,
            linewidth=self.linewidth,
            color=self.color,
            alpha=self.alpha,
        )


def _make_array_maybe(data):
    if data is None:
        return None
    else:
        return np.array(data, ndmin=1, copy=False)


class Plot(object):
    def __init__(self, xlabel=None, ylabel=None):

        self._xlabel = xlabel
        self._ylabel = ylabel
        self.objlist = []
        self._reset_fig()

    def add(self, *args):
        self.objlist += args
        self._reset_fig()

    def write(self, fname, dpi=None):
        fig, ax = self.get_fig()
        self.fig.savefig(fname, dpi=dpi)

    def show(self, dpi=None):
        with tempfile.TemporaryDirectory() as dir:
            fname = tempfile.mktemp(dir=dir, suffix='.png')
            self.write(fname, dpi=dpi)

            command = ['feh', fname]
            subprocess.check_call(command)

    def get_fig(self):
        if self.fig is None or self.ax is None:
            self._render_fig()
        return self.fig, self.ax

    def _render_fig(self):
        self.fig, self.ax = plt.subplots()

        if self._xlabel is not None:
            self.ax.set_xlabel(self._xlabel)
        if self._ylabel is not None:
            self.ax.set_ylabel(self._ylabel)


        for obj in self.objlist:
            obj._add_to_axes(self.ax)

    """
    def _add_points(self, pts):

        linestyle = 'none' if pts.linestyle is None else pts.linestyle
        self.ax.errorbar(
            pts.x, pts.y,
            xerr=pts.xerr,
            yerr=pts.yerr,
            marker=pts.marker,
            markersize=pts.size,
            linestyle=linestyle,
            linewidth=pts.linewidth,
            color=pts.color,
            capsize=pts.capsize,
            ecolor=pts.color,
            elinewidth=pts.linewidth,
            alpha=pts.alpha,
        )
        self.ax.scatter(
            pts.x, pts.y,
            s=pts.size,
            c=pts.color,
            edgecolors=pts.edgecolor,
            linewidths=pts.edgewidth,
            marker=pts.type,
            alpha=pts.alpha,
        )

        if pts.xerr is not None or pts.yerr is not None:
            self.ax.errorbar(
                pts.x, pts.y,
                xerr=pts.xerr,
                yerr=pts.yerr,
                fmt='none',  # no points just error bars
                ecolor=pts.color,
                alpha=pts.alpha,
            )
        """

    def _reset_fig(self):
        self.fig = None
        self.ax = None


def plot(x, y, xerr=None, yerr=None,
         marker='o',
         linestyle=None,
         linewidth=None,
         color=None,
         edgecolor=None,
         edgewidth=None,
         size=None,
         alpha=None,
         capsize=2,
         xlabel=None, ylabel=None,
         plt=None, show=False, dpi=None):

    if plt is None:
        plt = Plot(
            xlabel=xlabel,
            ylabel=ylabel,
        )
    assert isinstance(plt, Plot)

    pts = Points(
        x, y, xerr=xerr, yerr=yerr,
        marker=marker,
        size=size,
        linestyle=linestyle,
        linewidth=linewidth,
        color=color,
        edgecolor=edgecolor,
        alpha=alpha,
        capsize=capsize,
    )
    plt.add(pts)

    if show:
        plt.show(dpi=dpi)

    return plt
