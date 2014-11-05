from __future__ import print_function
import numpy

LABELS={}
LABELS['rproj'] = r'$R$ [$h^{-1}$ Mpc]'
LABELS['dsig'] = r'$\Delta\Sigma ~[M_{sun} pc^{-2}]$'
LABELS['osig'] = r'$\Delta\Sigma_\times ~ [M_{sun} pc^{-2}]$'

DEFAULT_MINVAL=1.0e-3


def plot_dsig(r, dsig, dsigerr, **kw):
    """
    plot delta sigma, potentially on a grid if more than one bin

    parameters
    ----------
    r: array
        array with shape [nbin,nrad], e.g from a bin struct
    dsig: array
        array with shape [nbin,nrad], e.g from a bin struct
    dsigerr: array
        array with shape [nbin,nrad], e.g from a bin struct

    kw: keywords
        extra keywords
    """
    import biggles


    nbin = r.shape[0]

    grid, arr = get_framed_array(nbin)

    visible=kw.get('visible',True)
    color=kw.get('color','black')
    type=kw.get('type','filled circle')
    is_ortho=kw.get('is_ortho',False)
    xlog=kw.get('xlog',True)
    ylog=kw.get('ylog',True)
    labels=kw.get('labels',None)

    ylabel=kw.get('ylabel',None)

    aspect_ratio=kw.get('aspect_ratio',None)
    if aspect_ratio is None:
        aspect_ratio=float(grid.nrow)/grid.ncol
    arr.aspect_ratio = aspect_ratio 

    arr.uniform_limits=True
    arr.xlog, arr.ylog = xlog, ylog

    arr.xlabel = LABELS['rproj']
    if ylabel is None:
        if is_ortho:
            ylabel = LABELS['osig']
        else:
            ylabel = LABELS['dsig']
    arr.ylabel=ylabel

    for i in xrange(nbin):
        row, col = grid.get_rowcol(i)

        ri=r[i,:]
        dsigi=dsig[i,:]
        dsigerri=dsigerr[i,:]

        plt=arr[row, col]
        xrng, yrng = _add_dsig_to_plot(plt,
                                       ri,
                                       dsigi,
                                       dsigerri,
                                       **kw)
        plt.xrange=xrng
        plt.yrange=yrng

        if labels is not None:
            lab=biggles.PlotLabel(0.9, 0.9, labels[i], halign='right')
            plt.add(lab)

    if visible:
        arr.show()

    return arr

def plot_dsig_one(r, dsig, dsigerr, **kw):
    """
    plot delta sigma

    useful if adding to an existing plot

    parameters
    ----------
    r: array
        radius
    dsig: array
        delta sigma
    dsigerr: array
        error on delta sigma
    """

    from biggles import FramedPlot

    visible=kw.get('visible',True)
    xlog=kw.get('xlog',True)
    ylog=kw.get('ylog',True)
    aspect_ratio=kw.get('aspect_ratio',1)

    is_ortho=kw.get('is_ortho',False)

    plt=kw.get('plt',None)

    if plt is None:
        plt=FramedPlot()
        plt.aspect_ratio=aspect_ratio
        plt.xlog=xlog
        plt.ylog=ylog

        plt.xlabel = LABELS['rproj']
        if is_ortho:
            plt.ylabel = LABELS['osig']
        else:
            plt.ylabel = LABELS['dsig']

    xrng, yrng = _add_dsig_to_plot(plt, r, dsig, dsigerr, **kw)
    plt.xrange=xrng
    plt.yrange=yrng

    if visible:
        plt.show()

    return plt


def _add_dsig_to_plot(plt, r, dsig, dsigerr, **kw):
    from biggles import Points, Curve, SymmetricErrorBarsY
    import esutil as eu

    xlog=kw.get('xlog',True)
    ylog=kw.get('ylog',True)
    minval=kw.get('minval',DEFAULT_MINVAL)
    color=kw.get('color','black')
    type=kw.get('type','filled circle')

    lineval=kw.get('lineval',0.0)

    if ylog:
        xrng,yrng=add_to_log_plot(plt, r, dsig, dsigerr, 
                                  color=color, 
                                  type=type,
                                  minval=minval)
    else:
        opts=Curve(r, dsig*0 + lineval)
        plt.add(opts)

        pts=Points(r, dsig, type=type, color=color)

        plt.add(pts)
        epts=SymmetricErrorBarsY(r, dsig, dsigerr, color=color)
        plt.add(epts)

        yrng=kw.get('yrange',None)
        xrng=kw.get('xrange',None)

        if xrng is None:
            if xlog:
                xrng=eu.plotting.get_log_plot_range(r)

        if yrng is None:
            nrad=dsig.ravel().size
            ytot=numpy.zeros(2*nrad)
            ytot[0:nrad] = dsig-dsigerr
            ytot[nrad:]  = dsig+dsigerr

            wt=1.0/dsigerr**2
            wm,we,wstd=eu.stat.wmom(dsig, wt, sdev=True)

            std=ytot.std()
            yrng=[wm-3.0*std, wm+3.0*std]
            '''

            if yrng[0] < 0:
                yrng[0] *= 1.1
            else:
                yrng[0] *= 0.9

            if yrng[1] < 0:
                yrng[1] *= 0.9
            else:
                yrng[1] *= 1.1
            '''

    return xrng, yrng

def add_to_log_plot(plt, x, y, yerr, 
                    color='black', type='filled circle',
                    minval=DEFAULT_MINVAL):
    import esutil as eu
    from biggles import Points, ErrorBarsY
    w, = numpy.where(y > minval)
    if w.size > 0:
        p = Points(x[w], y[w],type=type,color=color)
        plt.add(p)
    else:
        p=None


    ehigh = y + yerr
    elow  = y - yerr
    
    # only show errors where y+yerr is greater than the min value.
    w,=numpy.where(ehigh > minval)
    if w.size > 0:
        elow=numpy.where(elow < minval, minval, elow)
        pe = ErrorBarsY(x[w], elow[w], ehigh[w],color=color)
        plt.add(pe)
    else:
        pe=None

    odict={}
    odict['p'] = p
    odict['pe'] = p

    xrng = eu.plotting.get_log_plot_range(x)
    yrng = eu.plotting.get_log_plot_range(y,yerr)

    return xrng, yrng

class Grid(object):
    """
    represent a grid of plots
    """
    def __init__(self, nplot):
        self.set_grid(nplot)
    
    def set_grid(self, nplot):
        """
        set the grid given the number of plots
        """
        from math import sqrt

        self.nplot=nplot
        sq=int(sqrt(nplot))
        if nplot==sq*sq:
            self.nrow, self.ncol = sq,sq
        elif nplot <= sq*(sq+1):
            self.nrow, self.ncol = sq,sq+1
        else:
            self.nrow, self.ncol = sq+1,sq+1

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

        for i in xrange(nplot):
            row,col=grid.get_rowcol(nplot, i)
            arr[row,col].add( ... )
        """

        imax=self.nplot-1
        if index > imax:
            raise ValueError("index too large %d > %d" % (index,imax))

        row = index/self.ncol
        col = index % self.ncol

        return row,col

def get_framed_array(nplot):
    """
    get a grid and biggles.FramedArray with the right grid shape

    parameters
    ----------
    nplot: int
        total number of plots

    returns
    -------
    grid, arr: tuple
        element 0: Grid
        element 1: a biggles FramedArray
    """
    import biggles

    grid=Grid(nplot)
    arr=biggles.FramedArray(grid.nrow, grid.ncol)
    return grid, arr


