from __future__ import print_function
try:
    xrange
except:
    xrange=range

import numpy
import esutil as eu

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
    import images

    dolog=kw.pop('log',True)

    kw['more']=True
    hdict = eu.stat.histogram2d(x, y, **kw)

    hist = hdict['hist']
    if dolog:
        hist = numpy.log10( hist.clip(min=0.1) )

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

    if len(data.shape) != 2:
        raise ValueError("data should have shape [npoints,ndim]")

    ndim=data.shape[1]

    labels=kw.pop('labels',None)
    if labels is not None:
        nl=len(labels)
        assert len(labels)==ndim,"len(labels) = %d != %d" % (nl,ndim)

    grid=Grid(ndim)

    tab = kw.pop('plt',None)
    if tab is not None:
        add_to_existing_plots=True
    else:
        add_to_existing_plots=False
        tab=biggles.Table(grid.nrow, grid.ncol)
        tab.aspect_ratio=kw.pop('aspect_ratio',None)

        if tab.cols != grid.ncol or tab.rows != grid.nrow:
            m="input table has wrong dims.  Expected %s got %s"
            tup = ((grid.nrow,grid.ncol),(tab.rows,tab.cols))
            raise ValueError(m % tup)

    for dim in xrange(ndim):

        ddata = data[:,dim]

        mn, std = eu.stat.sigma_clip(ddata)

        binsize=binfac*std

        hc = biggles.make_histc(
            ddata,
            binsize=binsize,
            **kw
        )
        
        row,col=grid(dim)
        if add_to_existing_plots:
            plt = tab[row,col]
            plt.add(hc)
        else:
            plt = biggles.FramedPlot(aspect_ratio=1, **kw)
            plt.add(hc)
            tab[row,col]=plt

        
        if labels is not None:
            lab=labels[dim]
        else:
            lab='dim %d' % dim
        tab[row,col].xlabel=lab
    return tab
        
class Grid(object):
    """
    represent plots in a grid.  The grid is chosen
    based on the number of plots
    
    example
    -------
    grid=Grid(n)

    for i in xrange(n):
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

        self.nplot=nplot

        # first check some special cases
        if nplot==8:
            self.nrow, self.ncol = 2,4
        else:

            sq=int(sqrt(nplot))
            if nplot==sq*sq:
                self.nrow, self.ncol = sq,sq
            elif nplot <= sq*(sq+1):
                self.nrow, self.ncol = sq,sq+1
            else:
                self.nrow, self.ncol = sq+1,sq+1

        self.nplot_tot=self.nrow*self.ncol

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

        imax=self.nplot_tot-1
        if index > imax:
            raise ValueError("index too large %d > %d" % (index,imax))

        row = index//self.ncol
        col = index % self.ncol

        return row,col

    def __call__(self, index):
        return self.get_rowcol(index)


