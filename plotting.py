from __future__ import print_function
try:
    xrange
except:
    xrange=range

import esutil as eu

def multihist(data, binfac=0.1, **kw):
    """
    plot a histogram for each dimension of the data

    parameters
    ----------
    data: array
        array with shape [npoints, ndim]
    binfac: float
        The binsize for each dimension will be chosen as binfac*std(dimdata)
    """
    import biggles

    if len(data.shape) != 2:
        raise ValueError("data should have shape [npoints,ndim]")

    ndim=data.shape[1]

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

        
        tab[row,col].xlabel='dim %d' % dim
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


