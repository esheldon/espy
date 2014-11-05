from __future__ import print_function
import numpy

labels={}
labels['rproj'] = r'$R$ [$h^{-1}$ Mpc]'
labels['dsig'] = r'$\Delta\Sigma ~[M_{sun} pc^{-2}]$'
labels['osig'] = r'$\Delta\Sigma_\times ~ [M_{sun} pc^{-2}]$'


def plot_dsig(data, **keys):
    """
    This one stands alone. 
    """

    from biggles import FramedPlot, Points, SymmetricErrorBarsY
    import esutil as eu

    color=keys.get('color','black')
    type=keys.get('type','filled circle')
    nolabel=keys.get('nolabel',False)
    show=keys.get('show',True)
    minval=keys.get('minval',1.e-3)
    xlog=keys.get('xlog',True)
    ylog=keys.get('ylog',True)
    aspect_ratio=keys.get('aspect_ratio',1)
    plt=keys.get('plt',None)

    r=data['r'].ravel()
    dsig=data['dsig'].ravel()
    dsigerr=data['dsigerr'].ravel()

    if plt is None:
        plt=FramedPlot()
        plt.aspect_ratio=aspect_ratio
        plt.xlog=xlog
        plt.ylog=ylog

        if not nolabel:
            plt.xlabel = labels['rproj']
            plt.ylabel = labels['dsig']

    if ylog:
        xrng,yrng=add_to_log_plot(plt, r, dsig, dsigerr, 
                                      color=color, 
                                      type=type,
                                      minval=minval)
        plt.xrange = xrng
        plt.yrange = yrng

    else:
        zpts=Curve(r, dsig*0)
        plt.add(zpts)

        pts=Points(r, dsig, type=type, color=color)

        plt.add(pts)
        epts=SymmetricErrorBarsY(r, dsig, dsigerr, color=color)
        plt.add(epts)

        yrng=keys.get('yrange',None)
        xrng=keys.get('xrange',None)
        if yrng:
            plt.yrange=yrng
        if xrng:
            plt.xrange=xrng
        else:
            if xlog:
                plt.xrange=eu.plotting.get_log_plot_range(r)

    if show:
        plt.show()

    return plt

def add_to_log_plot(plt, x, y, yerr, 
                    color='black', type='filled circle',
                    minval=1.e-5):
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


