import numpy
from numpy import where
import esutil as eu
from esutil.numpy_util import where1
import copy

try:
    import biggles
    from biggles import FramedArray, FramedPlot, Points, \
            ErrorBarsY as ErrY, ErrorBarsX as ErrX, \
            SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX, \
            PlotKey, PlotLabel, Table, Curve
except:
    pass


labels={}
labels['rproj'] = r'$R$ [$h^{-1}$ Mpc]'
labels['dsig'] = r'$\Delta\Sigma ~[M_{sun} pc^{-2}]$'
labels['osig'] = r'$\Delta\Sigma_\times ~ [M_{sun} pc^{-2}]$'


def plot2dsig(r1, dsig1, dsig1err, r2, dsig2, dsig2err, **keys):
    """
    Plot delta sigma and a second delta sigma in two plots the second of which
    is linear.

    Parameters
    ----------
    show: bool, optional
        Show plot in a window, default True
    yrange1: [min,max]
        The y range of the delta sigma plot
    yrange2: [min,max]
        The y range of the ortho-delta sigma plot
    range4var: [min,max]
        The x range over which to calculate a osig variance
        and determine a plot range.  This is overridden by
        range2

    plot_label: string, optional
        a label for the top plot

    # label1,label2 go in a key in the bottom plot
    label1: string, optional
        a label for the first dsig
    label2: string, optional
        a label for the second dsig


    """

    color2 = 'red'
    ptype1='filled circle'
    size1=1.5
    ptype2='filled circle'
    size2=1.5

    xmul = keys.get('xmul',1.)

    show = keys.get('show',True)
    yrange1 = keys.get('yrange1',None)
    yrange2 = keys.get('yrange2',None)

    isortho = keys.get('ortho',False)

    # this is a y-log plot, use more powerful range determination
    print dsig1
    print dsig1err
    yrange1 = eu.plotting.get_log_plot_range(dsig1, err=dsig1err, 
                                             input_range=yrange1)

    rr=numpy.concatenate((r1,r2))
    xrng = eu.plotting.get_log_plot_range(rr)

    # this over-rides
    range4var = keys.get('range4var',None)
    if yrange2 is None:
        if range4var is not None:
            w=where1( (r2 >= range4var[0]) & (r2 <= range4var[1]))
            if w.size == 0:
                raise ValueError("no points in range [%d,%d]" % tuple(range4var))
            sdev = dsig2[w].std()
        else:
            sdev = dsig2.std()

        yrange2 = [-3.5*sdev, 3.5*sdev]


    label  = keys.get('plot_label',None)
    label1 = keys.get('label1',None)
    label2 = keys.get('label2',None)

    # The points and zero curve
    dsig1_p = Points(r1, dsig1, color='black', type=ptype1, size=size1)
    dsig1err_p = SymErrY(r1, dsig1, dsig1err, color='black')
    dsig1_p.label=label1

    dsig2_p = Points(r2*xmul, dsig2, color=color2, type=ptype2, size=size2)
    dsig2err_p = SymErrY(r2*xmul, dsig2, dsig2err, color=color2)
    dsig2_p.label=label2

    c=Curve([1.e-5,1.e5],[0,0], type='solid')


    arr = FramedArray(2,1)
    arr.cellspacing=1
    arr.aspect_ratio=2
    arr.xlabel = labels['rproj']
    if isortho:
        arr.ylabel = labels['osig']
    else:
        arr.ylabel = labels['dsig']
    arr.xrange = xrng


    arr[0,0].yrange = yrange1
    arr[0,0].xlog=True
    arr[0,0].ylog=True

    # biggles chokes if you give it negative data for a log plot
    arr[0,0].add(dsig1_p, dsig2_p)
    eu.plotting.add_log_error_bars(arr[0,0],'y',r1,dsig1,dsig1err,yrange1)
    eu.plotting.add_log_error_bars(arr[0,0],'y',r2,dsig2,dsig2err,yrange1,color=color2)

    if label is not None:
        arr[0,0].add(PlotLabel(0.9,0.9,label,halign='right'))


    arr[1,0].yrange = yrange2
    arr[1,0].xlog=True
    arr[1,0].add(c)
    arr[1,0].add(dsig1_p, dsig1err_p, dsig2_p, dsig2err_p)

    if label1 is not None or label2 is not None:
        key = PlotKey(0.9,0.15, [dsig1_p,dsig2_p], halign='right')
        arr[1,0].add(key)


    if show:
        arr.show()
    return arr


def plot2dsig_over(r1,dsig1,dsig1err,r2,dsig2, dsig2err, **keys):
    ptype1=keys.get('ptype1','filled circle')
    ptype2=keys.get('ptype2','filled circle')
    size1=keys.get('size1',1)
    size2=keys.get('size2',1)
    color1=keys.get('color1','red')
    color2=keys.get('color2','blue')
    label  = keys.get('label',None)
    label1 = keys.get('label1',None)
    label2 = keys.get('label2',None)
    xrng = keys.get('xrange',None)
    yrng = keys.get('yrange',None)
    show = keys.get('show',True)

    ylog = keys.get('ylog',True)
   
    plt=keys.get('plt',None)

    yall=numpy.concatenate((dsig1, dsig2))
    yerrall=numpy.concatenate((dsig1err, dsig2err))

    if yrng is None:
        if ylog:
            yrng = eu.plotting.get_log_plot_range(yall, err=yerrall, 
                                                  input_range=yrng)

    rr=numpy.concatenate((r1,r2))
    if xrng is None:
        xrng = eu.plotting.get_log_plot_range(rr)

    if plt is None:
        plt=FramedPlot()
    plt.xlog=True
    plt.ylog=ylog
    plt.xrange=xrng
    plt.yrange=yrng
    plt.xlabel = labels['rproj']
    plt.ylabel = labels['dsig']

    dsig1_p = Points(r1, dsig1, color=color1, type=ptype1, size=size1)
    dsig1err_p = SymErrY(r1, dsig1, dsig1err, color=color2)
    dsig1_p.label=label1

    dsig2_p = Points(r2, dsig2, color=color2, type=ptype2, size=size2)
    dsig2err_p = SymErrY(r2, dsig2, dsig2err, color=color2)
    dsig2_p.label=label2

    plt.add(dsig1_p, dsig2_p)

    if ylog:
        # biggles chokes if you give it negative data for a log plot
        eu.plotting.add_log_error_bars(plt,'y',r1,dsig1,dsig1err,yrng,
                                       color=color1)
        eu.plotting.add_log_error_bars(plt,'y',r2,dsig2,dsig2err,yrng,
                                       color=color2)
    else:
        err1 = biggles.SymmetricErrorBarsY(r1,dsig1,dsig1err,color=color1)
        err2 = biggles.SymmetricErrorBarsY(r2,dsig2,dsig2err,color=color2)
        plt.add(err1,err2)

        zerop = biggles.Curve(xrng, [0,0])
        plt.add(zerop)

    if label is not None:
        plt.add(PlotLabel(0.9,0.9,label,halign='right'))

    if label1 is not None or label2 is not None:
        key = PlotKey(0.9,0.15, [dsig1_p,dsig2_p], halign='right')
        plt.add(key)

    if show:
        plt.show()
    return plt

def plot2dsig_old(r, dsig1, dsig1err, dsig2, dsig2err, **keys):
    """
    Plot delta sigma and a second delta sigma in two plots the second of which
    is linear.

    Parameters
    ----------
    show: bool, optional
        Show plot in a window, default True
    range1: [min,max]
        The y range of the delta sigma plot
    range2: [min,max]
        The y range of the ortho-delta sigma plot
    range2var: [min,max]
        The x range over which to calculate a osig variance
        and determine a plot range.  This is overridden by
        range2

    plot_label: string, optional
        a label for the top plot

    # label1,label2 go in a key in the bottom plot
    label1: string, optional
        a label for the first dsig
    label2: string, optional
        a label for the second dsig


    """

    show = keys.get('show',True)
    range1 = keys.get('range1',None)
    range2 = keys.get('range2',None)

    # this over-rides
    range4var = keys.get('range4var',None)
    if range2 is None:
        if range4var is not None:
            w=where1( (r >= range4var[0]) & (r <= range4var[1]))
            if w.size == 0:
                raise ValueError("no points in range [%d,%d]" % tuple(range4var))
            sdev = dsig2[w].std()
        else:
            sdev = dsig2.std()

        range2 = [-3.5*sdev, 3.5*sdev]
        print 'range2:',range2


    label  = keys.get('plot_label',None)
    label1 = keys.get('label1',None)
    label2 = keys.get('label2',None)

    # for overplotting
    dsig1_p = Points(r, dsig1, color='red')
    dsig1err_p = SymErrY(r, dsig1, dsig1err, color='red')
    dsig1_p.label=label1

    dsig2_p = Points(r, dsig1, color='red')
    dsig2err_p = SymErrY(r, dsig1, dsig1err, color='black')
    dsig2_p.label=label1

    arr = FramedArray(2,1)
    arr.aspect_ratio=2
    arr.xlabel = labels['rproj']
    arr.ylabel = labels['dsig']

    eu.plotting.bscatter(r, dsig1, yerr=dsig1err,
                         xlog=True, ylog=True, 
                         yrange=range1,
                         show=False, 
                         plt=arr[0,0],
                         **keys)

    if label is not None:
        arr[0,0].add(PlotLabel(0.9,0.9,label,halign='right'))


    pdict=eu.plotting.bscatter(r, dsig2, yerr=dsig2err,
                               xlog=True, 
                               yrange=range2,
                               label=label2,
                               show=False,
                               dict=True,
                               plt=arr[1,0],
                               **keys)
    c=Curve([1.e-5,1.e5],[0,0], type='solid')
    arr[1,0].add(c)
    arr[1,0].add(dsig1_p,dsig1err_p)

    if label1 is not None or label2 is not None:
        key = PlotKey(0.9,0.9, [dsig1_p,pdict['p']], halign='right')
        arr[1,0].add(key)


    if show:
        arr.show()
    return arr




 
def plot_dsig(**keys):
    """
    This one stands alone. 
    """

    comb=keys.get('comb',None)
    r=keys.get('r',None)
    dsig=keys.get('dsig',None)
    dsigerr=keys.get('dsigerr',None)
    color=keys.get('color','black')
    type=keys.get('type','filled circle')
    nolabel=keys.get('nolabel',False)
    show=keys.get('show',True)
    minval=keys.get('minval',1.e-3)
    xlog=keys.get('xlog',True)
    ylog=keys.get('ylog',True)
    aspect_ratio=keys.get('aspect_ratio',1)
    plt=keys.get('plt',None)

    label=keys.get('label',None)

    if comb is not None:
        r=comb['r']
        dsig=comb['dsig']
        dsigerr=comb['dsigerr']
    else:
        if r is None or dsig is None or dsigerr is None:
            raise ValueError("Send a combined struct or r,dsig,dsigerr")

    if plt is None:
        plt=FramedPlot()
        plt.aspect_ratio=aspect_ratio
        plt.xlog=xlog
        plt.ylog=ylog

        if not nolabel:
            plt.xlabel = labels['rproj']
            plt.ylabel = labels['dsig']

    if ylog:
        od=add_to_log_plot(plt, r, dsig, dsigerr, 
                           color=color, 
                           type=type,
                           minval=minval)
        plt.xrange = od['xrange']
        plt.yrange = od['yrange']

        if label:
            od['p'].label=label
    else:
        zpts=Curve(r, dsig*0)
        plt.add(zpts)

        pts=Points(r, dsig, type=type, color=color)

        if label:
            pts.label=label

        plt.add(pts)
        if dsigerr is not None:
            epts=SymErrY(r, dsig, dsigerr, color=color)
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

    if ylog:
        od['plt'] = plt
    else:
        return plt

def plot_drho(comb=None, r=None, drho=None, drhoerr=None, 
              color='black',type='filled circle',
              nolabel=False, noshow=False, minval=1.e-5,
              aspect_ratio=1):
    """
    This one stands alone. 
    """

    if comb is not None:
        r=comb['rdrho']
        drho=comb['drho']
        drhoerr=comb['drhoerr']
    else:
        if r is None or drho is None or drhoerr is None:
            raise ValueError("Send a combined struct or r,drho,drhoerr")


    plt=FramedPlot()
    plt.aspect_ratio=aspect_ratio
    plt.xlog=True
    plt.ylog=True

    if not nolabel:
        plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        plt.ylabel = r'$\delta\rho ~ [M_{sun} pc^{-3}]$'


    od=add_to_log_plot(plt, r, drho, drhoerr, 
                       color=color, 
                       type=type,
                       minval=minval)

    # for drho we need even broader yrange
    plt.xrange = od['xrange']

    yr=od['yrange']
    plt.yrange = [0.5*yr[0], 3*yr[1]]

    if not noshow:
        plt.show()
    od['plt'] = plt
    return od

def plot_mass(comb=None, r=None, mass=None, masserr=None, 
              color='black',type='filled circle',
              nolabel=False, noshow=False, minval=1.e11,
              aspect_ratio=1):

    if comb is not None:
        r=comb['rmass']
        mass=comb['mass']
        masserr=comb['masserr']
    else:
        if r is None or mass is None or masserr is None:
            raise ValueError("Send a combined struct or r,mass,masserr")


    plt=FramedPlot()
    plt.aspect_ratio=aspect_ratio
    plt.xlog=True
    plt.ylog=True

    if not nolabel:
        plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        plt.ylabel = r'$M(<r) ~ [h^{-1} M_{sun}]$'


    od=add_to_log_plot(plt, r, mass, masserr, 
                       color=color, 
                       type=type,
                       minval=minval)

    plt.xrange = od['xrange']
    plt.yrange = od['yrange']

    if not noshow:
        plt.show()
    od['plt'] = plt
    return od






def add_to_log_plot(plt, x, y, yerr, 
                    color='black', type='filled circle',
                    minval=1.e-5):

    w = where1(y > minval)
    if w.size > 0:
        p = Points(x[w], y[w],type=type,color=color)
        plt.add(p)
    else:
        p=None


    ehigh = y + yerr
    elow  = y - yerr

    
    # only show errors where y+yerr is greater than the min value.
    w=where1(ehigh > minval)
    if w.size > 0:
        elow = where(elow < minval, minval, elow)
        pe = ErrY(x[w], elow[w], ehigh[w],color=color)
        plt.add(pe)
    else:
        pe=None

    odict={}
    odict['p'] = p
    odict['pe'] = p
    #odict['xrange'] = [0.3*x.min(), 2.25*x.max()]
    #odict['yrange'] = [0.3*elow.min(), 2.25*ehigh.max()]
    odict['xrange'] = eu.plotting.get_log_plot_range(x)
    odict['yrange'] = eu.plotting.get_log_plot_range(y,yerr)
    return odict


