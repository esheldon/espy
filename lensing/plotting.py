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

def plot_dsig_osig_tab(comb, show=True, dsigrange=None, osigrange=None, **keywords):

    tab = Table(2,1)
    tab.aspect_ratio=1.8
    tab.cellspacing=0
    tab.cellpadding=0

    # for overplotting
    dsig_p = Points(comb['r'], comb['dsig'], color='red')
    dsigerr_p = SymErrY(comb['r'], comb['dsig'], comb['dsigerr'], color='red')

    #arr.aspect_ratio=2
    #arr.xlabel = labels['rproj']

    #top_plt=FramedPlot()
    #top_plt.aspect_ratio=1
    #top_plt.xlog=True
    #top_plt.ylog=True


    top_plt = eu.plotting.bscatter(comb['r'], comb['dsig'], yerr=comb['dsigerr'],
                         xlog=True, ylog=True, 
                         ylabel=labels['dsig'],
                         yrange=dsigrange,
                         show=False, 
                         **keywords)
    top_plt.x1.draw_ticklabels=0

    bot_plt = eu.plotting.bscatter(comb['r'], comb['osig'], yerr=comb['dsigerr'],
                                   xlog=True, 
                                   xlabel=labels['rproj'],
                                   ylabel=labels['osig'],
                                   yrange=osigrange,
                                   show=False,
                                   **keywords)

    c=Curve([1.e-5,1.e5],[0,0], type='solid')
    bot_plt.add(c)

    bot_plt.add(dsig_p,dsigerr_p)

    #top_plt.aspect_ratio=1
    #bot_plt.aspect_ratio=0.2

    tab[0,0] = top_plt
    tab[1,0] = bot_plt

    if show:
        tab.show()
    return tab




def plot_dsig_osig_arr(comb, show=True, dsigrange=None, osigrange=None, **keywords):

    # for overplotting
    dsig_p = Points(comb['r'], comb['dsig'], color='red')
    dsigerr_p = SymErrY(comb['r'], comb['dsig'], comb['dsigerr'], color='red')
    dsig_p.label=r'$\Delta\Sigma_+$'

    arr = FramedArray(2,1)
    arr.aspect_ratio=2
    arr.xlabel = labels['rproj']

    eu.plotting.bscatter(comb['r'], comb['dsig'], yerr=comb['dsigerr'],
                         xlog=True, ylog=True, 
                         ylabel=labels['dsig'],
                         yrange=dsigrange,
                         show=False, 
                         plt=arr[0,0],
                         **keywords)


    pdict=eu.plotting.bscatter(comb['r'], comb['osig'], yerr=comb['dsigerr'],
                               xlog=True, 
                               xlabel=labels['rproj'],
                               ylabel=labels['osig'],
                               yrange=osigrange,
                               label=r'$\Delta\Sigma_\times$',
                               show=False,
                               dict=True,
                               plt=arr[1,0],
                               **keywords)
    c=Curve([1.e-5,1.e5],[0,0], type='solid')
    arr[1,0].add(c)

    arr[1,0].add(dsig_p,dsigerr_p)

    key = PlotKey(0.9,0.9, [dsig_p,pdict['p']], halign='right')
    arr[1,0].add(key)

    arr[0,0].ylabel = labels['dsig']
    arr[1,0].ylabel = labels['osig']
    arr.ylabel = labels['dsig']

    if show:
        arr.show()
    return arr



 
def plot_dsig(comb=None, r=None, dsig=None, dsigerr=None, 
              color='black',type='filled circle',
              nolabel=False, noshow=False, minval=1.e-3,
              aspect_ratio=1):
    """
    This one stands alone. 
    """


    if comb is not None:
        r=comb['r']
        dsig=comb['dsig']
        dsigerr=comb['dsigerr']
    else:
        if r is None or dsig is None or dsigerr is None:
            raise ValueError("Send a combined struct or r,dsig,dsigerr")

    plt=FramedPlot()
    plt.aspect_ratio=aspect_ratio
    plt.xlog=True
    plt.ylog=True

    if not nolabel:
        plt.xlabel = labels['rproj']
        plt.ylabel = labels['dsig']

    od=add_to_log_plot(plt, r, dsig, dsigerr, 
                       color=color, 
                       type=type,
                       minval=minval)

    plt.xrange = od['xrange']
    plt.yrange = od['yrange']

    if not noshow:
        plt.show()

    od['plt'] = plt
    return od

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


    ehigh = y + yerr
    elow  = y - yerr

    
    # only show errors where y+yerr is greater than the min value.
    w=where1(ehigh > minval)
    if w.size > 0:
        elow = where(elow < minval, minval, elow)
        pe = ErrY(x[w], elow[w], ehigh[w],color=color)
        plt.add(pe)

    odict={}
    odict['p'] = p
    odict['pe'] = p
    odict['xrange'] = [0.3*x.min(), 2.25*x.max()]
    odict['yrange'] = [0.3*elow.min(), 2.25*ehigh.max()]
    return odict


