import numpy
from numpy import where
from esutil.numpy_util import where1

try:
    import biggles
    from biggles import FramedArray, FramedPlot, Points, \
            ErrorBarsY as ErrY, ErrorBarsX as ErrX, \
            SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX, \
            PlotKey, PlotLabel, Table, Curve
except:
    pass


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
        plt.xlabel = r'$R$ [$h^{-1}$ Mpc]'
        plt.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

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


