import numpy
from . import files

_good_ranges={}
_good_ranges['exp'] = {'log10_flux':[-1.2, 1.4],
                       's2n_rat':[0.1, 0.7]}
_good_ranges['dev'] = {'log10_flux':[-1.2, 1.6],
                       's2n_rat':[0.1, 0.6]}
_good_ranges['bdf'] = {'log10_flux':[-1.2, 1.6],
                       's2n_rat':[0.01, 0.7]}

def select_by_s2n_flux(data, model, good=True):
    """
    Select a "good" range.  This is somewhat arbitrary
    based on find_good_s2n
    """

    flux, flux_err, T, T_err = get_flux_T(data, model)

    flux_s2n=flux/flux_err
    T_s2n=T/T_err

    logflux=numpy.log10(flux)

    s2nrat = T_s2n/flux_s2n

    logflux_range=_good_ranges[model]['log10_flux']
    s2n_rat_range=_good_ranges[model]['s2n_rat']

    if good:
        w,=numpy.where(   (logflux > logflux_range[0])
                        & (logflux < logflux_range[1])
                        & (s2nrat  > s2n_rat_range[0])
                        & (s2nrat  < s2n_rat_range[1]) )
    else:
        w,=numpy.where(   (logflux < logflux_range[0])
                        | (logflux > logflux_range[1])
                        | (s2nrat  < s2n_rat_range[0])
                        | (s2nrat  > s2n_rat_range[1]) )

    return w


def get_flux_T(data, model):
    pars_name='%s_pars' % model
    pcov_name='%s_pars_cov' % model
    flux=data[pars_name][:, 5:].sum(axis=1)
    flux_var=data[pcov_name][:,5:,5:].sum(axis=1).sum(axis=1)
    flux_err=numpy.sqrt(flux_var)
    T=data[pars_name][:,4]
    T_err=numpy.sqrt(data[pcov_name][:,4,4])

    return flux, flux_err, T, T_err

def find_good_s2n(version, model):
    """
    Make a plot of T_s2n/flux_s2n vs log10(flux)

    Indicates good ranges are
        exp
            log10 flux between [-1.2, 1.4]
            T_s2n/flux_s2n between [0.1,0.7]
        dev
            log10 flux between [-1.2, 1.6]
            T_s2n/flux_s2n between [0.1,0.6]

    """
    import biggles

    biggles.configure('screen','width',1200)
    biggles.configure('screen','height',1200)

    t=files.read_output(version)

    flux, flux_err, T, T_err = get_flux_T(t, model)
    flux_s2n = flux/flux_err
    T_s2n = T/T_err

    s2n_rat = T_s2n/flux_s2n

    logflux=numpy.log10(flux)

    wbad=select_by_s2n_flux(t, model, good=False)

    xlabel=r'$log_{10}(flux)$'
    ylabel=r'$(S/N)_T / (S/N)_{flux}$'

    xrange=[-2.6,2.2]
    yrange=[0.0,1.0]

    plt=biggles.FramedPlot()
    plt.title=model
    plt.xlabel=xlabel
    plt.ylabel=ylabel
    plt.xrange=xrange
    plt.yrange=yrange

    pts = biggles.Points(logflux, s2n_rat,
                         type='filled circle', size=0.3)
    bad_pts = biggles.Points(logflux[wbad], s2n_rat[wbad],
                             type='filled circle', size=0.3, color='red')
    
    plt.add(pts)
    plt.add(bad_pts)

    plt.show()

