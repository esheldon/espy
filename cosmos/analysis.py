import numpy
from . import files

_good_ranges={}
_good_ranges['exp'] = {'log10_flux':[-1.2, 1.4],
                       's2n_rat':[0.1, 0.7]}
_good_ranges['dev'] = {'log10_flux':[-1.2, 1.6],
                       's2n_rat':[0.1, 0.6]}

def select_by_s2n_flux(data, model, good=True):
    """
    Select a "good" range.  This is somewhat arbitrary
    based on find_good_s2n
    """
    pars_name='%s_pars' % model
    pcov_name='%s_pars_cov' % model

    flux_s2n=data[pars_name][:, 5]/numpy.sqrt(data[pcov_name][:,5,5])
    T_s2n=data[pars_name][:, 4]/numpy.sqrt(data[pcov_name][:,4,4])
    logflux=numpy.log10(data[pars_name][:,5])

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


def find_good_s2n(version):
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

    dev_flux_s2n=t['dev_pars'][:, 5]/numpy.sqrt(t['dev_pars_cov'][:,5,5])
    dev_T_s2n=t['dev_pars'][:, 4]/numpy.sqrt(t['dev_pars_cov'][:,4,4])
    dev_s2nrat = dev_T_s2n/dev_flux_s2n


    exp_flux_s2n=t['exp_pars'][:, 5]/numpy.sqrt(t['exp_pars_cov'][:,5,5])
    exp_T_s2n=t['exp_pars'][:, 4]/numpy.sqrt(t['exp_pars_cov'][:,4,4])
    exp_s2nrat = exp_T_s2n/exp_flux_s2n

    #tab = biggles.Table(2,1)
    exp_logflux=numpy.log10(t['exp_pars'][:,5])
    dev_logflux=numpy.log10(t['dev_pars'][:,5])

    wbad_exp=select_by_s2n_flux(t, 'exp', good=False)
    wbad_dev=select_by_s2n_flux(t, 'dev', good=False)

    xlabel=r'$log_{10}(flux)$'
    ylabel=r'$(S/N)_T / (S/N)_{flux}$'

    xrange=[-2.6,2.2]
    yrange=[0.0,1.0]

    exp_plt=biggles.FramedPlot()
    exp_plt.title='exp'
    exp_plt.xlabel=xlabel
    exp_plt.ylabel=ylabel
    exp_plt.xrange=xrange
    exp_plt.yrange=yrange

    exp_pts = biggles.Points(exp_logflux, exp_s2nrat,
                             type='filled circle', size=0.3)
    exp_bad_pts = biggles.Points(exp_logflux[wbad_exp], exp_s2nrat[wbad_exp],
                                 type='filled circle', size=0.3, color='red')
    
    exp_plt.add(exp_pts)
    exp_plt.add(exp_bad_pts)

    dev_plt=biggles.FramedPlot()
    dev_plt.title='dev'
    dev_plt.xlabel=xlabel
    dev_plt.ylabel=ylabel
    dev_plt.xrange=xrange
    dev_plt.yrange=yrange

    dev_pts = biggles.Points(dev_logflux, dev_s2nrat,
                             type='filled circle', size=0.3)
    dev_bad_pts = biggles.Points(dev_logflux[wbad_dev], dev_s2nrat[wbad_dev],
                                 type='filled circle', size=0.3, color='red')
    
    dev_plt.add(dev_pts, dev_bad_pts)

    exp_plt.show()
    dev_plt.show()


