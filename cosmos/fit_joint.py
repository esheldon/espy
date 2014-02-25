"""
Fit joint flux-size distribution
"""

import os
import numpy

from . import analysis
from . import files

NGAUSS=8
N_ITER=1000
MIN_COVAR=1.0e-6
COVARIANCE_TYPE='full'

def make_joint_gmm(means, covars, weights):
    """
    Make a GMM object from the inputs
    """
    from sklearn.mixture import GMM
    # we will over-ride values, pars here shouldn't matter except
    # for consistency
    gmm=GMM(n_components=NGAUSS,
            n_iter=N_ITER,
            min_covar=MIN_COVAR,
            covariance_type=COVARIANCE_TYPE)
    gmm.means_ = means.copy()
    gmm.covars_ = covars.copy()
    gmm.weights_ = weights.copy()

    return gmm

def fit_joint(version, model):

    ngauss=8
    data=files.read_output(version)

    w=select_by_flux(data, model)

    pname = '%s_pars' % model
    log_T_flux = numpy.log10( data[pname][w, 4:4+2] )

    gmm=fit_gmix2d(log_T_flux,ngauss)

    output=numpy.zeros(NGAUSS, dtype=[('means','f8',2),
                                      ('covars','f8',(2,2)),
                                      ('weights','f8')])
    output['means']=gmm.means_
    output['covars']=gmm.covars_
    output['weights']=gmm.weights_


    eps=files.get_dist_path(version,model,'joint-dist',ext='eps')
    gmmtest=make_joint_gmm(output['means'], output['covars'], output['weights'])

    log_flux_mode, log_T_near=plot_fits(log_T_flux, gmmtest, model, eps=eps)

    flux_mode = 10.0**log_flux_mode
    T_near = 10.0**log_T_near

    print
    print 'log_flux_mode:',log_flux_mode
    print 'log_T_near:   ',log_T_near
    print 'flux_mode:    ',flux_mode
    print 'T_near:       ',T_near

    h={'logfmode':log_flux_mode,
       'fmode':flux_mode,
       'logTnear':log_T_near,
       'Tnear':T_near}
    files.write_dist(version, model, 'joint-dist', output,header=h)

    return output, h


def fit_gmix2d(data, ngauss):
    """
    Send log10(T) and log10(flux).
    
    data is shape
        [npoints, ndim]
    so could be for instance
        data['pars'][:,4:4+2] 
    """
    from sklearn.mixture import GMM

    gmm=GMM(n_components=NGAUSS,
            n_iter=N_ITER,
            min_covar=MIN_COVAR,
            covariance_type=COVARIANCE_TYPE)

    gmm.fit(data)

    if not gmm.converged_:
        raise ValueError("did not converge")

    return gmm

_good_ranges={}
_good_ranges['exp'] = {'flux':[0.0, 100.0]}
_good_ranges['dev'] = {'flux':[0.0, 100.0]}
 
def select_by_flux(data, model):
    """
    Very loose cuts
    """

    flux, flux_err, T, T_err = analysis.get_flux_T(data, model)

    rng = _good_ranges[model]['flux']
    w,=numpy.where( (flux > rng[0]) & (flux < rng[1]) )
    return w


def plot_fits(data, gmm, model, show=False, eps=None):
    """
    """
    import esutil as eu
    import biggles

    num=data.shape[0]

    samples=gmm.sample(num*100)

    # min=0 just to make things line up

    plt_flux, log_flux_mode = \
            _plot_single(data[:,1], samples[:,1], put_mode=True)
    plt_flux.xlabel=r'$%s log_{10}(flux)$' % model

    # get the mean log_T near the flux mode
    w_near,=numpy.where( numpy.abs( data[:,1]-log_flux_mode) < 0.01 )
    log_T_near = data[w_near,0].mean()


    log_sigma = log_T_to_log_sigma(data[:,0])
    log_sigma_samples = log_T_to_log_sigma(samples[:,0])
    log_sigma_near = log_T_to_log_sigma(log_T_near)

    plt_sigma = _plot_single(log_sigma, log_sigma_samples,
                             put_point=log_sigma_near)
    plt_sigma.xlabel=r"$%s log_{10}(\sigma'')$" % model


    tab=biggles.Table(1, 2)
    tab[0,0] = plt_sigma
    tab[0,1] = plt_flux
    tab.aspect_ratio=0.5

    if show:
        tab.show()

    if eps:
        print eps
        d=os.path.dirname(eps)
        if not os.path.exists(d):
            os.makedirs(d)
        tab.write_eps(eps)

    return log_flux_mode, log_T_near

def log_T_to_log_sigma(log_T):
    T = 10.0**log_T
    sigma = numpy.sqrt(T/2.0 )

    log_sigma = numpy.log10( sigma )

    return log_sigma

def _plot_single(data, samples, put_mode=False, put_point=None):
    import biggles

    valmin=data.min()
    valmax=data.max()

    std = data.std()
    binsize=0.1*std

    hdict = get_norm_hist(data, min=valmin, max=valmax, binsize=binsize)
    sample_hdict = get_norm_hist(samples, min=valmin, max=valmax, binsize=binsize)

    hist=hdict['hist_norm']
    sample_hist=sample_hdict['hist_norm']


    ph = biggles.Histogram(hist, x0=valmin, width=4, binsize=binsize)
    ph.label = 'cosmos'
    sample_ph = biggles.Histogram(sample_hist, x0=valmin, 
                                  width=1, color='red', binsize=binsize)
    sample_ph.label = 'joint fit'


    key = biggles.PlotKey(0.9, 0.9, [ph, sample_ph], halign='right')

    plt = biggles.FramedPlot()

    plt.add( ph, sample_ph, key )

    if put_point:
        yp = numpy.interp( put_point, hdict['center'], hist )
        point = biggles.Point(put_point, yp,
                              type='asterisk', size=4, color='magenta')
        plt.add(point)

        return plt

    elif put_mode:
        wmode=sample_hist.argmax()
        mode_loc=sample_hdict['center'][wmode]
        mode_height = sample_hist[wmode]

        mpt=biggles.Point(mode_loc, mode_height, 
                          type='asterisk', size=4, color='magenta')
        plt.add(mpt)

        return plt, mode_loc

    else:
        return plt


def get_norm_hist(data, min=None, max=None, binsize=1):
    import esutil as eu

    hdict = eu.stat.histogram(data, min=min, max=max, binsize=binsize, more=True)

    hist_norm = hdict['hist']/float(hdict['hist'].sum())
    hdict['hist_norm'] = hist_norm

    return hdict
