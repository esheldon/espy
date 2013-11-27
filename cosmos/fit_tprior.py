import numpy
from . import files
from . import analysis

def fit_Tdist(version,show=True,compare_ngmix=False):

    ngauss=2
    data=files.read_output(version)

    wexp=analysis.select_by_s2n_flux(data, 'exp')
    wdev=analysis.select_by_s2n_flux(data, 'dev')

    exp_T = data['exp_pars'][wexp,4]
    dev_T = data['dev_pars'][wdev,4]
    exp_logT = numpy.log10( exp_T )
    dev_logT = numpy.log10( dev_T )

    exp_gmm=fit_gmix1d(exp_logT,ngauss)
    dev_gmm=fit_gmix1d(dev_logT,ngauss)

    if show:
        if compare_ngmix:
            plot_fits_2gauss_ngmix(exp_T, 'exp')
            plot_fits_2gauss_ngmix(dev_T, 'dev')
        else:
            plot_fits_2gauss(exp_logT, exp_gmm, 'exp')
            plot_fits_2gauss(dev_logT, dev_gmm, 'dev')

    print 'exp'
    print '    means:  ',exp_gmm.means_[:,0]
    print '    sigmas: ',numpy.sqrt(exp_gmm.covars_[:,0])
    print '    weights:',exp_gmm.weights_
    print 'dev'
    print '    means:  ',dev_gmm.means_[:,0]
    print '    sigmas: ',numpy.sqrt(dev_gmm.covars_[:,0])
    print '    weights:',dev_gmm.weights_

def plot_fits_2gauss(logT, gmm, type):
    import esutil as eu
    num=logT.size
    samp=gmm.sample(num).ravel()

    means=gmm.means_[:,0]
    sigmas=numpy.sqrt( gmm.covars_[:,0] )
    weights=gmm.weights_

    samp1= means[0] + sigmas[0]*numpy.random.randn(weights[0]*num)
    samp2= means[1] + sigmas[1]*numpy.random.randn(weights[1]*num)

    xlabel=r'$log_{10}T$'
    binsize=0.1
    plt=eu.plotting.bhist(logT,
                          binsize=binsize,
                          xlabel=xlabel,
                          title=type,
                          show=False)

    width=2
    plt=eu.plotting.bhist(samp,
                          binsize=binsize,
                          plt=plt,
                          color='blue',
                          width=width,
                          show=False)

    plt=eu.plotting.bhist(samp1,
                          binsize=binsize,
                          plt=plt,
                          color='red',
                          width=width,
                          show=False)

    plt=eu.plotting.bhist(samp2,
                          binsize=binsize,
                          plt=plt,
                          color='darkgreen',
                          width=width,
                          show=False)
    plt.show()

def plot_fits_2gauss_ngmix(T, type):
    """
    compare to the class we added to ngmix.priors
    """
    import ngmix
    import esutil as eu
    import biggles


    if type=='exp':
        p=ngmix.priors.TPriorCosmosExp()
        w,=numpy.where(T < 10)
    else:
        p=ngmix.priors.TPriorCosmosDev()
        w,=numpy.where(T < 100)

    #xlabel=r'$log_{10}T$'
    xlabel='T'
    binsize=0.1
    plt,hobj=eu.plotting.bhist(T[w],
                               binsize=binsize,
                               xlabel=xlabel,
                               title=type,
                               show=False,
                               gethist=True)

    hsum=hobj['hist'].sum()
    delta=hobj['high'][0]-hobj['low'][0]
    pvals=p.get_prob_array(hobj['center'])
    #pvals *= ( hsum/pvals.sum() * delta )
    pvals *= ( hsum/pvals.sum() )

    width=2
    c=biggles.Curve(hobj['center'], pvals, color='blue', type='solid', width=width)
    plt.add(c)
    plt.show()

def fit_gmix1d(data, ngauss):
    from sklearn.mixture import GMM
    gmm=GMM(n_components=ngauss, n_iter=1000)
    gmm.fit(data.reshape(data.size,1))

    if not gmm.converged_:
        raise ValueError("did not converge")

    return gmm


