from __future__ import print_function
import numpy
from . import files

logf_binsize=0.1
logf_min, logf_max =-5.0, 2.0
logr_binsize=0.05
logr_min,logr_max=-0.2,2.5
#n_binsize=0.02
n_binsize=0.08
n_min,n_max=0.11,5.99
g_binsize=0.01
g_min=0.0
g_max=0.9

MIN_COVAR=1.e-6

class CosmosSersicN(object):
    def __init__(self, n_gauss=10):
        import ngmix
        self.data=files.read_my_lackner_sersicn_fits(n_gauss)
        self.n_gauss=n_gauss

        self.gmnd=ngmix.gmix.GMixND(self.data['weights'],
                                    self.data['means'],
                                    self.data['covars'])

        self._make_gmm()

    def get_prob_scalar(self, n):
        """
        Additional checks on bounds over GMixND.  Scalar input.
        """
        from ngmix.priors import LOWVAL
        if not self.check_bounds_scalar(n):
            return LOWVAL
        else:
            return self.gmnd.get_prob_scalar(n)

    def get_lnprob_scalar(self, n):
        """
        Additional checks on bounds over GMixND.  Scalar input.
        """
        from ngmix.priors import LOWVAL
        if not self.check_bounds_scalar(n):
            return LOWVAL
        else:
            return self.gmnd.get_lnprob_scalar(n)

    def get_prob_array(self, n):
        """
        Additional checks on bounds over GMixND.  Array input.
        """
        from ngmix.priors import LOWVAL
        p=numpy.zeros(n.size)
        w=self.check_bounds_array(n)
        if w.size > 0:
            p[w]=self.gmnd.get_prob_array(n[w])
        return p

    def get_lnprob_array(self, n):
        """
        Additional checks on bounds over GMixND.  Array input.
        """
        from ngmix.priors import LOWVAL
        lnp=numpy.zeros(n.size) + LOWVAL
        w=self.check_bounds_array(n)
        if w.size > 0:
            lnp[w]=self.gmnd.get_lnprob_array(n[w])
        return lnp

    def check_bounds_scalar(self, n):
        """
        Check bounds on scalar input
        """
        from ngmix.gexceptions import GMixRangeError
        if n < n_min or n > n_max:
            raise GMixRangeError("n out of range")
        return self.gmnd.get_prob_scalar(n)

    def check_bounds_array(self, n):
        """
        Check bounds on scalar input
        """
        w,=numpy.where( (n > n_min) & (n < n_max) )
        return w
 
    def sample(self, n=None):
        if n is None:
            n=1
            is_scalar=True
        else:
            is_scalar=False

        samples=numpy.zeros(n)

        nleft=n
        ngood=0
        while nleft > 0:

            tsamples=self.gmm.sample(nleft).ravel()
            w=self.check_bounds_array(tsamples)

            if w.size > 0:
                first=ngood
                last=ngood+w.size

                samples[first:last] = tsamples[w]

                ngood += w.size
                nleft -= w.size

        if is_scalar:
            samples=samples[0]

        return samples


    def _make_gmm(self):
        """
        Make a GMM object for sampling
        """
        from sklearn.mixture import GMM

        # these numbers are not used because we set the means, etc by hand
        gmm=GMM(n_components=self.n_gauss,
                n_iter=10000,
                min_covar=1.0e-6,
                covariance_type='full')
        gmm.means_ = self.gmnd.means.reshape( (self.n_gauss,1) )
        gmm.covars_ = self.gmnd.covars.reshape( (self.n_gauss, 1, 1) )
        gmm.weights_ = self.gmnd.weights

        self.gmm=gmm 


def fit_sersic(pars=None, show=False, eps=None, n_gauss=20, n_iter=5000, min_covar=MIN_COVAR):
    """
    Fit functions to sersic n and ellipticity
    """
    if pars is None:
        pars=get_pars()
    
    #gpfitter=fit_g(pars[:,3])

    n_vals=pars[:,2:2+1]
    n_gmm=fit_gmix(n_vals, n_gauss, n_iter, min_covar=min_covar)

    n_rand=n_gmm.sample(pars.shape[0]*100)
    n_rand=n_rand.ravel()


    tab=make_plots(pars,n_rand=n_rand,n_title=r'$N_{gauss}=%s' % n_gauss)

    if show:
        tab.show()

    if eps is not None:
        print(eps)
        tab.write_eps(eps)

    write_sersicn_fit(n_gmm, n_gauss)

    return pars

def write_sersicn_fit(n_gmm, n_gauss):
    import fitsio
    output=numpy.zeros(n_gauss, dtype=[('means','f8'),
                                       ('covars','f8'),
                                       ('icovars','f8'),
                                       ('weights','f8')])

    output['means']=n_gmm.means_.ravel()
    output['covars']=n_gmm.covars_.ravel()
    output['weights']=n_gmm.weights_.ravel()
    output['icovars'] = 1.0/output['covars']

    fname=files.get_my_lackner_sersicn_fits(n_gauss)
    print(fname)
    with fitsio.FITS(fname,'rw',clobber=True) as fobj:
        fobj.write(output)

def fix_n(n):
    """
    Throw out some random points to fix odd spike
    """

    keep=numpy.ones(n.size)

    w,=numpy.where( (n > 1.049) & (n < 1.051) )
    n_region=w.size
    n_keep = int(n_region*2.0/60.0)
    print("keeping %d/%d from n region" % (n_keep,n_region))

    rvals = numpy.random.random(w.size)
    s=rvals.argsort()
    wrand = s[0:n_keep]

    keep[w] = 0
    keep[w[wrand]] = 1

    wkeep,=numpy.where(keep==1)

    print("keeping %d/%d overall" % (wkeep.size,n.size))
    return wkeep

def get_pars():
    """
    read lackner fits and remove flagged objects
    """
    data=files.read_fits_cat()

    sersicn=data['sersicfit'][:,2]
    w=fix_n(sersicn)
    data=data[w]

    logf=numpy.log10( data['sersicfit'][:,0] )
    logr=numpy.log10( data['sersicfit'][:,1] )
    q=data['sersicfit'][:,3]
    sersicn=data['sersicfit'][:,2]
    g = (1.0-q)/(1.0+q)


    w, = numpy.where(  (data['fit_status'][:,4] != 0)
                     & (data['fit_status'][:,4] != 5)
                     & (logf > logf_min)
                     & (logf < logf_max)
                     & (logr > logr_min)
                     & (logr < logr_max)
                     & (sersicn > n_min)
                     & (sersicn < n_max)
                     & (g >= g_min)
                     & (g < g_max) )

    pars=numpy.zeros( (w.size, 4) )
    pars[:,0] = logf[w]
    pars[:,1] = logr[w]
    pars[:,2] = sersicn[w]
    pars[:,3] = g[w]

    return pars


def make_plots(pars, n_rand=None, n_title=None):
    import esutil as eu
    import biggles

    tab=biggles.Table(2,2)
    logf_plt=eu.plotting.bhist(pars[:,0], binsize=logf_binsize, min=logf_min, max=logf_max,
                               xlabel=r'$log_{10}(F)$',
                               norm=1,
                               show=False)
    #eu.plotting.bhist(rpars[:,0], binsize=logf_binsize, min=logf_min, max=logf_max,
    #                  norm=1,
    #                  color='red',
    #                  show=False,
    #                  plt=logf_plt)

    logr_plt=eu.plotting.bhist(pars[:,1], binsize=logr_binsize, min=logr_min, max=logr_max,
                               xlabel=r'$log_{10}(r_{1/2})$',
                               norm=1,
                               show=False)
    #eu.plotting.bhist(rpars[:,1], binsize=logr_binsize, min=logr_min, max=logr_max,
    #                  norm=1,
    #                  color='red',
    #                  show=False,
    #                  plt=logr_plt)

    n_plt=eu.plotting.bhist(pars[:,2], binsize=n_binsize, min=n_min, max=n_max,
                            xlabel='n',
                            norm=1,
                            show=False)
    if n_rand is not None:
        eu.plotting.bhist(n_rand, binsize=n_binsize, min=n_min, max=n_max,
                          norm=1,
                          color='red',
                          show=False,
                          plt=n_plt)

    n_plt.title=n_title

    g_plt=eu.plotting.bhist(pars[:,3], binsize=g_binsize, min=g_min, max=g_max,
                            xlabel='|g|',
                            norm=1,
                            show=False)
    #eu.plotting.bhist(rpars[:,3], binsize=g_binsize, min=g_min, max=g_max,
    #                  norm=1,
    #                  color='red',
    #                  show=False,
    #                  plt=g_plt)


    tab[0,0]=logf_plt
    tab[0,1]=logr_plt
    tab[1,0]=n_plt
    tab[1,1]=g_plt

    return tab

def fit_gmix(data, n_gauss, n_iter, min_covar=MIN_COVAR):
    """
        gtot, T, flux
        etatot, log10(T), log10(flux)
    
    data is shape
        [npoints, ndim]
    """
    from sklearn.mixture import GMM

    gmm=GMM(n_components=n_gauss,
            n_iter=n_iter,
            min_covar=MIN_COVAR,
            covariance_type='full')

    gmm.fit(data)

    if not gmm.converged_:
        print("DID NOT CONVERGE")

    return gmm

def get_norm_hist(data, min=None, max=None, binsize=1):
    import esutil as eu

    hdict = eu.stat.histogram(data, min=min, max=max, binsize=binsize, more=True)

    hsum=float(hdict['hist'].sum())

    hist_norm = hdict['hist']/hsum
    hist_norm_err = sqrt(hdict['hist'])/hsum

    hdict['hist_norm'] = hist_norm
    hdict['hist_norm_err'] = hist_norm_err

    return hdict

def fit_g(g, gmax=0.9):
    """
    Fit only the g prior
    """
    import great3
    import esutil as eu
    import biggles
    import mcmc

    binsize=0.01

    hdict=get_norm_hist(g, min=0, binsize=binsize)

    yrange=[0.0, 1.1*(hdict['hist_norm']+hdict['hist_norm_err']).max()]
    plt=eu.plotting.bscatter(hdict['center'],
                             hdict['hist_norm'],
                             yerr=hdict['hist_norm_err'],
                             yrange=yrange,
                             show=False)

    ivar = numpy.ones(hdict['center'].size)
    w,=numpy.where(hdict['hist_norm_err'] > 0.0)
    if w.size > 0:
        ivar[w] = 1.0/hdict['hist_norm_err'][w]**2

    gpfitter=great3.fit_prior.GPriorFitterExp(hdict['center'],
                                              hdict['hist_norm'],
                                              ivar,
                                             gmax=gmax)

    gpfitter.do_fit()
    gpfitter.print_result()
    res=gpfitter.get_result()

    gvals=numpy.linspace(0.0, 1.0)
    model=gpfitter.get_model_val(res['pars'], g=gvals)
    crv=biggles.Curve(gvals, model, color='red')
    plt.add(crv)

    mcmc.plot_results(gpfitter.trials)
    plt.show()

    return gpfitter

def get_norm_hist(data, min=None, max=None, binsize=1):
    import esutil as eu

    hdict = eu.stat.histogram(data, min=min, max=max, binsize=binsize, more=True)

    hsum=float(hdict['hist'].sum())

    hist_norm = hdict['hist']/hsum
    hist_norm_err = numpy.sqrt(hdict['hist'])/hsum

    hdict['hist_norm'] = hist_norm
    hdict['hist_norm_err'] = hist_norm_err

    return hdict



