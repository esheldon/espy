from __future__ import print_function, division
import os
import numpy
from numpy import sqrt, zeros, where
from . import analyze, files

GMAX_HIGHS2N=0.985

def fit_g_prior(run, type='great-des', show=False, gmax=1.0, **keys):
    """
    this gmax is not the max in the cuts, but in the prior and histogram
    """
    import esutil as eu
    import ngmix

    # higher s/n requirements than for flux/T fits
    g=cache_data(run, 'g',
                 min_s2n=100.0,
                 min_Ts2n=30.0,
                 **keys)

    hd=eu.stat.histogram(g, nbin=100, min=0, max=gmax,more=True)

    xdata=hd['center']
    ydata=hd['hist'].astype('f8')

    if type=='no-exp':
        p=ngmix.priors.GPriorGreatDESNoExp(gmax=gmax)
    elif type=='great-des':
        p=ngmix.priors.GPriorGreatDES(gmax=gmax)
    elif type=='great-des2':
        p=ngmix.priors.GPriorGreatDES2()
    elif type=='ba':
        p=ngmix.priors.GPriorBA(0.2)
    elif type=='merf':
        p=ngmix.priors.GPriorMErf()
    else:
        raise ValueError("bad type '%s'" % type)

    p.dofit(hd['center'], hd['hist'].astype('f8'), show=show)


NKEEP_DEFAULT=100000
def fit_log_TF_prior(run,
                     ngauss=10,
                     min_covar=1.0e-6,
                     n_iter=5000,
                     show=False,
                     nkeep=NKEEP_DEFAULT,
                     **keys):

    partype='log_TF'
    par_labels=['log(T)','log(F)']

    _fit_prior_generic(run,
                       partype,
                       par_labels,

                       ngauss=ngauss,
                       min_covar=min_covar,
                       n_iter=n_iter,
                       show=show,
                       nkeep=nkeep,
                       **keys)


def fit_log_F_fracdev_prior(run,
                            ngauss=20,
                            min_covar=1.0e-6,
                            n_iter=5000,
                            show=False,
                            nkeep=NKEEP_DEFAULT,
                            **keys):

    partype='log_F_fracdev'
    par_labels=['log(F)','fracdev']

    _fit_prior_generic(run,
                       partype,
                       par_labels,

                       ngauss=ngauss,
                       min_covar=min_covar,
                       n_iter=n_iter,
                       show=show,
                       nkeep=nkeep,
                       **keys)

def _fit_prior_generic(run,
                       partype,
                       par_labels,
                       ngauss=10,
                       min_covar=1.0e-6,
                       n_iter=5000,
                       show=False,

                       **keys): # extra keys for cache and selection

    import esutil as eu
    import ngmix

    prior_file, eps_file = get_output_files(run, partype)

    log_TF=cache_data(run, partype, **keys)
    
    gm=ngmix.gmix.GMixND()
    gm.fit(log_TF, ngauss, n_iter=n_iter, min_covar=min_covar)

    do_plot_fits(log_TF, gm, par_labels, eps_file, show=show)
    write_prior(prior_file, gm.weights, gm.means, gm.covars)


def do_plot_fits(data, gm, par_labels, eps_file, nrand=1000000, show=False):
    import gaussmix
    
    
    ngauss,ndim = gm.means.shape
    tgmm = gaussmix.MyGMM(n_components=ngauss,
                          n_iter=100,
                          min_covar=1.0e-12,
                          covariance_type='full')
    tgmm.weights_ = gm.weights
    tgmm.means_ = gm.means
    tgmm.covars_ = gm.covars

    rsamples, comps = tgmm.samples_and_comps(nrand)
    plot_fits(data,
              rsamples,
              comps,
              show=show,
              par_labels=par_labels,
              eps=eps_file)
 



def write_prior(fname, weights, means, covars):
    import fitsio
    output=make_ngauss_output(weights, means, covars)
    print("writing:",fname)
    fitsio.write(fname, output, clobber=True)

def make_ngauss_output(weights, means, covars):
    ngauss,ndim = means.shape
    output=zeros(ngauss, dtype=[('weights','f8'),
                                ('means','f8',ndim),
                                ('covars','f8',(ndim,ndim))])
    output['weights']=weights
    output['means']=means
    output['covars']=covars
    return output

def cache_data(run,
               type,
               nkeep=None,
               **keys):
    """
    type should be g or log_TF
    """
    import esutil as eu
    import fitsio

    keys['min_s2n']=keys.get('min_s2n',1.0)
    keys['min_Ts2n']=keys.get('minTs2n',10.0)
    keys['max_g']=keys.get('max_g',GMAX_HIGHS2N)

    keys['fracdev_err_max']=keys.get('fracdev_err_max',0.05)
    keys['cut_fracdev_exact']=keys.get('cut_fracdev_exact',True)
    keys['fracdev_range']=keys.get('fracdev_range',[-1.0,2.0])


    tmpdir=os.environ['TMPDIR']
    cachename=os.path.join(tmpdir,'%s-%s.fits' % (run,type))

    if os.path.exists(cachename):
        print("reading:",cachename)
        tmp=fitsio.read(cachename)
    else:

        conf=files.read_config(run=run)

        dlist=[]
        for i in xrange(conf['ng']):
            fname=files.get_collated_file(run=run, gnum=i)

            print("reading:",fname)
            data=fitsio.read(fname)

            print("    selecting")
            data=analyze.select_good(data, **keys)

            dlist.append(data)

        data=eu.numpy_util.combine_arrlist(dlist)

        if type=='g':
            tmp=numpy.zeros(data.size, dtype=[('g','f4')])
            g1 = data['g'][:,0] - data['g'][:,0].mean()
            g2 = data['g'][:,1] - data['g'][:,1].mean()

            #g=sqrt(data['g'][:,0]**2 + data['g'][:,1]**2)
            g=sqrt(g1**2 + g2**2)
            tmp['g'] = g
        elif type=='log_TF':
            tmp=numpy.zeros(data.size, dtype=[('log_TF','f4',2)])
            tmp['log_TF'][:,0] = data['log_T']
            tmp['log_TF'][:,1] = data['log_flux']
        elif type=='log_F_fracdev':
            tmp=numpy.zeros(data.size, dtype=[('log_F_fracdev','f4',2)])
            tmp['log_F_fracdev'][:,0] = data['log_flux']
            tmp['log_F_fracdev'][:,1] = data['fracdev']
        else:
            raise ValueError("bad type '%s'" % type)

        print("writing:",cachename)
        fitsio.write(cachename, tmp, clobber=True)

    data = tmp[type]
    ntot=data.shape[0]
    print("nobj:",data.shape[0])
    if nkeep is not None and nkeep < ntot:
        print("fitting random subset: %d  (%.2f%%)" % (nkeep, 100*float(nkeep)/ntot) )
        ri=eu.random.random_indices(ntot, nkeep)
        data=data[ri,:]

    return data



def get_output_files(run, partype):
    prior_file=files.get_prior_file(run=run,
                                    partype=partype,
                                    ext='fits')
    eps_file=files.get_prior_file(run=run,
                                  partype=partype,
                                  ext='eps')
    d=os.path.dirname(prior_file)
    if not os.path.exists(d):
        os.makedirs(d)

    return prior_file, eps_file

def plot_fits(pars, samples, comps, dolog=True, show=False, eps=None, par_labels=None):
    """
    """
    import esutil as eu
    import biggles
    import images

    biggles.configure('screen','width', 1400)
    biggles.configure('screen','height', 800)

    num=pars.shape[0]
    ndim=pars.shape[1]

    nrow,ncol = images.get_grid(ndim) 
    tab=biggles.Table(nrow,ncol)


    for dim in xrange(ndim):
        plt = _plot_single(pars[:,dim], samples[:,dim], comps, do_ylog=True)
        if par_labels is not None:
            plt.xlabel=par_labels[dim]
        else:
            plt.xlabel=r'$P_%s$' % dim

        row=(dim)/ncol
        col=(dim) % ncol

        tab[row,col] = plt

    tab.aspect_ratio=nrow/float(ncol)

    if eps:
        import converter
        print(eps)
        d=os.path.dirname(eps)
        if not os.path.exists(d):
            os.makedirs(d)
        tab.write_eps(eps)
        converter.convert(eps, verbose=True, dpi=200)

    if show:
        tab.show()

def _plot_single(data, samples, comps, do_ylog=False):
    import biggles
    import esutil as eu
    import pcolors

    valmin=data.min()
    valmax=data.max()

    std = data.std()
    binsize=0.05*std

    ph,be,harr = biggles.make_histc(data, min=valmin, max=valmax,
                                    binsize=binsize,
                                    ylog=do_ylog, norm=1,
                                    get_hdata=True)
    sample_ph,sbe,sharr= biggles.make_histc(samples, min=valmin, max=valmax,
                                            binsize=binsize, color='red', ylog=do_ylog, norm=1,
                                            get_hdata=True)

    ph.label='data'
    sample_ph.label='fit'

    key = biggles.PlotKey(0.1, 0.9, [ph, sample_ph], halign='left')
    plt = biggles.FramedPlot()
    plt.add( ph, sample_ph, key )

    w,=where( (harr > 0) & (sharr > 0) )
    yrange=[min(harr[w].min(), sharr[w].min()),
            max(harr[w].max(), sharr[w].max())]

    if do_ylog:
        plt.ylog=True

    # now add the components
    h,rev=eu.stat.histogram(comps, rev=True)
    print(h)
    w,=where(h > 0)
    ncolors = w.size

    colors=pcolors.rainbow(ncolors)

    icolor=0
    for i in xrange(h.size):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]

            frac=float(w.size)/comps.size

            ph = biggles.make_histc(samples[w],
                                    #min=valmin, max=valmax,
                                    binsize=binsize,
                                    color=colors[icolor],
                                    ylog=do_ylog,
                                    norm=frac)
            plt.add(ph)
            icolor += 1

    plt.yrange=yrange
    return plt


