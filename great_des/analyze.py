from __future__ import print_function
import numpy
from numpy import diag, sqrt, newaxis, where
from . import files

MIN_ARATE=0.3
MAX_ARATE=0.6
MIN_S2N=10.0
MIN_TS2N=2.0

SN=0.16

def load_data(run, select=True, **keys):
    """
    load all g collated files into a list of structs
    """
    import fitsio
    conf=files.read_config(run=run)

    dlist=[]
    for i in xrange(conf['ng']):
        fname=files.get_collated_file(run=run, gnum=i)

        print("reading:",fname)
        data=fitsio.read(fname)
        data=add_T(data)

        if select:
            print("    selecting")
            data=select_good(data, **keys)

        dlist.append(data)
    return dlist

def add_T(data):
    import esutil as eu

    if 'T' not in data.dtype.names:
        print("adding T")
        T = data['pars'][:,4]
        Terr = sqrt(data['pars_cov'][:,4,4])
        ndata=eu.numpy_util.add_fields(data,[('T','f8'),('T_err','f8')])
        ndata['T']=T
        ndata['T_err']=Terr
        data=ndata

    return data

def select_good(data,
                min_arate=MIN_ARATE,
                max_arate=MAX_ARATE,
                min_s2n=MIN_S2N,
                min_ts2n=MIN_TS2N):
    """
    apply standard selection
    """

    Ts2n=data['T']/data['T_err']

    logic = (  (data['flags']==0)
             & (data['fit_flags']==0)
             & (data['s2n_w'] > min_s2n)
             & (Ts2n  > min_ts2n) )

    if 'arate' in data.dtype.names:
        logic = logic & (data['arate'] > min_arate) & (data['arate'] < max_arate)

    w,=where(logic)
    print("    kept %d/%d" % (w.size, data.size))

    data=data[w]
    return data

def get_weights(data):
    """
    a good set of weights
    """
    csum=data['g_cov'][:,0,0] + 2*data['g_cov'][:,0,1] + data['g_cov'][:,1,1]
    wts=1.0/(2*SN**2 + csum)
    return wts

def fit_m_c(dlist):
    """
    get m and c
    """
    import fitting

    ng = len(dlist)
    gtrue = numpy.zeros( (ng,2) )
    gdiff = gtrue.copy()
    gdiff_err = gtrue.copy()
    for i,data in enumerate(dlist):
        gtrue[i,:],gmean,gcov=calc_gmean(data)
        gdiff[i,:] = gmean - gtrue[i,:]
        gdiff_err[i,:] = sqrt(diag(gcov))

    lf1=fitting.fit_line(gtrue[:,0], gdiff[:,0], yerr=gdiff_err[:,0])
    lf2=fitting.fit_line(gtrue[:,1], gdiff[:,1], yerr=gdiff_err[:,1])
    fitters=[lf1,lf2]

    plt=plot_gdiff_vs_gtrue(gtrue, gdiff, gdiff_err, fitters=[lf1,lf2])
    #plt.show()

    return lf1,lf2

def plot_gdiff_vs_gtrue(gtrue, gdiff, gdiff_err, fitters=None):
    """
    plot gdiff vs gtrue, both components, return the plot
    """
    import biggles

    plt=biggles.FramedPlot()
    plt.xlabel='g_true'
    plt.ylabel=r'$\Delta g$'

    color1='blue'
    color2='red'
    pts1=biggles.Points(gtrue[:,0], gdiff[:,0], type='filled circle', color=color1)
    perr1=biggles.SymmetricErrorBarsY(gtrue[:,0], gdiff[:,0], gdiff_err[:,0], color=color1)
    pts2=biggles.Points(gtrue[:,1], gdiff[:,1], type='filled diamond', color=color2)
    perr2=biggles.SymmetricErrorBarsY(gtrue[:,1], gdiff[:,1], gdiff_err[:,1], color=color2)

    pts1.label=r'$g_1$'
    pts2.label=r'$g_2$'

    key=biggles.PlotKey(0.1, 0.9, [pts1,pts2], halign='left')
    
    plt.add( pts1, perr1, pts2, perr2, key )

    if fitters is not None:
        lf1,lf2=fitters

        ply1=lf1.get_poly()
        ply2=lf2.get_poly()
        c1=biggles.Curve(gtrue[:,0], ply1(gtrue[:,0]), color=color1)
        c2=biggles.Curve(gtrue[:,1], ply2(gtrue[:,1]), color=color2)
        plt.add(c1, c2)

    return plt

def calc_gmean(data):
    """
    get gtrue, gmeas, gcov
    """
    import jackknife
    gtrue = data['shear_true'].mean(axis=0)
    gmeas = numpy.zeros(2)

    wts=get_weights(data)

    '''
    g1sum = (data['g'][:,0]*wts).sum()
    g2sum = (data['g'][:,1]*wts).sum()

    g1sensum = (data['g_sens'][:,0]*wts).sum()
    g2sensum = (data['g_sens'][:,1]*wts).sum()

    gmeas[0]=g1sum/g1sensum
    gmeas[1]=g2sum/g2sensum
    print("gmeas: ",gmeas)
    '''

    wa=wts[:,newaxis]
    jdsum=data['g']*wa
    if 'g_sens' in data.dtype.names:
        jwsum=data['g_sens']*wa
    else:
        jwsum=numpy.ones( data['g'].shape )*wa
    #print(jdsum.shape)

    gmeas,gcov=jackknife.wjackknife(vsum=jdsum, wsum=jwsum)

    return gtrue, gmeas, gcov

class AnalyzerS2N(dict):
    """
    analyze m and c vs various parameters
    """
    def __init__(self, nbin=12, min_s2n=10., max_s2n=200.):
        self.nbin=nbin
        self.min_s2n=min_s2n
        self.max_s2n=max_s2n

    def go(self, dlist):
        """
        dlist is data for each g value

        data should already be trimmed
        """

        self._calc_m_c(dlist)

    def doplot(self, show=False):
        """
        plot m,c vs s2n
        """
        import biggles

        xrng=[0.5*self.s2n.min(), 1.5*self.s2n.max()]
        mplt=biggles.FramedPlot()
        cplt=biggles.FramedPlot()

        mplt.xlabel='S/N'
        mplt.ylabel='m'
        mplt.xlog=True
        mplt.xrange=xrng

        cplt.xlabel='S/N'
        cplt.ylabel='c'
        cplt.xlog=True
        cplt.xrange=xrng

        color1='blue'
        color2='red'

        mcurve1=biggles.Curve(self.s2n, self.m[:,0], type='solid',
                             color=color1)
        merr1=biggles.SymmetricErrorBarsY(self.s2n, self.m[:,0], self.merr[:,0],
                                          color=color1)
        mcurve2=biggles.Curve(self.s2n, self.m[:,1], type='dashed',
                             color=color2)
        merr2=biggles.SymmetricErrorBarsY(self.s2n, self.m[:,1], self.merr[:,1],
                                          color=color2)
        
        ccurve1=biggles.Curve(self.s2n, self.c[:,0], type='solid',
                             color=color1)
        cerr1=biggles.SymmetricErrorBarsY(self.s2n, self.c[:,0], self.cerr[:,0],
                                          color=color1)
        ccurve2=biggles.Curve(self.s2n, self.c[:,1], type='dashed',
                             color=color2)
        cerr2=biggles.SymmetricErrorBarsY(self.s2n, self.c[:,1], self.cerr[:,1],
                                          color=color2)
 
        key=biggles.PlotKey(0.1,0.9,[mcurve1,mcurve2],halign='left')

        mcurve1.label=r'$g_1$'
        mcurve2.label=r'$g_2$'
        ccurve1.label=r'$g_1$'
        ccurve2.label=r'$g_2$'

        zc=biggles.Curve( self.s2n, self.s2n*0 )

        mplt.add( mcurve1, merr1, mcurve2, merr2, zc, key )
        cplt.add( ccurve1, cerr1, ccurve2, cerr2, zc, key )

        if show:
            mplt.show()
            cplt.show()
        return mplt, cplt

    def _calc_m_c(self, dlist):
        """
        calculate m and c in bins of s/n

        returns
        -------
        s2n, m, merr, c, cerr
        """

        # get reverse indices for our binning
        revlist=[self._do_hist1(d) for d in dlist]

        m=numpy.zeros( (self.nbin,2) )
        merr=numpy.zeros( (self.nbin,2) )
        c=numpy.zeros( (self.nbin,2) )
        cerr=numpy.zeros( (self.nbin,2) )
        s2n=numpy.zeros(self.nbin)

        for i in xrange(self.nbin):
            wlist=[ rev[ rev[i]:rev[i+1] ] for rev in revlist ]

            s2n_sum=0.0
            wsum=0.0
            cut_dlist=[]
            for d,w in zip(dlist,wlist):
                td = d[w]

                wts = get_weights(td)
                s2n_sum += (wts*td['s2n_w']).sum()
                wsum += wts.sum()

                cut_dlist.append(td)

            lf1,lf2 = fit_m_c(cut_dlist)
            m[i,0],c[i,0] = lf1.pars
            m[i,1],c[i,1] = lf2.pars
            merr[i,0],cerr[i,0] = lf1.perr
            merr[i,1],cerr[i,1] = lf2.perr
            s2n[i] = s2n_sum/wsum

            print("s2n:",s2n[i])
            print("m1: %g +/- %g" % (m[i,0],merr[i,0]))
            print("m2: %g +/- %g" % (m[i,1],merr[i,1]))
            print("c1: %g +/- %g" % (c[i,0],cerr[i,0]))
            print("c2: %g +/- %g" % (c[i,1],cerr[i,1]))

        self.s2n=s2n
        self.m=m
        self.merr=merr
        self.c=c
        self.cerr = cerr

    def _do_hist1(self, data):
        import esutil as eu
        log_s2n = numpy.log10( data['s2n_w'] )
        minl = numpy.log10( self.min_s2n )
        maxl = numpy.log10( self.max_s2n )
        h,rev=eu.stat.histogram(log_s2n,
                                min=minl,
                                max=maxl,
                                nbin=self.nbin,
                                rev=True)

        return rev

def plot_e_vs_sigma(data,
                    use_true_fwhm=False,
                    yrange=None,
                    nbin=None,
                    nperbin=None,
                    min_lsigma=-1.,
                    max_lsigma=-0.1,
                    title=None,
                    show=False):
    import biggles
    import esutil as eu

    shear_true=data['shear_true'].mean(axis=0)

    colors=['blue','red']
    types=['filled circle','filled diamond']

    if use_true_fwhm:
        log_sigma = numpy.log10( data['fwhm']/2.35*0.27 )
    else:
        log_sigma = numpy.log10( numpy.sqrt( 0.5*data['T'] ) )

    weights = get_weights(data)

    bs1=eu.stat.Binner(log_sigma, data['g'][:,0], weights=weights)
    bs2=eu.stat.Binner(log_sigma, data['g'][:,1], weights=weights)

    bs1.dohist(min=min_lsigma,
               max=max_lsigma,
               nbin=nbin,
               nperbin=nperbin)
    bs2.dohist(min=min_lsigma,
               max=max_lsigma,
               nbin=nbin,
               nperbin=nperbin)

    bs1.calc_stats()
    bs2.calc_stats()

    plt=biggles.FramedPlot()
    plt.xrange=[0.9*min_lsigma, 1.1*max_lsigma]
    plt.aspect_ratio=1
    plt.title=title
    plt.xlabel=r'$log_{10}( \sigma [arcsec] )$'
    plt.ylabel='<e>'


    pts1=biggles.Points(bs1['xmean'],bs1['ymean'],type=types[0],color=colors[0])
    err1=biggles.SymmetricErrorBarsY(bs1['xmean'],bs1['ymean'],bs1['yerr'],
                                     color=colors[0])

    pts2=biggles.Points(bs2['xmean'],bs2['ymean'],type=types[1],color=colors[1])
    err2=biggles.SymmetricErrorBarsY(bs2['xmean'],bs2['ymean'],bs2['yerr'],
                                     color=colors[1])

    sc1 = biggles.Curve(bs1['xmean'], bs1['xmean']*0 + shear_true[0])
    sc2 = biggles.Curve(bs1['xmean'], bs1['xmean']*0 + shear_true[1])

    pts1.label=r'$e_1$'
    pts2.label=r'$e_2$'

    key=biggles.PlotKey(0.1,0.9,[pts1,pts2])

    plt.add(sc1,sc2,pts1,err1,pts2,err2,key)

    if yrange is not None:
        plt.yrange=yrange

    if show:
        plt.show()
    return plt
