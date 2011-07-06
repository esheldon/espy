import numpy
from  numpy import sqrt, log10, linspace, interp
import lensing


class NFWBiasFitter:
    def __init__(self, omega_m, z, r, **kwds):
        """
        Class:
            NFWBiasFitter
        Construction:
            fitter=NFWBiasFitter(omega_m, z, r, **linear_kwds)
        Inputs:
            omega_m, z: z is the redshift of the object.
        Keywords:
            withlin: True.  Include the linear bias term.
            Keywords for the Linear() class. See the Linear() class for the
            defaults
                omega_b
                sigma_8
                h
                ns

        Usage:
            fitter=NFWBiasFitter(omega_m, z, r, **linear_kwds)
            res = fitter.fit(dsig, dsigerr, guess, more=False)

        See the fit() method for a description of the outputs.

        """

        self.r_       = r.copy()

        self.z = z

        self.withlin = kwds.get('withlin',True)

        self.nfw = lensing.nfw.NFW(omega_m,z)

        if self.withlin:
            self.lin = lensing.linear.Linear(**kwds)
            print 'pre-computing linear dsig at %s points' % r.size
            self.lin_dsig_ = self.lin.dsig(self.r_)

    def fit(self, dsig, dsigerr, guess, more=False):
        """
        Class:
            NFWBiasFitter
        Method:
            fit
        Usage:
            fitter.fit(dsig, dsigerr, guess, more=False)

        Inputs:
            dsig,dsigerr: Must correspond to the radii given
                on construction.
            guess: [r200,c,B] or [r200,c] if withlin=False

        Outputs:
            if more=False:
                p,cov:  The parameters array and the covariance arra
            if more=True:
                a dict with  'p' and 'cov' along with
                    r200,r200_err
                    c, c_err
                    B, B_err  (-9999 if withlin=False)
                    m200, m200_err
            plus their errors

        """
        from scipy.optimize import curve_fit

        if dsig.size != self.r_.size or dsigerr.size != self.r_.size:
            raise ValueError("dsig,dsigerr must be same size as r")

        if self.withlin:
            npar = 3
            func = self.nfw_lin_dsig
        else:
            func = self.nfw.dsig
            npar = 2
        if len(guess) != npar:
            raise ValueError("parameter guess must have %s elements" % npar)

        print 'running curve_fit'
        p, cov = curve_fit(func, self.r_, dsig, sigma=dsigerr,
                           p0=guess)
        if not more:
            return p, cov
        else:
            return self.pars2more(p,cov)

    def nfw_lin_dsig(self, r, r200, c, B):
        if self.withlin:
            return self.nfw.dsig(r, r200, c) + self.lin_dsig(r,B)
        else:
            return self.nfw.dsig(r, r200, c)

    def lin_dsig(self, r, B):
        """
        r values must coinciced with linear dsig values
        """
        if r.size != self.lin_dsig_.size:
            raise ValueError("r and dsig must be same size")
        return B*self.lin_dsig_

    def pars2more(self,p,cov):
        r200 = p[0]
        r200_err = sqrt(cov[0,0])
        c = p[1]
        c_err = sqrt(cov[1,1])
        if self.withlin:
            B = p[2]
            B_err = sqrt(cov[2,2])
        else:
            B = -9999.0
            B_err = 9999.0

        m200 = self.nfw.m200(r200)
        m200_err = 3*m200*r200_err/r200

        res={'p':p,
             'cov':cov,
             'r200':r200,
             'r200_err':r200_err,
             'c':c,
             'c_err':c_err,
             'B':B,
             'B_err':B_err,
             'm200':m200,
             'm200_err':m200_err}
        return res




def fit_nfw_lin_dsig(omega_m, z, r, ds, dserr, guess, 
                     more=False, **kwds):
    """
    Name:
        fit_nfw_line_dsig
    Calling Sequence:
        fit_nfw_lin_dsig(omega_m, z, r, ds, dserr, guess, more=False,
                         **kwds)
    Inputs:
        omega_m:
        z: the redshift
        r: radius in Mpc
        ds: delta sigma in pc^2/Msun
        dserr: error
        guess: guesses for [r200,c,B]

    Keywords:
        withlin: True.  Include the linear bias term.
        more: Controls the output.  See Outputs section below.
        
        Keywords for the Linear() class. See the Linear() class for the
        defaults
            omega_b
            sigma_8
            h
            ns
    Outputs:
        if more=False:
            p,cov:  The parameters array and the covariance arra
        if more=True:
            a dict with  'p' and 'cov' along with
                r200,r200_err
                c, c_err
                B, B_err
                m200, m200_err
        plus their errors

    """
    fitter = NFWBiasFitter(omega_m, z, r, **kwds)
    return fitter.fit(ds, dserr, guess, more=more)

def fit_nfw_dsig(omega_m, z, r, ds, dserr, guess,rhofac=180):
    """
    Name:
        fit_nfw_dsig
    Calling Sequence:
        fit_nfw_dsig(omega_m, z, r, ds, dserr, guess)
    Inputs:
        omega_m: 
        z: the redshift
        r: radius in Mpc
        ds: delta sigma in pc^2/Msun
        dserr: error
        guess: guesses for [r200,c]

        rhofac=180 for calculation mass relative to
            rhofac times mean

    Outputs:
        This returns a dict with:

            r200
            c
            m200 

        plus errors and a covariance between r200 and c

        A mass relative to rhofac*rhmean is also computed.

    """
    from scipy.optimize import curve_fit

    n = lensing.nfw.NFW(omega_m, z)
    p, cov = curve_fit(n.dsig, r, ds, sigma=dserr, p0=guess)
    r200 = p[0]
    r200_err = sqrt(cov[0,0])
    c = p[1]
    c_err = sqrt(cov[1,1])

    rhoc = n.rhocrit

    m200 = n.m200(r200)
    m200_err = 3*m200*r200_err/r200

    # get rhofac times mean
    rm = n.r_fmean(r200, c, rhofac)
    rm_err = r200_err*(rm/r200)
    mm = n.m(rm, r200, c)
    mm_err = m200_err*(mm/m200)

    rtag = 'r%sm' % rhofac
    mtag = 'm%sm' % rhofac

    res={'r200':r200,
         'r200_err':r200_err,
         'c':c,
         'c_err':c_err,
         'cov':cov,
         'm200':m200,
         'm200_err':m200_err,
         rtag:rm,
         rtag+'_err':rm_err,
         mtag:mm,
         mtag+'_err':mm_err}
    return res




#
# fits byrun: these call the above functions
#

def fit_nfw_dsig_byrun(run, name, rrange=None, rhofac=180):
    """
    Fit an nfw profile to all bins
    """

    conf = lensing.files.read_config(run)
    din = lensing.files.lensbin_read(run,name)
    omega_m = conf['omega_m']

    rtag='r%01im' % rhofac
    mtag='m%01im' % rhofac
    newdt = [('rrange','f8',2),
             ('r200_fit','f8'),
             ('r200_fit_err','f8'),
             ('m200_fit','f8'),
             ('m200_fit_err','f8'),
             ('c_fit','f8'),
             ('c_fit_err','f8'),
             ('rc_fit_cov','f8',(2,2)),
             (rtag+'_fit','f8'),
             (rtag+'_fit_err','f8'),
             (mtag+'_fit','f8'),
             (mtag+'_fit_err','f8')]
    d = eu.numpy_util.add_fields(din, newdt)

    if rrange is None:
        rrange = [0,1.e6]

    r200guess = 1.0 # Mpc
    cguess = 5.0
    guess = numpy.array([r200guess,cguess],dtype='f8')
    for i in xrange(d.size):
        z = d['z_mean'][i]
        print 'omega_m:',omega_m
        print '      z:',z

        d['rrange'][i] = rrange
        r = d['r'][i]
        w=where1( (r > rrange[0]) & (r < rrange[1]) )

        r=r[w]
        ds = d['dsig'][i][w]
        dserr = d['dsigerr'][i][w]

        res = lensing.nfw.fit_nfw_dsig(omega_m, z, r, ds, dserr, guess,
                                       rhofac=rhofac)

        d['r200_fit'][i] = res['r200']
        d['r200_fit_err'][i] = res['r200_err']
        d['m200_fit'][i] = res['m200']
        d['m200_fit_err'][i] = res['m200_err']
        d['c_fit'][i] = res['c']
        d['c_fit_err'][i] = res['c_err']
        d['rc_fit_cov'][i] = res['cov']

        d[rtag+'_fit'][i] = res[rtag]
        d[rtag+'_fit_err'][i] = res[rtag+'_err']
        d[mtag+'_fit'][i] = res[mtag]
        d[mtag+'_fit_err'][i] = res[mtag+'_err']


        print '       c: %f +/- %f' % (d['c_fit'][i],d['c_fit_err'][i])
        print '    r200: %f +/- %f' % (d['r200_fit'][i],d['r200_fit_err'][i])
        print '    m200: %e +/- %e' % (d['m200_fit'][i],d['m200_fit_err'][i])
        print '   '+mtag+': %e +/- %e' % (d[mtag+'_fit'][i],d[mtag+'_fit_err'][i])
        if 'm200_mean' in d.dtype.names:
            m=d['m200_mean'][i]
            e=d['m200_err'][i]
            print '    m200 true: %e +/- %e' % (m,e)
        
    lensing.files.lensfit_write(d, run, name)



def fit_nfw_lin_dsig_byrun(run, name, withlin=True, rmax_from_true=False,
                           rmin=None, rmax=None):
    """
    Fit an nfw profile to all bins
    """

    conf = lensing.files.read_config(run)
    din = lensing.files.lensbin_read(run,name)
    omega_m = conf['omega_m']

    if withlin:
        npar = 3
        ex='lin'
    else:
        npar = 2
        ex=None
    newdt = [('rrange','f8',2),
             ('r200_fit','f8'),
             ('r200_fit_err','f8'),
             ('m200_fit','f8'),
             ('m200_fit_err','f8'),
             ('c_fit','f8'),
             ('c_fit_err','f8'),
             ('B_fit','f8'),
             ('B_fit_err','f8'),
             ('fit_cov','f8',(npar,npar))]
    d = eu.numpy_util.add_fields(din, newdt)

    rall = d['r'].ravel()
    if rmin is None:
        rmin = rall.min()
    if rmax is None:
        rmax = rall.max()
    rrange = [rmin,rmax]

    r200guess = 1.0 # Mpc
    cguess = 5.0
    Bguess = 5.0
    if withlin:
        guess = numpy.array([r200guess,cguess,Bguess],dtype='f8')
    else:
        guess = numpy.array([r200guess,cguess],dtype='f8')
    for i in xrange(d.size):
        z = d['z_mean'][i]
        print 'omega_m:',omega_m
        print '      z:',z

        r = d['r'][i]

        if rmax_from_true:
            # determine the rmax r200
            rrange = [rmin, 2*d['r200_mean'][i]]
            print '    using rrange:',rrange
        d['rrange'][i] = rrange
        w=where1( (r > rrange[0]) & (r < rrange[1]) )

        r=r[w]
        ds = d['dsig'][i][w]
        dserr = d['dsigerr'][i][w]

        res = lensing.fit.fit_nfw_lin_dsig(omega_m, z, r, ds, dserr, guess,
                                           withlin=withlin, more=True)

        d['r200_fit'][i] = res['r200']
        d['r200_fit_err'][i] = res['r200_err']
        d['m200_fit'][i] = res['m200']
        d['m200_fit_err'][i] = res['m200_err']
        d['c_fit'][i] = res['c']
        d['c_fit_err'][i] = res['c_err']
        d['B_fit'][i] = res['B']
        d['B_fit_err'][i] = res['B_err']
        d['fit_cov'][i] = res['cov']

        print '       c: %f +/- %f' % (d['c_fit'][i],d['c_fit_err'][i])
        print '       B: %f +/- %f' % (d['B_fit'][i],d['B_fit_err'][i])
        print '    r200: %f +/- %f' % (d['r200_fit'][i],d['r200_fit_err'][i])
        print '    m200: %e +/- %e' % (d['m200_fit'][i],d['m200_fit_err'][i])
        if 'm200_mean' in d.dtype.names:
            m=d['m200_mean'][i]
            e=d['m200_err'][i]
            print '    m200 true: %e +/- %e' % (m,e)
        
    lensing.files.lensfit_write(d, run, name, extra=ex)




def plot_nfw_lin_fits_byrun(run, name, npts=100, prompt=False, 
                            withlin=True,
                            ymin=0.01, ymax=2000.0):
    conf = lensing.files.read_config(run)
    if withlin:
        ex='lin'
        nex='lin'
    else:
        nex=''
        ex=None
    d = lensing.files.lensfit_read(run,name,extra=ex)
    omega_m = conf['omega_m']

    rravel = d['r'].ravel()
    xrange = [0.5*rravel.min(), 1.5*rravel.max()]

    #for i in xrange(d.size):
    i=0
    for dd in d:

        zrange = dd['z_range']
        mrange = dd['m200_range']

        if dd['rrange'][0] > 0:
            log_rmin = log10(dd['rrange'][0])
            log_rmax = log10(dd['rrange'][1])
        else:
            log_rmin = log10(dd['r'][0])
            log_rmax = log10(dd['r'][-1])
        rvals = 10.0**linspace(log_rmin,log_rmax,npts)

        plt = FramedPlot()  
        lensing.plotting.add_to_log_plot(plt, dd['r'],dd['dsig'],dd['dsigerr'])

        z = dd['z_mean']
        fitter = lensing.fit.NFWBiasFitter(omega_m,z,rvals,withlin=withlin)

        if withlin:
            yfit = fitter.nfw_lin_dsig(rvals, dd['r200_fit'],dd['c_fit'],dd['B_fit'])
            yfit_nfw = fitter.nfw.dsig(rvals,dd['r200_fit'],dd['c_fit'])
            yfit_lin = fitter.lin_dsig(rvals,dd['B_fit'])

            yfit = where(yfit < 1.e-5, 1.e-5, yfit)
            yfit_lin = where(yfit_lin < 1.e-5, 1.e-5, yfit_lin)

            cyfit = Curve(rvals,yfit,color='blue')
            cyfit_nfw = Curve(rvals,yfit_nfw,color='red')
            cyfit_lin = Curve(rvals,yfit_lin,color='orange')

            cyfit.label = 'Best Fit'
            cyfit_nfw.label = 'NFW'
            cyfit_lin.label = 'linear'

            key=PlotKey(0.1,0.3,[cyfit,cyfit_nfw,cyfit_lin])
            plt.add(cyfit,cyfit_nfw,cyfit_lin,key)
        else:
            yfit_nfw = fitter.nfw.dsig(rvals,dd['r200_fit'],dd['c_fit'])
            cyfit_nfw = Curve(rvals,yfit_nfw,color='blue')
            plt.add(cyfit_nfw)

        zlab='%0.2f < z < %0.2f' % (zrange[0],zrange[1])
        plt.add(PlotLabel(0.7,0.8,zlab))
        ll = (log10(mrange[0]),log10(mrange[1]))
        mlab = r'$%0.2f < logM_{200} < %0.2f$' % ll
        plt.add(PlotLabel(0.7,0.9,mlab))

        #yrange = [ymin,(dd['dsig']+dd['dsigerr']).max()*1.5]
        yrange = [ymin,ymax]
        plt.xrange = xrange
        plt.yrange = yrange
        plt.xlog=True
        plt.ylog=True
        plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        plt.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'
        plt.aspect_ratio=1
        if prompt:
            plt.show()
            rinput = raw_input('hit a key: ')
            if rinput == 'q':
                return
        else:
            d = lensing.files.lensbin_plot_dir(run,name)
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile=path_join(d,'desmocks-nfw%s-fit-%02i.eps' % (nex,i))
            print 'Writing epsfile:',epsfile
            plt.write_eps(epsfile)
        i += 1


    

def plot_nfwfits_byrun(run, name, prompt=False):
    conf = lensing.files.read_config(run)
    d = lensing.files.lensfit_read(run,name)
    omega_m = conf['omega_m']


    rvals = numpy.linspace(d['r'].min(), d['r'].max(),1000)
    for i in xrange(d.size):
        plt = FramedPlot()  
        lensing.plotting.add_to_log_plot(plt, 
                                          d['r'][i],
                                          d['dsig'][i],
                                          d['dsigerr'][i])

        z = d['z_mean'][i]
        n = lensing.nfw.NFW(omega_m, z)
        yfit = n.dsig(rvals, d['r200_fit'][i],d['c_fit'][i])
        plt.add(Curve(rvals,yfit,color='blue'))
        plt.xlog=True
        plt.ylog=True
        plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        plt.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'
        if prompt:
            plt.show()
            raw_input('hit a key: ')
        else:
            epsfile='/home/esheldon/tmp/plots/desmocks-nfwfit-%02i.eps' % i
            print 'Writing epsfile:',epsfile
            plt.write_eps(epsfile)







def test_fit_nfw_lin_dsig(rmin=0.01):

    from biggles import FramedPlot,Points,SymmetricErrorBarsY,Curve,PlotKey
    omega_m=0.25
    z=0.25


    r200 = 1.0
    c = 5.0
    B=10.0

    rmax = 50.0
    log_rmin = log10(rmin)
    log_rmax = log10(rmax)
    npts = 30
    r = 10.0**linspace(log_rmin,log_rmax,npts)

    fitter = NFWBiasFitter(omega_m,z,r)
    ds = fitter.dsig(r, r200, c, B)
    # 10% errors
    dserr = 0.1*ds
    ds += dserr*numpy.random.standard_normal(ds.size)

    guess = numpy.array([r200,c,B],dtype='f8')
    # add 10% error to the guess
    guess += 0.1*guess*numpy.random.standard_normal(guess.size)

    res = fitter.fit(ds,dserr,guess, more=True)

    r200_fit=res['r200']
    r200_err = res['r200_err']
    c_fit=res['c']
    c_err = res['c_err']
    B_fit=res['B']
    B_err = res['B_err']


    print 'Truth:'
    print '    r200: %f' % r200
    print '       c: %f' % c
    print '       B: %f' % B
    print 'r200_fit: %f +/- %f' % (r200_fit,r200_err)
    print '   c_fit: %f +/- %f' % (c_fit,c_err)
    print '   B_fit: %f +/- %f' % (B_fit,B_err)
    print 'Cov:'
    print res['cov']

 
    rfine = 10.0**linspace(log_rmin,log_rmax,100)
    fitter2 = NFWBiasFitter(omega_m,z,rfine)

    yfit = fitter2.dsig(rfine, r200_fit, c_fit, B_fit)
    yfit_nfw = fitter2.nfw.dsig(rfine, r200_fit, c_fit)
    yfit_lin = fitter2.lin_dsig(rfine,B_fit)

    plt=FramedPlot()
    plt.add(Points(r,ds,type='filled circle'))
    plt.add(SymmetricErrorBarsY(r,ds,dserr))

    cyfit = Curve(rfine,yfit,color='blue')
    cyfit_nfw = Curve(rfine,yfit_nfw,color='red')
    cyfit_lin = Curve(rfine,yfit_lin,color='orange')

    cyfit.label = 'Best Fit'
    cyfit_nfw.label = 'NFW'
    cyfit_lin.label = 'linear'

    key=PlotKey(0.1,0.3,[cyfit,cyfit_nfw,cyfit_lin])
    plt.add(cyfit,cyfit_nfw,cyfit_lin,key)

    plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
    plt.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

    plt.xrange = [0.5*rmin, 1.5*rmax]
    plt.yrange = [0.5*(ds-dserr).min(), 1.5*(ds+dserr).max()]

    plt.xlog=True
    plt.ylog=True
    plt.show()







def test_fit_nfw_dsig(rmin=0.01):

    from biggles import FramedPlot,Points,SymmetricErrorBarsY,Curve
    omega_m=0.25
    z=0.25
    n = lensing.nfw.NFW(omega_m, z)

    r200 = 1.0
    c = 5.0

    rmax = 5.0
    log_rmin = log10(rmin)
    log_rmax = log10(rmax)
    npts = 25
    logr = numpy.linspace(log_rmin,log_rmax,npts)
    r = 10.0**logr

    ds = n.dsig(r, r200, c)
    # 10% errors
    dserr = 0.1*ds
    ds += dserr*numpy.random.standard_normal(ds.size)

    guess = numpy.array([r200,c],dtype='f8')
    # add 10% error to the guess
    guess += 0.1*guess*numpy.random.standard_normal(guess.size)

    res = fit_nfw_dsig(omega_m, z, r, ds, dserr, guess)

    r200_fit = res['r200']
    r200_err = res['r200_err']

    c_fit = res['c']
    c_err = res['c_err']

    print 'Truth:'
    print '    r200: %f' % r200
    print '       c: %f' % c
    print 'r200_fit: %f +/- %f' % (r200_fit,r200_err)
    print '   c_fit: %f +/- %f' % (c_fit,c_err)
    print 'Cov:'
    print res['cov']

    logr = numpy.linspace(log_rmin,log_rmax,1000)
    rlots = 10.0**logr
    yfit = n.dsig(rlots,r200_fit,c_fit)

    plt=FramedPlot()
    plt.add(Points(r,ds,type='filled circle'))
    plt.add(SymmetricErrorBarsY(r,ds,dserr))
    plt.add(Curve(rlots,yfit,color='blue'))

    plt.xlabel = r'$r$ [$h^{-1}$ Mpc]'
    plt.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

    plt.xrange = [0.5*rmin, 1.5*rmax]
    plt.yrange = [0.5*(ds-dserr).min(), 1.5*(ds+dserr).max()]

    plt.xlog=True
    plt.ylog=True
    plt.show()



