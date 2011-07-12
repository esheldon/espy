import os
import sys
from sys import stdout
import numpy
from numpy import log10,sqrt,linspace,where
import lensing
import esutil as eu
from esutil.ostools import path_join
from esutil.stat import histogram
from esutil.numpy_util import where1

try:
    import biggles
    from biggles import FramedArray, FramedPlot, Points, \
            ErrorBarsY, ErrorBarsX, \
            SymmetricErrorBarsY, SymmetricErrorBarsX, \
            PlotKey, PlotLabel, Table, Curve
except:
    pass


def plot_mzbin(run, nmass, nz, dops=False, addfits=False):
    """
    Run is a run of the lensing code.

    This one bins by m200 and plots the average delta sigma
    in each bin.  We must generalize this

    For z binning, the different z bins are overplotted
    """

    conf = lensing.files.json_read(run)

    name = 'm%02iz%01i' % (nmass,nz)
    data=lensing.files.lensbin_read(run,name)

    biggles.configure('screen','width', 1140)
    biggles.configure('screen','height', 1140)
    biggles.configure('fontsize_min',1.0)

    if nmass == 12:
        nrow = 3
        ncol = 4
    else:
        raise ValueError("Unsupported nmass: %s" % nmass)

    pa = FramedArray( nrow, ncol)
    #pa.aspect_ratio = 1.0/1.61803399
    pa.aspect_ratio = 1.0/1.5

    pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
    pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'
    pa.xrange = [0.01,60.0]
    pa.yrange = [1.e-2,8000]
    pa.xlog = True
    pa.ylog = True


    fmt = "%-2i %-2i %-2i %8i %11.5f %11.5f\n"
    row = -1
    colors = ['blue','magenta','red']

    i = 0
    for mi in xrange(nmass):
        col = i % ncol
        if col == 0:
            row += 1

        mrange = data['m200_range'][i]

        pplist = []
        for zi in xrange(nz):
            zr = data['z_range'][i]
            pp = lensing.plotting.add_to_log_plot(pa[row,col],
                                                   data['r'][i],
                                                   data['dsig'][i],
                                                   data['dsigerr'][i],
                                                   color=colors[zi])
            pp.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            pplist.append(pp)
            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.2,pplist, fontsize=0.1)
            pa[row,col].add(key)


        ll = (log10(mrange[0]),log10(mrange[1]))
        lab = r'$%0.2f < logM_{200} < %0.2f$' % ll
        pl = PlotLabel(.5, .85, lab)
        pa[row,col].add(pl)


    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        epsfile = path_join(d, 'lensbin-%s-%s-allplot.eps' % (run,name))
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
    else:
        pa.show()

def add_ranges_to_plot(plt, data, tag, trange, zrange=None, 
                       color='black', minval=1.e-5):

    ehigh = comb['dsig'] + comb['dsigerr']
    elow  = comb['dsig'] - comb['dsigerr']
    elow = numpy.where(elow < minval, minval, elow)

    wg = where1(comb['dsig'] > minval)
    p = Points(comb['r'][wg], 
               comb['dsig'][wg],type='filled circle',color=color)
    plt.add(p)

    wg = where1(ehigh > minval)
    pe = ErrorBarsY(comb['r'][wg], elow[wg], ehigh[wg],color=color)
    plt.add(pe)

    return p

def mzbin_name(nmass,nz):
    name='m%02iz%01i' % (nmass,nz)
    return name
def mzbin_byrun(run, nmass, nz):
    """
    Do the binning and write out a file
    """
    name='m%02iz%01i' % (nmass,nz)
    d = lensing.files.lensout_collate(run)
    res = mzbin(d, nmass, nz)
    lensing.files.lensbin_write(res,run,name) 

def mzbin(data, nmass, nz):

    mlow,mhigh = m200bins(nmass)
    zlow,zhigh = zbins(nz)

    nrbin = data['rsum'][0].size
    dt = lensbin_dtype(nrbin, ['m200','r200','z'])

    bs = numpy.zeros(nmass*nz, dtype=dt)

    i=0
    for ml,mh in zip(mlow,mhigh):
        mrange = [ml,mh]

        ll = (log10(mrange[0]),log10(mrange[1]))
        print '%0.2f < log M200 < %0.2f' % ll
        for zl,zh in zip(zlow,zhigh):
            zrange = (zl,zh)
            print '    %0.2f < z < %0.2f' % zrange
            comb,w = combine_lensum_from_ranges(data,'m200',mrange,zrange=zrange,
                                                getind=True)
    
            bs['r'][i] = comb['r']
            bs['dsig'][i] = comb['dsig']
            bs['dsigerr'][i] = comb['dsigerr']
            bs['osig'][i] = comb['osig']
            bs['npair'][i] = comb['npair']

            mn,err,sdev = lensave(data,'z',ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = zrange

            mn,err,sdev = lensave(data,'m200',ind=w, sdev=True)
            bs['m200_mean'][i] = mn
            bs['m200_err'][i] = err
            bs['m200_sdev'][i] = sdev
            bs['m200_range'][i] = mrange

            mn,err,sdev = lensave(data,'r200',ind=w, sdev=True)
            bs['r200_mean'][i] = mn
            bs['r200_err'][i] = err
            bs['r200_sdev'][i] = sdev
            bs['r200_range'][i] = mrange



            i+=1

    return bs


def plot_mzbin_invmass_byrun(run, nmass, nz, type='',
                             residual=False,
                             yrange=None,
                             dops=False):

    name=mzbin_name(nmass,nz)
    conf = lensing.files.json_read(run)
    d = lensing.files.lensinv_read(run,name)

    plt = FramedPlot()
    colors = ['blue','magenta','red']
    i=0

    tag = 'm200'+type+'_inv'
    for mi in xrange(nmass):
        plist=[]
        for zi in xrange(nz):
            zr = d['z_range'][i]

            mt = d['m200_mean'][i]
            mterr = d['m200_err'][i]
            minv = d[tag][i]
            minverr = d[tag+'_err'][i]

            if residual:
                if True:
                    yp = (mt - minv)/mt
                    yperr = sqrt( (mterr/mt)**2 + (minverr/minv)**2 )
                else:
                    yp = mt - minv
                    yperr = sqrt(mterr**2 + minverr**2)
            else:
                yp = minv 
                yperr = minverr

            p=Points([mt],[yp],type='filled circle',color=colors[zi])
            p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            plist.append(p)

            plt.add(p)
            plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
            plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.9,plist)
            plt.add(key)

    plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
    if not residual:
        plt.yrange = [0.5*d[tag].min(), 1.5*d[tag].max()]
        plt.ylabel = r'$M_{200}^{inv} [h^{-1} M_{sun}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-inv)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
    plt.xlog = True
    if not residual:
        plt.ylog = True

    plt.aspect_ratio=1

    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'lensinv-m200'+type+'-residual-vs-true-%s-%s.eps' % (run,name)
        else:
            f = 'lensinv-m200'+type+'-vs-true-%s-%s.eps' % (run,name)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)

    else:
        plt.show()



def plot_mzbin_mass_byrun(run, nmass, nz, 
                          withlin=True,
                          residual=False,
                          noerr=False,
                          fudge=None,tag='m200',
                          yrange=None,
                          dops=False):
    if withlin:
        ex='lin'
        nex='-lin'
    else:
        nex=''
        ex=None

    name=mzbin_name(nmass,nz)
    conf = lensing.files.json_read(run)
    d = lensing.files.lensfit_read(run,name,extra=ex)

    plt = FramedPlot()
    colors = ['blue','magenta','red']
    i=0

    if fudge is not None:
        d[tag+'_fit'] *= fudge
        d[tag+'_fit_err'] *= fudge
    for mi in xrange(nmass):
        plist=[]
        for zi in xrange(nz):
            zr = d['z_range'][i]

            mt = d['m200_mean'][i]
            mterr = d['m200_err'][i]
            mfit = d[tag+'_fit'][i]
            mfiterr = d[tag+'_fit_err'][i]

            if residual:
                if True:
                    yp = (mt - mfit)/mt
                    yperr = sqrt( (mterr/mt)**2 + (mfiterr/mfit)**2 )
                else:
                    yp = mt - mfit
                    yperr = sqrt(mterr**2 + mfiterr**2)
            else:
                yp = mfit 
                yperr = mfiterr

            p=Points([mt],[yp],type='filled circle',color=colors[zi])
            p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            plist.append(p)

            plt.add(p)
            if not noerr:
                plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
                plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.9,plist)
            plt.add(key)

    plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
    if not residual:
        plt.yrange = [0.5*d[tag+'_fit'].min(), 1.5*d[tag+'_fit'].max()]
        plt.ylabel = r'$M_{200}^{fit} [h^{-1} M_{sun}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        #plt.ylabel = r'$M_{true}-M_{fit} [h^{-1} M_{sun}]$'
        #plt.ylabel = r'$(M_{200}^{true}-M_{200}^{fit})/M_{200}^{true}$'
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-fit)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
    plt.xlog = True
    if not residual:
        plt.ylog = True

    plt.aspect_ratio=1

    if not withlin:
        plt.add(PlotLabel(0.7,0.1,"no linear bias"))

    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'lensfit-m200-residual-vs-true-%s-%s%s.eps' % (run,name,nex)
        else:
            f = 'lensfit-m200-vs-true-%s-%s%s.eps' % (run,name,nex)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)

    else:
        plt.show()


def invert_byrun(run, name, show=False, rmax=None, massin=False, massout=False):
    """
    Invert for drho and mass

    For desmocks we want massin=True
    """

    conf = lensing.files.json_read(run)
    din = lensing.files.lensbin_read(run,name)
    omega_m = conf['omega_m']

    nr = din['r'][0].size
    ndrho = nr-1

    newdt = [('rdrho','f8',ndrho),
             ('drho','f8',ndrho),
             ('drhoerr','f8',ndrho),
             ('drhocov','f8',(ndrho,ndrho)),
             ('drho_noendc','f8',ndrho),
             ('rmass','f8',ndrho),
             ('mass','f8',ndrho),
             ('masserr','f8',ndrho),
             ('masscov','f8',(ndrho,ndrho)),
             ('massin','f8',ndrho),
             ('massinerr','f8',ndrho),
             ('massincov','f8',(ndrho,ndrho)),
             ('massout','f8',ndrho),
             ('massouterr','f8',ndrho),
             ('massoutcov','f8',(ndrho,ndrho)),
             ('r200_inv','f8'),
             ('r200_inv_err','f8'),
             ('m200_inv','f8'),
             ('m200_inv_err','f8'),
             ('r200in_inv','f8'),
             ('r200in_inv_err','f8'),
             ('m200in_inv','f8'),
             ('m200in_inv_err','f8'),
             ('r200out_inv','f8'),
             ('r200out_inv_err','f8'),
             ('m200out_inv','f8'),
             ('m200out_inv_err','f8')]
    d = eu.numpy_util.add_fields(din, newdt)

    pdir = lensing.files.lensbin_plot_dir(run,name)
    if not os.path.exists(pdir):
        os.makedirs(pdir)

    inv = lensing.invert.Inverter()
    for i in xrange(d.size):
        z = d['z_mean'][i]
        print 'omega_m:',omega_m
        print '      z:',z

        r = d['r'][i]
        dsig =d['dsig'][i]
        dsigerr =d['dsigerr'][i]

        #dsig = where(dsig < 1.e-3, 1.e-3, dsig)
        #w=where1(dsig < 1.e-5)
        #if w.size > 0:
        # replace with interpolation


        res = inv.invert(r, dsig, dserr=dsigerr)
        d['rdrho'][i]       = res['rdrho']
        d['drho'][i]        = res['drho']
        d['drhoerr'][i]     = res['drhoerr']
        d['drhocov'][i]     = res['drhocov']
        d['drho_noendc'][i] = res['drho_noendc']

        d['rmass'][i]       = res['rmass']
        d['mass'][i]        = res['mass']
        d['masserr'][i]     = res['masserr']
        d['masscov'][i]     = res['masscov']

        d['massin'][i]      = res['massin']
        d['massinerr'][i]      = res['massinerr']
        d['massincov'][i]   = res['massincov']
        d['massout'][i]     = res['massout']
        d['massouterr'][i]     = res['massouterr']
        d['massoutcov'][i]  = res['massoutcov']

        epsfile=path_join(pdir,'invert-%s-%s-%02i.eps' % (run,name,i))
        inv.plot(res, epsfile=epsfile, show=show)
        if show:
            k = raw_input('hit a key (q to quit): ')
            if k == 'q':
                return

        rmass   = d['rmass'][i]
        for mt in ['','in','out']:
            print "Getting '%s' virial mass" % mt
            mass = d['mass'+mt][i]
            masserr = d['mass'+mt+'err'][i]
            try:
                mvf = lensing.invert.MvirFinder(omega_m, z, 200)



                r200,r200_err,m200,m200_err=mvf.find_mvir(rmass,mass,masserr)

                d['r200'+mt+'_inv'][i] = r200
                d['r200'+mt+'_inv_err'][i] = r200_err
                d['m200'+mt+'_inv'][i] = m200
                d['m200'+mt+'_inv_err'][i] = m200_err

                print '    r200: %f +/- %f' % (r200,r200_err)
                print '    m200: %e +/- %e' % (m200,m200_err)
                if 'm200_mean' in d.dtype.names:
                    m=d['m200_mean'][i]
                    e=d['m200_err'][i]
                    print '    m200 true: %e +/- %e' % (m,e)
            except:
                d['r200'+mt+'_inv'][i] = -9999
                d['r200'+mt+'_inv_err'][i] = 9999
                d['m200'+mt+'_inv'][i] = -9999
                d['m200'+mt+'_inv_err'][i] = 9999


                print "Error getting '%s' virial mass:" % mt
                print sys.exc_info()

            
    lensing.files.lensinv_write(d, run, name)





def plot_nfw_lin_fits_byrun(run, name, npts=100, prompt=False, 
                            withlin=True,
                            ymin=0.01, ymax=2000.0):
    conf = lensing.files.json_read(run)
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


    
def plot_nfw_lin_fits_byrun_old(run, name, npts=100, prompt=False, 
                            withlin=True,
                            ymin=0.01, ymax=2000.0):
    conf = lensing.files.json_read(run)
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
    conf = lensing.files.json_read(run)
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

def nfw_lin_fits_byrun(run, name, withlin=True, rmax_from_true=False,
                       rmin=None, rmax=None):
    """
    Fit an nfw profile to all bins
    """

    conf = lensing.files.json_read(run)
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


def nfwfits_byrun(run, name, rrange=None, rhofac=180):
    """
    Fit an nfw profile to all bins
    """

    conf = lensing.files.json_read(run)
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

def combine_lensum_from_ranges(data, tag, trange, zrange=None,getind=False):
    logic = (data[tag] >= trange[0]) & (data[tag] < trange[1])
    if zrange is not None:
        logic = logic & \
            (data['z'] >= zrange[0]) & (data['z'] < zrange[1])
    w=where1(logic)

    comb = combine_lensout(data[w])

    if getind:
        return comb, w
    else:
        return comb
    


def zbins(nbin):
    if nbin == 3:
        low  = [0.0,  0.34, 0.53]
        high = [0.34, 0.53, 1.4]
    elif nbin == 2:
        low  = [0.0, 0.44]
        high = [0.44, 1.4]
    else:
        raise ValueError("Unsupported nbin: %s\n" % nbin)
    
    return low,high

def m200bins(nbin):
    if nbin == 12:
        # I ran test_logbin with nbin=14 and combined last three
        low = [5.03e+12,
               7.28568896657e+12,
               1.05369194771e+13,
               1.52390079477e+13,
               2.20393981122e+13,
               3.18744547422e+13,
               4.60983943365e+13,
               6.66697509837e+13,
               9.64210524077e+13,
               1.39448838645e+14,
               2.01677726116e+14,
               2.91676184662e+14]

        high = [7.28568896657e+12,
                1.05369194771e+13,
                1.52390079477e+13,
                2.20393981122e+13,
                3.18744547422e+13,
                4.60983943365e+13,
                6.66697509837e+13,
                9.64210524077e+13,
                1.39448838645e+14,
                2.01677726116e+14,
                2.91676184662e+14,
                8.83e+14]
    else:
        raise ValueError("Unsupported nbin: %s\n" % nbin)

    low=numpy.array(low,dtype='f8')
    high=numpy.array(high,dtype='f8')

    return low,high

def test_logbin(data, nbin, tag, nzbin=0):
    """
    Do a log binning.  For binning by mass estimators, 
    probably want to combine the last 2 or 3 bins because
    of the exponential cutoff
    """

    log_data = numpy.log10(data[tag])

    # we can't just take the returned means because we want
    # the mean in the linear mass
    hdict = histogram(log_data, nbin=nbin, more=True, rev=True)

    rev = hdict['rev']
    mean_data = numpy.zeros(nbin,dtype='f8')
    stdout.write("Getting mean of '%s'\n" % tag)
    for i in xrange(nbin):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]
            mean_data[i],err = lensave(data, tag, ind=w)

    h = hdict['hist']
    low=10.0**hdict['low']
    high=10.0**hdict['high']
    fmt = "%-2i %8i %11.5e %11.5e %11.5e\n"
    for i in xrange(nbin):
        stdout.write(fmt % (i+1,h[i],low[i],high[i],mean_data[i]))

    print eu.numpy_util.arr2str(low)
    print eu.numpy_util.arr2str(high)
    
def lensave(data, tag, ind=None, sdev=False):
    if ind is None:
        wts = data['weight']
        tdata = data[tag]
    else:
        wts = data['weight'][ind]
        tdata = data[tag][ind]

    return eu.stat.wmom(tdata, wts, calcerr=True, sdev=sdev)


def combine_lensout(lout):

    nlens = lout.size
    nbin = lout['rsum'][0].size
    dt = lenscomb_dtype()

    comb=numpy.zeros(nbin,dtype=dt)

    for i in xrange(nbin):
        npair = lout['npair'][:,i].sum()
        rsum  = lout['rsum'][:,i].sum()

        comb['npair'][i] = npair
        comb['r'][i] = rsum/npair

        wsum = lout['wsum'][:,i].sum()
        dsum = lout['dsum'][:,i].sum()
        osum = lout['osum'][:,i].sum()

        comb['dsig'][i] = dsum/wsum
        comb['osig'][i] = osum/wsum

        comb['dsigerr'][i] = numpy.sqrt(1.0/wsum)

    return comb

def lenscomb_dtype():
    dt=[('r','f8'),
        ('dsig','f8'),
        ('dsigerr','f8'),
        ('osig','f8'),
        ('npair','i8')]
    return numpy.dtype(dt)

def lensbin_dtype(nrbin, bintags):
    if not isinstance(bintags,list):
        bintags = [bintags]

    dt = []
    for bt in bintags:
        tn = bt+'_range'
        dt.append( (tn,'f8',2) )
        tn = bt+'_mean'
        dt.append( (tn,'f8') )
        tn = bt+'_err'
        dt.append( (tn,'f8') )
        tn = bt+'_sdev'
        dt.append( (tn,'f8') )

    nrbin = int(nrbin)
    dt += [('r','f8',nrbin),
           ('dsig','f8',nrbin),
           ('dsigerr','f8',nrbin),
           ('osig','f8',nrbin),
           ('npair','i8',nrbin)]
    return numpy.dtype(dt)


