"""
Routines to deproject \Delta\Sigma to the 3-d \delta \rho and mass.

Ported from Dave Johnston's IDL code.
"""

import os
import sys
from sys import stdout
import lensing
import numpy
from numpy import log, log10, exp, sqrt, linspace, logspace, where, \
        arange, zeros, \
        isfinite, \
        pi as PI, \
        dot
import esutil as eu
from esutil.numpy_util import where1
from esutil.misc import colprint
from esutil.ostools import path_join


try:
    import biggles
    from biggles import FramedArray, FramedPlot, Points, \
            ErrorBarsY as ErrY, ErrorBarsX as ErrX, \
            SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX, \
            PlotKey, PlotLabel, Table, Curve
except:
    pass

try:
    import scipy.integrate
    from scipy.integrate import romberg
except:
    pass


def invert_byrun(run, name, show=False, rmax=None):
    """
    Invert for drho and mass

    For desmocks we want to try the massin since there are problems at large
    scales.  the shape of the mass profile is wrong on small scales, but total
    mass might be OK.

    """

    conf = lensing.files.cascade_config(run)
    din = lensing.files.lensbin_read(run,name)
    omega_m = conf['cosmo_config']['omega_m']

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

            
    lensing.files.sample_write(d, 'invert', run, name=name)





def test_invert(relerr=0.2,dowrite=False, epsfile=None):
    inv = Inverter()
    inv.test(relerr=relerr, dowrite=dowrite, epsfile=epsfile)

def invert(r, ds, dserr=None, dscov=None, 
           nsigo=1, 
           endslope=1.0, 
           endslope_npts=3,
           noendc=False,
           nocov=False):

    """
    
    Convenience function. Loads an Inverter() instance and performs the
    inversion to drho and to mass. Uses the special formulas for getting
    mass instead of integrating drho, as the latter has a hole in the
    middle.

    """

    inv=Inverter(nsigo=nsigo,
                 endslope=endslope,
                 endslope_npts=endslope_npts,
                 noendc=noendc)
    return inv.invert(r, ds, dserr=dserr, dscov=dscov, nocov=nocov)


class Inverter:
    def __init__(self, 
                 nsigo=1, 
                 endslope=1.0, 
                 endslope_npts=3,
                 noendc=False):
        self.nsigo=nsigo
        self.endslope=endslope
        self.endslope_npts=endslope_npts
        self.noendc=noendc

    def invert(self, r, ds, dserr=None, dscov=None, nocov=False):
        """

        Perform the inversion to drho and to mass. Uses the special formulas
        for getting mass instead of integrating drho, as the latter has a
        hole in the middle.

        """

        dserr,dscov=geterrcov(dserr,dscov)

        # first get drho and then mass comes from a combination
        # of delta sigma and delta rho
        print 'Getting drho'
        res = self.ds2drho(r, ds, dscov=dscov, nocov=nocov)

        rdrho = res['rdrho']
        drho = res['drho']
        drhocov = res['drhocov']
        drhoerr = res['drhoerr']

        print 'Getting mass'
        res2 = self.dsdrho2mass(r,ds,dscov,rdrho,drho,drhocov,nocov=nocov)
        for key in res2:
            res[key] = res2[key]
        return res

    def ds2drho(self, r, ds, dserr=None, dscov=None, nocov=False):
        """
        """

        #if r.max() > 100:
        #    raise ValueError("I have problems with large radius")
        nr = r.size

        dserr,dscov=geterrcov(dserr,dscov)

        # has to do with where to do final analytic integration start this
        # number times rmax, go to infinity. A tiny correction.
        nrmax=30.0     

        # get parameters for hybrid interpolation.
        slope = zeros(nr,dtype='f8')
        amp   = zeros(nr,dtype='f8')
        off   = zeros(nr,dtype='f8')

        for i in xrange(nr-1):
            amp[i],slope[i],off[i] = interp_hybrid_pars(r[i],r[i+1],
                                                        ds[i],ds[i+1],
                                                        dserr[i],dserr[i+1],
                                                        nsig=self.nsigo)
            #if not nocov: 
            #    print r[i],ds[i],ds[i+1],amp[i],slope[i],off[i]
            if not isfinite(slope[i]):
                raise ValueError("not finite slope found at i=%d" % i)

        # Dave kept the negative of the slope
        slope = -slope
        slope[-1] = self.endslope

        # get the last few points to fix the endpoint powerlaw amplitude
        amp[-1] = get_endpoint_amp(r,ds,dserr,
                                   self.endslope_npts,self.endslope)

        # now we know the local PL representation of DS including endpoint 
        # extrapolation model. Perfrom the numerical
        # integrations for each bin and add them up to get drho.
        # Use Johnston et al. (2004) formulas.

        drho_noendc = zeros(nr-1)
        endc = zeros(nr-1)

        for j in xrange(nr-1):
            for i in xrange(j, nr-1):
                # numerical integration is in here

                integ1 = xi_pl_int(r[i], r[i+1], r[j], slope[i]+1)

                integ2 = xi_pl_int(r[i], r[i+1], r[j], 0.0+1)

                drho_1  =  amp[i]*(2.0-slope[i])*integ1
                drho_2  = -off[i]*(2.0-0.0)*integ2
                drho_noendc[j] = drho_noendc[j] + drho_1+drho_2

            if not self.noendc:
                rmax = r.max()
                rrmax=nrmax*rmax
                # numerical integration except for last bit to infinity

                integ=xi_pl_int(rmax,rrmax,r[j],slope[-1]+1)

                endc[j]=amp[-1]*(2.0-slope[-1])*integ

                # last bit to infinity. tiny correction for typical
                # slopes

                endc[j]=endc[j]+\
                    amp[-1]*(2.0-slope[1])/PI * rrmax**(-(slope[-1]+1))/(slope[-1]+1)

        drho = drho_noendc + endc

        if not nocov:

            print 'Getting drhocov'
            drhocov = self.calc_drhocov(r, ds, dscov)

            drhoerr = sqrt( numpy.diag(drhocov) )
        else:
            drhocov=None
            drhoerr=None


        out={}
        out['ds_slope'] = slope
        out['ds_amp'] = amp
        out['ds_offset'] = off
        out['r'] = r.copy()
        out['dsig'] = ds.copy()
        out['dsigerr'] = dserr.copy()
        out['dsigcov'] = dscov.copy()
        out['rdrho'] = r[0:nr-1].copy()
        out['drho_noendc'] = drho_noendc
        out['drho'] = drho
        out['drhoerr'] = drhoerr
        out['drhocov'] = drhocov
        out['endc'] = endc
        return out

    def calc_drhocov(self, r, ds, dscov):

        dserr = sqrt( numpy.diagonal(dscov) )
        # the jacobian
        nr = ds.size
        mat = zeros( (nr,nr-1) )
        dsp = ds.copy()
        dsm = ds.copy()

        # fill in the jacobian by symmetric numerical differentiation.
        for i in xrange(nr):
            # a small number
            h = abs( 0.005*( abs(ds[i])+abs(dserr[i]) ) )

            dsp[i] += h
            dsm[i] -= h

            resp = self.ds2drho(r, dsp, dscov=dscov, nocov=True)
            resm = self.ds2drho(r, dsm, dscov=dscov, nocov=True)

            # the derivative as a difference.
            drhop = resp['drho']
            drhom = resm['drho']
            der=(drhop-drhom)/(2*h)
            mat[i,:]=der

            # puts back the offset values to originals
            dsp[i]=ds[i]
            dsm[i]=ds[i]

        # now propagate the covar
        # in IDL: drhocov=mat##(dscov##transpose(mat))

        tmp = dot(dscov,mat)
        drhocov = dot(tmp.T, mat)

        drhoerr2 = numpy.diag(drhocov)
        w=where1(drhoerr2 <= 0)
        if w.size > 0:
            raise ValueError("Some drho diagonals are negative")
        return drhocov



    def dsdrho2mass(self, r, ds, dscov, rdrho, drho, drhocov, nocov=False):

        """
        You usually won't call this directly.
        """

        ndrho = rdrho.size
        drhoerr = sqrt(numpy.diag(drhocov))

        # get parameters for hybrid interpolation over drho.
        slope = zeros(ndrho,dtype='f8')
        amp   = zeros(ndrho,dtype='f8')
        off   = zeros(ndrho,dtype='f8')

        for i in xrange(ndrho-1):
            amp[i],slope[i],off[i] = interp_hybrid_pars(rdrho[i],rdrho[i+1],
                                                        drho[i],drho[i+1],
                                                        drhoerr[i],drhoerr[i+1],
                                                        nsig=self.nsigo)
            if not isfinite(slope[i]):
                raise ValueError("not finite slope found at i=%d" % i)

        # Dave kept the negative of the slope
        slope = -slope

        # drho is steeper than delta sigma: slope+1
        slope[-1] = self.endslope +1

        # get the last few points to fix the endpoint powerlaw amplitude
        amp[-1] = get_endpoint_amp(rdrho,drho,drhoerr,
                                   self.endslope_npts,(self.endslope+1))

        mass_in  = zeros(ndrho,dtype='f8')
        mass_out = mass_in.copy()

        mass_in_c  = mass_in.copy()
        mass_out_c = mass_out.copy()


        # do numerical integrations for each interval and add them up
        # using Johnston et al. (2004) Formulas.

        imin = 0
        Rs = rdrho[imin]

        # mass_in with no correction yet
        for j in xrange(1,ndrho):
            # this will only use up to j-1, so should be fine
            # to subscript rdrho and amp,slope,off with i+1
            for i in xrange(j):
                x0=sqrt((rdrho[i]/Rs)**2-1)
                x1=sqrt((rdrho[i+1]/Rs)**2-1)
                integ1 = mass_pl_int(x0,x1,slope[i],'in')

                # second intergral can be done analytically 
                integ2=(x1-x0)+ 2*(x1**3-x0**3)/3.0
                mass_in[j] += 2*amp[i]*integ1*Rs**(-slope[i]) - 2*off[i]*integ2

        # mass_out and corrections
        for j in xrange(ndrho):
            for i in xrange(j,ndrho-1):
                x0 = sqrt( (rdrho[i]/rdrho[j])**2 - 1)
                x1 = sqrt( (rdrho[i+1]/rdrho[j])**2 - 1)

                integ1=mass_pl_int(x0,x1,slope[i],'out')
                integ2=mass_pl_int(x0,x1,0,'out')

                x0_c=sqrt( (rdrho[i]/Rs)**2 - 1)
                x1_c=sqrt( (rdrho[i+1]/Rs)**2 - 1)

                integ1_c = mass_pl_int(x0_c,x1_c,slope[i],'out')
                # second part analytically
                integ2_c = \
                    (x1_c-x0_c)  \
                    + 2*(x1_c**3-x0_c**3)/3.0 \
                    - 2*((1+x1_c**2)**(3/2.0) \
                    - (1+x0_c**2)**(3/2.0))/3.0

                mass_out[j] += \
                    2*amp[i]*integ1*rdrho[j]**(-slope[i]) - 2*off[i]*integ2
                mass_in_c[j] += \
                    2*amp[i]*integ1_c*Rs**(-slope[i]) - 2*off[i]*integ2_c


            # the infinity correction.
            # corr for mass_out
            x0=sqrt((rdrho[-1]/rdrho[j])**2-1)
            intoc=mass_pl_int(x0,(x0+1)*20.,slope[i],'out')
            mass_out_c[j]=2*amp[i]*intoc*rdrho[j]**(-slope[i])

            # corr for mass_in
            x0_c=sqrt((rdrho[-1]/Rs)**2-1)
            int_c=mass_pl_int(x0_c,(x1_c+1)*20.,slope[i],'out')
            mass_in_c[j] += 2*amp[i]*int_c*Rs**(-slope[i])

        mass_in =  (mass_in+ds[imin]/Rs)*PI*Rs**3
        mass_out =  (mass_out+mass_out_c+ds[0:ndrho]/rdrho)*PI*rdrho**3
        mass_in_c = mass_in_c*PI*Rs**3
        mass_in =   mass_in+mass_in_c

        fac = 1.e12
        mass_in  *= fac
        mass_out *= fac

        mass = (mass_in + mass_out)/2.0

        if not nocov:
            print 'Getting masscov'
            mincov, moutcov = self.calc_masscov(r, ds, dscov, 
                                                rdrho, drho, drhocov)
            masscov = (mincov+moutcov)/2.0
            masserr = sqrt( numpy.diag(masscov) )

        else:
            mincov=None
            moutcov=None
            masscov=None
            masserr=None

        out={}
        out['drho_slope'] = slope
        out['drho_amp'] = amp
        out['drho_offset'] = off
        out['rmass'] = rdrho.copy()
        out['massin'] = mass_in
        if mincov is not None:
            out['massinerr'] = sqrt( numpy.diag(mincov) )
        out['massincov'] = mincov
        out['massout'] = mass_out
        if moutcov is not None:
            out['massouterr'] = sqrt( numpy.diag(moutcov) )
        out['massoutcov'] = moutcov
        out['mass'] = mass
        out['masscov'] = masscov
        out['masserr'] = masserr

        return out

    def calc_masscov(self, 
                     r, ds, dscov, 
                     rdrho, drho, drhocov):

        dserr = sqrt( numpy.diagonal(dscov) )

        # the jacobians
        nr = ds.size
        mati = zeros( (nr,nr-1) )
        mato = zeros( (nr,nr-1) )
        dsp = ds.copy()
        dsm = ds.copy()

        # fill in the jacobians by symmetric numerical differentiation.
        print
        for i in xrange(nr):
            stdout.write('.');stdout.flush()
            # a small number
            h = abs( 0.005*( abs(ds[i])+abs(dserr[i]) ) )

            dsp[i] += h
            dsm[i] -= h

            # first get drho with altered values
            rhoresp = self.ds2drho(r, dsp, dscov=dscov, nocov=True)
            rhoresm = self.ds2drho(r, dsm, dscov=dscov, nocov=True)


            # now get mass with altered values
            drhop = rhoresp['drho']
            drhom = rhoresm['drho']
            mresp = self.dsdrho2mass(r,dsp,dscov,rdrho,drhop,drhocov,
                                     nocov=True)
            mresm = self.dsdrho2mass(r,dsm,dscov,rdrho,drhom,drhocov,
                                     nocov=True)

            deri = (mresp['massin'] - mresm['massin'])/(2*h)
            dero = (mresp['massout']- mresm['massout'])/(2*h)
            mati[i,:] = deri
            mato[i,:] = dero

            # puts back the offset values to originals
            dsp[i]=ds[i]
            dsm[i]=ds[i]
        print


        tmp = dot(dscov,mati)
        mincov = dot(tmp.T, mati)
        minerr2 = numpy.diag( mincov )
        w=where1(minerr2 <= 0)
        if w.size > 0:
            raise ValueError("mass_in_cov has negative diagonals")

        tmp = dot(dscov,mato)
        moutcov = dot(tmp.T, mato)
        mouterr2 = numpy.diag( moutcov )
        w=where1(mouterr2 <= 0)
        if w.size > 0:
            raise ValueError("mass_out_cov has negative diagonals")

        return mincov, moutcov

    def plot(self, odict, epsfile=None, show=True, overdrho=True):
        """
        odict is the output of invert()
        """

        tab = Table(2,2)
        od = lensing.plotting.plot_dsig(odict, noshow=True)
        # show some log-interpolated points of dsig
        nper = 10
        r = odict['r']
        for i in xrange(r.size-1):
            amp = odict['ds_amp'][i]
            slope = odict['ds_slope'][i]
            off = odict['ds_offset'][i]
            x = logspace(log10(r[i]),log10(r[i+1]),nper)
            y = amp*x**(-slope) - off
            od['plt'].add(Curve(x,y,color='blue'))

        tab[0,0] = od['plt']

        # plot uncorrected with no error bar
        rod = lensing.plotting.plot_drho(r=odict['rdrho'],
                                         drho=odict['drho_noendc'],
                                         drhoerr=odict['drhoerr']*0,
                                         noshow=True)
        rod2 = lensing.plotting.add_to_log_plot(rod['plt'],
                                                odict['rdrho'],
                                                odict['drho'],
                                                odict['drhoerr'],
                                                color='blue')
        pnocorr = rod['p']
        p = rod2['p']
        pnocorr.label = 'no end correction'
        p.label = 'end correction'

        rod['plt'].add(PlotKey(0.1,0.2,[pnocorr,p]))
        tab[0,1] = rod['plt']


        # for the mass stuff
        # first show log-interpolated points of drho to see how
        # it looks

        rod_interp = lensing.plotting.plot_drho(odict, noshow=True)
        nper = 10
        for i in xrange(odict['drho'].size-1):
            amp = odict['drho_amp'][i]
            slope = odict['drho_slope'][i]
            off = odict['drho_offset'][i]
            x = logspace(log10(r[i]),log10(r[i+1]),nper)
            y = amp*x**(-slope) - off
            y = where(y < 1.e-5, 1.e-5, y)
            rod_interp['plt'].add(Curve(x,y,color='blue'))

        tab[1,0] = rod_interp['plt']
       

        mod = lensing.plotting.plot_mass(odict, noshow=True)
        mp = mod['p']
        mp.label = 'Mass Avg'

        odin=lensing.plotting.add_to_log_plot(mod['plt'], 
                                              odict['rmass'], 
                                              odict['massin'], 
                                              odict['massinerr'],
                                              color='green', type='filled diamond')
        pin = odin['p']
        pin.label = 'Mass In'
        mod['plt'].add(pin)

        odout=lensing.plotting.add_to_log_plot(mod['plt'], 
                                               odict['rmass'], 
                                               odict['massout'], 
                                               odict['massouterr'],
                                               color='magenta', type='filled diamond')
        
        pout = odout['p']
        pout.label = 'Mass Out'

        rhoc0 = 0.27752
        mdel = (4./3.)*PI*200*rhoc0*1.e12*odict['rmass']**3

        pdel = Curve(odict['rmass'],mdel,color='blue')
        pdel.label = r'$200 \rho_{c} V$'
        mod['plt'].add(pdel)
        mod['plt'].add(PlotKey(0.1,0.9,[mp,pin,pout,pdel]))

        tab[1,1] = mod['plt']


        if epsfile is not None:
            print 'Writing to epsfile:',epsfile
            tab.write_eps(epsfile)
        if show:
            tab.show()

 
    def test(self, relerr=0.2,dowrite=False, epsfile=None):
        omega_m=0.25
        z = 0.3
        nfw = lensing.nfw.NFW(omega_m, z)
        lin = lensing.linear.Linear(omega_m=omega_m)

        r200=0.5
        c=5.0
        B=1.0

        print 'generating fake data with %i%% errors' % (relerr*100,)
        nr = 21
        logrmin=-2.0
        logrmax=1.5
        r = logspace(logrmin,logrmax,nr)
        ds = nfw.dsig(r,r200,c) + B*lin.dsig(r)
        dserr = ds*relerr
        ds += dserr*numpy.random.standard_normal(nr)

        if dowrite:
            f='~/tmp/fake-dsig.dat'
            print 'Writing fake data to:',f
            colprint(r,ds,dserr,format='%0.8f', file=f)

        print 'generating fine grid truth values'
        rtrue = logspace(logrmin,logrmax,nr*5)
        dstrue  = nfw.dsig(rtrue,r200,c) + B*lin.dsig(rtrue)
        drhotrue = nfw.rho(rtrue,r200,c)  + B*lin.drho(rtrue)
        # add linear term
        mnfw = nfw.m(rtrue, r200,c)
        mlin = lin.m(rtrue, B)
        mtrue = mnfw + mlin

        print 'Doing inversion'
        odict = self.invert(r,ds,dserr=dserr)

        tab = Table(2,2)

        od = lensing.plotting.plot_dsig(odict, noshow=True)

        # show some log-interpolated points of dsig
        nper = 10
        for i in xrange(nr-1):
            amp = odict['ds_amp'][i]
            slope = odict['ds_slope'][i]
            off = odict['ds_offset'][i]
            x = logspace(log10(r[i]),log10(r[i+1]),nper)
            y = amp*x**(-slope) - off
            od['plt'].add(Curve(x,y,color='blue'))

        od['plt'].add(Curve(rtrue,dstrue,color='red'))
        tab[0,0] = od['plt']

        # plot uncorrected with no error bar
        rod = lensing.plotting.plot_drho(r=odict['rdrho'],
                                         drho=odict['drho_noendc'],
                                         drhoerr=odict['drhoerr']*0,
                                         noshow=True)
        rod2 = lensing.plotting.add_to_log_plot(rod['plt'],
                                                odict['rdrho'],
                                                odict['drho'],
                                                odict['drhoerr'],
                                                color='blue')
        
        pnocorr = rod['p']
        p = rod2['p']
        pnocorr.label = 'no end correction'
        p.label = 'end correction'

        plist = [pnocorr,p]
        drhotrue_curve = Curve(rtrue,drhotrue,color='red')
        drhotrue_curve.label = 'Truth'
        rod['plt'].add(drhotrue_curve)
        plist.append(drhotrue_curve)

        rod['plt'].add(PlotKey(0.1,0.2,plist))


        tab[0,1] = rod['plt']



        # for the mass stuff
        # first show log-interpolated points of drho to see how
        # it looks

        rod_interp = lensing.plotting.plot_drho(odict, noshow=True)
        nper = 10
        for i in xrange(odict['drho'].size-1):
            amp = odict['drho_amp'][i]
            slope = odict['drho_slope'][i]
            off = odict['drho_offset'][i]
            x = logspace(log10(r[i]),log10(r[i+1]),nper)
            y = amp*x**(-slope) - off
            y = where(y < 1.e-5, 1.e-5, y)
            rod_interp['plt'].add(Curve(x,y,color='blue'))

        tab[1,0] = rod_interp['plt']
       
        mod = lensing.plotting.plot_mass(odict, noshow=True)
        
        mp = mod['p']
        mp.label = 'Recovered Mass'

        nfw_curve = Curve(rtrue,mnfw,color='green')
        nfw_curve.label = 'NFW'
        mod['plt'].add(nfw_curve)

        lin_curve = Curve(rtrue,mlin,color='orange')
        lin_curve.label = 'linear'
        mod['plt'].add(lin_curve)

        true_curve = Curve(rtrue,mtrue,color='magenta')
        true_curve.label = 'True'
        mod['plt'].add(true_curve)

        #m200 = nfw.m200(r200)
        #m200p = Points([r200],[m200],type='cross',color='blue')
        #m200p.label = r'$M_{200} true'

        #mod['plt'].add(m200p)

        # now find r200 from the mass profile
        #mvf = MvirFinder(omega_m, z, 200)
        #rvir,rvir_err,mvir,mvir_err = mvf.find_mvir(odict['rmass'],
        #                                            odict['mass'],
        #                                            odict['masserr'])

        #mvp = Points([rvir],[mvir],type='filled square',color='red')
        #mvexp = SymErrX([rvir],[mvir],[rvir_err],color='red')
        #mveyp = SymErrY([rvir],[mvir],[mvir_err],color='red')
        #mvp.label = r'$M_{200}'

        #mod['plt'].add(mvp,mvexp,mveyp)

        #mod['plt'].add(PlotKey(0.1,0.9,[mp,true_curve,nfw_curve,lin_curve,
        #                                m200p,mvp]))
        mod['plt'].add(PlotKey(0.1,0.9,[mp,true_curve,nfw_curve,lin_curve]))


        tab[1,1] = mod['plt']


        if dowrite:
            f='~/tmp/fake-dsig-invert.dat'
            print 'Writing inverted fake data to:',f
            colprint(odict['rdrho'], odict['drho'],odict['mass'],odict['masserr'],
                     format='%20.8e',file=f)

        if epsfile is not None:
            print 'Writing eps file:',epsfile
            tab.write_eps(epsfile)
        else:
            tab.show()


    def test_masserr(self, 
                     nrad=35, 
                     nr200=8, 
                     nb=8, 
                     c=5.0):
        """

        This is a bit of a MISNOMER, as the mass is the mass whatever
        you call it. If in the sim you count particles, you will get
        part of this "two halo term" bias there too.
        
        But the original purpose was
            
            Create profiles without noise and see how biased our answer
            is as a function of the bias.

            Use lots of points to avoid interpolation error

        """

        outfile='~/tmp/invert-masserr-r%02i-b%02i.rec' % (nr200,nb)
        print 'Will write to output:',outfile

        omega_m=0.25
        z=0.3
        nfw = lensing.nfw.NFW(omega_m, z)
        lin = lensing.linear.Linear(omega_m=omega_m)

        r200vals = linspace(0.3,1.3,nr200)
        Bvals = linspace(0.5,1.4,nb)

        dtype = [('r200','f8'),
                 ('m200','f8'),
                 ('B','f8',nb),
                 ('m200meas','f8',nb)]
        out = zeros(nr200, dtype=dtype)

        # lots of r values so interpolation isn't the
        # source of error
        logrmin=-2.0
        logrmax=1.5
        r = logspace(logrmin,logrmax,nrad)
        mvf = MvirFinder(omega_m, z, 200)

        i=1
        for ri in xrange(nr200):
            r200 = r200vals[ri]
            m200 = nfw.m200(r200)

            out['r200'][ri] = r200
            out['m200'][ri] = m200

            for bi in xrange(Bvals.size):
                print '%s/%s' % (i, nr200*nb)
                i+=1
                B = Bvals[bi]
                out['B'][ri] = B

                ds = nfw.dsig(r,r200,c) + B*lin.dsig(r)
                # fake this just so the code will work
                dserr = ds*0.01

                odict = self.invert(r,ds,dserr=dserr)
                rvir,rvir_err,mvir,mvir_err = mvf.find_mvir(odict['rmass'],
                                                            odict['mass'],
                                                            odict['masserr'])
                out['m200meas'][ri][bi] = mvir
                print 'r200:',r200,'  B:',B,'  mtrue:',m200,'  m200 recover:',mvir

        print 'writing data:',outfile
        hdr={'omega_m':omega_m,
             'z':z}
        eu.io.write(outfile,out,header=hdr)

    def plot_masserr(self, filename):
        import pcolors
        data=eu.io.read(filename)
        nr200 = data.size
        nb = data['B'][0].size

        plt = FramedPlot()
        colors = pcolors.rainbow(nr200, 'hex')
        clist = []
        for ri in xrange(nr200):
            r200 = data['r200'][ri]
            m200 = data['m200'][ri]

            p = Points(data['B'][ri], (m200-data['m200meas'][ri])/m200, 
                       type='filled circle', color=colors[ri])
            crv = Curve(data['B'][ri], (m200-data['m200meas'][ri])/m200, 
                        color=colors[ri])
            crv.label = 'r200: %0.2f' % r200
            clist.append(crv)
            plt.add(p,crv)

        key = PlotKey(0.1,0.4,clist)
        plt.add(key)
        plt.xlabel=r'$\Omega_m \sigma_8^2 D(z)^2 b(M,z)$'
        plt.ylabel = r'$(M_{200}^{true}-M)/M_{200}^{true}$'

        epsfile = filename.replace('.rec','.eps')
        print 'writing eps:',epsfile
        plt.write_eps(eu.ostools.expand_path(epsfile))
 
class MvirFinder:
    def __init__(self, omega_m, z, delrho, mean=False):
        self.omega_m=omega_m
        self.z=z
        self.delrho = delrho

        self.rhoc0 = 0.27752

        if mean:
            self.rho_std = self.rhoc0*omega_m*(1+z)**3
        else:
            Ez = omega_m*(1+z)**3 + (1.-omega_m)
            self.rho_std = self.rhoc0*Ez

        self.rmin = 0.1
        self.rmax = 3.0

    def find_mvir(self, r, mass, masserr):
        """
        Look where the mass curve crosses the mass from the 
        background times delrho

        r in Mpc, mass in Msun

        """
        rap, massap, massaperr, msap = self.get_aperture(r,mass,masserr)

        wl=where1(msap <  massap)
        wg=where1(msap >= massap)

        i1 = rap[wl].argmax()
        i2 = rap[wg].argmin()
        i1 = wl[i1]
        i2 = wg[i2]
        r1 = rap[i1]
        r2 = rap[i2]

        if i2-i1 != 1:
            raise ValueError("indices i1,i2 are not consecutive")

        
        al=(log(massap[i2])-log(massap[i1]))/(log(rap[i2])-log(rap[i1]))
        A=massap[i1]/(rap[i1]**al)

        B=self.mfactor()
        rvir=(A/B)**(1./(3.-al))

        mvir=B*rvir**3
        mvir_err=(massaperr[i1]+massaperr[i2])/2.0
        rvir_err=rvir*(mvir_err/mvir)/3.0

        return rvir, rvir_err, mvir, mvir_err

    def get_aperture(self, r, mass, masserr):
        """
        get r,mass,masserr and the mass from background
        within our aperture (rmin,rmax)
        """


        w=where1((r >= self.rmin) & (r <= self.rmax))

        rap=r[w]
        massap=mass[w]
        massaperr=masserr[w]

        # mass from background (standard rho)
        msap = self.mstd(r[w])

        wl = where1(msap < massap)
        wg = where1(msap >= massap)
        if wl.size == w.size or wg.size == w.size:
            raise ValueError("lines do not cross in "
                             "aperture: %s,%s" % (self.rmin,self.rmax))

        if rap[wl].max() > rap[wg].min():
            raise ValueError("mass is not monatonic")

        return rap, massap, massaperr, msap


    def mstd(self, r):
        """
        get volume*delrho*rho_std, the mass within aperture
        assuming standard density
        """
        return self.mfactor()*r**3

    def mfactor(self):
        """
        to get mass in aperture due to background, multiply by r**3
        """
        return (4./3.)*PI*self.delrho*self.rho_std*1.e12

    def test(self):
        """
        Generate an NFW mass and get it's r200
        """
        pass



def get_endpoint_amp(x,y,err,npts,endslope):
    wuse=arange(npts)+x.size-npts
    xuse=x[wuse]**(-endslope)
    yuse=y[wuse]
    vuse=err[wuse]**2

    # inverse variance weighted and get the Linear Least Squares best fit
    arg1 = yuse*xuse/vuse
    arg2 = xuse*xuse/vuse
    amp_end = arg1.sum()/arg2.sum()
    return amp_end

def geterrcov(err,cov):
    if cov is None:
        if err is None:
            raise RuntimeError("Send one of err or cov")
        cov = numpy.diagflat(err**2)
    else:
        err2 = numpy.diagonal(cov)
        wbad=where1(err2 <= 0)
        if wbad.size > 0:
            raise ValueError("cov must be positive definite")
        err = sqrt(err2)
    return err,cov


def xi_pl_int(R0, R1, r, b):
    """
    perform the numerical integral 

        \\frac{1}{\\pi} \\int_R0^R1 \\frac{dR -W'(R)}{\\sqrt{R^2-r_j^2} 

    when -W'(R) is a powerlaw R^(-b) with b positive
    transforms variables first
    """

    x0=sqrt((R0/r)**2-1)
    x1=sqrt((R1/r)**2-1)

    p = -(b+1.0)/2.0
    #ival = romberg(xi_int_kernal, x0, x1, args=(p,))
    ival = intfunc(xi_int_kernal, x0, x1, (p,))
    return ival*r**(-b)/PI
    
def xi_int_kernal(x, p):
    return (1.0 +x**2)**p

def mass_pl_int(x0,x1,b,type):
    """
    perform the numerical integral 

        \\int_X0^X1 dx (1+x^2)^(-b/2) f(x)

    where 
    
        f(x)=1+2*x^2                 if type=='in'
        f(x)=1+2*x^2 -2x*sqrt(1+x^2) if type == 'out'

    """


    p = -b/2.0
    if type == 'in':
        func = mass_int_kernal_in
    else:
        func = mass_int_kernal_out

    #return romberg(func, x0, x1, args=(p,))
    return intfunc(func, x0, x1, (p,))



def mass_int_kernal_in(x, p):
    return (1.0 + 2*x**2)*(1.0 + x**2)**p

def mass_int_kernal_out(x, p):
    # use taylor expansion for large numbers

    if x < 10.0: 
        return (1. +x**2)**p * (1. +2*x*(x-sqrt(1. +x**2)))
    elif x >= 10.0 and x < 100.0:
        return (1. +x**2)**p * (1./(4.*x**2)-1./(8.*x**4)+5./(16.*x**6))
    else:
        return (1. +x**2)**p * (1./(4.*x**2)-1./(8.*x**4))


def intfunc(func, x0, x1, pars_tuple):
    #return romberg(func, x0, x1, args=pars_tuple)
    ret,err = scipy.integrate.quad(func, x0, x1, args=pars_tuple)
    return ret


def interp_hybrid(x0, x1, y0, y1, e0, e1, x, nsig=1, more=False):
    # this functions takes an interval x0,x1 and fuction vales at those
    # points y0 y1 and errors at those points e1 e2 and figures out an
    # appropriate interpolation for points x. 
    #
    # The hybrid interpolation is given by y=Amp* x^(slope) -off. This has
    # better stability than trying to fit powerlaw to points close to
    # zero or even negative points (impossible). The offset is determined as the amount
    # to add to y to make the point nsig*sigma above zero.

    amp, slope, off = interp_hybrid_pars(x0, x1, y0, y1, e0, e1, nsig=1)
    y = amp*x**slope - off
    if more:
        return y, amp, slope, off
    else:
        return y

def interp_hybrid_pars(x0, x1, y0, y1, e0, e1, nsig=1):
    off = (-min([y0-nsig*e0,y1-nsig*e1]))
    if off < 0:
        off=0.0

    dx = log(x1)-log(x0)
    slope=(log(y1+off)-log(y0+off))/dx
    amp=exp((log(x1)*log(y0+off)-log(x0)*log(y1+off))/dx)

    return amp, slope, off

def test_interp_hybrid():
    """
    Send y0,y1 with one or both less than zero to test the
    hybrid offset scheme
    """
    slope = -2.0
    #xvals = 0.1+linspace(0.0,8.0,9)
    xvals = 10.0**linspace(0.0,1.0,10)
    yvals = xvals**slope
    #yerr = 0.5*yvals
    #yerr = sqrt(yvals)
    yerr = yvals.copy()
    yerr[:] = 0.05

    #xfine = 0.1+linspace(0.0,8.0,1000)
    xfine = 10.0**linspace(0.0,1.0,1000)
    yfine = xfine**slope 

    #yerr = yvals.copy()
    #yerr[:] = 2

    plt=FramedPlot()
    plt.xrange = [0.5*xvals.min(), 1.5*xvals.max()]
    plt.xlog=True
    #plt.ylog=True
    plt.add(Points(xvals,yvals,type='filled circle',size=1))
    plt.add(Curve(xfine,yfine,color='blue'))
    #w=where1( (yvals-yerr) > 1.e-5 )
    #plt.add(SymmetricErrorBarsY(xvals[w],yvals[w],yerr[w]))
    plt.add(SymmetricErrorBarsY(xvals,yvals,yerr))

    # make points in between consecutive xvals,yvals so we
    # can use hybrid 2-point function
    xi = numpy.zeros(xvals.size-1,dtype='f8')
    yi = numpy.zeros(xi.size,dtype='f8')
    for i in xrange(xi.size):

        logx = (log10(xvals[i+1])+log10(xvals[i]))/2.0
        xi[i] = 10.0**logx
        yi[i],amp,slope,off = interp_hybrid(xvals[i], xvals[i+1], 
                                            yvals[i], yvals[i+1],
                                            yerr[i], yerr[i+1], 
                                            xi[i],more=True)
        
        print 'amp:',amp
        print 'slope:',slope
        print 'off:',off

    print xvals
    print xi
    print yi
    plt.add( Points(xi, yi, type='filled circle', size=1, color='red'))

    plt.show()


