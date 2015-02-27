from __future__ import print_function
import os
import numpy
import esutil
from esutil import cosmology as cosmo
from esutil.ostools import path_join
from esutil.numpy_util import where1
from math import pi as PI

class ScinvCalculator:
    def __init__(self,
                 zlmin,
                 zlmax,
                 nzl,
                 zsmin,
                 zsmax,
                 npts=100, 
                 omega_m=0.3, 
                 omega_l=0.7,
                 omega_k=0.0,
                 H0=100.0, 
                 flat=True):
        """
        Specialized calculator for integrating over a p(zs) to get a mean
        inverse critical density on a pre-determined grid of zlens.
            
        parameters
        ----------
        zlmin: float
            Minimum zlens point
        zlmax: float
            Maximum zlens point
        nzl: int
            Number of points in zlens grid
        zsmin: float
            Min val over which to integrate p(zs)
        zsmax: float
            Max val over which to integrate p(zs)
        npts: int, optional
            Number of points for integration over p(zs)
        omega_m: float, optional
            default 0.3
        omega_l: float, optional
            default 0.7
        omega_k: float, optional
            default 0.0
        flat: bool, optional
            default True, flat universe

        NOTE: npts for the distance calculations is always 5, which is
            good to 1.e-8 to redshift 1

        usage:
            nzl=65
            zlmin=0.09
            zlmax=0.95
            zsmin=0.0
            zsmax=1.5
            npts=100 # integration of p(z)
            scalc=ScinvCalculator(zlmin, zlmax, nzl, zsmin, zsmax, npts=npts)

            mean_scinv = scalc.calc_mean_scinv(zs, pzs)
        """

        self.omega_m=omega_m
        self.omega_l=omega_l
        self.omega_k=omega_k
        self.flat=flat
        self.H0=H0

        self.setup_cosmology()
        self.setup_gauleg(zsmin, zsmax, npts)
        self.setup_zl(zlmin, zlmax, nzl)
        self.setup_scinv_grid()

       
    def calc_mean_scinv(self, zs, pz):
        """
        pz must correspond exactly to the zsvals input

        pz will be interpolated onto the gauss-legendra
        abcissa

        """

        if pz.size != zs.size:
            raise ValueError("pz(%d) and zs(%d) must be same "
                             "size" % (pz.size, zs.size))

        mean_scinv = numpy.zeros(self.nzl, dtype='f8')

        # get p(z) interpolated to the gauleg integration points
        pzinterp = self.interpolate_pofz(zs, pz)

        for i in range(self.nzl):
            # we've pre-computed scinv at zl,zs locations of relevance
            # now just multiply by the interpolated p(z) and do the
            # integral
            numerator = self.f1*self.scinv[i,:]*pzinterp*self.wii
            denom = self.f1*pzinterp*self.wii

            mean_scinv[i] = numerator.sum()/denom.sum()

        return mean_scinv

    def interpolate_pofz(self, z, pz):
        """
        interpolate p(z) onto the points used for gauleg integration
        """

        if self.zsmin < z[0] or self.zsmax > z[-1]:
            tup=(z[0],z[-1],self.zsmin,self.zsmax)
            raise ValueError("attempt to interpolate outside of range "
                             "[%g,%g] <-> [%g,%g] " % tup)

        pzvals = esutil.stat.interplin(pz, z, self.zsvals_int)

        pzvals.clip(min=0.0, out=pzvals)
        return pzvals

    def setup_cosmology(self):
        """
        Create the cosmology object used for sigmacrit calculations
        """
        self.cosmo = cosmo.Cosmo(omega_m=self.omega_m, 
                                 omega_l=self.omega_l, 
                                 omega_k=self.omega_k, 
                                 H0=self.H0,
                                 flat=self.flat)

    def setup_scinv_grid(self):
        """
        Set up the pre-computed grid of scinv(zl) over which we will interpolate
        """
        #print("Precomputing scinv on a grid of "
        #      "zlmin: %g zlmax: %g nzl: %d "
        #      "npts: %d" % (self.zlmin, self.zlmax, self.nzl, self.npts))

        self.scinv = numpy.zeros((self.nzl, self.npts),dtype='f8')

        c=self.cosmo
        for i in range(self.nzl):
            zl = self.zlvals[i]
            self.scinv[i,:] = c.sigmacritinv(zl,self.zsvals_int)

    def setup_zl(self, zlmin, zlmax, nzl):
        """
        the points where inverse sigmacrit will be evaluated
        """
        self.zlmin=zlmin
        self.zlmax=zlmax
        self.nzl = nzl
        self.zlvals = numpy.linspace(zlmin, zlmax, nzl)

    def setup_gauleg(self, zsmin, zsmax, npts):
        """
        set up the gauss-legendre weights and x vals used for integration over
        zs
        """

        self.zsmin=zsmin
        self.zsmax=zsmax
        self.npts = npts

        self.xii, self.wii = esutil.integrate.gauleg(-1.0, 1.0, npts)

        self.f1 = (zsmax - zsmin)/2.0
        self.f2 = (zsmax + zsmin)/2.0

        # the gauss-legendre integration points: must interpolate the 
        # input p(zs) to this grid
        self.zsvals_int = self.xii*self.f1 + self.f2



def generate_scinv_pbs():
    """
    Generate pbs files for scinv calculations using the cas photozs
    """
    import zphot
    from glob import glob

    pbsdir = '~/pbs/add-scinv'
    pbsdir=os.path.expanduser(pbsdir)
    if not os.path.exists(pbsdir):
        os.mkdir(pbsdir) 

    zobj = zphot.cas.CasZphot()
    dir = zobj.output_dir()
    pattern = os.path.join(dir, '*.dat')
    files = glob(pattern)
    
    for f in files:

        base = os.path.basename(f)
        shortname = base.replace('pofz.','')
        shortname = shortname.replace('.dat','')

        pbsname = base.replace('.dat','.pbs')
        pbsname = os.path.join(pbsdir, pbsname)
        
        pbslog = pbsname.replace('.pbs','.pbslog')
        logf=pbsname.replace('.pbs','.pbs.log')

        pbs="""

#!/bin/bash
#PBS -S /bin/bash
#PBS -N %s
#PBS -j oe
#PBS -o %s
#PBS -m a
#PBS -V
#PBS -r n
#PBS -W umask=0022
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1

source /global/data/products/eups/bin/setups.sh

setup python
setup esutil -r ~/exports/esutil-work

logf="%s"
python ~esheldon/python/zphot/bin/add-scinv.py %s &> "$logf"
        """ % (shortname,pbslog,logf,f)

        print(pbsname)
        fobj=open(pbsname,'w')
        fobj.write(pbs)
        fobj.close()


def load_test_data(pzrun):
    chunk=0
    z_file = zphot.weighting.z_file(pzrun,chunk=chunk)
    zsfull = zphot.weighting.read_z(z_file)
    zs = (zsfull['zmin']+zsfull['zmax'])/2
    data = zphot.weighting.read_pofz(pzrun,chunk)

    return zs, data

def test_with_pzrun(pzrun, beg, end):
    zs, data=load_test_data(pzrun)

    tester=Tester(zs, data, label=pzrun)
    tester.test_scinv_dz(beg, end)

def test_des(version, type, beg, end):
    import des
    dpofz=des.pz.DESPofz(version, type)

    data=dpofz[beg:end]

    label='%s-%s' % (version,type)

    tester=Tester(dpofz.zvals, data, label=label)
    tester.test_scinv_dz(beg, end, yrange=[0,3e-4])


class Tester(object):
    def __init__(self, zs, data, label='test'):
        self.data=data
        self.zs=zs
        self.label=label

        d=self.plot_dir()
        if not os.path.exists(d):
            os.makedirs(d)

    def test_scinv_dz(self, beg, end, yrange=[0,2.1e-4], show=False, reload=False, type='png'):
        """

        Test accuracy of interpolating scinv as a function of dzl, the
        lens redshift spacing.

        """
        import biggles
        from biggles import Points,FramedPlot,PlotKey, Table,Histogram,Curve
        from time import time
        import lensing
        import pcolors

        biggles.configure('default','fontface','HersheySans')
        biggles.configure('default','fontsize_min',1.3)

        zsmin=self.zs[0]
        zsmax=self.zs[-1]

        zlmin = 0.00
        zlmax = 1.0

        #dzl_vals = numpy.linspace(0.001,0.015,10)
        dzl_vals = numpy.linspace(0.001,0.015,4)
        nzl_vals = ( (zlmax-zlmin)/dzl_vals ).astype('i8')

        numcheck = len(dzl_vals)
        colors=pcolors.rainbow(numcheck, 'hex')
        scalc = []
        for nzl in nzl_vals:
            s=ScinvCalculator(zlmin, zlmax, nzl, 
                              zsmin, zsmax, npts=100)
            scalc.append(s)

        times = numpy.zeros(numcheck, dtype='f8')
         
        # we'll fill this in
        #scinv_all = numpy.zeros( (numcheck, scalc[0].zlvals.size) )

        xlim = [0, zsmax]
        for i in xrange(beg,end):
            scinv_all = []
            pz = self.data['pofz'][i]

            #print(pz)

            for j in xrange(numcheck):
                dzl = dzl_vals[j]
                nzl = nzl_vals[j]
                print("    nzl: %s dzl: %g" % (nzl, dzl))
                tm0=time()
                #scinv_all[j,:] = scalc[j].calc_mean_scinv(pz)
                sc=scalc[j].calc_mean_scinv(self.zs, pz)
                #sc=sc.clip(min=0.0)
                #print("sc",j,sc)
                scinv_all.append(sc)
                times[j] += time()-tm0

            print("\nplotting")


            # plot the p(z)
            tab = Table(3,1)

            binsize=self.zs[1]-self.zs[0]
            pzh = Histogram(pz, x0=self.zs[0], binsize=binsize)
            plt_pzh = FramedPlot()
            plt_pzh.xrange = xlim

            plt_pzh.xtitle=r'$z_s$'
            plt_pzh.ytitle=r'$P(z_s)$'
            plt_pzh.add(pzh)
            tab[0,0] = plt_pzh
            #plt_pzh.show()

            # plot scinv for each dzl
            plt_scinv = FramedPlot()
            plt_scinv.xrange = xlim

            scinv_plots=[]
            for j in xrange(numcheck):
                dzl = dzl_vals[j]
                nzl = nzl_vals[j]
                p = Curve(scalc[j].zlvals, scinv_all[j], type='solid',
                          color=colors[j])
                p.label = r'$nz_{lens}: %s dz_{lens}: %0.3f$' % (nzl,dzl)
                
                plt_scinv.add(p)
                scinv_plots.append(p)

            scinv_key = PlotKey(0.95,0.9,scinv_plots,halign='right')
            plt_scinv.add(scinv_key)

            plt_scinv.ylabel=r'$\Sigma_{crit}^{-1}(z_{lens})$'
            plt_scinv.xlabel=r'$z_{lens}$'
            plt_scinv.yrange = yrange

            #plt_scinv.show()
            tab[1,0] = plt_scinv

            # %diff to best dz

            plt_pdiff=FramedPlot()
            plt_pdiff.xrange = xlim
            plt_pdiff.yrange = [-0.05,0.05]
            pdiff_plots=[]

            zl_interp = numpy.linspace(zlmin,zlmax,1000)
            scinv_interp_best = esutil.stat.interplin(scinv_all[0], scalc[0].zlvals, zl_interp)

            w=where1(scinv_interp_best > 0)
            for j in xrange(numcheck):
                dzl = dzl_vals[j]
                nzl = nzl_vals[j]

                scinv_interp = esutil.stat.interplin(scinv_all[j], scalc[j].zlvals, zl_interp)

                if w.size > 0:
                    pdiff = scinv_interp[w]/scinv_interp_best[w] - 1.0
                    p = Curve(zl_interp[w], pdiff, type='solid',color=colors[j])
                else:
                    pdiff = numpy.ones(scinv_interp.size)
                    p = Curve(zl_interp, pdiff, type='solid',color=colors[j])

                p.label = r'$nz_{lens}: %s dz_{lens}: %0.3f$' % (nzl,dzl)

                plt_pdiff.add(p)
                pdiff_plots.append(p)

            key = PlotKey(0.95,0.9,pdiff_plots,halign='right')
            plt_pdiff.add(key)

            plt_pdiff.ylabel=r'$\Sigma_{crit}^{-1} /  \Sigma_{crit}^{-1}_{best} - 1$'
            plt_pdiff.xlabel=r'$z_{lens}$'

            tab[2,0] = plt_pdiff

            if show:
                tab.show()

            plotfile=self.dzl_plot_file(i, type)
            print("writing to file:",plotfile)
            if type == 'png':
                tab.write_img(1000,1000,plotfile)
            else:
                tab.write_eps(plotfile)

        for j in xrange(numcheck):
            dzl = dzl_vals[j]
            print("time dzl=%s: %s" % (dzl,times[j]))



    def test_scinv_npts(self, nplot, show=False, reload=False, type='png'):
        """

        Test accuracy as a function of the numer of points used in the
        integration.


        """

        dzl = 0.015
        zlmin = 0.02
        zlmax = 0.6

        from biggles import Points,FramedPlot,PlotKey, Table,Histogram,Curve
        from time import time
        import lensing
        import pcolors

        if self.data is None or reload:
            self.load_example_data()


        # this is old ScinvCalculator, need to make work
        # with new one
        scalc1000 = ScinvCalculator(self.zs, dzl, zlmin, zlmax, npts=1000)


        nptsvals = [100,200,300,400,500,600,700,800,900]
        numcheck = len(nptsvals)
        #colors=['black','magenta','blue','green','orange','red']
        colors=pcolors.rainbow(len(nptsvals), 'hex')
        scalc = []
        for npts in nptsvals:
            scalc.append(ScinvCalculator(self.zs, dzl, zlmin, zlmax, npts=npts))

        times = numpy.zeros(numcheck, dtype='f8')
        time1000 = 0.0
         
        # we'll fill this in
        scinv_all = numpy.zeros( (numcheck, scalc1000.zlvals.size) )


        xlim = [0, scalc1000.zsvals.max()]
        for i in xrange(nplot):
            pz = self.data['pofz'][i]

            print("Doing 1000...",end='')
            
            tm0 = time()
            scinv1000 = scalc1000.calc_mean_scinv(pz)
            time1000 += time()-tm0

            print("done")


            for j in xrange(numcheck):
                npts = nptsvals[j]
                print("%d " % npts,end='')
                tm0=time()
                scinv_all[j,:] = scalc[j].calc_mean_scinv(pz)
                times[j] += time()-tm0

            print("\nplotting")


            # plot the p(z)
            tab = Table(3,1)

            binsize=scalc1000.zsvals[1]-scalc1000.zsvals[0]
            pzh = Histogram(pz, x0=scalc1000.zsvals[0], binsize=binsize)
            plt_pzh = FramedPlot()
            plt_pzh.xrange = xlim

            plt_pzh.xtitle=r'$z_s$'
            plt_pzh.ytitle=r'$P(z_s)$'
            plt_pzh.add(pzh)
            tab[0,0] = plt_pzh

            # plot scinv for each npts value
            plt_scinv = FramedPlot()
            plt_scinv.xrange = xlim

            scinv_plots=[]
            for j in xrange(numcheck):
                npts = nptsvals[j]
                p = Curve(scalc[j].zlvals, scinv_all[j,:], type='solid',
                          color=colors[j])
                p.label = 'npts: %d' % npts
                
                plt_scinv.add(p)
                scinv_plots.append(p)

            scinv_key = PlotKey(0.95,0.9,scinv_plots,halign='right')
            plt_scinv.add(scinv_key)

            plt_scinv.ylabel=r'$\langle \Sigma_{crit}^{-1}(z_{lens}) \rangle$'
            plt_scinv.xlabel=r'$z_{lens}$'
            plt_scinv.yrange = [0,2.1e-4]

            tab[1,0] = plt_scinv

            # ratio to 1000 points

            plt_rat=FramedPlot()
            plt_rat.xrange = xlim
            plt_rat.yrange = [1-1.e-2,1+1.e-2]
            rat_plots=[]
            for j in xrange(numcheck):
                npts = nptsvals[j]
                w=where1(scinv1000 > 0)
                ratio = scinv_all[j,w]/scinv1000[w]
                #ratio=scinv_all[j,:]/scinv1000[:]

                p = Curve(scalc[j].zlvals[w], ratio,type='solid',color=colors[j])
                p.label = 'npts: %d' % npts

                plt_rat.add(p)
                rat_plots.append(p)

            key = PlotKey(0.95,0.9,rat_plots,halign='right')
            plt_rat.add(key)

            plt_rat.ylabel=r'$\langle \Sigma_{crit}^{-1} \rangle / \langle \Sigma_{crit}^{-1} \rangle_{1000}$'
            plt_rat.xlabel=r'$z_{lens}$'

            tab[2,0] = plt_rat

            if show:
                tab.show()


            plotfile=self.npts_plot_file(i, type)
            print("writing to file:",plotfile)
            if type == 'png':
                tab.write_img(1000,1000,plotfile)
            else:
                tab.write_eps(plotfile)

        print("time npts=1000:",time1000)
        for j in xrange(numcheck):
            npts=nptsvals[j]
            print("time npts=%s: %s" % (npts,times[j]))

    def plot_dir(self):
        d=os.environ['LENSDIR']
        d = path_join(d,'sigmacrit-tests',self.label)
        return d

    def npts_plot_file(self, index, type='png'):
        dir=self.plot_dir()
        f='sigmacrit-%s-test-npts-%06i.png' % (self.label,index)
        f=path_join(dir,f)
        return f

    def dzl_plot_file(self, index, type='png'):
        dir=self.plot_dir()
        f='sigmacrit-%s-test-dzl-%06i.png' % (self.label,index)
        f=path_join(dir,f)
        return f


