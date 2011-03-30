from __future__ import print_function
import os
import numpy
import esutil
from esutil import cosmology as cosmo
from esutil.ostools import path_join
from esutil.numpy_util import where1
import zphot
from math import pi as PI

def make_zlvals(dzl, zlmin, zlmax):
    """

    Given the input, generate the zlvals for the <scinv> as a function of zl.
    Make sure the last point is >= zlmax

    """

    n_zlens = int( (zlmax-zlmin)/dzl ) + 1
    while True:
        zlvals = zlmin + dzl*numpy.arange(n_zlens, dtype='f8')
        if zlvals[-1] >= zlmax:
            return zlvals
        else:
            n_zlens += 1


class ScinvCalculator:
    def __init__(self, zsvals, dzl, zlmin, zlmax,
                 npts=100, 
                 omega_m=0.3, 
                 omega_l=0.7,
                 omega_k=0.0,
                 h=1.0,
                 flat=True):
        """

        Specialized calculator for integrating over a p(zs) to get a mean
        inverse critical density on a pre-determined grid of zlens.
            
            NOTE: npts for the distance calculations is always 5, which is
            good to 1.e-8 to redshift 1

        The input zsvals are used for the integration.  Each call to
        calc_mean_scinv must send a p(z) on that input grid. 

        The dzl,zlmin,zlmax are used to create a grid in zlens, and it is on
        this grid that the final mean inverse critical density will be
        computed.  It seems a dzl of 0.015 works pretty well.  This results
        in 19 lens z vals between 0.02 and 0.30

        usage:
            # initialize for 
            #  dzl=0.015
            #  zlmin=0.02
            #  zlmax=0.6
            #  npts=100 # integration of p(z)
            scalc=ScinvCalculator(zs, 0.015, 0.02, 0.6, npts=100)

            # now for a given function p(zsource) evaluated at points
            # zs, return the scinv(zlens).

            mean_scinv = scalc.calc_mean_scinv(pzsource)
        """


        self.zsvals = zsvals
        zsmax=zsvals.max()
        zsmin=zsvals.min()

        self.npts = npts

        self.zlvals = make_zlvals(dzl, zlmin, zlmax)
        self.n_zlens = self.zlvals.size

        # now gauss-legendre weights and x vals used for integration
        # over zs
        self.xii, self.wii = esutil.integrate.gauleg(-1.0, 1.0, npts)

        self.f1 = (zsmax - zsmin)/2.0
        self.f2 = (zsmax + zsmin)/2.0

        self.zsvals_int = self.xii*self.f1 + self.f2


        print("Precomputing scinv on a grid of dzl: %0.3f nzl: %d npts: %d... " % \
              (dzl,self.n_zlens,npts),end='')
        self.scinv = numpy.zeros((self.n_zlens, npts),dtype='f8')

        c = cosmo.Cosmo(omega_m=omega_m, 
                        omega_l=omega_l, 
                        omega_k=omega_k, 
                        h=h, flat=flat)

        for i in range(self.n_zlens):
            zl = self.zlvals[i]
            self.scinv[i,:] = c.sigmacritinv(zl,self.zsvals_int)
        print("done")


    def interpolate_pofz(self, z, pz):
        """
        pz must be evaluated at the same locations as self.zsvals
        """
        pzvals = esutil.stat.interplin(pz, z, self.zsvals_int)
        return pzvals
       
    def calc_mean_scinv(self, pz):
        """
        pz must correspond exactly to the zsvals input

        pz will be interpolated onto the gauss-legendra
        abcissa

        """
        if pz.size != self.zsvals.size:
            raise ValueError("pz must be same size as zsvals")

        mean_scinv = numpy.zeros(self.zlvals.size, dtype='f8')

        # get p(z) interpolated to the zsvals_int points
        pzinterp = self.interpolate_pofz(self.zsvals, pz)

        for i in range(self.zlvals.size):
            # we've pre-computed scinv at zl,zs locations of relevance
            # now just multiply by the interpolated p(z) and do the
            # integral
            numerator = self.f1*self.scinv[i,:]*pzinterp*self.wii
            denom = self.f1*pzinterp*self.wii

            mean_scinv[i] = numerator.sum()/denom.sum()

        return mean_scinv





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



class Tester:
    def __init__(self, pzrun):
        self.pzrun = pzrun
        self.wo = zphot.weighting.WeightedOutputs()

        # get the zs values
        # the z values for the histogram
        z_file = self.wo.z_file(pzrun,chunk=0)
        zsfull = zphot.weighting.read_z(z_file)
        self.zs = (zsfull['zmin']+zsfull['zmax'])/2


        self.data=None

    def load_example_data(self, chunk=0):
        self.data = zphot.weighting.read_pofz_byrun(self.pzrun,chunk)

    def test_scinv_dz(self, nplot, show=False, reload=False, type='png'):
        """

        Test accuracy of interpolating scinv as a function of dzl, the
        lens redshift spacing.


        """

        zlmin = 0.02
        zlmax = 0.6

        from biggles import Points,FramedPlot,PlotKey, Table,Histogram,Curve
        from time import time
        import lensing
        import pcolors

        if self.data is None or reload:
            self.load_example_data()


        dzl_vals = numpy.linspace(0.001,0.015,10)
        numcheck = len(dzl_vals)
        colors=pcolors.rainbow(numcheck, 'hex')
        scalc = []
        for dzl in dzl_vals:
            scalc.append(ScinvCalculator(self.zs, dzl, zlmin, zlmax, npts=100))

        times = numpy.zeros(numcheck, dtype='f8')
         
        # we'll fill this in
        #scinv_all = numpy.zeros( (numcheck, scalc[0].zlvals.size) )


        xlim = [0, scalc[0].zsvals.max()]
        for i in xrange(nplot):
            scinv_all = []
            pz = self.data['pofz'][i]

            for j in xrange(numcheck):
                dzl = dzl_vals[j]
                print("%f " % dzl,end='')
                tm0=time()
                #scinv_all[j,:] = scalc[j].calc_mean_scinv(pz)
                scinv_all.append( scalc[j].calc_mean_scinv(pz) )
                times[j] += time()-tm0

            print("\nplotting")


            # plot the p(z)
            tab = Table(3,1)

            binsize=scalc[0].zsvals[1]-scalc[0].zsvals[0]
            pzh = Histogram(pz, x0=scalc[0].zsvals[0], binsize=binsize)
            plt_pzh = FramedPlot()
            plt_pzh.xrange = xlim

            plt_pzh.xtitle=r'$z_s$'
            plt_pzh.ytitle=r'$P(z_s)$'
            plt_pzh.add(pzh)
            tab[0,0] = plt_pzh

            # plot scinv for each dzl
            plt_scinv = FramedPlot()
            plt_scinv.xrange = xlim

            scinv_plots=[]
            for j in xrange(numcheck):
                dzl = dzl_vals[j]
                p = Curve(scalc[j].zlvals, scinv_all[j], type='solid',
                          color=colors[j])
                p.label = r'$dz_{lens}: %0.3f$' % dzl
                
                plt_scinv.add(p)
                scinv_plots.append(p)

            scinv_key = PlotKey(0.95,0.9,scinv_plots,halign='right')
            plt_scinv.add(scinv_key)

            plt_scinv.ylabel=r'$\langle \Sigma_{crit}^{-1}(z_{lens})\rangle$'
            plt_scinv.xlabel=r'$z_{lens}$'
            plt_scinv.yrange = [0,2.1e-4]

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

                scinv_interp = esutil.stat.interplin(scinv_all[j], scalc[j].zlvals, zl_interp)

                pdiff = scinv_interp[w]/scinv_interp_best[w] - 1.0

                p = Curve(zl_interp[w], pdiff, type='solid',color=colors[j])
                p.label = r'$dz_{lens}: %0.3f$' % dzl

                plt_pdiff.add(p)
                pdiff_plots.append(p)

            key = PlotKey(0.95,0.9,pdiff_plots,halign='right')
            plt_pdiff.add(key)

            plt_pdiff.ylabel=r'$\langle \Sigma_{crit}^{-1} \rangle / \langle \Sigma_{crit}^{-1} \rangle_{best} - 1$'
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
        d = path_join(d,'sigmacrit-tests')
        return d

    def npts_plot_file(self, index, type='png'):
        dir=self.plot_dir()
        f='sigmacrit-%s-test-npts-%06i.png' % (self.pzrun,index)
        f=path_join(dir,f)
        return f

    def dzl_plot_file(self, index, type='png'):
        dir=self.plot_dir()
        f='sigmacrit-%s-test-dzl-%06i.png' % (self.pzrun,index)
        f=path_join(dir,f)
        return f


