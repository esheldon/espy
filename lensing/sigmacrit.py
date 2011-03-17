from sys import stdout
import numpy
import esutil
from esutil import cosmology as cosmo
import zphot
from math import pi as PI

class ScinvCalculator:
    def __init__(self, zsmin, zsmax,
                 npts=100, 
                 n_zlens=20, zlmin=0.02, zlmax=0.6,
                 omega_m=0.3, 
                 omega_l=0.7,
                 omega_k=0.0,
                 h=1.0,
                 flat=True):
        """

        Specialized calculator for integrating over a p(zs) to get a mean
        inverse critical density on a grid of zlens.
            
            NOTE: npts for the distance calculations is always 5, which is
            good to 1.e-8 to redshift 1

        The input zsmin,zsmax are used to pre-compute gauss-legendre xvals and
        weights.   The zs and p(zs) input later should be on a subset of this
        range.

        The zlmin,zlmax are used to create a grid in zlens, and it is on this
        grid that the final mean inverse critical density will be computed.

        usage:
            # initialize for the given range in zlens and zsource
            zsmin=0.03
            zsmax=1.47
            scalc=ScinvCalculator(zsmin, zsmax, 
                                  zlmin=0.02, zlmax=0.6, npts=100)

            # now for a given function p(zsource) evaluated at points
            # zsource, return the scinv(zlens).

            mean_scinv = scalc.calc_mean_sigmacritinv(zsource, pzsource)
        """



        self.npts = npts
        self.n_zlens = n_zlens

        zlvals = numpy.arange(n_zlens, dtype='f8')
        self.zlvals = esutil.numpy_util.arrscl(zlvals, zlmin, zlmax)

        # now gauss-legendre weights and x vals used for integration
        # over zs
        self.xii, self.wii = esutil.math_util.gauleg(-1.0, 1.0, npts)

        self.f1 = (zsmax - zsmin)/2.0
        self.f2 = (zsmax + zsmin)/2.0

        self.zsvals_int = self.xii*self.f1 + self.f2


        stdout.write("Precomputing scinv on a grid of zl,npts=%d... " % npts)
        stdout.flush()
        # precompute
        self.scinv = numpy.zeros((n_zlens, npts),dtype='f8')

        c = cosmo.Cosmo(omega_m=omega_m, 
                        omega_l=omega_l, 
                        omega_k=omega_k, 
                        h=h, flat=flat)

        for i in range(n_zlens):
            zl = self.zlvals[i]
            self.scinv[i,:] = c.sigmacritinv(zl,self.zsvals_int)
        stdout.write("done\n")


    def interpolate_pofz(self, z, pz):
        """
        pz must be evaluated at the same locations as self.zsvals
        """
        pzvals = esutil.stat.interplin(pz, z, self.zsvals_int)
        return pzvals
       
    def calc_mean_sigmacritinv(self, z, pz):
        """
        pz will be interpolated to the self.zsvals_int locations
        """
        mean_scinv = numpy.zeros(self.zlvals.size, dtype='f8')

        pzinterp = self.interpolate_pofz(z, pz)

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

        stdout.write("%s\n" % pbsname)
        fobj=open(pbsname,'w')
        fobj.write(pbs)
        fobj.close()



class Tester:
    def __init__(self):
        self.data=None
        self.zphot = zphot.cas.CasZphot()

    def get_example_data(self):
        dir = self.zphot.output_dir()
        self.example_fname = 'pofz.ra4.2h.dat'
        fname=os.path.join(dir, self.example_fname)
        stdout.write("Reading sample data: %s\n" % fname)
        self.data = self.zphot.read_raw(fname)

    def plot_dir(self):
        return os.path.join(self.zphot.basedir,'test_sigmacritinv')

    def plot_file(self, num):
        dir=self.plot_dir()
        f=self.example_fname.replace('.dat','')
        f = f+'-%d.png' % num
        f = os.path.join(dir, f)
        return f


    def test_sigmacritinv(self, nplot, reload=False):

        from time import time
        import lensing

        plot_dir = self.plot_dir()
        import lensing
        if self.data is None or reload:
            self.get_example_data()


        # this is old ScinvCalculator, need to make work
        # with new one
        scalc1000 = lensing.ScinvCalculator(1000)


        nptsvals = [10,20,40,60,80,100]
        numcheck = len(nptsvals)
        colors=['black','magenta','blue','green','orange','red']
        scalc = []
        for npts in nptsvals:
            scalc.append(ScinvCalculator(npts))

        times = numpy.zeros(numcheck, dtype='f8')
        time1000 = 0.0
         
        # we'll fill this in
        scinv_all = numpy.zeros( (numcheck, scalc1000.zlvals.size) )

        figsize=(6.8,8.8)
        plt = esutil.plotting.setuplot(backend='Agg',
                                       params={'figure.figsize':figsize})
        #fig1 = plt.subplot(2,1,1)
        #fig2 = plt.subplot(2,1,2)

        xlim = [0, scalc1000.zsvals.max()]
        for i in range(nplot):
            pz = self.data['pz'][i]

            stdout.write("Doing 1000..."); stdout.flush()
            
            tm0 = time()
            scinv1000 = scalc1000.calc_mean_sigmacritinv(pz)
            time1000 += time()-tm0

            stdout.write("done\n"); stdout.flush()


            for j in range(numcheck):
                npts = nptsvals[j]
                stdout.write("%d " % npts)
                stdout.flush()
                #scinv_old = lensing.tools.mean_sigmacritinv(zl,self.zsvals,pz,
                #                                        npts=npts)
                tm0=time()
                scinv_all[j,:] = scalc[j].calc_mean_sigmacritinv(pz)
                times[j] += time()-tm0

            stdout.write("\nplotting\n"); stdout.flush()

            

            # plot the p(z)
            plt.clf()

            plt.subplot(3,1,1)
            plt.plot(scalc1000.zsvals, pz)
            plt.ylabel(r'$P(z_{source})$')
            plt.xlabel(r'$z_{source}$')
            plt.xlim(xlim)




            # now plot the inverse critical densities
            plt.subplot(3,1,2)

            for j in range(numcheck):
                npts = nptsvals[j]
                plt.plot(scalc[j].zlvals, scinv_all[j,:],
                         color=colors[j], label='npts=%d' % npts)

            leg = plt.legend()
            leg.draw_frame(False)
            plt.ylabel(r'$\langle \Sigma_{crit}(z_{lens}) \rangle$')
            plt.xlabel(r'$z_{lens}$')
            plt.xlim(xlim)


            # now ratios w.r.g. npts=1000
            plt.subplot(3,1,3)
            for j in range(numcheck):
                npts = nptsvals[j]
                ratio=scinv_all[j,:]/scinv1000[:]

                plt.plot(scalc[j].zlvals, ratio, 
                         color=colors[j], label='npts=%d' % npts)

            leg = plt.legend()
            leg.draw_frame(False)
            plt.ylabel(r'$\langle \Sigma_{crit} \rangle / \langle \Sigma_{crit} \rangle_{1000}$')
            plt.xlabel(r'$z_{lens}$')
            plt.xlim(xlim)
            plt.ylim([0.95,1.05])

            stdout.write("done\n")
            stdout.flush()

            plotfile=self.plot_file(i)
            stdout.write("writing to file: %s  ..." % plotfile)
            plt.savefig(plotfile)
            stdout.write("done\n")

        stdout.write("time npts=1000: %s\n" % time1000)
        for j in range(numcheck):
            npts=nptsvals[j]
            stdout.write("time npts=%s: %s\n" % (npts,times[j]))


