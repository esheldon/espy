import os
from . import shapesim
import esutil as eu
from esutil.misc import wlog
from esutil.numpy_util import where1
import numpy
from numpy import median

from lensing.util import shear_fracdiff, e2gamma, gamma2e

class SimPlotter(dict):
    def __init__(self, run):
        c = shapesim.read_config(run)
        for k,v in c.iteritems():
            self[k] = v

        self.simc = shapesim.read_config(c['sim'])

        d = shapesim.get_plot_dir(run)
        if not os.path.exists(d):
            os.makedirs(d)
        
        self._data=None

    def doplots(self, 
                s2meas=False,
                type='diff',
                s2max=None, 
                yrange=None, 
                reduce_key=False, 
                show=True):
        import biggles
        import pcolors
        import converter

        data = self.read_data(s2meas=s2meas, s2max=s2max)

        epsfile = shapesim.get_plot_file(self['run'],type,
                                         s2max=s2max,
                                         s2meas=s2meas,
                                         yrange=yrange)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        plt = biggles.FramedPlot()
        plt.aspect_ratio=1
        #plt.xlabel=r'$\gamma$'
        plt.xlabel=r'ellipticity'

        if type == 'diff':
            plt.ylabel=r'$\Delta \gamma$'
        elif type == 'fdiff':
            plt.ylabel=r'$\Delta \gamma/\gamma$'
        else:
            raise ValueError("type should be 'diff or 'fdiff'")
 
        allplots=[]
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            #s2 = median(st['s2_meas'])
            if s2meas:
                s2 = median(st['s2_meas'])
            else:
                s2 = median(st['s2'])

            s = st['etrue'].argsort()

            #if 'e_chol' in st.dtype.names:
            if False:
                print 'using chol'
                e_meas = st['e_chol'][s]
            elif 'e_meas' not in st.dtype.names:
                e_meas = gamma2e(st['gamma_meas'][s])
            else:
                e_meas = st['e_meas'][s]

            etrue = st['etrue'][s]
            #etrue = st['e_uw'][s]

            wbad=where1(e_meas != e_meas)

            if wbad.size != 0:
                wlog("found bad:",e_meas[wbad])

            fdiff = shear_fracdiff(etrue,e_meas)
            # straight diff
            gammadiff = fdiff*st['gamma'][s]

            if type == 'diff':
                yplot=gammadiff
            elif type == 'fdiff':
                yplot=fdiff
            else:
                raise ValueError("type should be 'diff or 'fdiff'")

            fdiff_abs = numpy.abs(fdiff)
            wf = fdiff_abs.argmax()
            if fdiff_abs[wf] > 0.05:
                wlog("  large fdiff %g for s2: "
                     "%g(%d) e: %g:" % (fdiff[wf],s2,i,etrue[wf]))
            label = r'%0.3f' % s2
            cr = biggles.Curve(etrue, yplot, color=colors[i])
            cr.label = label

            plt.add(cr)
            allplots.append(cr)

        fsize=1.5
        if not reduce_key:
            key = biggles.PlotKey(0.9,0.9, allplots, halign='right', 
                                  fontsize=fsize)
        else:
            # pick a few
            nplot=len(allplots)
            tplots = [allplots[0], 
                      allplots[nplot*1/4], 
                      allplots[nplot/2], 
                      allplots[nplot*3/4], 
                      allplots[-1]]
            key = biggles.PlotKey(0.9,0.9, tplots, halign='right', fontsize=fsize)

        plt.add(key)

        klabtext=r'$<\sigma^2_{psf}/\sigma^2_{gal}>$'
        if self['s2n'] > 0:
            klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.95,0.95,klabtext,
                                 fontsize=1.5,halign='right')
        plt.add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.1,0.9, plab, halign='left')
        plt.add(l)

        if self.simc['psfmodel'] == 'turb':
            siglab=r'$FWHM_{PSF}: %.1f$ pix' % self.simc['psf_fwhm']
        else:
            psf_sigma = self.simc['psf_sigma']
            siglab=r'$\sigma_{PSF}: %.1f$ pix' % psf_sigma

        s2n=self['s2n']
        if s2n > 0:
            if s2n > 1000:
                ls2n = numpy.log10(s2n)
                ls2n = r'$10^{%.1f}$' % ls2n
            else:
                ls2n = '%.0f' % s2n

            #siglab+=r'$ S/N: %(s2n)d N_{trial}: %(ntrial)d$' % self
            siglab+=' S/N: %s' % ls2n
        #else:
        #    siglab+=r'$  N_{trial}: %(ntrial)d$' % self

        sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                               fontsize=2.5)
        plt.add(sl)

        if not reduce_key:
            plt.xrange = [0,1.4]
        if yrange is not None:
            plt.yrange = yrange

        wlog("Writing plot file:",epsfile)
        if show:
            plt.show()
        plt.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)




    def read_data(self, s2meas=False, s2max=None):
        if self._data is None:
            wlog("reading data")
            self._data = shapesim.read_all_outputs(self['run'],
                                                   average=True,verbose=True)
        
        alldata = self._data
        ntot=len(alldata)
        nkeep=ntot
        if s2max is not None:
            if s2meas:
                tag='s2_meas'
            else:
                tag='s2'
            keepdata = []
            for st in alldata:
                med_s2 = median(st[tag])
                wlog("med s2:",med_s2," max:",s2max)
                if med_s2 < s2max:
                    wlog("keeping")
                    keepdata.append(st.copy())
                else:
                    wlog("not keeping")
            nkeep = len(keepdata)
            wlog("kept %d/%d with s2 < %.3g" % (nkeep,ntot,s2max)) 
        else:
            keepdata = alldata

        return keepdata
