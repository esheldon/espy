import os
from . import shapesim
import lensing
import esutil as eu
from esutil.misc import wlog
from esutil.numpy_util import where1
import numpy
from numpy import median

from lensing.util import shear_fracdiff, e2gamma, gamma2e, g1g2_to_e1e2
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

    def doplots_vs_e(self, 
                s2meas=False,
                type='diff',
                s2max=None, 
                yrange=None, 
                show=True):
        import biggles
        import pcolors
        import converter

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data(s2meas=s2meas, s2max=s2max)

        epsfile = shapesim.get_plot_file(self['run'],type,
                                         s2max=s2max,
                                         s2meas=s2meas,
                                         yrange=yrange)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1
        arr.xlabel=r'ellipticity'
        arr.ylabel = r'$\Delta \gamma$'

 
        plots1=[]
        plots2=[]
        allplots=[]
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            if s2meas:
                s2 = median(st['s2_meas'])
            else:
                s2 = median(st['s2'])

            s = st['etrue'].argsort()

            etrue = st['etrue'][s]
            
            shear1diff = st['shear1'][s] - shear_true.g1
            shear2diff = st['shear2'][s] - shear_true.g2

            label = r'%0.3f' % s2
            cr1 = biggles.Curve(etrue, shear1diff, color=colors[i])
            cr2 = biggles.Curve(etrue, shear2diff, color=colors[i])
            cr1.label = label
            cr2.label = label

            arr[0,0].add(cr1)
            arr[1,0].add(cr2)
            if i < 15:
                plots1.append(cr1)
            else:
                plots2.append(cr1)

        fsize=1.5
        key1 = biggles.PlotKey(0.9,0.85, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.9,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$<\sigma^2_{psf}/\sigma^2_{gal}>$'
        if self['s2n'] > 0:
            klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.95,0.92,klabtext,
                                 fontsize=1.5,halign='right')
        arr[0,0].add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.9,0.1, plab, halign='right')
        arr[1,0].add(l)

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
        arr[1,0].add(sl)


        g1lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_1$ = %.2g' % shear_true.g1, halign='left')
        g2lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_2$ = %.2g' % shear_true.g2, halign='left')
        arr[0,0].add(g1lab)
        arr[1,0].add(g2lab)

        arr.xrange = [0,1.4]
        if yrange is not None:
            arr.yrange = yrange

        wlog("Writing plot file:",epsfile)
        if show:
            arr.show()
        arr.write_eps(epsfile)
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
