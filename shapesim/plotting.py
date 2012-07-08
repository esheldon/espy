import os
from . import shapesim
import lensing
import esutil as eu
from esutil.misc import wlog
from esutil.numpy_util import where1
import numpy
from numpy import median, zeros

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


    def plots_shear_vs_s2n(self, 
                           skip1=[],
                           skip2=[],
                           s2meas=False,
                           type='diff',
                           s2max=None, 
                           yrange=None, 
                           show=True):
        import biggles
        import pcolors
        import converter

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data(s2meas=s2meas, s2max=s2max,
                              skip1=skip1,skip2=skip2)

        epsfile = shapesim.get_plot_file(self['run'],type,yrange=yrange)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1

 
        plots1=[]
        plots2=[]
        allplots=[]
        max_s2n_len=0
        for st in data:
            if st['s2'].size > max_s2n_len:
                max_s2n_len = st['s2'].size
        avg=zeros(max_s2n_len, dtype=[('s2n','f8'),('shear1','f8'),('shear2','f8'),('n','i8')])

        # looping over s2
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            if s2meas:
                s2 = median(st['s2_meas'])
            else:
                s2 = median(st['s2'])

            if 's2n_matched' in st.dtype.names:
                s2n_name='s2n_matched'
                xlabel = r'$S/N_{matched}$'
            elif 's2n_uw' in st.dtype.names:
                s2n_name='s2n_uw'
                xlabel = r'$S/N_{uw}$'
            else:
                # assuming matched
                s2n_name='s2n'
                xlabel = r'$S/N_{matched}$'
            s = st[s2n_name].argsort()

            s2n = st[s2n_name][s]

            if self['run'][0:5] == 'deswl':
                tag1='gamma1_meas'
                tag2='gamma2_meas'
                # convention
                st[tag1] = -st[tag1]
            else:
                tag1='shear1'
                tag2='shear2'
           
            if type == 'diff':
                yvals1 = st[tag1][s] - shear_true.g1
                yvals2 = st[tag2][s] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1][s]
                yvals2 = st[tag2][s]

            else:
                raise ValueError("bad plot type: '%s'" % type)


            label = r'%0.3f' % s2
            pr1 = biggles.Curve(s2n, yvals1, color=colors[i])
            pr2 = biggles.Curve(s2n, yvals2, color=colors[i])
            pr1.label = label
            pr2.label = label


            arr[0,0].add(pr1)
            arr[1,0].add(pr2)
            #arr[0,0].add(err1)
            #arr[1,0].add(err2)
            if i < 15:
                plots1.append(pr1)
            else:
                plots2.append(pr1)

            if st['s2'].size == avg['n'].size:
                avg['n'] += 1
                avg['s2n'] += s2n
                avg['shear1'] += yvals1
                avg['shear2'] += yvals2

        avg['s2n'] = avg['s2n']/avg['n']
        avg['shear1'] = avg['shear1']/avg['n']
        avg['shear2'] = avg['shear2']/avg['n']
        avg1 = biggles.Points(avg['s2n'],avg['shear1'],type='filled circle',size=2)
        avg2 = biggles.Points(avg['s2n'],avg['shear2'],type='filled circle',size=2)
        avg1.label = 'average'

        if len(plots2) > 0:
            plots2.append(avg1)
        else:
            plots1.append(avg1)

        fsize=2
        key1 = biggles.PlotKey(0.9,0.85, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.9,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$<\sigma^2_{psf}/\sigma^2_{gal}>$'
        klab = biggles.PlotLabel(0.89,0.92,klabtext,
                                 fontsize=2.5,halign='right')
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
        siglab += r'$ e_{tot}: %.2f$' % st['etrue'].mean()

        sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                               fontsize=2.5)
        arr[1,0].add(sl)

        arr.xlabel=xlabel


        # might be a diff
        arr[0,0].add(avg1)
        arr[1,0].add(avg2)

        if type == 'val':
            arr.ylabel = r'$\gamma$'
            expect1 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()],
                                    [shear_true.g1,shear_true.g1])
            expect1.label = r'$\gamma_1$ = %.2g' % shear_true.g1
            expect2 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()],
                                    [shear_true.g2,shear_true.g2])
            expect2.label = r'$\gamma_2$ = %.2g' % shear_true.g2

            ekey1 = biggles.PlotKey(0.1,0.9, [expect1], halign='left', 
                                   fontsize=3)
            ekey2 = biggles.PlotKey(0.1,0.9, [expect2], halign='left', 
                                   fontsize=3)

            arr[0,0].add(expect1,ekey1)
            arr[1,0].add(expect2,ekey2)

        else:
            arr.ylabel = r'$\Delta \gamma$'
            g1lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_1$ = %.2g' % shear_true.g1, halign='left')
            g2lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_2$ = %.2g' % shear_true.g2, halign='left')

            arr[0,0].add(g1lab)
            arr[1,0].add(g2lab)

            expect1 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()], [0,0])
            expect2 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()], [0,0])

            arr[0,0].add(expect1)
            arr[1,0].add(expect2)



        arr.xrange = [0.1*s2n.min(),s2n.max()*1.4]
        if yrange is not None:
            arr.yrange = yrange

        wlog("Writing plot file:",epsfile)
        if show:
            arr.show()
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)




    def plot_shear_vs_e(self, 
                        skip1=[],
                        skip2=[],
                        s2meas=False,
                        type='diff',
                        s2max=None, 
                        yrange=None, 
                        show=True):
        import biggles
        import pcolors
        import converter

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data(s2meas=s2meas, s2max=s2max,
                              skip1=skip1,skip2=skip2)

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

        fsize=2
        key1 = biggles.PlotKey(0.9,0.85, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.9,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$<\sigma^2_{psf}/\sigma^2_{gal}>$'
        #if self['s2n'] > 0:
        #    klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.9,0.92,klabtext,
                                 fontsize=fsize,halign='right')
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


    def plot_ediff_vs_e(self, 
                        skip1=[],
                        skip2=[],
                        s2meas=False,
                        type='diff',
                        s2max=None, 
                        yrange=None, 
                        show=True):
        """
        Plot the measured ellipticity minus true vs true
        """
        import biggles
        import pcolors
        import converter

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        data = self.read_data(s2meas=s2meas, s2max=s2max,
                              skip1=skip1,skip2=skip2)

        epsfile = shapesim.get_plot_file(self['run'],'etot-'+type,
                                         s2max=s2max,
                                         s2meas=s2meas,
                                         yrange=yrange)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        plt=biggles.FramedPlot()
        plt.xlabel=r'e_{tot}'
        plt.ylabel = r'$\Delta e_{tot}$'

 
        plots=[]
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            if s2meas:
                s2 = median(st['s2_meas'])
            else:
                s2 = median(st['s2'])

            s = st['etrue'].argsort()

            etrue = st['etrue'][s]
            ediff = st['e_meas'][s] - etrue

            label = r'%0.3f' % s2
            cr = biggles.Curve(etrue, ediff, color=colors[i])
            cr.label = label

            plt.add(cr)
            plots.append(cr)

        fsize=2
        key = biggles.PlotKey(0.9,0.87, plots, halign='right', 
                              fontsize=fsize)
        plt.add(key)

        klabtext=r'$<\sigma^2_{psf}/\sigma^2_{gal}>$'
        #if self['s2n'] > 0:
        #    klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.86,0.94,klabtext,
                                 fontsize=fsize,halign='right')
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
        plt.aspect_ratio=1

        plt.xrange = [0,1.4]
        if yrange is not None:
            plt.yrange = yrange

        wlog("Writing plot file:",epsfile)
        if show:
            plt.show()
        plt.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)




    def read_data(self, s2meas=False, s2max=None, skip1=[], skip2=[]):
        if self._data is None:
            wlog("reading data")
            self._data = shapesim.read_all_outputs(self['run'],
                                                   skip1=skip1,skip2=skip2,
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


