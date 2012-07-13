import os
from . import shapesim
import lensing
import esutil as eu
from esutil.misc import wlog
from esutil.numpy_util import where1
from esutil.stat import sigma_clip
import numpy
from numpy import median, zeros, sqrt
import copy

from lensing.util import shear_fracdiff, e2gamma, gamma2e, g1g2_to_e1e2

class SimPlotter(dict):
    def __init__(self, run, **keys):
        c = shapesim.read_config(run)
        for k,v in c.iteritems():
            self[k] = v

        self.simc = shapesim.read_config(c['sim'])
        self.psf_e1 = self.simc.get('psf_e1',None)
        self.psf_e2 = self.simc.get('psf_e2',None)

        self.psf_estring=''
        if self.psf_e1 is not None:
            self.psf_estring = r'$e_{psf}$: %.2f,%.2f' % (self.psf_e1,self.psf_e2)

        d = shapesim.get_plot_dir(run)
        if not os.path.exists(d):
            os.makedirs(d)
        
        self._data=None

        title=copy.deepcopy(self['run'])
        title = title.replace('gmix-fit','gmix')
        self.title = title.replace('-',' ')
        self.maketitle=keys.get('maketitle',False)

        self.skip1 = keys.get('skip1',[])
        self.skip2 = keys.get('skip2',[])

        self.s2min = keys.get('s2min',None)

    def get_title(self, title=None):
        if title is None and self.maketitle:
            title=self.title
        return title

    def plots_shear_vs_s2n(self, 
                           type='diff',
                           xrng=None,
                           yrng=None, 
                           doavg=False,
                           docum=False,
                           title=None,
                           show=True):
        import biggles
        import pcolors
        import converter

        if doavg:
            extra='-avg'
        else:
            extra=''

        puterr=False
        if docum:
            puterr=True
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data(docum=docum)

        epsfile = shapesim.get_plot_file(self['run'],type+extra,yrng=yrng)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1
        
        title=self.get_title(title=title)
        if title:
            arr.title=title
 
        plots1=[]
        plots2=[]
        allplots=[]

        if doavg:
            max_s2n_len=0
            for st in data:
                if st['s2'].size > max_s2n_len:
                    max_s2n_len = st['s2'].size
            avg=zeros(max_s2n_len, dtype=[('s2n','f8'),('shear1','f8'),('shear2','f8'),('n','i8')])

        # looping over s2
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

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

            s2n = st[s2n_name]

            if self['run'][0:5] == 'deswl':
                if docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='gamma1_meas'
                    tag2='gamma2_meas'
                # convention
                st[tag1] = -st[tag1]
            else:
                if docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='shear1'
                    tag2='shear2'

            if docum and i == (len(data)-1):
                # save for later
                pts1_0 = st[tag1].copy()
                pts2_0 = st[tag2].copy()
                if type == 'diff':
                    pts1_0 -= shear_true.g1
                    pts2_0 -= shear_true.g2
                #err1_0 = [pts1_0.std()]*pts1_0.size
                #err2_0 = [pts2_0.std()]*pts2_0.size

                #n=pts1_0.size
                #start=1
                #err1_0 = pts1_0[start:].std()
                #err2_0 = pts2_0[start:].std()
                #err1_0 = [min(err1_0,err2_0)]*n
                #err2_0 = err1_0

                err1_0 = st['shear1cum_alt_err']
                err2_0 = st['shear2cum_alt_err']
                print err1_0
          
            if type == 'diff':
                yvals1 = st[tag1] - shear_true.g1
                yvals2 = st[tag2] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1]
                yvals2 = st[tag2]

            else:
                raise ValueError("bad plot type: '%s'" % type)

            if docum:
                label = r'< %0.3f' % s2
            else:
                label = r'%0.3f' % s2
            pr1 = biggles.Curve(s2n, yvals1, color=colors[i])
            pr2 = biggles.Curve(s2n, yvals2, color=colors[i])
            pr1.label = label
            pr2.label = label


            arr[0,0].add(pr1)
            arr[1,0].add(pr2)
            
            if not docum and i == (len(data)-1):
                terr1 = [yvals1.std()]*yvals1.size
                terr2 = [yvals2.std()]*yvals2.size
                err1p = biggles.SymmetricErrorBarsY(s2n, yvals1, terr1,width=4)
                err2p = biggles.SymmetricErrorBarsY(s2n, yvals2, terr2,width=4)
                arr[0,0].add(err1p)
                arr[1,0].add(err2p)
            if i < 15:
                plots1.append(pr1)
            else:
                plots2.append(pr1)

            if doavg:
                if st['s2'].size == avg['n'].size:
                    avg['n'] += 1
                    avg['s2n'] += s2n
                    avg['shear1'] += yvals1
                    avg['shear2'] += yvals2

        if doavg:
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

        if docum:
            err1p = biggles.SymmetricErrorBarsY(s2n, pts1_0, err1_0)
            err2p = biggles.SymmetricErrorBarsY(s2n, pts2_0, err2_0)
            arr[0,0].add(err1p)
            arr[1,0].add(err2p)

        fsize=2
        key1 = biggles.PlotKey(0.9,0.92, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.9,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$'
        klab = biggles.PlotLabel(0.76,0.92,klabtext,
                                 fontsize=2.5,halign='right')
        arr[0,0].add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.9,0.1, plab, halign='right')
        arr[1,0].add(l)

        if self.simc['psfmodel'] == 'turb':
            siglab = r'$FWHM: %.1f$ pix' % self.simc['psf_fwhm']
        else:
            psf_sigma = self.simc['psf_sigma']
            siglab = r'$\sigma: %.1f$ pix' % psf_sigma
        siglab += ' '+self.psf_estring
        siglab += r'$ e_{gal}^{tot}: %.2f$' % st['etrue'].mean()

        sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                               fontsize=2.5)
        arr[1,0].add(sl)

        arr.xlabel=xlabel


        # might be a diff
        if doavg:
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



        if xrng is None:
            xrng = [0.1*s2n.min(),s2n.max()*1.4]
        arr.xrange = xrng

        if yrng is not None:
            arr.yrange = yrng

        wlog("Writing plot file:",epsfile)
        if show:
            arr.show()
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)




    def plot_shear_vs_e(self, 
                        type='diff',
                        yrng=None, 
                        doavg=False,
                        docum=False,
                        title=None,
                        show=True):
        import biggles
        import pcolors
        import converter

        if doavg:
            extra='-avg'
        else:
            extra=''
        if docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data(docum=docum)

        epsfile = shapesim.get_plot_file(self['run'],type+extra,
                                         s2min=self.s2min,
                                         yrng=yrng)
        wlog("will plot to:",epsfile)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1
        arr.xlabel=r'ellipticity'
        arr.ylabel = r'$\Delta \gamma$'

        title=self.get_title(title=title)
        if title:
            arr.title=title

        plots1=[]
        plots2=[]
        allplots=[]

        if doavg:
            max_e_len=0
            for st in data:
                if st['etrue'].size > max_e_len:
                    max_e_len = st['s2'].size
            avg=zeros(max_e_len, dtype=[('e','f8'),('shear1','f8'),('shear2','f8'),('n','i8')])

        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            s2 = median(st['s2'])

            etrue = st['etrue']
            
            if self['run'][0:5] == 'deswl':
                if docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='gamma1_meas'
                    tag2='gamma2_meas'
                # convention
                st[tag1] = -st[tag1]
            else:
                if docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='shear1'
                    tag2='shear2'

            if docum and i == (len(data)-1):
                # save for later
                pts1_0 = st[tag1].copy()
                pts2_0 = st[tag2].copy()
                if type == 'diff':
                    pts1_0 -= shear_true.g1
                    pts2_0 -= shear_true.g2
                n=pts1_0.size
                start=1
                err1_0 = pts1_0[start:].std()
                err2_0 = pts2_0[start:].std()
                err1_0 = [min(err1_0,err2_0)]*n
                err2_0 = err1_0
 
            if type == 'diff':
                yvals1 = st[tag1] - shear_true.g1
                yvals2 = st[tag2] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1]
                yvals2 = st[tag2]

            else:
                raise ValueError("bad plot type: '%s'" % type)


            if docum:
                label = r'< %0.3f' % s2
            else:
                label = r'%0.3f' % s2

            cr1 = biggles.Curve(etrue, yvals1, color=colors[i])
            cr2 = biggles.Curve(etrue, yvals2, color=colors[i])
            cr1.label = label
            cr2.label = label

            arr[0,0].add(cr1)
            arr[1,0].add(cr2)
            if i < 15:
                plots1.append(cr1)
            else:
                plots2.append(cr1)

            if doavg:
                if st['s2'].size == avg['n'].size:
                    avg['n'] += 1
                    avg['e'] += etrue
                    avg['shear1'] += yvals1
                    avg['shear2'] += yvals2

        if doavg:
            avg['e'] = avg['e']/avg['n']
            avg['shear1'] = avg['shear1']/avg['n']
            avg['shear2'] = avg['shear2']/avg['n']
            avg1 = biggles.Points(avg['e'],avg['shear1'],type='filled circle',size=2)
            avg2 = biggles.Points(avg['e'],avg['shear2'],type='filled circle',size=2)
            avg1.label = 'average'

            if len(plots2) > 0:
                plots2.append(avg1)
            else:
                plots1.append(avg1)

        if docum:
            err1p = biggles.SymmetricErrorBarsY(etrue, pts1_0, err1_0)
            err2p = biggles.SymmetricErrorBarsY(etrue, pts2_0, err2_0)
            arr[0,0].add(err1p)
            arr[1,0].add(err2p)


        fsize=2
        key1 = biggles.PlotKey(0.9,0.92, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.9,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$'
        #if self['s2n'] > 0:
        #    klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.76,0.92,klabtext,
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
        siglab += ' '+self.psf_estring

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


        #g1lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_1$ = %.2g' % shear_true.g1, halign='left')
        #g2lab = biggles.PlotLabel(0.1,0.9, r'$\gamma_2$ = %.2g' % shear_true.g2, halign='left')
        #arr[0,0].add(g1lab)
        #arr[1,0].add(g2lab)

        # might be a diff
        if doavg:
            arr[0,0].add(avg1)
            arr[1,0].add(avg2)


        if type == 'val':
            arr.ylabel = r'$\gamma$'
            expect1 = biggles.Curve([0,1],
                                    [shear_true.g1,shear_true.g1])
            expect1.label = r'$\gamma_1$ = %.2g' % shear_true.g1
            expect2 = biggles.Curve([0,1],
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

            expect1 = biggles.Curve([0,1], [0,0])
            expect2 = biggles.Curve([0,1], [0,0])

            arr[0,0].add(expect1)
            arr[1,0].add(expect2)





        arr.xrange = [0,1.4]
        if yrng is not None:
            arr.yrange = yrng

        wlog("Writing plot file:",epsfile)
        if show:
            arr.show()
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)


    def plot_ediff_Rshear_vs_e(self, 
                               docum=False,
                               s2max=None, 
                               yrng=None, 
                               yrng2=None,
                               title=None,
                               show=True):
        """
        Plot the measured ellipticity minus true vs true
        """
        import biggles
        import pcolors
        import converter

        extra=''
        if docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        data = self.read_data(docum=docum)

        doR=True
        if 'Rshear' not in data[0].dtype.names:
            doR=False

        epsfile_etot = shapesim.get_plot_file(self['run'],'etot'+extra,
                                              yrng=yrng)
        epsfile_R = shapesim.get_plot_file(self['run'],'Rshear'+extra,
                                              yrng=yrng2)
        wlog("will plot to:",epsfile_etot)
        wlog("will plot to:",epsfile_R)

        colors=pcolors.rainbow(len(data), 'hex')

        biggles.configure('PlotKey','key_vsep',1.0)
        biggles.configure('screen','width',1200)

        #tab=biggles.Table(1,2)

        plt1=biggles.FramedPlot()
        plt1.xlabel=r'e'
        plt1.ylabel = r'$\Delta e$'
        #plt1.x1.draw_ticklabels = 0

        plt2=biggles.FramedPlot()
        plt2.xlabel=r'e'
        plt2.ylabel = r'$\Delta R/R$'

        plt1.aspect_ratio=1
        plt2.aspect_ratio=1
 
        title=self.get_title(title=title)
        if title:
            plt1.title=title
            plt2.title=title

        plots1=[]
        plots2=[]
        for i,st in enumerate(reversed(data)):
            wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            s2 = median(st['s2'])

            etrue = st['etrue']
            emeas = st['e_meas']
            ediff = emeas  - etrue

            label = r'%0.3f' % s2
            cr = biggles.Curve(etrue, ediff, color=colors[i])
            cr.label = label
            plt1.add(cr)
            plots1.append(cr)

            if doR:
                #R_true = 1-.5*etrue**2
                R_true = st['Rshear_true']
                R_meas = st['Rshear']
                R_fdiff = R_meas/R_true-1

                Rp = biggles.Curve(etrue, R_fdiff, color=colors[i])

                plt2.add(Rp)
            #if i < 15:
            #    plots1.append(cr)
            #else:
            #    plots2.append(cr)


        fsize=2
        key1 = biggles.PlotKey(0.9,0.92, plots1, halign='right', 
                              fontsize=fsize)
        plt1.add(key1)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$'
        #if self['s2n'] > 0:
        #    klabtext += ' (S/N)'
        klab = biggles.PlotLabel(0.72,0.92,klabtext,
                                 fontsize=fsize,halign='right')
        plt1.add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.1,0.9, plab, halign='left')
        plt1.add(l)



        if self.simc['psfmodel'] == 'turb':
            siglab=r'$FWHM_{PSF}: %.1f$ pix' % self.simc['psf_fwhm']
        else:
            psf_sigma = self.simc['psf_sigma']
            siglab=r'$\sigma_{PSF}: %.1f$ pix' % psf_sigma
        siglab += ' '+self.psf_estring

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
        plt1.add(sl)

        plt1.xrange = [0,1.4]
        if yrng is not None:
            plt1.yrange = yrng

        wlog("plot to:",epsfile_etot)
        if show:
            plt1.show()
        plt1.write_eps(epsfile_etot)
        converter.convert(epsfile_etot,dpi=100,verbose=True)

        if doR:
            wlog("plot to:",epsfile_R)
            plt2.add(key1)
            plt2.add(klab)
            plt2.add(l)
            plt2.add(sl)
            plt2.xrange = [0,1.4]

            if yrng2 is not None:
                plt2.yrange = yrng2

            if show:
                plt2.show()
            plt2.write_eps(epsfile_R)

            converter.convert(epsfile_R,dpi=100,verbose=True)




    def read_data(self, docum=False):
        if self._data is None:
            wlog("reading data")
            """
            self._data = shapesim.read_all_outputs(self['run'],
                                                   skip1=skip1,skip2=skip2,
                                                   average=True,
                                                   docum=docum,
                                                   verbose=True)
            """
            self._data = shapesim.read_averaged_outputs(self['run'], 
                                                        docum=docum, 
                                                        skip1=self.skip1) 
        if self.s2min is not None:
            keepdata = self.limit_s2(self._data, self.s2min)
            return keepdata
        else:
            return self._data

    def limit_s2(self, datalist, s2min):
        out=[]
        for i,d in enumerate(datalist):
            s2 = median(d['s2'])
            if s2 > s2min:
                out.append(d)

        return out
