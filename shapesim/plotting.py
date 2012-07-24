import os
from . import shapesim
import lensing
import esutil as eu
from esutil.misc import wlog
from esutil.numpy_util import where1
from esutil.stat import sigma_clip
import numpy
from numpy import median, zeros, sqrt, array
import copy
import fitting

from lensing.util import shear_fracdiff, e2gamma, gamma2e, g1g2_to_e1e2

GRATIO=1.61803399
class MultiPlotterBase(dict):
    """
    Make plots for multiple runs.
    """
    def __init__(self, set, **keys):
        self.set = set
        self.set_runs()

        # for now
        keys['docum'] = True

        for k,v in keys.iteritems():
            self[k] = v

        self.plotters=[]
        for run in self.runs:
            self.plotters.append(SimPlotter(run,**keys))

        d = shapesim.get_plot_dir(set)
        if not os.path.exists(d):
            os.makedirs(d)

    def set_runs(self):
        if self.set == 'set-s2n-edg01':
            # elliptial psfs
            runs = ['gmix-fit-edg03r02',
                    'gmix-fit-edg04r01',
                    'gmix-fit-edg05r01',
                    'gmix-fit-edg06r01',
                    'gmix-fit-edg07r01',
                    'gmix-fit-edg08r01']
        elif self.set == 'set-s2n-edg02':
            # round psfs
            runs = ['gmix-fit-edg09r01',
                    'gmix-fit-edg02r02',
                    'gmix-fit-edg10r01',
                    'gmix-fit-edg11r01',
                    'gmix-fit-edg12r01',
                    'gmix-fit-edg13r01']

        elif self.set == 'set-e-gg01':
            runs = ['gmix-fit-gg04r09',
                    'gmix-fit-gg04r04',
                    'gmix-fit-gg04r05',
                    'gmix-fit-gg04r06',
                    'gmix-fit-gg04r07',
                    'gmix-fit-gg04r08']
        else:
            raise ValueError("don't know about set %s" % self.set)
        self.runs = runs

    def get_title(self):
        title=self.get('title',None)
        if title:
            title=title.replace('-',' ')
        return title

class MultiPlotterVsE(MultiPlotterBase):
    """
    Make delta shear plots as a function of true e; each
    panel will be from a different run with different S/N.
    """
    def __init__(self, set, **keys):
        super(MultiPlotterVsE,self).__init__(set, **keys)

    def doplots(self):
        import biggles
        import converter
        import pcolors
        
        #biggles.configure("screen","width",1100)
        biggles.configure("default","fontsize_min",1.5)
        #biggles.configure('_HalfAxis','ticklabels_style',{'fontsize':2.0})
        nrun=len(self.runs)
        
        nrow=2
        ncol=3
        arr = biggles.FramedArray(nrow,ncol)
        is2list = [5,11,17,19]
        n_is2 = len(is2list)
        #colors=pcolors.rainbow(n_is2, 'hex')
        colors = ['blue','magenta','green','red']

        if self['docum']:
            tag1='shear1cum'
            tag2='shear2cum'
            errtag1='shear1cum_err'
            errtag2='shear2cum_err'
        else:
            tag1='shear1'
            tag2='shear2'
            errtag1='shear1err'
            errtag2='shear2err'


        n_s2n = len(self.runs)
        td = self.plotters[0].read_data()
        n_s2 = len(td)
        ne = td[0].size

        #plt=biggles.FramedPlot()

        dt=[('g1true','f8',nrun),
            ('g2true','f8',nrun),
            ('g1meas','f8',nrun),
            ('g1err','f8',nrun),
            ('g2meas','f8',nrun),
            ('g2err','f8',nrun)]
        
        kplots=[]
        for i_s2n,plotter in enumerate(self.plotters):
            irow = i_s2n / ncol
            icol = i_s2n % ncol
            td = plotter.read_data()
            etrue = td[0]['etrue']

            shear_true = plotter.get_shear_true()

            s2n = plotter['s2n']

            for i_is2 in xrange(n_is2):
                st = td[is2list[i_is2]]
                s2 = st['s2'].mean()

                diff1 = st[tag1] - shear_true.g1
                diff2 = st[tag2] - shear_true.g2

                p1 = biggles.Points(etrue,diff1,
                                    type='filled circle', 
                                    color=colors[i_is2])
                c1 = biggles.Curve(etrue,diff1,color=colors[i_is2])

                arr[irow,icol].add(p1,c1) 

                if i_s2n == 0:
                    label = '%.2f' % s2
                    if self['docum']:
                        label = '< '+label
                    c1.label = label
                    kplots.append(c1)

                if i_is2 == 0:
                    g1err  = st[errtag1]
                    g2err  = st[errtag2]
                    perr1 = biggles.SymmetricErrorBarsY(etrue,
                                                        diff1,g1err,
                                                        color=colors[i_is2])
                    arr[irow,icol].add(perr1)

                    z1=biggles.Curve(etrue,zeros(etrue.size))
                    arr[irow,icol].add(z1)

            if i_s2n == 0:
                g1labs = r'$\gamma_1: %.2g$' % shear_true.g1
                g2labs = r'$\gamma_1: %.2g$' % shear_true.g2
                g1lab = biggles.PlotLabel(0.1,0.9,g1labs,halign='left')
                arr[irow,icol].add(g1lab)

            if s2n > 1000:
                ls2n = numpy.log10(s2n)
                ls2n = r'$10^{%.1f}$' % ls2n
            else:
                ls2n = '%.0f' % s2n

            s2nlab = biggles.PlotLabel(0.9,0.9,'S/N: %s' % ls2n,
                                     fontsize=2.5,halign='right')
            arr[irow,icol].add(s2nlab)


        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$'
        klab = biggles.PlotLabel(0.5,0.3,klabtext,
                                 halign='right')

        key=biggles.PlotKey(0.9,0.3,kplots,halign='right')
        arr[0,0].add(klab,key)


        yrng=self.get('yrange',None)
        if yrng:
            arr.yrange = yrng
        arr.xrange=array([0.01,0.85])
        arr.uniform_limits=1
        arr.xlabel = r'$e_{true}$'
        arr.ylabel = r'$\Delta \gamma_1$'
        #arr.aspect_ratio=1/1.61803399
        arr.aspect_ratio=1/1.4

        title=self.get_title()
        if title:
            arr.title=title
        arr.show()

        epsfile = shapesim.get_plot_file(self.set,'vs-e',yrng=yrng)
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)


class MultiPlotterVsShear(MultiPlotterBase):
    """
    Make delta shear plots as a function of true shear.
    uou should feed this run with all the same models for
    object, and same values of S/N, with the only difference
    the shear

    For 
        \gamma_{meas} = \gamma_{true}*m + b

    then the diff gives
       
        \Delta \gamma = \gamma_{true} (m-1) + b
    """
    def __init__(self, set, **keys):
        super(MultiPlotterVsShear,self).__init__(set, **keys)

    def doplots(self):
        import biggles
        import converter
        import pcolors
        
        scale=.01
        #biggles.configure("screen","width",2000)
        #biggles.configure("screen","height",1100)
        biggles.configure("default","fontsize_min",1.2)
        #biggles.configure('_HalfAxis','ticklabels_style',{'fontsize':2.0})
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        biggles.configure('PlotKey','key_vsep',3)
        biggles.configure('PlotKey','key_width',15)
        biggles.configure('_ErrorBar','barsize',1)
        nrun=len(self.runs)
        
        td = self.plotters[0].read_data()
        n_s2 = len(td)
        n_s2n = td[0].size
        if n_s2 != 4:
            raise ValueError("adapt for n_s2 != 4")

        nrow=5
        ncol=3
        arr = biggles.FramedArray(nrow,ncol)
        # colors of each s2 bin
        colors = ['blue','magenta','green','red']
        #linetypes=['dotted','dashed','dotdashed','solid']
        linetypes=['solid','dotdashed','dashed','dotted']

        s2n_name='s2n_matched'
        tag1='shear1'
        tag2='shear2'
        errtag1='shear1err'
        errtag2='shear2err'



        #plt=biggles.FramedPlot()

        dt=[('g1true','f8',nrun),
            ('g2true','f8',nrun),
            ('g1meas','f8',nrun),
            ('g1err','f8',nrun),
            ('g2meas','f8',nrun),
            ('g2err','f8',nrun)]
        
        #for i_s2n in xrange(n_s2n):

        for i in xrange(nrow*ncol):
            irow = i / ncol
            icol = i % ncol
            arr[irow,icol].add(biggles.Curve([-1000]*2,[-1000]*2))
            #arr[irow,icol].yrange = [-0.0025,0.0025]
            arr[irow,icol].yrange = [-0.035,0.035]

        fdtype=[('s2','f8',n_s2n), ('s2n','f8',n_s2n),
                ('c','f8',n_s2n), ('cerr','f8',n_s2n),
                ('m','f8',n_s2n),('merr','f8',n_s2n)]
        fits1 = zeros(n_s2, dtype=fdtype)
        fits2 = zeros(n_s2, dtype=fdtype)

        fcurves=[]
        nplot=0
        is2ns = list(reversed(xrange(n_s2n)))
        #for i_s2n in xrange(n_s2n):
        for i_s2n in is2ns:
            #irow = i_s2n / ncol
            #icol = i_s2n % ncol
            irow = nplot / ncol
            icol = nplot % ncol
            td = self.plotters[0].read_data()
            s2n = td[0][s2n_name][i_s2n]

            data = zeros(n_s2,dtype=dt)
            for irun,plotter in enumerate(self.plotters):
                d = plotter.read_data()
                shear = plotter.get_shear_true()

                # redundant
                
                s2vals=[]
                for is2 in xrange(n_s2):

                    s2,ellip = shapesim.get_s2_e(plotter.simc, is2, 0)
                    s2vals.append(s2)

                    data['g1true'][is2,irun] = shear.g1
                    data['g2true'][is2,irun] = shear.g2
                    data['g1meas'][is2,irun] = d[is2][tag1][i_s2n]
                    data['g2meas'][is2,irun] = d[is2][tag2][i_s2n]
                    data['g1err'][is2,irun] = d[is2][errtag1][i_s2n]
                    data['g2err'][is2,irun] = d[is2][errtag2][i_s2n]


            for is2 in xrange(n_s2):
                s2=s2vals[is2]

                g1true = data['g1true'][is2,:]
                g2true = data['g2true'][is2,:]
                diff1  = data['g1meas'][is2,:] - g1true
                diff2  = data['g2meas'][is2,:] - g2true
                p1 = biggles.Points(g1true/scale,diff1,
                                    type='filled circle', color=colors[is2])
                c1 = biggles.Curve(g1true/scale,diff1,color=colors[is2],
                                   type=linetypes[is2])
                if i_s2n == (n_s2n-1):
                    cfake = biggles.Curve(g1true-1000,diff1,color=colors[is2],
                                          type=linetypes[is2])
                    label = '%.2f' % s2
                    cfake.label = label
                    fcurves.append(cfake)
                arr[irow,icol].add(p1,c1) 

                g1err  = data['g1err'][is2,:]
                g2err  = data['g2err'][is2,:]
                perr1 = biggles.SymmetricErrorBarsY(g1true/scale,
                                                    diff1,g1err,
                                                    color=colors[is2])
                arr[irow,icol].add(perr1)
                if is2 == 0:
                    z1=biggles.Curve([-50,50], [0,0])
                    arr[irow,icol].add(z1)


                linfit1 = fitting.LineFitter(g1true, diff1, g1err)
                linfit2 = fitting.LineFitter(g1true, diff2, g2err)
                fits1['s2'][is2,i_s2n] = s2
                fits1['s2n'][is2,i_s2n] = s2n
                fits1['m'][is2,i_s2n] = linfit1.pars[0] + 1
                fits1['merr'][is2,i_s2n] = linfit1.perr[0]
                fits1['c'][is2,i_s2n] = linfit1.pars[1]
                fits1['cerr'][is2,i_s2n] = linfit1.perr[1]
                fits2['s2'][is2,i_s2n] = s2
                fits2['s2n'][is2,i_s2n] = s2n
                fits2['m'][is2,i_s2n] = linfit2.pars[0] +1
                fits2['merr'][is2,i_s2n] = linfit2.perr[0]
                fits2['c'][is2,i_s2n] = linfit2.pars[1]
                fits2['cerr'][is2,i_s2n] = linfit2.perr[1]

                if is2 == 0 and nplot==0:
                    print '%15s %15s %15s +/- %15s %15s +/- %15s' % ('s2','s2n','m','err','c','err')
                print '%15s %15s %15g +/- %15g %15g +/- %15g' % (s2,s2n,
                                                                 linfit1.pars[0]+1,linfit1.perr[0],
                                                                 linfit1.pars[1],linfit1.perr[1])

            s2nlab = biggles.PlotLabel(0.9,0.9,'S/N: %d' % s2n,
                                     fontsize=2.5,halign='right')
            arr[irow,icol].add(s2nlab)
            if irow == 3:
                arr[irow,icol].yrange = [-0.012,0.012]
            elif irow == 2:
                arr[irow,icol].yrange = [-0.005,0.005]
            elif irow == 1:
                arr[irow,icol].yrange = [-0.0035,0.0035]
            elif irow == 0:
                arr[irow,icol].yrange = [-0.0035,0.0035]

            nplot+=1

        """
        nleft = nrow*ncol-nplot
        for i in xrange(nplot,nleft):
            irow = i / ncol
            icol = i % ncol
            arr[irow,icol].add(biggles.Curve([-1000]*2,[-1000]*2))
            arr[irow,icol].yrange = [-0.0025,0.0025]
        """
        fsize=2
        key=biggles.PlotKey(0.85,0.9,fcurves,halign='right',fontsize=fsize)
        arr[nrow-1,ncol-1].add(key, *fcurves)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$: '
        klab = biggles.PlotLabel(0.6,0.9,klabtext,
                                 fontsize=fsize,halign='right')
        arr[nrow-1,ncol-1].add(klab)

        simc = self.plotters[0].simc
        objmodel = simc['objmodel']
        psfmodel = simc['psfmodel']
        plab='%s %s' % (objmodel,psfmodel)

        lowest=0.15
        l = biggles.PlotLabel(0.075,lowest+0, plab, halign='left')
        arr[nrow-1,ncol-1].add(l)

        if simc['psfmodel'] == 'turb':
            siglab = r'$FWHM: %.1f$ pix' % simc['psf_fwhm']
        else:
            psf_sigma = simc['psf_sigma']
            siglab = r'$\sigma: %.1f$ pix' % psf_sigma
        #siglab += ' '+self.plotters[0].psf_estring
        elab = r'$e_{gal}^{tot}: %.2f$' % td[0]['etrue'].mean()

        sl = biggles.PlotLabel(0.075,lowest+0.45, siglab, halign='left', 
                               fontsize=2.5)
        psf_estring=self.plotters[0].psf_estring
        if psf_estring:
            psfel = biggles.PlotLabel(0.075,lowest+0.3, psf_estring, 
                                      halign='left', 
                                      fontsize=2.5)
            arr[nrow-1,ncol-1].add(psfel)

        el = biggles.PlotLabel(0.075,lowest+0.15, elab, halign='left', 
                               fontsize=2.5)
        arr[nrow-1,ncol-1].add(sl,el)



        #yrng=self.get('yrange',None)
        #if yrng:
        #    arr.yrange = yrng
        yrng=None
        arr.xrange=array([-0.005,0.059])/scale
        #arr.uniform_limits=1
        arr.xlabel = r'$\gamma_{true}/%.2g$' % scale
        arr.ylabel = r'$\Delta \gamma$'
        #arr.aspect_ratio=1/1.61803399
        #arr.aspect_ratio=1/1.4
        arr.aspect_ratio=1.2

        title=self.get_title()
        if title:
            arr.title=title
        arr.show()

        epsfile = shapesim.get_plot_file(self.set,'vs-shear',yrng=yrng)
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)

        marr = biggles.FramedArray(1,2)
        carr = biggles.FramedArray(1,2)

        kplots=[]
        for is2 in xrange(n_s2):
            mpts1=biggles.Points(fits1['s2n'][is2,:], fits1['m'][is2,:],
                                 color=colors[is2])
            cpts1=biggles.Points(fits1['s2n'][is2,:], fits1['c'][is2,:],
                                 color=colors[is2])
            mcur1=biggles.Curve(fits1['s2n'][is2,:], fits1['m'][is2,:],
                                color=colors[is2],type=linetypes[is2])
            ccur1=biggles.Curve(fits1['s2n'][is2,:], fits1['c'][is2,:],
                                color=colors[is2],type=linetypes[is2])

            merr1=biggles.SymmetricErrorBarsY(fits1['s2n'][is2,:], 
                                              fits1['m'][is2,:],
                                              fits1['merr'][is2,:],
                                              color=colors[is2])
            cerr1=biggles.SymmetricErrorBarsY(fits1['s2n'][is2,:], 
                                              fits1['c'][is2,:],
                                              fits1['cerr'][is2,:],
                                              color=colors[is2])
            mpts2=biggles.Points(fits2['s2n'][is2,:], fits2['m'][is2,:],
                                 color=colors[is2])
            cpts2=biggles.Points(fits2['s2n'][is2,:], fits2['c'][is2,:],
                                 color=colors[is2])
            mcur2=biggles.Curve(fits2['s2n'][is2,:], fits2['m'][is2,:],
                                color=colors[is2],type=linetypes[is2])
            ccur2=biggles.Curve(fits2['s2n'][is2,:], fits2['c'][is2,:],
                                color=colors[is2],type=linetypes[is2])



            merr2=biggles.SymmetricErrorBarsY(fits2['s2n'][is2,:], 
                                              fits2['m'][is2,:],
                                              fits2['merr'][is2,:],
                                              color=colors[is2])
            cerr2=biggles.SymmetricErrorBarsY(fits2['s2n'][is2,:], 
                                              fits2['c'][is2,:],
                                              fits2['cerr'][is2,:],
                                              color=colors[is2])


            label = '%.2f' % s2vals[is2]
            ccur1.label = label
            ccur2.label = label
            kplots.append(ccur1)

            marr[0,0].add(mpts1,merr1,mcur1)
            marr[0,1].add(mpts2,merr2,mcur2)

            carr[0,0].add(cpts1,cerr1,ccur1)
            carr[0,1].add(cpts2,cerr2,ccur2)

        lab1 = biggles.PlotLabel(0.9,0.1,r'$\gamma_1$',halign='right')
        lab2 = biggles.PlotLabel(0.9,0.1,r'$\gamma_2$',halign='right')
        marr[0,0].add(lab1)
        marr[0,1].add(lab2)
        carr[0,0].add(lab1)
        carr[0,1].add(lab2)

        xlabel = r'$S/N_{matched}$'
        marr.xlabel = xlabel
        marr.ylabel = 'm (calibration bias)'
        marr.aspect_ratio=1/GRATIO
        carr.xlabel = xlabel
        carr.ylabel = 'c (additive bias)'
        carr.aspect_ratio=1/GRATIO

        key=biggles.PlotKey(0.85,0.9,kplots,halign='right',fontsize=fsize)
        carr[0,1].add(key)
        marr[0,1].add(key)

        carr[0,0].add(biggles.Curve([-50,500],[0,0]))
        carr[0,1].add(biggles.Curve([-50,500],[0,0]))
        marr[0,0].add(biggles.Curve([-50,500],[1,1]))
        marr[0,1].add(biggles.Curve([-50,500],[1,1]))

        marr.yrange=[.8,1.2]

        xrng=[-1,110]
        carr.xrange=xrng
        marr.xrange=xrng

        marr.show()
        carr.show()
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
            self.psf_estring = r'$e_{psf}$: %.2g,%.2g' % (self.psf_e1,self.psf_e2)

        d = shapesim.get_plot_dir(run)
        if not os.path.exists(d):
            os.makedirs(d)
        
        self.noerr = keys.get('noerr',False)
        self._data=None

        title=copy.deepcopy(self['run'])
        title = title.replace('gmix-fit','gmix')
        self.title = title.replace('-',' ')
        self.maketitle=keys.get('maketitle',False)

        self.skip1 = keys.get('skip1',[])
        self.skip2 = keys.get('skip2',[])

        self.s2min = keys.get('s2min',None)

        self.docum = keys.get('docum',False)

    def get_title(self, title=None):
        if title is None and self.maketitle:
            title=self.title
        return title

    def get_shear_true(self):
        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])
        return shear_true

    def plots_shear_vs_s2n(self, 
                           type='diff',
                           xrng=None,
                           yrng=None, 
                           doavg=False,
                           title=None,
                           show=True):
        import biggles
        import pcolors
        import converter

        biggles.configure("default","fontsize_min",2)

        if doavg:
            extra='-avg'
        else:
            extra=''

        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        shear_true = self.get_shear_true()

        data = self.read_data()

        epsfile = shapesim.get_plot_file(self['run'],type+extra,yrng=yrng)
        wlog("will plot to:",epsfile)

        if len(data) == 4:
            colors=['red','forestgreen','NavajoWhite3','blue']
            types=['dotted','dashed','dotdashed','solid']
        else:
            colors=pcolors.rainbow(len(data), 'hex')
            types=['solid']*len(data)

        biggles.configure('PlotKey','key_vsep',1.0)
        biggles.configure("default","fontsize_min",1.5)
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
                if self.docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='gamma1_meas'
                    tag2='gamma2_meas'
                # convention
                st[tag1] = -st[tag1]
            else:
                if self.docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='shear1'
                    tag2='shear2'

            if self.docum and i == (len(data)-1):
                # save for later
                pts1_0 = st[tag1].copy()
                pts2_0 = st[tag2].copy()
                if type == 'diff':
                    pts1_0 -= shear_true.g1
                    pts2_0 -= shear_true.g2

                err1_0 = st['shear1cum_err']
                err2_0 = st['shear2cum_err']
          
            if type == 'diff':
                yvals1 = st[tag1] - shear_true.g1
                yvals2 = st[tag2] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1]
                yvals2 = st[tag2]

            else:
                raise ValueError("bad plot type: '%s'" % type)

            if self.docum:
                label = r'< %0.3f' % s2
            else:
                label = r'%0.3f' % s2
            pr1 = biggles.Points(s2n, yvals1, color=colors[i],type='filled circle',size=1.5)
            pr2 = biggles.Points(s2n, yvals2, color=colors[i],type='filled circle',size=1.5)
            cr1 = biggles.Curve(s2n, yvals1, color=colors[i],type=types[i],width=2.5)
            cr2 = biggles.Curve(s2n, yvals2, color=colors[i],type=types[i],width=2.5)
            cr1.label = label
            cr2.label = label


            arr[0,0].add(cr1,pr1)
            arr[1,0].add(cr2,pr2)
            
            #if not self.docum and i == (len(data)-1):
            if True:
                g1err = [st['shear1err'].max()]*st['shear1err'].size
                g2err = [st['shear2err'].max()]*st['shear2err'].size
                err1p = biggles.SymmetricErrorBarsY(s2n, yvals1, g1err,
                                                    color=colors[i])
                err2p = biggles.SymmetricErrorBarsY(s2n, yvals2, g2err,
                                                    color=colors[i])
                if not self.noerr:
                    arr[0,0].add(err1p)
                    arr[1,0].add(err2p)

            if i < 15:
                plots1.append(cr1)
            else:
                plots2.append(cr1)

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
            avg1 = biggles.Points(avg['s2n'],avg['shear1'],
                                  type='filled circle',size=2)
            avg2 = biggles.Points(avg['s2n'],avg['shear2'],
                                  type='filled circle',size=2)
            avg1.label = 'average'

            if len(plots2) > 0:
                plots2.append(avg1)
            else:
                plots1.append(avg1)

        if self.docum and not self.noerr:
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
        klab = biggles.PlotLabel(0.715,0.92,klabtext,
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
            xrng = [0,s2n.max()*1.4]
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
                        title=None,
                        show=True):
        import biggles
        import pcolors
        import converter

        if doavg:
            extra='-avg'
        else:
            extra=''
        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        slist=self.simc['shear']
        shear_true = lensing.shear.Shear(g1=slist[0],g2=slist[1])

        data = self.read_data()

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
                if self.docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='gamma1_meas'
                    tag2='gamma2_meas'
                # convention
                st[tag1] = -st[tag1]
            else:
                if self.docum:
                    tag1='shear1cum'
                    tag2='shear2cum'
                else:
                    tag1='shear1'
                    tag2='shear2'

            if self.docum and i == (len(data)-1):
                # save for later
                pts1_0 = st[tag1].copy()
                pts2_0 = st[tag2].copy()
                if type == 'diff':
                    pts1_0 -= shear_true.g1
                    pts2_0 -= shear_true.g2
                err1_0 = st['shear1cum_err']
                err2_0 = st['shear2cum_err']
 
            if type == 'diff':
                yvals1 = st[tag1] - shear_true.g1
                yvals2 = st[tag2] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1]
                yvals2 = st[tag2]

            else:
                raise ValueError("bad plot type: '%s'" % type)


            if self.docum:
                label = r'< %0.3f' % s2
            else:
                label = r'%0.3f' % s2

            cr1 = biggles.Curve(etrue, yvals1, color=colors[i])
            cr2 = biggles.Curve(etrue, yvals2, color=colors[i])
            cr1.label = label
            cr2.label = label

            arr[0,0].add(cr1)
            arr[1,0].add(cr2)

            if not self.docum and i == (len(data)-1) and not self.noerr:
                err1p = biggles.SymmetricErrorBarsY(etrue, yvals1, st['shear1err'], 
                                                    width=4)
                err2p = biggles.SymmetricErrorBarsY(etrue, yvals2, st['shear2err'],
                                                    width=4)
                arr[0,0].add(err1p)
                arr[1,0].add(err2p)



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
            avg1 = biggles.Points(avg['e'],avg['shear1'],
                                  type='filled circle',size=2)
            avg2 = biggles.Points(avg['e'],avg['shear2'],
                                  type='filled circle',size=2)
            avg1.label = 'average'

            if len(plots2) > 0:
                plots2.append(avg1)
            else:
                plots1.append(avg1)

        if self.docum and not self.noerr:
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
        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'byellip':
            raise ValueError("Can only make plots vs e for 'byellip' runs")

        data = self.read_data()

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




    def read_data(self):
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
                                                        docum=self.docum, 
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
