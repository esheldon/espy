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
import fitsio
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

        for k,v in keys.iteritems():
            self[k] = v

        if 'show' not in self:
            self['show'] = True

        self.fill_color='grey80'

        self.plotters=[]
        for run in self.runs:
            self.plotters.append(SimPlotter(run,**keys))

        d = shapesim.get_plot_dir(set)
        if not os.path.exists(d):
            os.makedirs(d)

    def set_runs(self):
        if self.set == 'set-s2n-edg01':
            raise ValueError("don't do this plot")
            # elliptial psfs
            runs = ['gmix-fit-edg03r02',
                    'gmix-fit-edg04r01',
                    'gmix-fit-edg05r01',
                    'gmix-fit-edg06r01',
                    'gmix-fit-edg07r01',
                    'gmix-fit-edg08r01']
        elif self.set == 'set-s2n-edg02':
            # vs s2n but also calibration error
            # round psfs
            runs = ['gmix-fit-edg09r01',
                    'gmix-fit-edg02r02',
                    'gmix-fit-edg10r01',
                    'gmix-fit-edg11r01',
                    'gmix-fit-edg12r01',
                    'gmix-fit-edg13r01']
        elif self.set == 'set-s2n-edg03':
            # vs s2n but also calibration error
            # round psfs
            # admom s/n and delta function
            runs = ['gmix-fit-edg09r02',
                    'gmix-fit-edg02r03',
                    'gmix-fit-edg10r02',
                    'gmix-fit-edg11r02',
                    'gmix-fit-edg12r02',
                    'gmix-fit-edg13r02']
        elif self.set == 'set-s2n-edg04':
            # same as set edg03 but more
            # runs averaged in
            runs=['gmix-fit-edg09r02r03',
                  'gmix-fit-edg02r03r04',
                  'gmix-fit-edg10r02r03',
                  'gmix-fit-edg11r02r03',
                  'gmix-fit-edg12r02r03',
                  'gmix-fit-edg13r02r03']

        elif self.set == 'set-s2n-et01':
            # also calibration error
            runs = ['gmix-fit-et05r01',
                    'gmix-fit-et06r01',
                    'gmix-fit-et07r01',
                    'gmix-fit-et08r01',
                    'gmix-fit-et09r01',
                    'gmix-fit-et10r01']
        elif self.set == 'set-s2n-et02':
            # also calibration error
            # this one admom s/n and delta func
            runs = ['gmix-fit-et05r02',
                    'gmix-fit-et06r02',
                    'gmix-fit-et07r02',
                    'gmix-fit-et08r02',
                    'gmix-fit-et09r02',
                    'gmix-fit-et10r02']
        elif self.set == 'set-s2n-et03':
            # also calibration error
            # this one admom s/n and delta func
            runs = ['gmix-fit-et05r02r03',
                    'gmix-fit-et06r02r03',
                    'gmix-fit-et07r02r03',
                    'gmix-fit-et08r02r03',
                    'gmix-fit-et09r02r03',
                    'gmix-fit-et10r02r03']

        elif self.set == 'set-e-gg01':
            # vs galaxy ellipticity
            runs = ['gmix-fit-gg05r02',
                    'gmix-fit-gg05r03',
                    'gmix-fit-gg05r04',
                    'gmix-fit-gg05r05',
                    'gmix-fit-gg05r06',
                    'gmix-fit-gg05r07']
        elif self.set == 'set-e-gg02':
            # vs galaxy ellipticity
            # admom S/N
            runs = ['gmix-fit-gg05r08',
                    'gmix-fit-gg05r09',
                    'gmix-fit-gg05r10',
                    'gmix-fit-gg05r11',
                    'gmix-fit-gg05r12',
                    'gmix-fit-gg05r13']
        elif self.set == 'set-e-edg01':
            runs = ['gmix-fit-edg07r02',
                    'gmix-fit-edg07r03',
                    'gmix-fit-edg07r04',
                    'gmix-fit-edg07r05',
                    'gmix-fit-edg07r06',
                    'gmix-fit-edg07r07']

        elif self.set == 'set-epsf-edg01':
            # elliptical psfs, matched S/N
            runs = ['gmix-fit-edg14r01',
                    'gmix-fit-edg03r02',
                    'gmix-fit-edg04r01',
                    'gmix-fit-edg05r01',
                    'gmix-fit-edg06r01',
                    'gmix-fit-edg07r01',
                    'gmix-fit-edg08r01']
        elif self.set == 'set-epsf-edg02':
            # elliptical psfs, admom S/N, delta func
            runs = ['gmix-fit-edg14r02',
                    'gmix-fit-edg03r05',
                    'gmix-fit-edg04r02',
                    'gmix-fit-edg05r02',
                    'gmix-fit-edg06r02',
                    'gmix-fit-edg07r08',
                    'gmix-fit-edg08r02']

        elif self.set == 'set-nbias01':
            # make sure these are all admom S/N, delta
            # function
            runs = ['gmix-fit-et05r02r03',
                    'gmix-fit-et06r02r03',
                    'gmix-fit-et07r02r03',
                    'gmix-fit-et08r02r03',
                    'gmix-fit-et09r02r03',
                    'gmix-fit-et10r02r03',
                    'gmix-fit-dt03r05r06r07r08']
        else:
            raise ValueError("don't know about set %s" % self.set)

        self.runs = runs

    def get_title(self):
        title=self.get('title',None)
        if title:
            title=title.replace('-',' ')
        return title

class MultiPlotterNoiseBias(MultiPlotterBase):
    def __init__(self, set, **keys):
        super(MultiPlotterNoiseBias,self).__init__(set, **keys)
        self.do_setup()

    def do_setup(self):
        td = self.plotters[0].read_data()
        self.n_s2 = len(td)
        self.n_s2n = td[0].size
        if self.n_s2 != 4:
            raise ValueError("adapt for n_s2 != 4")
        #self.colors = ['blue','magenta','green','red']
        self.colors=list(reversed(['red','forestgreen','NavajoWhite3','blue']))
        self.linetypes=['solid','dotdashed','dashed','dotted']
        self.point_types=['filled circle','filled diamond','filled square','filled triangle']

    def doplots(self):
        import biggles
        import pcolors
        s2n_name='s2n_meas'
        s2n_lab=r'$(S/N)_{uw}$'
        ylabel = r'$\Delta \gamma$'
        yrange=[-0.005,0.005]
        fsize=2

        bias_file=shapesim.get_bias_file('set-s2n-et03', 'bias')
        fits_g1=eu.io.read(bias_file,ext=1)
        n_s2 = fits_g1.size
        n_s2n = fits_g1['s2n_meas'].size

        if n_s2 != 4:
            raise ValueError("expected n_s2==4 for fits")

        """
        dt=[('s2n','f8',n_s2n),
            ('s2n_meas','f8',n_s2n),
            ('g1corr','f8',n_s2n),
            ('g2corr','f8',n_s2n)]
 
        
        corrdata=zeros(n_s2, dtype=dt)
        """

        colors=pcolors.rainbow(len(self.runs), 'hex')


        for i_s2 in xrange(n_s2):
            plt=biggles.FramedPlot()
            s2data={'s2n_meas':[],
                    'shear1corr':[],
                    'shear1diff':[],
                    'shear1err':[]}

            # we will only interpolate within this range
            s2n_fit = fits_g1[s2n_name][i_s2]
            m_fit= 1+fits_g1['m'][i_s2]

            s2 = fits_g1['s2'][i_s2][0]

            runlist=[]
            for irun,plotter in enumerate(self.plotters):
                run=plotter['run']
                runlist.append(run)

                d = plotter.read_data()
                if len(d) != 4:
                    raise ValueError("expected n_s2==4 for data")

                shear = plotter.get_shear_true()


                s2n = d[i_s2][s2n_name]
                print 's2:',s2
                print 's2vals:',d[i_s2]['s2']
                s2n_min = s2n_fit.min() - 0.2
                s2n_max = s2n_fit.max() + 0.2
                w=where1(  (s2n >= s2n_min) & (s2n <= s2n_max) )
                if w.size > 0:
                    #m = 1+eu.stat.interplin(m_fit,s2n_fit,s2n[w])
                    m = numpy.interp(s2n[w], s2n_fit, m_fit)
                    print 'mvals:',m
                    shear1corr = d[i_s2]['shear1'][w]/m
                    shear1err = d[i_s2]['shear1err'][w]
                    shear1diff = shear1corr-shear.g1

                    pts = biggles.Points(s2n[w],shear1diff,
                                         type='filled circle',
                                         color=colors[irun])
                    perr=biggles.SymmetricErrorBarsY(s2n[w],shear1diff,shear1err,
                                                     type='filled circle',
                                                     color=colors[irun])
 
                    plt.add(pts,perr)

                    s2data['s2n_meas'] += list(s2n[w])
                    s2data['shear1corr'] += list(shear1corr)
                    s2data['shear1diff'] += list(shear1diff)
                    s2data['shear1err'] += list(shear1err)

            if len(s2data['shear1corr']) == 0:
                wlog(" ==> interpolated no values!!")
            else:

                s2n_pmax = 1.05*array(s2data['s2n_meas']).max()

                """
                pts = biggles.Points(s2data['s2n_meas'],
                                     s2data['shear1diff'],
                                     type='filled circle')
                perr=biggles.SymmetricErrorBarsY(s2data['s2n_meas'],
                                                 s2data['shear1diff'],
                                                 s2data['shear1err'])
                plt.add(pts,perr)
                """
                plt.xlabel = s2n_lab
                plt.ylabel = ylabel
                plt.yrange=yrange

                plt.add(biggles.Curve([0,s2n_pmax],[0,0]))

                klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$:'
                plab=biggles.PlotLabel(0.9,0.9,klabtext + ': %.2f' % s2,
                                      halign='right',fontsize=fsize)

                plt.add(plab)
                if self['show']:
                    plt.show()


class MultiPlotterVsE(MultiPlotterBase):
    """
    Make delta shear plots as a function of true e; each
    panel will be from a different run with different S/N.
    """
    def __init__(self, set, **keys):
        super(MultiPlotterVsE,self).__init__(set, **keys)
        self.do_setup()
        self['use_rb'] = keys.get('use_rb',True)

    def do_setup(self):
        import pcolors
        td = self.plotters[0].read_data()
        self.n_s2n = len(self.runs)
        self.n_s2 = len(td)
        self.ne = td[0].size
        if self.n_s2 != 4:
            raise ValueError("adapt for n_s2 != 4")
        #self.colors = ['blue','magenta','green','red']
        self.colors=list(reversed(['red','forestgreen','NavajoWhite3','blue']))
        self.linetypes=['solid','dotdashed','dashed','dotted']
        self.point_types=['filled circle','filled diamond','filled square','filled triangle']


    def doplots(self):
        import biggles
        import converter
        

        #biggles.configure("screen","width",1100)
        biggles.configure("default","fontsize_min",1.5)
        biggles.configure('PlotKey','key_vsep',3.0)
        biggles.configure('PlotKey','key_width',15)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        #biggles.configure('_HalfAxis','ticklabels_style',{'fontsize':2.0})
        nrun=len(self.runs)
        
        nrow=2
        ncol=3
        arr = biggles.FramedArray(nrow,ncol)

        tag1='shear1'
        tag2='shear2'
        errtag1='shear1err'
        errtag2='shear2err'


        td = self.plotters[0].read_data()

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

            for is2 in xrange(self.n_s2):
                st = td[is2]
                s2 = st['s2'].mean()

                diff1 = st[tag1] - shear_true.g1
                diff2 = st[tag2] - shear_true.g2

                p1 = biggles.Points(etrue,diff1,
                                    type=self.point_types[is2],
                                    color=self.colors[is2],
                                    size=3)
                c1 = biggles.Curve(etrue,diff1,
                                   type=self.linetypes[is2],
                                   color=self.colors[is2],
                                   width=5.)

                arr[irow,icol].add(p1,c1) 

                if i_s2n == 0:
                    if self['use_rb']:
                        label = '%.2f' % sqrt(1/s2)
                    else:
                        label = '%.2f' % s2

                    c1.label = label
                    kplots.append(c1)

                g1err  = st[errtag1]
                g2err  = st[errtag2]
                perr1 = biggles.SymmetricErrorBarsY(etrue,
                                                    diff1,g1err,
                                                    color=self.colors[is2])
                arr[irow,icol].add(perr1)

                z1=biggles.Curve([-50,50],[0,0])
                arr[irow,icol].add(z1)

            #if i_s2n == 0:
            #    g1labs = r'$\gamma_1: %.2g$' % shear_true.g1
            #    g2labs = r'$\gamma_2: %.2g$' % shear_true.g2
            #    g1lab = biggles.PlotLabel(0.1,0.9,g1labs,halign='left')
            #    arr[irow,icol].add(g1lab)
            if i_s2n==0:
                psf_e1=self.plotters[0].simc.get('psf_e1',None)
                psf_e2=self.plotters[0].simc.get('psf_e2',None)
                if psf_e1 is not None:
                    psfetxt = r'$e_{PSF}: %.2g,%.2g$' % (psf_e1,psf_e2)
                    psfelab = biggles.PlotLabel(.9,.8,psfetxt,halign='right')
                    arr[irow,icol].add(psfelab)

            if s2n > 1000:
                ls2n = numpy.log10(s2n)
                ls2n = r'$10^{%.1f}$' % ls2n
            else:
                ls2n = '%.0f' % s2n

            s2nlab = biggles.PlotLabel(0.9,0.9,'S/N: %s' % ls2n,
                                     fontsize=2.5,halign='right')
            arr[irow,icol].add(s2nlab)

        if self['use_rb']:
            klabtext=r'$\sigma_{gal}/\sigma_{psf}$:'
        else:
            klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$:'

        klab = biggles.PlotLabel(0.45,0.3,klabtext,
                                 halign='right')

        key=biggles.PlotKey(0.8,0.3,kplots,halign='right')
        arr[0,0].add(klab,key)


        yrng=self.get('yrange',None)
        if yrng:
            arr.yrange = yrng
        arr.xrange=array([0.01,0.85])
        arr.uniform_limits=1
        arr.xlabel = r'$e_{true}$'
        arr.ylabel = r'$\Delta \gamma$'
        #arr.aspect_ratio=1/1.61803399
        arr.aspect_ratio=1/1.4

        title=self.get_title()
        if title:
            arr.title=title
        if self['show']:
            arr.show()

        extra=''
        if self['use_rb']:
            extra='-rb'
        epsfile = shapesim.get_plot_file(self.set,'vs-e'+extra,yrng=yrng)
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
        self.do_setup()

        self['use_rb'] = keys.get('use_rb',True)

    def do_setup(self):
        import pcolors

        try:
            ii=self.runs.index('gmix-fit-et05r01')
            newdata = []
            data=self.plotters[ii].read_data()
            for d in data:
                d=d[1:]
                newdata.append(d)
            self.plotters[ii]._data=newdata
        except:
            pass

        td = self.plotters[0].read_data()
        self.n_s2 = len(td)
        self.n_s2n = td[0].size
        if self.n_s2 != 4:
            raise ValueError("adapt for n_s2 != 4")
        #self.colors = ['blue','magenta','green','red']
        self.colors=list(reversed(['red','forestgreen','NavajoWhite3','blue']))
        self.linetypes=['solid','dotdashed','dashed','dotted']
        self.point_types=['filled circle','filled diamond','filled square','filled triangle']



    def doplots(self):
        import biggles
        import converter
        import pcolors
        
        scale=.01
        #biggles.configure("screen","width",2000)
        #biggles.configure("screen","height",1100)
        biggles.configure("default","fontsize_min",1.2)
        #biggles.configure('_HalfAxis','ticklabels_style',{'fontsize':2.0})
        biggles.configure('_HalfAxis','ticks_size',3.)
        biggles.configure('_HalfAxis','subticks_size',3./2)
        biggles.configure('PlotKey','key_vsep',3)
        biggles.configure('PlotKey','key_width',15)
        biggles.configure('_ErrorBar','barsize',1)
        nrun=len(self.runs)
        
        td = self.plotters[0].read_data()
        n_s2 = self.n_s2
        n_s2n = self.n_s2n

        if n_s2n == 12:
            nrow=4
            ncol=3
        else:
            nrow=5
            ncol=3
        arr = biggles.FramedArray(nrow,ncol)
        # colors of each s2 bin
        colors = self.colors
        linetypes=self.linetypes

        s2n_name='s2n_admom'
        #s2n_name='s2n_meas'
        tag1='shear1'
        tag2='shear2'
        errtag1='shear1err'
        errtag2='shear2err'



        #plt=biggles.FramedPlot()

        dt=[('s2n','f8',nrun),
            ('s2n_meas','f8',nrun),
            ('g1true','f8',nrun),
            ('g2true','f8',nrun),
            ('g1meas','f8',nrun),
            ('g1err','f8',nrun),
            ('g2meas','f8',nrun),
            ('g2err','f8',nrun)]
        
        for i in xrange(nrow*ncol):
            irow = i / ncol
            icol = i % ncol
            arr[irow,icol].add(biggles.Curve([-1000]*2,[-1000]*2))
            #arr[irow,icol].yrange = [-0.0025,0.0025]
            arr[irow,icol].yrange = [-0.035,0.009]

        fdtype=[('s2','f8',n_s2n), ('s2n','f8',n_s2n),('s2n_meas','f8',n_s2n),
                ('c','f8',n_s2n), ('cerr','f8',n_s2n),
                ('m','f8',n_s2n),('merr','f8',n_s2n)]
        fits1 = zeros(n_s2, dtype=fdtype)
        fits2 = zeros(n_s2, dtype=fdtype)

        fcurves=[]
        nplot=0
        is2ns = list(reversed(xrange(n_s2n)))
        for i_s2n in is2ns:
            irow = nplot / ncol
            icol = nplot % ncol
            td = self.plotters[0].read_data()
            s2n = td[0][s2n_name][i_s2n]

            data = zeros(n_s2,dtype=dt)

            # runs correspond to different applied shear values
            for irun,plotter in enumerate(self.plotters):
                d = plotter.read_data()
                shear = plotter.get_shear_true()

                # redundant
                
                s2vals=[]
                for is2 in xrange(n_s2):

                    s2,ellip = shapesim.get_s2_e(plotter.simc, is2, 0)
                    s2vals.append(s2)

                    data['s2n'][is2,irun] = d[is2][s2n_name][i_s2n]
                    data['s2n_meas'][is2,irun] = d[is2]['s2n_meas'][i_s2n]

                    data['g1true'][is2,irun] = shear.g1
                    data['g2true'][is2,irun] = shear.g2
                    data['g1meas'][is2,irun] = d[is2][tag1][i_s2n]
                    data['g2meas'][is2,irun] = d[is2][tag2][i_s2n]
                    data['g1err'][is2,irun] = d[is2][errtag1][i_s2n]
                    data['g2err'][is2,irun] = d[is2][errtag2][i_s2n]


            for is2 in xrange(n_s2):
                s2=s2vals[is2]

                s2n_meas = data['s2n_meas'][is2,:]
                g1true = data['g1true'][is2,:]
                g2true = data['g2true'][is2,:]
                g1err  = data['g1err'][is2,:]
                g2err  = data['g2err'][is2,:]
                diff1  = data['g1meas'][is2,:] - g1true
                diff2  = data['g2meas'][is2,:] - g2true

                linfit1 = fitting.LineFitter(g1true, diff1, g1err)
                linfit2 = fitting.LineFitter(g1true, diff2, g2err)

                fits1['s2'][is2,i_s2n] = s2
                fits1['s2n'][is2,i_s2n] = s2n
                fits1['s2n_meas'][is2,i_s2n] = s2n_meas.mean()

                fits1['m'][is2,i_s2n] = linfit1.pars[0]
                fits1['merr'][is2,i_s2n] = linfit1.perr[0]
                fits1['c'][is2,i_s2n] = linfit1.pars[1]
                fits1['cerr'][is2,i_s2n] = linfit1.perr[1]

                fits2['s2'][is2,i_s2n] = s2
                fits2['s2n'][is2,i_s2n] = s2n
                fits2['s2n_meas'][is2,i_s2n] = s2n_meas.mean()
                fits2['m'][is2,i_s2n] = linfit2.pars[0]
                fits2['merr'][is2,i_s2n] = linfit2.perr[0]
                fits2['c'][is2,i_s2n] = linfit2.pars[1]
                fits2['cerr'][is2,i_s2n] = linfit2.perr[1]


                p1 = biggles.Points(g1true/scale,diff1,
                                    type=self.point_types[is2], 
                                    color=colors[is2],
                                    size=4)
                #c1 = biggles.Curve(g1true/scale,diff1,color=colors[is2],
                #                   type=linetypes[is2])
                lwidth=5.
                c1 = biggles.Curve(g1true/scale,linfit1(g1true/scale),
                                   color=colors[is2],
                                   type=linetypes[is2],
                                   width=lwidth)
                if i_s2n == (n_s2n-1):
                    cfake = biggles.Curve(g1true-1000,diff1,color=colors[is2],
                                          type=linetypes[is2],
                                          width=lwidth)
                    if self['use_rb']:
                        label = '%.2f' % sqrt(1/s2)
                    else:
                        label = '%.2f' % s2
                    cfake.label = label
                    fcurves.append(cfake)
                arr[irow,icol].add(p1,c1) 

                perr1 = biggles.SymmetricErrorBarsY(g1true/scale,
                                                    diff1,g1err,
                                                    color=colors[is2])
                arr[irow,icol].add(perr1)
                if is2 == 0:
                    z1=biggles.Curve([-50,50], [0,0])
                    arr[irow,icol].add(z1)



                if is2 == 0 and nplot==0:
                    print '%15s %15s %15s +/- %15s %15s +/- %15s' % ('s2','s2n','m','err','c','err')
                print '%15s %15s %15g +/- %15g %15g +/- %15g' % (s2,s2n,
                                                                 linfit1.pars[0]+1,linfit1.perr[0],
                                                                 linfit1.pars[1],linfit1.perr[1])

            s2nlab = biggles.PlotLabel(0.1,0.9,'S/N: %d' % s2n,
                                     fontsize=2.5,halign='left')
            arr[irow,icol].add(s2nlab)
            if irow == 3:
                arr[irow,icol].yrange = [-0.007,0.01]
            elif irow == 2:
                arr[irow,icol].yrange = [-0.0024,0.0039]
            elif irow == 1:
                #arr[irow,icol].yrange = [-0.0029,0.0029]
                arr[irow,icol].yrange = [-0.0009,0.0029]
            elif irow == 0:
                #arr[irow,icol].yrange = [-0.0019,0.0019]
                arr[irow,icol].yrange = [-0.0009,0.0019]

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
        if self.n_s2n == 12:
            keyrow=nrow-1
            keycol=0
            pos=(.80,.30)
        else:
            keyrow=nrow-1
            keycol=ncol-1
            pos=(.85,.85)

        key=biggles.PlotKey(pos[0],pos[1],fcurves,halign='right',fontsize=fsize)
        arr[keyrow,keycol].add(key, *fcurves)

        if self['use_rb']:
            klabtext=r'$\sigma_{gal}/\sigma_{psf}$: '
        else:
            klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$: '
        klab = biggles.PlotLabel(0.5,pos[1],klabtext,
                                 fontsize=fsize,halign='right')
        arr[keyrow,keycol].add(klab)

        simc = self.plotters[0].simc
        objmodel = simc['objmodel']
        psfmodel = simc['psfmodel']
        plab='%s %s' % (objmodel,psfmodel)

        lowest=0.15
        l = biggles.PlotLabel(0.075,lowest+0, plab, halign='left')
        arr[nrow-1,ncol-1].add(l)

        """
        if simc['psfmodel'] == 'turb':
            siglab = r'$FWHM: %.1f$ pix' % simc['psf_fwhm']
        else:
            psf_sigma = simc['psf_sigma']
            siglab = r'$\sigma: %.1f$ pix' % psf_sigma
        sl = biggles.PlotLabel(0.075,lowest+0.45, siglab, halign='left', 
                               fontsize=2.5)
        arr[nrow-1,ncol-1].add(sl)
        """
        elab = r'$e_{gal}^{tot}: %.2f$' % td[0]['etrue'].mean()

        psf_estring=self.plotters[0].psf_estring
        if psf_estring:
            psfel = biggles.PlotLabel(0.075,lowest+0.3, psf_estring, 
                                      halign='left', 
                                      fontsize=2.5)
            arr[nrow-1,ncol-1].add(psfel)

        el = biggles.PlotLabel(0.075,lowest+0.15, elab, halign='left', 
                               fontsize=2.5)
        arr[nrow-1,ncol-1].add(el)



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
        if self['show']:
            arr.show()

        extra=''
        if self['use_rb']:
            extra='-rb'
        epsfile = shapesim.get_plot_file(self.set,'vs-shear'+extra,yrng=yrng)
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)

        self.do_plot_fits(fits1, fits2)

    def do_plot_fits(self, fits1, fits2):
        """
        gamma2 plots look exactly the same, now only to gamma1
        """
        import biggles
        import converter
        mplt = biggles.FramedPlot()
        cplt = biggles.FramedPlot()
        biggles.configure("default","fontsize_min",1.5)
        biggles.configure('PlotKey','key_vsep',1.5)
        biggles.configure('PlotKey','key_width',12)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)

        #mupper = biggles.Curve([-500,500],[0.004,0.004],color='grey')
        #mlower = biggles.Curve([-500,500],[-0.004,-0.004],color='grey')
        #mplt.add(mupper)
        #mplt.add(mlower)
        mplt.add( biggles.FillBetween([-500,500], [0.004,0.004], 
                                      [-500,500], [-0.004,-0.004],
                                      color=self.fill_color))
        cplt.add( biggles.FillBetween([-500,500], [0.0004,0.0004], 
                                      [-500,500], [-0.0004,-0.0004],
                                      color=self.fill_color))

        n_s2 = self.n_s2
        n_s2n = self.n_s2n
        colors = self.colors
        linetypes=self.linetypes

        kplots=[]
        for is2 in xrange(n_s2):

            s2,ellip = shapesim.get_s2_e(self.plotters[is2].simc, is2, 0)

            mpts1=biggles.Points(fits1['s2n'][is2,:], fits1['m'][is2,:],
                                 type=self.point_types[is2], 
                                 color=colors[is2])
            cpts1=biggles.Points(fits1['s2n'][is2,:], fits1['c'][is2,:],
                                 type=self.point_types[is2], 
                                 color=colors[is2])
            width=3.
            mcur1=biggles.Curve(fits1['s2n'][is2,:], fits1['m'][is2,:],
                                color=colors[is2],type=linetypes[is2],
                                width=width)
            ccur1=biggles.Curve(fits1['s2n'][is2,:], fits1['c'][is2,:],
                                color=colors[is2],type=linetypes[is2],
                                width=width)

            merr1=biggles.SymmetricErrorBarsY(fits1['s2n'][is2,:], 
                                              fits1['m'][is2,:],
                                              fits1['merr'][is2,:],
                                              color=colors[is2])
            cerr1=biggles.SymmetricErrorBarsY(fits1['s2n'][is2,:], 
                                              fits1['c'][is2,:],
                                              fits1['cerr'][is2,:],
                                              color=colors[is2])
            mpts2=biggles.Points(fits2['s2n'][is2,:], fits2['m'][is2,:],
                                 type=self.point_types[is2], 
                                 color=colors[is2])
            cpts2=biggles.Points(fits2['s2n'][is2,:], fits2['c'][is2,:],
                                 type=self.point_types[is2], 
                                 color=colors[is2])
            mcur2=biggles.Curve(fits2['s2n'][is2,:], fits2['m'][is2,:],
                                color=colors[is2],type=linetypes[is2],
                                width=width)
            ccur2=biggles.Curve(fits2['s2n'][is2,:], fits2['c'][is2,:],
                                color=colors[is2],type=linetypes[is2],
                                width=width)



            merr2=biggles.SymmetricErrorBarsY(fits2['s2n'][is2,:], 
                                              fits2['m'][is2,:],
                                              fits2['merr'][is2,:],
                                              color=colors[is2])
            cerr2=biggles.SymmetricErrorBarsY(fits2['s2n'][is2,:], 
                                              fits2['c'][is2,:],
                                              fits2['cerr'][is2,:],
                                              color=colors[is2])

            
            if self['use_rb']:
                label = '%.2f' % sqrt(1/s2)
            else:
                label = '%.2f' % s2
            ccur1.label = label
            ccur2.label = label
            kplots.append(ccur1)

            mplt.add(mpts1,merr1,mcur1)
            cplt.add(cpts1,cerr1,ccur1)

        """
        lab1 = biggles.PlotLabel(0.9,0.1,r'$\gamma_1$',halign='right')
        lab2 = biggles.PlotLabel(0.9,0.1,r'$\gamma_2$',halign='right')
        marr[0,0].add(lab1)
        marr[0,1].add(lab2)
        carr[0,0].add(lab1)
        carr[0,1].add(lab2)
        """

        xlabel = r'$S/N$'
        mplt.xlabel = xlabel
        mplt.ylabel = 'Calibration error (m-1)'
        mplt.aspect_ratio=1 #/GRATIO
        cplt.xlabel = xlabel
        cplt.ylabel = 'Additive error c'
        cplt.aspect_ratio=1 #/GRATIO

        if self['use_rb']:
            klabtext=r'$\sigma_{gal}/\sigma_{psf}$: '
        else:
            klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$: '
        fsize=3
        labypos=.9
        klab = biggles.PlotLabel(.65,labypos,klabtext,
                                 fontsize=fsize,halign='right')
        key=biggles.PlotKey(.85,labypos,kplots,halign='right',fontsize=fsize)
        cplt.add(key)
        mplt.add(key)
        cplt.add(klab)
        mplt.add(klab)

        cplt.add(biggles.Curve([-50,500],[0,0]))
        mplt.add(biggles.Curve([-50,500],[0,0]))

        mplt.yrange=[-0.02,.1]

        xrng=[-1,110]
        cplt.xrange=xrng
        cplt.yrange=[-.0029,.0029]
        mplt.xrange=xrng

        if self['show']:
            mplt.show()
            cplt.show()

        extra=''
        if self['use_rb']:
            extra='-rb'

        c_epsfile = shapesim.get_plot_file(self.set,'c-vs-shear'+extra)
        m_epsfile = shapesim.get_plot_file(self.set,'m-vs-shear'+extra)

        cplt.write_eps(c_epsfile)
        mplt.write_eps(m_epsfile)
        converter.convert(c_epsfile,dpi=100,verbose=True)
        converter.convert(m_epsfile,dpi=100,verbose=True)

        bias_file=shapesim.get_bias_file(self.set, 'bias')
        wlog("writing fits to file:",bias_file)
        with fitsio.FITS(bias_file,mode='rw',clobber=True) as fobj:
            fobj.write(fits1)
            fobj.write(fits2)

class MultiPlotterVsEpsf(MultiPlotterBase):
    """
    Additive shear bias as a function of PSF ellipticity
    in bins of S/N
    """
    def __init__(self, set, **keys):
        super(MultiPlotterVsEpsf,self).__init__(set, **keys)
        self.do_setup()

        self['use_rb'] = keys.get('use_rb',True)

    def do_setup(self):
        import pcolors
        td = self.plotters[0].read_data()
        self.n_s2 = len(td)
        self.n_s2n = td[0].size
        if self.n_s2 != 4:
            raise ValueError("adapt for n_s2 != 4")
        #self.colors = ['blue','magenta','green','red']
        self.colors=list(reversed(['red','forestgreen','NavajoWhite3','blue']))
        self.linetypes=['solid','dotdashed','dashed','dotted']
        self.point_types=['filled circle','filled diamond','filled square','filled triangle']



    def doplots(self):
        import biggles
        import converter
        import pcolors
        
        scale=.01
        #biggles.configure("screen","width",2000)
        #biggles.configure("screen","height",1100)
        biggles.configure("default","fontsize_min",1.5)
        #biggles.configure('_HalfAxis','ticklabels_style',{'fontsize':2.0})
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        biggles.configure('PlotKey','key_vsep',3)
        biggles.configure('PlotKey','key_width',15)
        biggles.configure('_ErrorBar','barsize',1)
        nrun=len(self.runs)
        
        td = self.plotters[0].read_data()
        n_s2 = self.n_s2
        n_s2n = self.n_s2n

        nrow=2
        ncol=3
        arr = biggles.FramedArray(nrow,ncol)
        # colors of each s2 bin
        colors = self.colors
        linetypes=self.linetypes

        s2n_name='s2n_admom'
        tag1='shear1'
        tag2='shear2'
        errtag1='shear1err'
        errtag2='shear2err'



        #plt=biggles.FramedPlot()

        dt=[('psf_e1','f8',nrun),
            ('psf_e2','f8',nrun),
            ('g1true','f8',nrun),
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

            arr[irow,icol].add( biggles.FillBetween([-500,500], [0.0004,0.0004], 
                                                    [-500,500], [-0.0004,-0.0004], 
                                                    color=self.fill_color))

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

                psf_e1=plotter.simc['psf_e1']
                psf_e2=plotter.simc['psf_e2']
                psfe = lensing.shear.Shear(e1=psf_e1,e2=psf_e2)

                # redundant
                
                s2vals=[]
                for is2 in xrange(n_s2):

                    s2,ellip = shapesim.get_s2_e(plotter.simc, is2, 0)
                    s2vals.append(s2)

                    data['psf_e1'][is2,irun] = psfe.e1
                    data['psf_e2'][is2,irun] = psfe.e2

                    data['g1true'][is2,irun] = shear.g1
                    data['g2true'][is2,irun] = shear.g2
                    data['g1meas'][is2,irun] = d[is2][tag1][i_s2n]
                    data['g2meas'][is2,irun] = d[is2][tag2][i_s2n]
                    data['g1err'][is2,irun] = d[is2][errtag1][i_s2n]
                    data['g2err'][is2,irun] = d[is2][errtag2][i_s2n]


            for is2 in xrange(n_s2):
                s2=s2vals[is2]

                psf_e1 = data['psf_e1'][is2,:]
                psf_e2 = data['psf_e2'][is2,:]


                diff1  = data['g1meas'][is2,:] - data['g1true'][is2,:]
                diff2  = data['g2meas'][is2,:] - data['g2true'][is2,:]
                p1 = biggles.Points(psf_e1/scale,diff1,
                                    type=self.point_types[is2], 
                                    color=colors[is2],
                                    size=3)
                width=5
                c1 = biggles.Curve(psf_e1/scale,diff1,color=colors[is2],
                                   type=linetypes[is2],
                                   width=width)
                if i_s2n == (n_s2n-1):
                    cfake = biggles.Curve(psf_e1-1000,diff1,color=colors[is2],
                                          type=linetypes[is2],
                                          width=width)
                    if self['use_rb']:
                        label = '%.2f' % sqrt(1/s2)
                    else:
                        label = '%.2f' % s2
                    cfake.label = label
                    fcurves.append(cfake)
                arr[irow,icol].add(p1,c1) 

                g1err  = data['g1err'][is2,:]
                g2err  = data['g2err'][is2,:]
                perr1 = biggles.SymmetricErrorBarsY(psf_e1/scale,
                                                    diff1,g1err,
                                                    color=colors[is2])
                arr[irow,icol].add(perr1)
                if is2 == 0:
                    z1=biggles.Curve([-50,50], [0,0])
                    arr[irow,icol].add(z1)


                """
                linfit1 = fitting.LineFitter(psf_e1, diff1, g1err)
                linfit2 = fitting.LineFitter(psf_e2, diff2, g2err)
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
                """

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
        xmin = 0.75
        key=biggles.PlotKey(xmin,.3,fcurves,halign='right',fontsize=fsize)
        arr[0,0].add(key, *fcurves)

        if self['use_rb']:
            klabtext=r'$\sigma_{gal}/\sigma_{psf}$: '
        else:
            klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}$: '
        klab = biggles.PlotLabel(xmin-0.3,.3,klabtext,
                                 fontsize=fsize,halign='right')
        arr[0,0].add(klab)

        simc = self.plotters[0].simc
        objmodel = simc['objmodel']
        psfmodel = simc['psfmodel']
        plab='%s %s' % (objmodel,psfmodel)

        lowest=0.8
        step=0.10
        l = biggles.PlotLabel(0.075,lowest+0, plab, halign='left')
        arr[0,0].add(l)

        elab = r'$e_{gal}^{tot}: %.2f$' % td[0]['etrue'].mean()

        el = biggles.PlotLabel(0.075,lowest+step, elab, halign='left', 
                               fontsize=2.5)
        arr[0,0].add(el)



        yrng=self.get('yrange',None)
        if yrng:
            arr.yrange = yrng
        arr.xrange=array([-0.005,0.069])/scale
        #arr.uniform_limits=1
        arr.xlabel = r'$e_{PSF}/%.2g$' % scale
        arr.ylabel = r'Additive error c'
        #arr.aspect_ratio=1/1.61803399
        #arr.aspect_ratio=1/1.4
        arr.aspect_ratio=1/1.45

        title=self.get_title()
        if title:
            arr.title=title
        if self['show']:
            arr.show()

        extra=''
        if self['use_rb']:
            extra='-rb'
        epsfile = shapesim.get_plot_file(self.set,'vs-epsf'+extra,yrng=yrng)
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)



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

        self.fill_color='grey80'

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
                           s2n_name=None,
                           xrng=None,
                           yrng=None, 
                           title=None,
                           show=True):
        import biggles
        import pcolors
        import converter

        biggles.configure("default","fontsize_min",2)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        biggles.configure('PlotKey','key_width',13)
        biggles.configure('PlotKey','key_vsep',1.0)

        if s2n_name is None:
            s2n_name='s2n_admom'

        extra=''

        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        is_shearmag=False
        if 'shearmag' in self.simc:
            is_shearmag=True
            shearmag = self.simc['shearmag']
        else:
            shear_true = self.get_shear_true()

        data = self.read_data()

        epsfile = shapesim.get_plot_file(self['run'],type+extra,yrng=yrng)
        wlog("will plot to:",epsfile)

        if len(data) == 4:
            colors=['red','forestgreen','NavajoWhite3','blue']
            linetypes=['dotted','dashed','dotdashed','solid']
        else:
            colors=pcolors.rainbow(len(data), 'hex')
            linetypes=['solid']*len(data)
        point_types=['filled circle','filled diamond','filled square','filled triangle']

        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1
        
        title=self.get_title(title=title)
        if title:
            arr.title=title
 

        if self['run'][0:5] == 'deswl':
            tag1='gamma1_meas'
            tag2='gamma2_meas'
            # convention
            st[tag1] = -st[tag1]
        else:
            if is_shearmag:
                tag1='sheardiff'
                tag2='osheardiff'
                errtag1='sheardifferr'
                errtag2='osheardifferr'
                lab1=r'$\gamma_+$'
                lab2=r'$\gamma_\times$'
            else:
                tag1='shear1'
                tag2='shear2'
                errtag1='shear1err'
                errtag2='shear2err'
                lab1=r'$\gamma_1$'
                lab2=r'$\gamma_2$'


        plots1=[]
        plots2=[]
        allplots=[]


        # looping over s2
        for i,st in enumerate(reversed(data)):
            #wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            s2 = median(st['s2'])

            #s2n_name='s2n_admom'
            xlabel = 'S/N'

            s2n = st[s2n_name]

          
            if type == 'diff':
                if is_shearmag:
                    yvals1 = st[tag1]
                    yvals2 = st[tag2]
                else:
                    yvals1 = st[tag1] - shear_true.g1
                    yvals2 = st[tag2] - shear_true.g2
            elif type == 'val':
                yvals1 = st[tag1]
                yvals2 = st[tag2]

            else:
                raise ValueError("bad plot type: '%s'" % type)

            label = r'%0.3f' % s2

            cwidth=5.
            pr1 = biggles.Points(s2n, yvals1, color=colors[i],
                                 size=2,
                                 type=point_types[i])
            pr2 = biggles.Points(s2n, yvals2, color=colors[i],
                                 size=2,
                                 type=point_types[i])
            cr1 = biggles.Curve(s2n, yvals1, color=colors[i],type=linetypes[i],width=cwidth)
            cr2 = biggles.Curve(s2n, yvals2, color=colors[i],type=linetypes[i],width=cwidth)
            cr1.label = label
            cr2.label = label


            arr[0,0].add(cr1,pr1)
            arr[1,0].add(cr2,pr2)
            
            if True:
                g1err = st[errtag1]
                g2err = st[errtag2]
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


        fsize=2
        key1 = biggles.PlotKey(0.85,0.92, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.85,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}: $'
        klab = biggles.PlotLabel(.65,.92,klabtext,
                                 fontsize=2.5,halign='right')
        arr[0,0].add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.9,0.1, plab, halign='right')
        arr[1,0].add(l)

        if 'Tpsf' not in self.simc:
            if self.simc['psfmodel'] == 'turb':
                siglab = r'$FWHM: %.1f$ pix' % self.simc['psf_fwhm']
            else:
                psf_sigma = self.simc['psf_sigma']
                siglab = r'$\sigma: %.1f$ pix' % psf_sigma
            siglab += ' '+self.psf_estring
            if 'etrue' in st.dtype.names:
                siglab += r'$ e_{gal}^{tot}: %.2f$' % st['etrue'].mean()

                sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                                       fontsize=2.5)
                arr[1,0].add(sl)

        arr.xlabel=xlabel


        arr.ylabel = r'$\Delta \gamma$'

        if is_shearmag:
            g1lab_txt = lab1 + ' = %.2g' % shearmag
            g2lab_txt = lab2
        else:
            g1lab_txt = lab1 + ' = %.2g' % shear_true.g1
            g2lab_txt = lab2 + ' = %.2g' % shear_true.g2

        g1lab = biggles.PlotLabel(0.1,0.9, g1lab_txt, halign='left')

        if 'dt' in self['run']:
            g2x = 0.7
        else:
            g2x = 0.1

        g2lab = biggles.PlotLabel(g2x,0.9, g2lab_txt, halign='left')

        arr[0,0].add(g1lab)
        arr[1,0].add(g2lab)

        expect1 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()], [0,0])
        expect2 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()], [0,0])

        arr[0,0].add(expect1)
        arr[1,0].add(expect2)



        """
        if xrng is None:
            xrng = [0,s2n.max()*1.4]
        arr.xrange = xrng
        """

        if yrng is not None:
            arr.yrange = yrng

        wlog("Writing plot file:",epsfile)
        if show:
            arr.show()
        arr.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)


    def plots_shear1_frac_vs_s2n(self, 
                                 xrng=None,
                                 yrng=None, 
                                 title=None,
                                 show=True):
        """
        special plotting just shear1 as a fraction of
        true
        """
        import biggles
        import pcolors
        import converter

        type='frac'

        biggles.configure("default","fontsize_min",2)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        biggles.configure('PlotKey','key_width',13)
        biggles.configure('PlotKey','key_vsep',1.0)

        extra=''

        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        is_shearmag=False
        if 'shearmag' in self.simc:
            is_shearmag=True
            shearmag = self.simc['shearmag']
        else:
            shear_true = self.get_shear_true()

        data = self.read_data()

        epsfile = shapesim.get_plot_file(self['run'],type+extra,yrng=yrng)
        wlog("will plot to:",epsfile)

        if len(data) == 4:
            colors=['red','forestgreen','NavajoWhite3','blue']
            linetypes=['dotted','dashed','dotdashed','solid']
        else:
            colors=pcolors.rainbow(len(data), 'hex')
            linetypes=['solid']*len(data)
        point_types=['filled circle','filled diamond','filled square','filled triangle']

        plt=biggles.FramedPlot()
        plt.add( biggles.FillBetween([-500,500], [0.004,0.004], 
                                     [-500,500], [-0.004,-0.004],
                                     color=self.fill_color))


        title=self.get_title(title=title)
        if title:
            plt.title=title
 

        tag1='shear1'
        tag2='shear2'
        errtag1='shear1err'
        errtag2='shear2err'
        lab1=r'$\gamma$'


        plots1=[]
        allplots=[]


        # looping over s2
        for i,st in enumerate(reversed(data)):
            #wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            s2 = median(st['s2'])

            s2n_name='s2n_admom'
            xlabel = 'S/N'

            s2n = st[s2n_name]

          
            yvals1 = (st[tag1] - shear_true.g1)/shear_true.g1

            label = r'%0.3f' % s2

            cwidth=5.
            pr1 = biggles.Points(s2n, yvals1, color=colors[i],
                                 size=2,
                                 type=point_types[i])
            cr1 = biggles.Curve(s2n, yvals1, color=colors[i],type=linetypes[i],width=cwidth)
            cr1.label = label

            plt.add(cr1,pr1)
            
            g1err = st[errtag1]/shear_true.g1
            err1p = biggles.SymmetricErrorBarsY(s2n, yvals1, g1err,
                                                color=colors[i])
            plt.add(err1p)

            plots1.append(cr1)


        fsize=2
        key1 = biggles.PlotKey(0.85,0.92, plots1, halign='right', 
                               fontsize=fsize)
        plt.add(key1)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}: $'
        klab = biggles.PlotLabel(.62,.92,klabtext,
                                 fontsize=2.5,halign='right')
        plt.add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.9,0.1, plab, halign='right')
        plt.add(l)

        plt.xlabel=xlabel


        plt.ylabel = r'$\Delta \gamma/\gamma$'

        g1lab_txt = lab1 + ' = %.2g' % shear_true.g1

        g1lab = biggles.PlotLabel(0.1,0.9, g1lab_txt, halign='left')

        plt.add(g1lab)

        erange=[0.2*s2n.min(),1.05*s2n.max()]
        expect1 = biggles.Curve([0.2*s2n.min(),1.05*s2n.max()], [0,0])

        plt.add(expect1)
        plt.xrange=[0.,erange[1]]

        if yrng is not None:
            plt.yrange = yrng

        wlog("Writing plot file:",epsfile)
        if show:
            plt.show()
        plt.write_eps(epsfile)
        converter.convert(epsfile,dpi=100,verbose=True)


    def plots_shear_vs_err(self, 
                           xrng=None,
                           yrng=None, 
                           title=None,
                           show=True):
        import biggles
        import pcolors
        import converter

        biggles.configure("default","fontsize_min",2)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',2.5/2)
        biggles.configure('PlotKey','key_width',13)
        biggles.configure('PlotKey','key_vsep',1.0)

        extra=''

        if self.docum:
            extra+='-cum'

        runtype = self.get('runtype','byellip')
        if runtype != 'bys2n':
            raise ValueError("Can only make plots vs s2n for 'bys2n' runs")

        is_shearmag=False
        if 'shearmag' in self.simc:
            is_shearmag=True
            shearmag = self.simc['shearmag']
        else:
            shear_true = self.get_shear_true()

        data = self.read_data()

        epsfile = shapesim.get_plot_file(self['run'],'ediff'+extra,yrng=yrng)
        wlog("will plot to:",epsfile)

        if len(data) == 4:
            colors=['red','forestgreen','NavajoWhite3','blue']
            linetypes=['dotted','dashed','dotdashed','solid']
        else:
            colors=pcolors.rainbow(len(data), 'hex')
            linetypes=['solid']*len(data)
        point_types=['filled circle','filled diamond','filled square','filled triangle']

        arr=biggles.FramedArray(2,1)
        #arr.aspect_ratio=1
        
        title=self.get_title(title=title)
        if title:
            arr.title=title
 

        if self['run'][0:5] == 'deswl':
            tag1='gamma1_meas'
            tag2='gamma2_meas'
            # convention
            st[tag1] = -st[tag1]
        else:
            if is_shearmag:
                tag1='sheardiff'
                tag2='osheardiff'
                errtag1='sheardifferr'
                errtag2='osheardifferr'
                lab1=r'$\gamma_+$'
                lab2=r'$\gamma_\times$'
            else:
                tag1='shear1'
                tag2='shear2'
                errtag1='shear1err'
                errtag2='shear2err'
                lab1=r'$\gamma_1$'
                lab2=r'$\gamma_2$'


        plots1=[]
        plots2=[]
        allplots=[]


        err_range=[1.e9,-1.e9]
        # looping over s2
        for i,st in enumerate(reversed(data)):
            #wlog("s2:",median(st['s2']),"s2_meas:",median(st['s2_meas']))

            s2 = median(st['s2'])

            xlabel = r'$\sigma(\gamma) per galaxy$'

            #g1sens_mean=st['g1sensum']/st['nsum']
            err = st[errtag1]*sqrt(st['nsum'])#*g1sens_mean

            minerr=err.min()
            maxerr=err.max()
            if minerr < err_range[0]:
                err_range[0]=minerr
            if maxerr > err_range[1]:
                err_range[1]=maxerr

          
            yvals1 = st[tag1] - shear_true.g1
            yvals2 = st[tag2] - shear_true.g2

            label = r'%0.3f' % s2

            cwidth=5.
            pr1 = biggles.Points(err, yvals1, color=colors[i],
                                 size=2,
                                 type=point_types[i])
            pr2 = biggles.Points(err, yvals2, color=colors[i],
                                 size=2,
                                 type=point_types[i])
            cr1 = biggles.Curve(err, yvals1, color=colors[i],type=linetypes[i],width=cwidth)
            cr2 = biggles.Curve(err, yvals2, color=colors[i],type=linetypes[i],width=cwidth)
            cr1.label = label
            cr2.label = label


            arr[0,0].add(cr1,pr1)
            arr[1,0].add(cr2,pr2)
            
            if True:
                g1err = st[errtag1]
                g2err = st[errtag2]
                err1p = biggles.SymmetricErrorBarsY(err, yvals1, g1err,
                                                    color=colors[i])
                err2p = biggles.SymmetricErrorBarsY(err, yvals2, g2err,
                                                    color=colors[i])
                if not self.noerr:
                    arr[0,0].add(err1p)
                    arr[1,0].add(err2p)

            if i < 15:
                plots1.append(cr1)
            else:
                plots2.append(cr1)


        fsize=2
        key1 = biggles.PlotKey(0.85,0.92, plots1, halign='right', 
                               fontsize=fsize)
        arr[0,0].add(key1)
        if len(plots2) > 0:
            key2 = biggles.PlotKey(0.85,0.92, plots2, halign='right', 
                                   fontsize=fsize)
            arr[1,0].add(key2)

        klabtext=r'$\sigma^2_{psf}/\sigma^2_{gal}: $'
        klab = biggles.PlotLabel(.65,.92,klabtext,
                                 fontsize=2.5,halign='right')
        arr[0,0].add(klab)
        objmodel = self.simc['objmodel']
        psfmodel = self.simc['psfmodel']


        plab='%s %s' % (objmodel,psfmodel)
        l = biggles.PlotLabel(0.9,0.1, plab, halign='right')
        arr[1,0].add(l)

        if 'Tpsf' not in self.simc:
            if self.simc['psfmodel'] == 'turb':
                siglab = r'$FWHM: %.1f$ pix' % self.simc['psf_fwhm']
            else:
                psf_sigma = self.simc['psf_sigma']
                siglab = r'$\sigma: %.1f$ pix' % psf_sigma
            siglab += ' '+self.psf_estring
            if 'etrue' in st.dtype.names:
                siglab += r'$ e_{gal}^{tot}: %.2f$' % st['etrue'].mean()

                sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                                       fontsize=2.5)
                arr[1,0].add(sl)

        arr.xlabel=xlabel


        arr.ylabel = r'$\Delta \gamma$'

        if is_shearmag:
            g1lab_txt = lab1 + ' = %.2g' % shearmag
            g2lab_txt = lab2
        else:
            g1lab_txt = lab1 + ' = %.2g' % shear_true.g1
            g2lab_txt = lab2 + ' = %.2g' % shear_true.g2

        g1lab = biggles.PlotLabel(0.1,0.9, g1lab_txt, halign='left')

        if 'dt' in self['run']:
            g2x = 0.7
        else:
            g2x = 0.1

        g2lab = biggles.PlotLabel(g2x,0.9, g2lab_txt, halign='left')

        arr[0,0].add(g1lab)
        arr[1,0].add(g2lab)

        expect1 = biggles.Curve([0.5*err_range[0],1.5*err_range[1]], [0,0])
        expect2 = biggles.Curve([0.5*err_range[0],1.5*err_range[1]], [0,0])

        arr[0,0].add(expect1)
        arr[1,0].add(expect2)


        if yrng is not None:
            arr.yrange = yrng

        arr.xrange=[.01,1.]
        arr.xlog=True

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

        if 'Tpsf' not in self.simc:
            if self.simc['psfmodel'] == 'turb':
                siglab=r'$FWHM_{PSF}: %.1f$ pix' % self.simc['psf_fwhm']
            else:
                psf_sigma = self.simc['psf_sigma']
                siglab=r'$\sigma_{PSF}: %.1f$ pix' % psf_sigma
            siglab += ' '+self.psf_estring

            sl = biggles.PlotLabel(0.075,0.1, siglab, halign='left', 
                                   fontsize=2.5)

        s2n=self['s2n']
        if s2n > 0:
            if s2n > 1000:
                ls2n = numpy.log10(s2n)
                ls2n = r'$10^{%.1f}$' % ls2n
            else:
                ls2n = '%.0f' % s2n

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
            data = shapesim.read_averaged_outputs(self['run'], 
                                                  docum=self.docum, 
                                                  skip1=self.skip1) 
            self._data=data
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
