from numpy import zeros, sqrt, linspace
import cluster_step
from cluster_step import files, stats


from esutil.numpy_util import aprint

import biggles
from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY, \
        PlotLabel, PlotKey

from fitting import LineFitter

class BiasFitter(object):
    def __init__(self, run, **keys):

        self.run     = run
        psfnums = keys.get('psfnums',None)
        self.psfnums = files.get_psfnums(psfnums)
        self.objtype = keys.get('objtype',None)

        self.s2n_min=keys.get('s2n_min',20.0)

        self.s2_max=keys.get('s2_max',0.5)

        self.s2n_field=keys.get('s2n_field','s2n_w')

        # starting point for labels
        self.lab1_loc=[1.-0.075, 0.1]
        self.lab1_halign='right'
        self.lab1_yshift=0.075

        self.lab2_loc=[0.075,0.075]
        self.lab2_halign='left'
        self.lab2_yshift=+0.075

        self.read_data()
        self.set_averages()
        aprint(self.avg, fancy=True)

        self.do_fits()

    def read_data(self):

        self.datadict={}
        for shnum in files.SHNUMS:
            print 'reading shnum:',shnum,
            data=files.read_output_set(self.run, 
                                       self.psfnums, 
                                       shnum,
                                       objtype=self.objtype,
                                       s2_max=self.s2_max,
                                       s2n_field=self.s2n_field,
                                       s2n_min=self.s2n_min)
            print data.size
            if data.size==0:
                raise ValueError("no data read")
            self.datadict[shnum]=data
     
    def set_averages(self):
        nsh=len(files.SHNUMS)
        dt=[('g1true','f8'),
            ('g2true','f8'),
            ('g1','f8'),('g1_err','f8'),
            ('g2','f8'),('g2_err','f8'),
            ('g1sens','f8'),('g1sens_err','f8'),
            ('g2sens','f8'),('g2sens_err','f8')]
        out=zeros(nsh, dtype=dt)
        for i,shnum in enumerate(files.SHNUMS):
            data=self.datadict[shnum]
            shres=stats.get_mean_shear_shstruct(data)

            out['g1true'][shnum-1] = cluster_step.sh1exp[shnum]
            out['g2true'][shnum-1] = cluster_step.sh2exp[shnum]

            for n in shres:
                out[n][shnum-1] = shres[n]

        self.avg=out

    def do_fits(self):

        g1diff=self.avg['g1']-self.avg['g1true']
        self.g1fit=LineFitter(self.avg['g1true'],
                              g1diff,
                              self.avg['g1_err'])

        g2diff=self.avg['g2']-self.avg['g2true']
        self.g2fit=LineFitter(self.avg['g2true'],
                              g2diff,
                              self.avg['g2_err'])

        print self.g1fit
        print self.g2fit


    def make_plots(self):
        arr=FramedArray(2,1)
        arr.uniform_limits=1
        arr.aspect_ratio=1.618
        arr.xlabel=r'$\gamma_{true}$'
        arr.ylabel=r'$\Delta \gamma$'

        arr.yrange=[-0.01,0.01]

        # the data
        g1diff=self.avg['g1']-self.avg['g1true']
        g1_pts=Points(self.avg['g1true'], g1diff,
                      type='filled circle')
        g1err_pts=SymmetricErrorBarsY(self.avg['g1true'], g1diff, self.avg['g1_err'])

        g2diff=self.avg['g2']-self.avg['g2true']
        g2_pts=Points(self.avg['g2true'], g2diff, type='filled circle')
        g2err_pts=SymmetricErrorBarsY(self.avg['g2true'], g2diff, self.avg['g2_err'])

        # the fits
        fitdata1=self.g1fit(self.avg['g1true'])
        fitdata2=self.g2fit(self.avg['g2true'])

        fit1c = Curve(self.avg['g1true'], fitdata1, color='blue')
        fit2c = Curve(self.avg['g2true'], fitdata2, color='blue')

        # plot labels
        m1lab,c1lab,m2lab,c2lab=self.get_mc_labels()
        psflab=self.get_psf_label()
        typelab=self.get_objtype_label()
        sizelab=self.get_size_label()
        s2nlab=self.get_s2n_label()
        sh1lab,sh2lab=self.get_shear_labels()

        arr[0,0].add(fit1c, g1_pts, g1err_pts,m1lab,c1lab,sh1lab,
                     psflab,s2nlab,sizelab)
        arr[1,0].add(fit2c, g2_pts, g2err_pts,m2lab,c2lab,sh2lab)

        if typelab:
            arr[0,0].add(typelab)

        self.plt=arr

    def get_shear_labels(self):
        
        sh1labs=r'$\gamma_1$'
        sh2labs=r'$\gamma_2$'
        #yshift=2*self.lab1_yshift
        xloc=self.lab2_loc[0]
        #yloc=self.lab1_loc[1] + yshift
        yloc=1-0.075

        sh1lab=PlotLabel(xloc, yloc, sh1labs, halign=self.lab1_halign)
        sh2lab=PlotLabel(xloc, yloc, sh2labs, halign=self.lab1_halign)

        return sh1lab,sh2lab

    def get_psf_label(self):

        if len(self.psfnums)==len(files.PSFNUMS) is None:
            labs = 'all'
        else:
            psfnums=[str(p) for p in self.psfnums]
            labs=','.join(psfnums)

        labs='psf: %s' % labs
        lab=PlotLabel(self.lab2_loc[0], self.lab2_loc[1], labs, 
                      halign=self.lab2_halign)
        return lab

    def get_s2n_label(self):
        labs=r'$%s < %.2f$' % (self.s2n_field,self.s2n_min)
        yshift=1*self.lab2_yshift
        lab=PlotLabel(self.lab2_loc[0],self.lab2_loc[1]+yshift,
                      labs,
                      halign=self.lab2_halign)
        return lab


    def get_size_label(self):
        sratio=sqrt(1/self.s2_max)
        #labs=r'$\sigma^2_{psf}/\sigma^2_{gal} < %.2f$' % self.s2_max
        labs=r'$\sigma_{gal}/\sigma_{psf} > %.2f$' % sratio
        yshift=2*self.lab2_yshift
        lab=PlotLabel(self.lab2_loc[0],self.lab2_loc[1]+yshift,
                      labs,
                      halign=self.lab2_halign)
        return lab

    def get_objtype_label(self):
        """
        objtype
        s2
        s/n (include field name)
        """
        if self.objtype is None:
            return None
        labs='objtype: %s' % self.objtype
        yshift=3*self.lab2_yshift
        lab=PlotLabel(self.lab2_loc[0],self.lab2_loc[1]+yshift,
                      labs,
                      halign=self.lab2_halign)
        return lab 

    def get_mc_labels(self):
        halign='right'

        xstart=self.lab1_loc[0]
        ystart=self.lab1_loc[1]
        yshift=self.lab1_yshift

        #m1=self.g1fit.pars[0]-1
        m1=self.g1fit.pars[0]
        m1err=self.g1fit.perr[0]
        c1=self.g1fit.pars[1]
        c1err=self.g1fit.perr[1]

        #m2=self.g2fit.pars[0]-1
        m2=self.g2fit.pars[0]
        m2err=self.g2fit.perr[0]
        c2=self.g2fit.pars[1]
        c2err=self.g2fit.perr[1]

        m1labs=r'$m: %.4f \pm\ %.4f$' % (m1,m1err)
        c1labs=r'$c: %.4f \pm\ %.4f$' % (c1,c1err)

        m2labs=r'$m: %.4f \pm\ %.4f$' % (m2,m2err)
        c2labs=r'$c: %.4f \pm\ %.4f$' % (c2,c2err)

        m1lab=PlotLabel(xstart,ystart, m1labs, halign=halign)
        c1lab=PlotLabel(xstart,ystart+yshift, c1labs, halign=halign)

        m2lab=PlotLabel(xstart,ystart, m2labs, halign=halign)
        c2lab=PlotLabel(xstart,ystart+yshift, c2labs, halign=halign)

        return m1lab,c1lab,m2lab,c2lab

    def show(self):
        if not hasattr(self,'plt'):
            self.make_plots()

        self.plt.show()

