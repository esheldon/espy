"""
    %prog [options]

shearnums 1-8
psfnums 1-6
"""

import sys
import os
import numpy
import lensing
from numpy import zeros, sqrt, where

import cluster_step
from cluster_step import files, stats
from cluster_step import sh1exp, sh2exp

import esutil as eu
from esutil.numpy_util import aprint
import biggles
from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY, \
        PlotLabel

from optparse import OptionParser
parser=OptionParser(__doc__)

parser.add_option('-r','--run',default=None,
                  help='The run id, required')
parser.add_option('-s','--shnum',default=None,
                  help='The shear number, required')

parser.add_option('-p','--psfnums',default='1,2,3,4,5,6',
                  help='restrict to these PSFs, comma separated')

parser.add_option('-f','--field',default='s2n_w',
                  help="bin by this field, default s2n_w")

parser.add_option('-n','--nbin',default=20,
                  help="number of logarithmic bins, default %default")

parser.add_option('--s2n',default='10,200', help="s/n range, %default")
parser.add_option('--sratio',default='1.0,1.e6',
                  help='sratio range, %default')

parser.add_option('--show',action='store_true',
                  help="show the plot on the screen")


parser.add_option('--frac',action='store_true',
                  help=("show the fraction relative to the"
                        "expected truth"))

parser.add_option('-P','--progress',action='store_true',
                  help="show the progress bar")
parser.add_option('--ignore-missing',action='store_true',
                  help="ignore missing files")


parser.add_option('-y','--yrange',default=None,
                  help="y range for plot")

class ShearPlotter(object):
    def __init__(self):
        biggles.configure( 'default', 'fontsize_min', 2)
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if options.run is None or options.shnum is None:
            parser.print_help()
            sys.exit(1)

        self.yrange=options.yrange
        if self.yrange is not None:
            self.yrange=[float(s) for s in options.yrange.split(',')]

        self.run=options.run
        self.shnum=int(options.shnum)

        self.psfnums=[int(s) for s in options.psfnums.split(',')]

        self.nbin=int(options.nbin)

        self.s2n_range=[float(s) for s in options.s2n.split(',')]
        self.sratio_range=[float(s) for s in options.sratio.split(',')]

        self.doshow = options.show

        self.bin_field=options.field

        reader=files.Reader(run=self.run, 
                            shnums=self.shnum,
                            psfnums=self.psfnums,
                            s2n_range=self.s2n_range,
                            sratio_range=self.sratio_range,
                            setname=None,
                            ignore_missing=options.ignore_missing,
                            progress=options.progress)

        self.data=reader.get_data()
        
        self.set_bindata()
        if self.options.frac:
            self.make_frac_plot()
        else:
            self.make_plot()


    def set_bindata(self):
        from cluster_step.pipe import CLUSTERSTEP_PSF_STAR,CLUSTERSTEP_GAL

        p=files.read_prior_original(type='nosplit')
        print p.dtype
        wm,=where((p['mag'] > 20) & (p['mag'] < 21))
        se1=numpy.zeros(wm.size)
        se2=numpy.zeros(wm.size)
        for i in xrange(wm.size):
            se1[i],se2[i] = lensing.util.g1g2_to_e1e2(p['g'][wm[i],0],p['g'][wm[i],1])
        etot2=se1**2 + se2**2
        mesq = etot2.mean()
        R = 1.0 - 0.5*mesq
        print R

        data=self.data
        wstar,=where( (data['am_flags']==0) 
                        & ((data['uw_type'] & CLUSTERSTEP_PSF_STAR) != 0) )

        wgal,=where( (data['am_flags']==0) 
                        & ((data['uw_type'] & CLUSTERSTEP_GAL) != 0)
                        & (data['am_s2n'] > 100) & (data['am_s2n'] < 200))

        am_star_irr_mean=data['am_irr'][wstar].mean()
        am_star_irc_mean=data['am_irc'][wstar].mean()
        am_star_icc_mean=data['am_icc'][wstar].mean()
        am_star_T_mean = am_star_irr_mean+am_star_icc_mean

        am_gal_irr=data['am_irr'][wgal]
        am_gal_irc=data['am_irc'][wgal]
        am_gal_icc=data['am_icc'][wgal]
        am_gal_T = am_gal_irr+am_gal_icc

        wgal2,=where(am_gal_T/am_star_T_mean > 2)
        wgal=wgal[wgal2]


        star_isum=data['uw_isum'][wstar].sum()
        star_irr = data['uw_irrsum'][wstar].sum()/star_isum
        star_irc = data['uw_ircsum'][wstar].sum()/star_isum
        star_icc = data['uw_iccsum'][wstar].sum()/star_isum

        gal_isum=data['uw_isum'][wgal].sum()
        gal_irr = data['uw_irrsum'][wgal].sum()/gal_isum
        gal_irc = data['uw_ircsum'][wgal].sum()/gal_isum
        gal_icc = data['uw_iccsum'][wgal].sum()/gal_isum

        irr=gal_irr-star_irr
        irc=gal_irc-star_irc
        icc=gal_icc-star_icc
        T = irr+icc
        e1=(icc-irr)/T
        e2=2*irc/T

        print 0.5*e1/R, 0.5*e2/R
        stop

        """
        self.bindata=stats.logbin_shear_data(self.data, self.bin_field, 
                                             nbin=self.nbin, 
                                             min=self.s2n_range[0],
                                             max=self.s2n_range[1])
        aprint(self.bindata, header=True, page=False, fancy=True)
        """


    def get_labels(self):

        run_mess='run: %s' % self.run
        set_mess='set: %s' % self.setname
        nobj_mess='nobj: %s' % self.data.size

        sratio_range=self.sratio_range
        if sratio_range[1] > 1000:
            srat_mess = r'$\sigma_{gal}/\sigma_{psf} > %.2f$' % sratio_range[0]
        else:
            srat_mess=r'$%.2f < \sigma_{gal}/\sigma_{psf} < %.2f$' % tuple(sratio_range)

        smess = 's: %s' % self.shnum

        pstr=[str(s) for s in self.psfnums]
        pstr=','.join(pstr)
        pmess="p: %s" % pstr

        halign='left'
        x=0.05
        #y=0.9
        y=0.7
        #inc=-0.075
        inc=-0.05

        runlab=PlotLabel(x, y, run_mess, halign=halign)
        y += inc
        snumlab=PlotLabel(x, y, smess, halign=halign)
        y += inc
        pnumlab=PlotLabel(x, y, pmess, halign=halign)
        y += inc
        setlab=PlotLabel(x, y, set_mess, halign=halign)
        y += inc
        sratlab=PlotLabel(x, y, srat_mess, halign=halign)
        y += inc
        nobjlab=PlotLabel(x, y, nobj_mess, halign=halign)
        
        return (runlab,snumlab,pnumlab,setlab,sratlab,nobjlab)


    def show(self):
        self.plt.show()

    def make_plot(self):

        bindata=self.bindata
        bin_field=self.bin_field
        plt=FramedPlot()
        
        plt.uniform_limits=1
        plt.xlog=True
        #plt.xrange=[0.5*bindata[bin_field].min(), 1.5*bindata[bin_field].max()]
        plt.xrange=self.s2n_range
        plt.xlabel=bin_field
        plt.ylabel = r'$\gamma$'
        if self.yrange is not None:
            plt.yrange=self.yrange

        xdata=bindata[bin_field]
        xerr=bindata[bin_field+'_err']

        if self.shnum in sh1exp:
            g1exp=zeros(xdata.size)+sh1exp[self.shnum]
            g2exp=zeros(xdata.size)+sh2exp[self.shnum]
            g1exp_plt=Curve(xdata, g1exp)
            g2exp_plt=Curve(xdata, g2exp)
            plt.add(g1exp_plt)
            plt.add(g2exp_plt)


        xerrpts1 = SymmetricErrorBarsX(xdata, bindata['g1'], xerr)
        xerrpts2 = SymmetricErrorBarsX(xdata, bindata['g2'], xerr)

        type='filled circle'
        g1color='blue'
        g2color='red'
        g1pts = Points(xdata, bindata['g1'], type=type, color=g1color)
        g1errpts = SymmetricErrorBarsY(xdata, bindata['g1'], bindata['g1_err'], color=g1color)
        g2pts = Points(xdata, bindata['g2'], type=type, color=g2color)
        g2errpts = SymmetricErrorBarsY(xdata, bindata['g2'], bindata['g2_err'], color=g2color)

        g1pts.label=r'$\gamma_1$'
        g2pts.label=r'$\gamma_2$'

        key=biggles.PlotKey(0.9,0.5,[g1pts,g2pts],halign='right')

        plt.add( xerrpts1, g1pts, g1errpts )
        plt.add( xerrpts2, g2pts, g2errpts )
        plt.add(key)

        labels=self.get_labels()

        plt.add(*labels)
        plt.aspect_ratio=1

        self.plt=plt

    def make_frac_plot(self):
        if self.shnum not in sh1exp:
            raise ValueError("you must know the expected value")

        if sh1exp[self.shnum] != 0:
            gfield='g1'
            gtrue=sh1exp[self.shnum]
        elif sh2exp[self.shnum] != 0:
            gfield='g2'
            gtrue=sh2exp[self.shnum]
        else:
            raise ValueError("all expected are listed as zero")

        bindata=self.bindata
        bin_field=self.bin_field
        plt=FramedPlot()
        
        plt.title=self.get_title()

        plt.xlog=True
        plt.xrange=[0.5*bindata[bin_field].min(), 1.5*bindata[bin_field].max()]
        plt.xlabel=bin_field
        ylabel=r'$\Delta \gamma/\gamma$'
        plt.ylabel = ylabel

        xdata=bindata[bin_field]

        zero=zeros(xdata.size)
        zero_plt=Curve(xdata, zero)
        plt.add(zero_plt)

        xfill=[xdata.min(), xdata.max()]

        plt.add( biggles.FillBetween(xfill, [0.004,0.004], 
                                     xfill, [-0.004,-0.004],
                                     color='grey80'))


        gfrac = bindata[gfield]/gtrue-1
        gfrac_err = bindata[gfield+'_err']/gtrue
        type='filled circle'
        color='blue'
        gpts = Points(xdata, gfrac, type=type, color=color)
        gerrpts = SymmetricErrorBarsY(xdata, gfrac, gfrac_err,color=color)

        plt.add( gpts, gerrpts )

        self.plt=plt

 
    
    def write(self):
        import converter

        pstr=[str(s) for s in self.psfnums]
        pstr='-'.join(pstr)

        extra=self.get_epsfile_extra()
        path=files.get_summary_plot_path(ftype='shear',
                                         run=self.run,
                                         shnum=self.shnum,
                                         psfnum=pstr,
                                         extra=extra)
        pngpath=path.replace('.eps','.png')
        dir=os.path.dirname(path)
        if not os.path.exists(dir):
            try:
                os.makedirs(dir)
            except:
                pass

        print 'writing:',path
        self.plt.write_eps(path)
        print 'writing:',path.replace('.eps','.png')
        converter.convert(path, dpi=90)
        #print 'writing:',pngpath
        #self.plt.write_img(800,800,pngpath)

    def get_epsfile_extra(self):
        extra=[self.setname]
        sr=self.sratio_range
        extra += ['srat%0.2f' % sr[0]]
        extra += ['v%s' % self.bin_field]

        if self.options.frac:
            extra += ['frac']
        extra='-'.join(extra)
        extra = extra.replace('_','-')
        return extra


def main():
    sp=ShearPlotter()
    if sp.doshow:
        sp.show()
    sp.write()

main()
