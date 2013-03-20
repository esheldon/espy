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

        reader=files.Reader(run=self.run, 
                            shnums=self.shnum,
                            psfnums=self.psfnums,
                            s2n_range=self.s2n_range,
                            sratio_range=self.sratio_range,
                            setname=None,
                            ignore_missing=options.ignore_missing,
                            verbose=True,
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
        wstar,=where( (data['uw_flags']==0) 
                        & ((data['uw_type'] & CLUSTERSTEP_PSF_STAR) != 0) )

        #s2n_min=700
        #s2n_max=1000
        s2n_min=70
        s2n_max=80
        Tadmom=data['am_irr']+data['am_icc']
        wgal,=where( (data['uw_flags']==0) 
                        & ((data['uw_type'] & CLUSTERSTEP_GAL) != 0)
                        & (Tadmom > 12)
                        & (data['am_s2n'] > s2n_min) 
                        & (data['am_s2n'] < s2n_max))

        print 'ngal:',wgal.size

        star_isum = data['uw_isum'][wstar].sum()
        star_irr  = data['uw_irrsum'][wstar].sum()/star_isum
        star_irc  = data['uw_ircsum'][wstar].sum()/star_isum
        star_icc  = data['uw_iccsum'][wstar].sum()/star_isum

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
        plt=FramedPlot()
        
        plt.uniform_limits=1
        plt.xlog=True
        plt.xrange=self.s2n_range
        plt.xlabel='S/N'
        plt.ylabel = r'$\gamma$'
        if self.yrange is not None:
            plt.yrange=self.yrange

        xdata=bindata['s2n']

        if self.shnum in sh1exp:
            g1exp=zeros(xdata.size)+sh1exp[self.shnum]
            g2exp=zeros(xdata.size)+sh2exp[self.shnum]
            g1exp_plt=Curve(xdata, g1exp)
            g2exp_plt=Curve(xdata, g2exp)
            plt.add(g1exp_plt)
            plt.add(g2exp_plt)



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

        plt.add( g1pts, g1errpts )
        plt.add( g2pts, g2errpts )
        plt.add(key)

        labels=self.get_labels()

        plt.add(*labels)
        plt.aspect_ratio=1

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
