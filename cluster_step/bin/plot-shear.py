"""
    %prog [options] run shearnum

shearnums 1-8
psfnums 1-6
"""

import sys
import os
from numpy import zeros, sqrt

import cluster_step
from cluster_step import files, stats

import esutil as eu
from esutil.numpy_util import aprint
import biggles
from biggles import FramedPlot, FramedArray, Points, \
        Curve, SymmetricErrorBarsX, SymmetricErrorBarsY

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option('-p','--psfnums',default=None,
                  help='restrict to these PSFs, comma separated')
parser.add_option('-f','--field',default='s2n_w',
                  help="bin by this field, default s2n_w")
parser.add_option('-n','--nperbin',default=4000,
                  help="bin by this field, default s2n_w")
parser.add_option('-t','--type',default=None,
                  help="limit to objects best fit by this model")
parser.add_option('-s','--show',action='store_true',
                  help="show the plot on the screen")

parser.add_option('--s2',default=None,
                  help='restrict s2 less than this value')

parser.add_option('--frac',action='store_true',
                  help=("show the fraction relative to the"
                        "expected truth"))

parser.add_option('-P','--progress',action='store_true',
                  help="show the progress bar")
sh1exp={8:0.15}
sh2exp={8:0}


class ShearPlotter(object):
    def __init__(self):
        options,args = parser.parse_args(sys.argv[1:])

        self.options=options

        if len(args) < 2:
            parser.print_help()
            sys.exit(45)

        self.run=args[0]
        self.shnum=int(args[1])
        self.nperbin=int(options.nperbin)
        self.objtype=options.type
        self.doshow = options.show

        if self.objtype:
            print 'selecting type:',self.objtype

        self.s2_max=options.s2
        if self.s2_max:
            self.s2_max=float(self.s2_max)

        self.bin_field=options.field

        self.data=files.read_output_set(self.run, 
                                        options.psfnums, 
                                        self.shnum, 
                                        objtype=self.objtype,
                                        s2_max=self.s2_max,
                                        progress=options.progress)
        
        self.set_bindata()
        self.set_psfnums_string()
        if self.options.frac:
            self.make_frac_plot()
        else:
            self.make_plot()

    def set_bindata(self):
        self.bindata=stats.bin_shear_data(self.data, self.bin_field, self.nperbin)
        aprint(self.bindata, header=True, page=False, fancy=True)

    def set_psfnums_string(self):
        self.psfnums_string=None
        if self.options.psfnums:
            pn='-'.join(self.options.psfnums.split(','))
            self.psfnums_string=pn

    def get_title(self):
        title=self.run


        if self.psfnums_string:
            title='%s-p%s' % (title,self.psfnums_string)

        title += '-s%s' % self.shnum

        if self.objtype:
            title = '%s %s' % (title,self.objtype)

        if self.options.s2:
            title = r'%s $\sigma^2_{psf}/\sigma^2_{gal} < %s$' % (title,self.options.s2)
            #sratio=sqrt(1/self.s2_max)
            #title = r'%s $\sigma_{gal}/\sigma_{psf} > %.2f$' % (title,sratio)

        return title


    def show(self):
        self.plt.show()

    def make_plot(self):

        bindata=self.bindata
        bin_field=self.bin_field
        arr=FramedArray(2,1)
        
        arr.title=self.get_title()

        arr.uniform_limits=1
        arr.xlog=True
        arr.xrange=[0.5*bindata[bin_field].min(), 1.5*bindata[bin_field].max()]
        arr.xlabel=bin_field
        arr.ylabel = r'$\gamma$'

        xdata=bindata[bin_field]
        xerr=bindata[bin_field+'_err']

        if self.shnum in sh1exp:
            g1exp=zeros(xdata.size)+sh1exp[self.shnum]
            g2exp=zeros(xdata.size)+sh2exp[self.shnum]
            g1exp_plt=Curve(xdata, g1exp)
            g2exp_plt=Curve(xdata, g2exp)
            arr[0,0].add(g1exp_plt)
            arr[1,0].add(g2exp_plt)


        xerrpts1 = SymmetricErrorBarsX(xdata, bindata['g1'], xerr)
        xerrpts2 = SymmetricErrorBarsX(xdata, bindata['g2'], xerr)

        g1pts = Points(xdata, bindata['g1'])
        g1errpts = SymmetricErrorBarsY(xdata, bindata['g1'], bindata['g1_err'])
        g2pts = Points(xdata, bindata['g2'])
        g2errpts = SymmetricErrorBarsY(xdata, bindata['g2'], bindata['g2_err'])

        arr[0,0].add( xerrpts1, g1pts, g1errpts )
        arr[1,0].add( xerrpts2, g2pts, g2errpts )


        self.plt=arr

    def make_frac_plot(self):
        if self.shnum not in sh1exp:
            raise ValueError("you must know the expected value")

        if sh1exp[self.shnum] != 0:
            gfield='g1'
            gtrue=sh1exp[self.shnum]
        elif sh1exp[self.shnum] != 0:
            gfield='g2'
            gtrue=sh2exp[self.shnum]
        else:
            raise ValueError("none of expected values are > 0")

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
        gpts = Points(xdata, gfrac)
        gerrpts = SymmetricErrorBarsY(xdata, gfrac, gfrac_err)

        plt.add( gpts, gerrpts )

        self.plt=plt

 
    
    def write(self):
        import converter

        biggles.configure( 'default', 'fontsize_min', 2)

        extra=self.get_epsfile_extra()
        path=files.get_summary_plot_path(ftype='shear',
                                         run=self.run,
                                         shnum=self.shnum,
                                         psfnum=self.psfnums_string,
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
        extra=[]
        if self.objtype:
            extra += [self.objtype]
        if self.s2_max:
            extra += ['s2%0.2f' % self.s2_max]
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
