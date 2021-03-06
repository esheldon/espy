"""
Binner classes
    LambdaBinner: bin by redmapper lambda
    N200Binner: bin by ngals_r200
    MZBinner: bin by true mass and redshift

    To use these:
        nb = N200Binner(12)
        nb.bin_byrun('05')

Standalone functions:
    reduce_from_ranges
    logbin_linear_edges
"""
from __future__ import print_function
import os
import sys
import copy
from sys import stdout
import numpy
from numpy import log10,sqrt,linspace,where, zeros
import esutil as eu
from esutil.ostools import path_join
from esutil.stat import histogram
from esutil.numpy_util import where1
from esutil.ostools import makedirs_fromfile

import converter

import lensing

from .util import lens_wmom

try:
    import biggles
    from biggles import FramedArray, FramedPlot, Points, \
            ErrorBarsY, ErrorBarsX, \
            SymmetricErrorBarsY, SymmetricErrorBarsX, \
            PlotKey, PlotLabel, Table, Curve
except:
    pass

def instantiate_binner(type, nbin=None):
    if nbin is None:
        # this is the bin type name, corresponds to a file
        # bin-{type}.yaml
        b = AnyBinner(type)
    else:
        if type == 'lambda':
            b = LambdaBinner(nbin)
        elif type=='ilum':
            b = ILumBinner(nbin)
        elif type == 'n200':
            b = N200Binner(nbin)
        elif type == 'mz':
            b = MZBinner(nbin)
        elif type == 'vz':
            b = VoidZBinner(nbin)
        else:
            raise ValueError("unsupported binner type: '%s'" % type)
    return b


def bin_lenses_byrun(run, type, nbin=None):
    b=instantiate_binner(type, nbin=nbin)
    b.bin_byrun(run)

class BinnerBase(dict):
    def __init__(self, nbin, fs='nfs', **keys):
        self['nbin'] = nbin
        self.update(keys)

        self.fs=fs
        self.set_bin_ranges()

        self.dpi=150

    def get_name(self):
        raise RuntimeError("override this method")

    def bin_byrun(self, run):
        """
        Do the binning and write out a file
        """

        name=self.get_name()
        d = lensing.files.sample_read(type='collated',sample=run,fs=self.fs)
        res = self.bin(d)

        outdir = lensing.files.sample_dir(type='binned',sample=run,name=name)
        if not os.path.exists(outdir):
            print("Making output dir:",outdir)
            os.makedirs(outdir)
        lensing.files.sample_write(data=res,
                                   type='binned',
                                   sample=run,
                                   name=name,
                                   clobber=True)

    def bin(self, data):
        raise RuntimeError("override this method")

    def set_bin_ranges(self, binnum=None):
        raise RuntimeError("override this method")

    def bin_ranges(self, binnum=None):
        if binnum is not None:
            if (binnum < 0) or (binnum > (self['nbin']-1)):
                raise ValueError("binnum out of bounds: [%d,%d]" % (0,self['nbin']-1))
            return self.lowlim[binnum], self.highlim[binnum]

        return self.lowlim, self.highlim

    def bin_label(self, binnum):
        raise RuntimeError("override this method")

    def get_bin_erf_weights_1d(self, binnum, x, sigma):
        """
        For a 1-d binning and given sigma, get a 2-sided erf
        weighting each point.  This is a probalistic way to 
        assign objects to a bin range.
        """
        low, high = self.bin_ranges(binnum=binnum)
        weights = two_sided_erf_weights(low, high, x, sigma)
        return weights

    def compare_random(self, lensrun, type, binnum, randrun, **keys):
        xrng=keys.get('xrange',None)
        yrng=keys.get('yrange',None)

        extra='randmatch-%s' % randrun
        data = lensing.files.sample_read(type=type, sample=lensrun, name=self.get_name())
        rand = lensing.files.sample_read(type='binned', sample=lensrun, name=self.get_name(), 
                                         extra=extra)

        zpts=biggles.Curve( [0.9*data['r'].min(), 1.1*data['r'].max()], [0,0])

        plt=biggles.FramedPlot()
        plt.aspect_ratio=1
        plt.xlog=True
        if xrng is not None:
            plt.xrange=xrng
        if yrng is not None:
            plt.yrange=yrng

        plt.xlabel=r'$r [h^{-1} Mpc]$'
        plt.ylabel=r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'

        pts=biggles.Points(data['r'][binnum], data['dsig'][binnum],
                           type='filled circle')
        epts=biggles.SymmetricErrorBarsY(data['r'][binnum], data['dsig'][binnum], 
                                         data['dsigerr'][binnum])

        rpts=biggles.Points(rand['r'][binnum], rand['dsig'][binnum],
                           type='filled circle',color='red')
        repts=biggles.SymmetricErrorBarsY(rand['r'][binnum], rand['dsig'][binnum], 
                                          rand['dsigerr'][binnum],color='red')

        pts.label=r'$\Delta\Sigma_+$'
        rpts.label=r'$\Delta\Sigma^{rand}_+$'

        key=biggles.PlotKey(0.9,0.9,[pts,rpts],halign='right')
        lab=biggles.PlotLabel(0.1,0.9,self.bin_label(binnum),halign='left')

        plt.add(zpts,pts,epts,rpts,repts,key,lab)
            
        epsfile=lensing.files.sample_file(type=type+'-plots',
                                          sample=lensrun,
                                          name=self.get_name(),
                                          extra='%s-%02i' % (randrun,binnum),
                                          ext='eps')
        plt.write_eps(epsfile)
        converter.convert(epsfile, dpi=self.dpi, verbose=True)

    def plot_dsig_osig_byrun_bin(self, run, type, binnum, **keys):
        """

        Plot a single bin, with dsig and osig.

        See lensing.plotting.plot2dsig
        """

        linear=keys.get('linear',False)

        name = self.get_name()
        data=lensing.files.sample_read(type=type,sample=run,name=name)

        max_binnum = self['nbin']-1
        if (binnum < 0) or (binnum > max_binnum):
            raise ValueError("binnum is out of range [0,%d]" % max_binnum)
        data = data[binnum]

        keys['plot_label'] = self.bin_label(binnum)
        keys['xmul'] = 1.1

        if linear:

            yrange=keys.get('yrange',None)
            xrng=keys.get('xrange',None)

            if xrng is None:
                xrng=eu.plotting.get_log_plot_range(data['r'])

            fac=1.1
            rfac=data['r']*fac
            opts=Points(rfac, data['osig'], color='red', type='filled diamond')
            oepts=SymmetricErrorBarsY(rfac, data['osig'], data['dsigerr'],color='red')
            oc=Curve(rfac, data['osig'], color='red')

            pts=Points(data['r'], data['dsig'], type='filled circle')
            c=Curve(data['r'], data['dsig'])
            epts=SymmetricErrorBarsY(data['r'], data['dsig'], data['dsigerr'])

            opts.label=r'$\Delta\Sigma_\times$'
            pts.label=r'$\Delta\Sigma_+$'

            zpts=Curve(xrng, [0,0])

            key=PlotKey(0.9,0.9,[pts,opts],halign='right')

            lab=PlotLabel(0.1,0.9,keys['plot_label'],halign='left')

            plt=FramedPlot()
            plt.xlog=True
            plt.aspect_ratio=1
            plt.yrange=yrange
            plt.xrange=xrng
            plt.xlabel=lensing.plotting.labels['rproj']
            plt.ylabel=lensing.plotting.labels['dsig']

            plt.add(zpts,opts,oepts,oc,pts,epts,c,key,lab)


            show=keys.get('show',False)
            if show:
                plt.show()
            
            """
            keys1=copy.deepcopy(keys)
            keys1['show']=False

            plt1=lensing.plotting.plot_dsig(r=data['r'], 
                                            dsig=data['osig'],
                                            dsigerr=data['dsigerr'],
                                            ylog=False,
                                            color='red',
                                            **keys1)

            plt=lensing.plotting.plot_dsig(comb=data, ylog=False, plt=plt1,**keys,
                                           label=r'$\Delta\Sigma_+$')
            """


        else:
            if 'yrange' in keys:
                keys['yrange2'] = keys['yrange']
            plt=lensing.plotting.plot2dsig(data['r'], 
                                           data['dsig'], data['dsigerr'],
                                           data['r'],
                                           data['osig'], data['dsigerr'],
                                           **keys)

        return plt

    def plot_dsig_osig_byrun(self, run, type, **keys):
        """

        Plot all bins dsig+osig but one at a time, not in
        an array

        """
        for binnum in xrange(self['nbin']):
            plt=self.plot_dsig_osig_byrun_bin(run, type, binnum, **keys)
            epsfile=lensing.files.sample_file(type=type+'-plots',
                                              sample=run,
                                              name=self.get_name(),
                                              extra='osig-comp-%02i' % binnum,
                                              ext='eps')
            makedirs_fromfile(epsfile)
            stdout.write("Plotting to file: %s\n" % epsfile)
            plt.write_eps(epsfile)
            converter.convert(epsfile, dpi=self.dpi, verbose=True)


    def plot_dsig_byrun_1var(self, run, type, yrnge=None, xrnge=None, show=False):
        """

        Make an array of plots with each plot a bin in one variable.  
        
        It is expected that methods like self.bin_label() are meaningfully
        over-ridden

        parameters
        ----------
        run: string
            The lensing run id
        type:
            The type to read, e.g. 'binned', 'corrected'

        """

        name = self.get_name()
        data=lensing.files.sample_read(type=type,sample=run,name=name)

        # this is for the screen: currently tuned for my big screen!
        biggles.configure('screen','width', 1140)
        biggles.configure('screen','height', 1140)
        biggles.configure('fontsize_min',1.0)
        biggles.configure('_HalfAxis','ticks_size',3)
        biggles.configure('_HalfAxis','subticks_size',1.5)
        biggles.configure('linewidth',1)

        if self['nbin'] == 12:
            nrow = 3
            ncol = 4
            aspect_ratio = 1.0/1.5
        elif self['nbin'] == 16:
            nrow = 4
            ncol = 4
            aspect_ratio = 1.0
        elif self['nbin']== 10:
            nrow=2
            ncol=5
            aspect_ratio=2./5.
        elif self['nbin']== 9:
            nrow=3
            ncol=3
            aspect_ratio=1.0
        elif self['nbin'] == 6:
            nrow=2
            ncol=3
            aspect_ratio=2./3.
        elif self['nbin'] == 4:
            nrow = 2
            ncol = 2
            aspect_ratio=1.0
        elif self['nbin'] == 2:
            nrow = 1
            ncol = 2
            aspect_ratio=float(nrow)/ncol
        elif self['nbin']==1:
            nrow=1
            ncol=1
            aspect_ratio=1.0
        else:
            raise ValueError("Unsupported nbin: %s" % self['nbin'])

        if self['nbin']==1:
            pa = FramedPlot()
        else:
            pa = FramedArray(nrow, ncol)

        pa.aspect_ratio = aspect_ratio

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'

        if xrnge is None:
            xrnge = [0.01,60.0]
        if yrnge is None:
            yrnge = [1.e-2,8000]
        pa.xrange = xrnge
        pa.yrange = yrnge
        pa.xlog = True
        pa.ylog = True


        row = -1

        i = 0
        for i in xrange(self['nbin']):
            col = i % ncol
            if col == 0:
                row += 1

            if self['nbin']==1:
                send_plt=pa
            else:
                send_plt=pa[row,col]
            eu.plotting.bscatter(data['r'][i],
                                 data['dsig'][i],
                                 yerr=data['dsigerr'][i],
                                 xrange=xrnge,
                                 yrange=yrnge,
                                 xlog=True,ylog=True,
                                 show=False, 
                                 plt=send_plt,
                                 size=2)
            label = self.bin_label(i)
            pl = PlotLabel(.85, .85, label, halign='right')

            send_plt.add(pl)


        if show:
            pa.show()

        epsfile=lensing.files.sample_file(type=type+'-plots',
                                          sample=run,
                                          name=name,
                                          extra='allplot',
                                          ext='eps')
        d = os.path.dirname(epsfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
        converter.convert(epsfile, dpi=self.dpi, verbose=True)


    def plot_osig_byrun_1var(self, run, type, show=False):
        """

        Make an array of orthotangential with each plot a bin in one variable.  
        
        It is expected that methods like self.bin_label() are meaningfully
        over-ridden

        parameters
        ----------
        run: string
            The lensing run id
        type:
            The type to read, e.g. 'binned', 'corrected'

        """

        name = self.get_name()
        data=lensing.files.sample_read(type=type,sample=run,name=name)

        # this is for the screen: currently tuned for my big screen!
        biggles.configure('screen','width', 1140)
        biggles.configure('screen','height', 1140)
        biggles.configure('fontsize_min',1.0)
        biggles.configure('_HalfAxis','ticks_size',3)
        biggles.configure('_HalfAxis','subticks_size',1.5)
        biggles.configure('linewidth',1)

        if self['nbin'] == 12:
            nrow = 3
            ncol = 4
            aspect_ratio = 1.0/1.5
        elif self['nbin'] == 16:
            nrow = 4
            ncol = 4
            aspect_ratio = 1.0
        elif self['nbin']== 10:
            nrow=2
            ncol=5
            aspect_ratio=2./5.
        elif self['nbin']== 9:
            nrow=3
            ncol=3
            aspect_ratio=1.
        elif self['nbin'] == 2:
            nrow=2
            ncol=1
            aspect_ratio=2.0
        else:
            raise ValueError("Unsupported nbin: %s" % self['nbin'])

        pa = FramedArray(nrow, ncol)
        pa.aspect_ratio = aspect_ratio

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma_{\times} [M_{\odot} pc^{-2}]$'

        xrnge = [0.01,60.0]
        yrnge = [-20,20]
        pa.xrange = xrnge
        pa.yrange = yrnge
        pa.xlog = True


        row = -1

        low, high = self.bin_ranges()

        i = 0
        for i in xrange(self['nbin']):
            col = i % ncol
            if col == 0:
                row += 1

            eu.plotting.bscatter(data['r'][i],
                                 data['osig'][i],
                                 yerr=data['dsigerr'][i],
                                 xrange=xrnge,
                                 yrange=yrnge,
                                 xlog=True,ylog=False,
                                 show=False, 
                                 plt=pa[row,col],
                                 size=2)
            label = self.bin_label(i)
            pl = PlotLabel(.85, .85, label, halign='right')

            pa[row,col].add(pl)

            line=biggles.Curve([xrnge[0],xrnge[1]],[0,0])
            pa[row,col].add(line)


        if show:
            pa.show()

        epsfile=lensing.files.sample_file(type=type+'-plots',sample=run,name=name,extra='osig-allplot',ext='eps')
        d = os.path.dirname(epsfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
        converter.convert(epsfile, dpi=self.dpi, verbose=True)

    def plot_dsig_2runs(self, run1, run2, type,show=False, **kw):
        """

        Make an array of plots with each plot a bin in one variable.  
        Overplot a second run.
        
        It is expected that methods like self.bin_label() are meaningfully
        over-ridden

        parameters
        ----------
        run: string
            The lensing run id
        type:
            The type to read, e.g. 'binned', 'corrected'

        """

        name = self.get_name()
        data1=lensing.files.sample_read(type=type,sample=run1,name=name)
        data2=lensing.files.sample_read(type=type,sample=run2,name=name)
        color1='blue'
        color2='red'
        sym1='filled circle'
        sym2='filled diamond'
        width=4
        size1=2
        size2=3

        # this is for the screen: currently tuned for my big screen!
        biggles.configure('screen','width', 1140)
        biggles.configure('screen','height', 1140)
        biggles.configure('fontsize_min',1.0)
        biggles.configure('_HalfAxis','ticks_size',3)
        biggles.configure('_HalfAxis','subticks_size',1.5)
        biggles.configure('linewidth',0.6)

        doarray=True
        if self['nbin'] == 12:
            nrow = 3
            ncol = 4
            aspect_ratio = 1.0/1.5
        elif self['nbin'] == 16:
            nrow = 4
            ncol = 4
            aspect_ratio = 1.0
        elif self['nbin']==1:
            doarray=False
        else:
            raise ValueError("Unsupported nbin: %s" % self['nbin'])

        if doarray:
            pa = FramedArray(nrow, ncol)
            pa.aspect_ratio = aspect_ratio
        else:
            pa=FramedPlot()
            pa.aspect_ratio=1

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'

        linear=kw.get('linear',False)
        xrnge=kw.get('xrange',None)
        yrnge=kw.get('yrange',None)

        if xrnge is None:
            xrnge = [0.01,60.0]
        if not linear:
            pa.ylog = True
            if yrnge is None:
                yrnge = [1.e-2,8000]
        else:
            if yrnge is None:
                ymin=min(data1['dsig'].min(), data2['dsig'].min())
                if ymin > 0:
                    ymin=0.9*ymin
                else:
                    ymin=1.1*ymin
                ymax=max(data1['dsig'].max(), data2['dsig'].max())
                if ymax > 0:
                    ymax=1.1*ymax
                else:
                    ymax=0.9*ymax


        pa.xrange = xrnge
        pa.yrange = yrnge
        pa.xlog = True



        row = -1

        i = 0
        for i in xrange(self['nbin']):
            if doarray:
                col = i % ncol
                if col == 0:
                    row += 1
                plt=pa[row,col]
            else:
                plt=pa

            pdict2=eu.plotting.bscatter(data2['r'][i],
                                        data2['dsig'][i],
                                        yerr=data2['dsigerr'][i],
                                        xrange=xrnge,
                                        yrange=yrnge,
                                        xlog=True,ylog=True,
                                        show=False, 
                                        plt=plt,
                                        color=color2,
                                        label=run2,
                                        type=sym2,
                                        size=size2,width=width,
                                        dict=True)

            pdict1=eu.plotting.bscatter(data1['r'][i],
                                        data1['dsig'][i],
                                        yerr=data1['dsigerr'][i],
                                        xrange=xrnge,
                                        yrange=yrnge,
                                        xlog=True,ylog=True,
                                        show=False, 
                                        plt=plt,
                                        color=color1,
                                        type=sym1,
                                        label=run1,
                                        size=size1,width=width,
                                        dict=True)

            if i == (self['nbin']-1):
                key=PlotKey(0.2,0.2,[pdict1['p'],pdict2['p']])
                plt.add(key)

            label = self.bin_label(i)
            pl = PlotLabel(.85, .85, label, halign='right')

            plt.add(pl)

            if linear:
                plt.add(biggles.Curve([xrnge[0],xrnge[1]],[0.0,0.0]))


        if show:
            pa.show()

        if linear:
            extra='linear-allplot-'+run2
        else:
            extra='allplot-'+run2
        epsfile=lensing.files.sample_file(type=type+'-plots',
                                          sample=run1,
                                          name=name,
                                          extra=extra,
                                          ext='eps')
        d = os.path.dirname(epsfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
        converter.convert(epsfile, dpi=self.dpi, verbose=True)



class AnyBinner(BinnerBase):
    """
    bin any
    """

    def __init__(self, bin_conf_name, fs='nfs', **keys):
        self.update(keys)
        self['name']=bin_conf_name
        self.bin_conf = lensing.files.read_config('bin',bin_conf_name)
        self.bin_info = self.bin_conf['bin_info']

        self['nbin'] = len(self.bin_info)

        self.fs=fs
        self.dpi=150

    def get_name(self):
        return self['name']
    def get_nbin(self):
        return self['nbin']

    def bin(self, data):
        """
        call also call base method bin_byrun
        """
        nrbin = data['rsum'][0].size

        bi0 = self.bin_info[0]['bins']
        bintags = [ri[0] for ri in bi0]
        print("bintags:",bintags)
        bs = lensbin_struct(nrbin, bintags=bintags, n=self['nbin'])

        print("len(data):",len(data))
        print("len(bin_info):",len(self.bin_info))
        for i,bin_info in enumerate(self.bin_info):
            range_info = bin_info['bins']

            comb,w = reduce_from_ranges_many(data,
                                             range_info,
                                             getind=True)
        
            bs['nlenses'][i] = w.size
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            for bi in range_info:
                field_name=bi[0]
                mn='%s_mean' % field_name
                en='%s_err' % field_name
                sn='%s_sdev' % field_name
                rn='%s_range' % field_name

                # now the things we are averaging by lens weight
                mean,err,sdev = lens_wmom(data, field_name, ind=w, sdev=True)
                bs[mn][i] = mean
                bs[en][i] = err
                bs[sn][i] = sdev
                bs[rn][i] = data[field_name][w].min(), data[field_name][w].max()

            i+=1

        return bs

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        range_info = self.bin_info[binnum]['bins']
        comb,w = reduce_from_ranges_many(data,
                                         range_info,
                                         getind=True)
 
        return w



    def bin_label(self, binnum):
        return self.bin_info[binnum]['label']


class VoidZBinner(BinnerBase):
    range_type='()'

    def get_name(self):
        return 'z-%02d' % self['nbin']

    def set_bin_ranges(self):
        if self['nbin'] == 4:
            lowlim  = [0.02,0.05, 0.098, 0.148]
            highlim = [0.05,0.098,0.148, 0.20]
        elif self['nbin'] == 2:
            lowlim  = [0.02,0.1]
            highlim = [0.1,0.2]
        elif self['nbin'] == 1:
            lowlim = [0.02]
            highlim= [0.20]

        self.lowlim = lowlim
        self.highlim = highlim

    def bin(self, data):
        """
        call also call base method bin_byrun
        """
        from math import ceil

        low, high = self.bin_ranges()

        nrbin = data['rsum'][0].size
        bs = lensbin_struct(nrbin, bintags=['z','radius'], n=self['nbin'])

        i=0
        for l,h in zip(low,high):
            zrange = [l,h]

            print("l,h:",l,h)
            print('%0.2f < z < %0.2f' % tuple(zrange))

            print("    reducing and jackknifing by lens")
            comb,w = reduce_from_ranges(data,
                                        'z',
                                        zrange, 
                                        range_type=self.range_type,
                                        getind=True)
        
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            # now the things we are averaging by lens weight
            mn,err,sdev = lens_wmom(data,'z',ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = data['z'][w].min(), data['z'][w].max()

            mn,err,sdev = lens_wmom(data,'radius',ind=w, sdev=True)
            bs['radius_mean'][i] = mn
            bs['radius_err'][i] = err
            bs['radius_sdev'][i] = sdev
            bs['radius_range'][i] = zrange
            bs['radius_minmax'][i] = \
                data['radius'][w].min(), data['radius'][w].max()

            # fix up last one
            #if i == (bs.size-1):
            #    bs['lambda_range'][i,1] = ceil(bs['lambda_minmax'][i,1])

            i+=1

        return bs

    def bin_label(self, binnum):
        zrange = self.bin_ranges(binnum)
        return r'$%0.2f < z < %0.2f$' % zrange

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        low, high = self.bin_ranges()
        logic = get_range_logic(data, 'z', [low[binnum], high[binnum]], self.range_type)
        return where1(logic)




def define_lambda_bins(sample, lastmin=58.):
    l=lensing.lcat.read_original(sample=sample)
    w=where1(l['z_lambda'] < 0.4)
    l=l[w]

    # 12 is unimportant here
    lb = instantiate_binner('lambda', 12)
    lb.define_bins(l['lambda_chisq'], lastmin=lastmin)

class LambdaBinner(BinnerBase):
    """
    RedMapper binner on lambda
    """
    lambda_field = 'lambda_chisq'
    range_type='()'

    # play with the range
    e_field = 'pe'
    e_range=[0.0, 0.2]
    e_range_type='[)'

    z_field = 'z_lambda'
    z_range_type='()'

    def get_name(self):
        return 'lambda-%02d' % self['nbin']

    def bin(self, data):
        """
        call also call base method bin_byrun
        """
        from math import ceil
        lambda_field = self.lambda_field
        z_field = self.z_field

        low, high = self.bin_ranges()

        nrbin = data['rsum'][0].size
        bs = lensbin_struct(nrbin, bintags=['lambda','z'], n=self['nbin'])

        i=0
        for l,h in zip(low,high):
            lamrange = [l,h]

            print("l,h:",l,h)
            if h is not None:
                print('%0.2f < lambda < %0.2f' % tuple(lamrange))
            else:
                print('lambda > %0.2f' % l)

            print("    reducing and jackknifing by lens")
            range_info=[(lambda_field, lamrange, self.range_type)]

            if 'zrange' in self and self['zrange'] is not None:
                zrange=self['zrange']
                zrange_info = (z_field, zrange, self.z_range_type)
                range_info.append(zrange_info)
            elif hasattr(self,'z_lowlim'):
                zrange=[self.z_lowlim[i], self.z_highlim[i]]
                zrange_info = (z_field, zrange, self.z_range_type)
                range_info.append(zrange_info)

            if self.e_field in data.dtype.names:
                e_range_info = (self.e_field, self.e_range, self.e_range_type)
                range_info.append(e_range_info)

            comb,w = reduce_from_ranges_many(data,
                                             range_info,
                                             getind=True)
        
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            # now the things we are averaging by lens weight
            mn,err,sdev = lens_wmom(data,z_field,ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = data[z_field][w].min(), data[z_field][w].max()

            mn,err,sdev = lens_wmom(data,lambda_field,ind=w, sdev=True)
            bs['lambda_mean'][i] = mn
            bs['lambda_err'][i] = err
            bs['lambda_sdev'][i] = sdev
            bs['lambda_range'][i] = lamrange
            bs['lambda_minmax'][i] = \
                data[lambda_field][w].min(), data[lambda_field][w].max()

            # fix up last one
            if i == (bs.size-1):
                bs['lambda_range'][i,1] = ceil(bs['lambda_minmax'][i,1])

            i+=1

        return bs

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        low, high = self.bin_ranges()
        logic = get_range_logic(data, self.lambda_field, [low[binnum], high[binnum]], self.range_type)
        return where1(logic)



    def bin_label(self, binnum):
        lrange = self.bin_ranges(binnum)
        if lrange[0] is None and lrange[1] is None:
            raise ValueError("expected at least one in range to be not None")
        elif lrange[1] == 1.e6:
            return r'$\lambda\ > %0.1f$' % lrange[0]
        elif lrange[0] is None:
            return r'$\lambda\ < %0.1f$' % lrange[1]
        else:
            return r'$%0.1f < \lambda\ < %0.1f$' % lrange


    def define_bins(self, lambda_input, lastmin=58., alpha=0.6666, show=False, prompt=False):
        """
        Data should be trimmed to a relevant sample, e.g. 
        in redshift and lambda.

        choose a lower bound for lambda on the last bin, then
        work downward keeping N*lambda^alpha = constant

        Note this will not produce a fixed number of bins; we may
        want to adjust.
        """
        
        lam = lambda_input.copy()
        lam.sort()

        # reverse sort, biggest to smallest
        lam = numpy.fromiter(reversed(lam), dtype='f8')

        ind = numpy.arange(lam.size,dtype='i8')

        w_last = where1(lam > lastmin)
        mlam_last = lam[w_last].sum()/w_last.size

        print("wlast.size:",w_last.size,"mlam_last:",mlam_last)

        cref = 0.5*log10(w_last.size) + alpha*log10(mlam_last)

        # now look at mean lambda*N for rest by using cumulative
        # sums


        lam = lam[w_last[-1]:]
        lam_last = lastmin
        binnum=0

        minbin=[]
        maxbin=[]
        while 1:
            nd = 1+numpy.arange(lam.size)

            mlam = lam.cumsum()/nd
            cval = 0.5*log10(nd) + alpha*log10(mlam)

            #wthis = where1(cval > cref)
            wthis = where1(cval < cref)

            lam_min = lam[wthis[-1]]
            lam_max = lam_last

            minbin.append(lam_min)
            maxbin.append(lam_max)
            if show:
                print("numleft:",lam.size,"wthis.size",wthis.size)

                plt=FramedPlot()
                curve = Curve(lam, cval)
                plt.add(curve)
                #points = Points(lam, cval, type='filled circle', size=0.5)
                #plt.add(points)

                oc = Curve([lam.min(), lam.max()],[cref,cref],color='blue')
                plt.add(oc)

                #p = biggles.Point(lam[wthis.min()], cval[wthis.min()], color='orange', type='filled circle')
                p = biggles.Point(lam[wthis[-1]], cval[wthis[-1]], color='orange', type='filled circle')
                plt.add(p)

                binnum -= 1
                blab=PlotLabel(0.9,0.9,'bin %d' % binnum, halign='right')

                rlab=PlotLabel(0.9,0.85,r'%0.2f $ < \lambda < $ %0.2f' % (lam_min,lam_max), halign='right')
                nlab=PlotLabel(0.9,0.80,'N: %d' % wthis.size, halign='right')

                plt.add(blab)
                plt.add(rlab)
                plt.add(nlab)
                plt.show()

                if prompt:
                    key=raw_input('hit a key: ')
                    if key == 'q':
                        return
            lam_last = lam_min
            lam=lam[wthis[-1]+1:]


            if len(lam) == 0:
                break
        
        minbin = list(reversed(minbin))
        maxbin = list(reversed(maxbin))
        minbin.append(maxbin[-1])
        maxbin.append(None)
        
        minstr=[]
        maxstr=[]
        for i in xrange(len(minbin)):
            if minbin[i] is not None:
                minstr.append('%0.1f' % minbin[i])
            else:
                minstr.append('None')
            if maxbin[i] is not None:
                maxstr.append('%0.1f' % maxbin[i])
            else:
                maxstr.append('None')

        minstr = '[' + ', '.join(minstr) +']'
        maxstr = '[' + ', '.join(maxstr) +']'

        #for i in xrange(len(minbin)):
        #    print('%i %0.1f %0.1f' % (i+1,minbin[i],maxbin[i]))

        print("nbin:",len(minbin))
        print(minstr)
        print(maxstr)

    def set_bin_ranges(self):
        if self['nbin'] == 12:
            # ran define_bins
            lowlim  = [10.0, 10.4, 11.7, 13.3, 15.1, 17.3, 19.8, 23.0, 27.2, 32.9, 41.6, 58.0]
            highlim = [10.4, 11.7, 13.3, 15.1, 17.3, 19.8, 23.0, 27.2, 32.9, 41.6, 58.0, 1.e6]
        elif self['nbin'] == 10:
            lowlim=[20.0, 21.5, 24.2, 27.0, 30.4, 34.4, 39.5, 46.3, 56.1, 75.0]
            highlim=[21.5, 24.2, 27.0, 30.4, 34.4, 39.5, 46.3, 56.1, 75.0, 1.e6]
        elif self['nbin']==9:
            lowlim=[20.0, 21.3, 24.4, 27.7, 31.7, 36.6, 43.3, 53.0, 71.0]
            highlim=[21.3, 24.4, 27.7, 31.7, 36.6, 43.3, 53.0, 71.0, 1.e6]
        elif self['nbin'] == 6:
            lowlim = [20.0, 20.8, 25.9, 32.2, 41.3, 58.0]
            highlim = [20.8, 25.9, 32.2, 41.3, 58.0, 1.e6]
        elif self['nbin'] == 1:
            lowlim=[20.0]
            highlim=[1.e6]
            self.z_lowlim=[0.1]
            self.z_highlim=[0.3]
        else:
            raise ValueError("Unsupported nbin: %d\n", self['nbin'])

        #self.lowlim = numpy.array(lowlim,dtype='f8')
        #self.highlim = numpy.array(highlim,dtype='f8')
        self.lowlim = lowlim
        self.highlim = highlim

def two_sided_erf_weights(low, high, x, sigma):
    """
    x and sigma should be arrays of the same length
    """
    from scipy.special import erf
    weights = zeros( x.size )

    mn=0.5*(low + high)

    wlow,=where(x < mn)
    whigh,=where(x >= mn)

    if wlow.size > 0:
        wts=0.5*( 1.0+erf((x[wlow]-low)/sigma[wlow]) )
        weights[wlow] = wts

    if whigh.size > 0:
        wts=0.5*( 1.0+erf((high-x[whigh])/sigma[whigh]) )
        weights[whigh] = wts

    return weights

class ILumBinner(BinnerBase):
    """
    lrg des binner on lambda
    """
    range_type='()'

    # in units of L* 
    ilum_field = 'ilum'
    z_field = 'zred2'
    def get_name(self):
        return 'ilum-%02d' % self['nbin']

    def bin(self, data):
        """
        call also call base method bin_byrun
        """
        from math import ceil
        ilum_field = self.ilum_field
        z_field = self.z_field

        low, high = self.bin_ranges()

        nrbin = data['rsum'][0].size
        bs = lensbin_struct(nrbin, bintags=['ilum','z'], n=self['nbin'])

        i=0
        for l,h in zip(low,high):
            lamrange = [l,h]

            print("l,h:",l,h)
            if h is not None:
                print('%0.2f < L/L_* < %0.2f' % tuple(lamrange))
            else:
                print('L/L_*> %0.2f' % l)

            print("    reducing and jackknifing by lens")
            comb,w = reduce_from_ranges(data,
                                        ilum_field,
                                        lamrange, 
                                        range_type=self.range_type,
                                        getind=True)
        
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            # now the things we are averaging by lens weight
            mn,err,sdev = lens_wmom(data,z_field,ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = data[z_field][w].min(), data[z_field][w].max()

            mn,err,sdev = lens_wmom(data,ilum_field,ind=w, sdev=True)
            bs['ilum_mean'][i] = mn
            bs['ilum_err'][i] = err
            bs['ilum_sdev'][i] = sdev
            bs['ilum_range'][i] = lamrange
            bs['ilum_minmax'][i] = \
                data[ilum_field][w].min(), data[ilum_field][w].max()

            # fix up last one
            if i == (bs.size-1):
                bs['ilum_range'][i,1] = ceil(bs['ilum_minmax'][i,1])

            i+=1

        return bs

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        low, high = self.bin_ranges()
        logic = get_range_logic(data, self.ilum_field, [low[binnum], high[binnum]], self.range_type)
        return where1(logic)


    def bin_label(self, binnum):
        lrange = self.bin_ranges(binnum)
        if lrange[0] is None and lrange[1] is None:
            raise ValueError("expected at least one in range to be not None")
        elif lrange[1] == 1.e6:
            return r'$L/L_*> %0.1f$' % lrange[0]
        elif lrange[0] is None:
            return r'$L/L_* < %0.1f$' % lrange[1]
        else:
            return r'$%0.1f < L/L_* < %0.1f$' % lrange


    def set_bin_ranges(self):
        if self['nbin'] == 2:
            lowlim = [0.1, 1.0]
            highlim = [1.0, 1.e6]
        elif self['nbin']==4:
            lowlim  = numpy.array([-0.92, -0.35, 0.00, 0.40])
            highlim = numpy.array([-0.35,  0.00, 0.40, 1.70])
            lowlim  = 10**lowlim
            highlim = 10**highlim

            highlim[-1] = 1.e6

        elif self['nbin']==12:
            # this doesn't really work for the gaussianish ilum distribution
            lowlim=numpy.array([-0.92, -0.71880262, -0.51880262, -0.31880262, -0.11880262,
                                0.08119738,  0.28119738,  0.48119738,  0.68119738,  0.88119738,
                                1.08119738,  1.28119738])
            highlim=numpy.array([-0.71880262, -0.51880262, -0.31880262, -0.11880262,  0.08119738,
                                 0.28119738,  0.48119738,  0.68119738,  0.88119738,  1.08119738,
                                 1.28119738,  1.7])
            lowlim  = 10**lowlim
            highlim = 10**highlim


        else:
            raise ValueError("Unsupported nbin: %d\n", self['nbin'])

        self.lowlim = lowlim
        self.highlim = highlim


class N200Binner(BinnerBase):
    """
    ngals_r200 binner for original MaxBCG
    """
    #def __init__(self, nbin):
    #    self['nbin'] = nbin
    #    self.set_bin_ranges()
    range_type='[]'
    def get_name(self):
        return 'n200-%02d' % self['nbin']


    def bin(self, data):

        nlow, nhigh = self.bin_ranges()

        nrbin = data['rsum'][0].size
        bs = lensbin_struct(nrbin, bintags=['ngals200','z'], n=self['nbin'])

        i=0
        for Nl,Nh in zip(nlow,nhigh):
            Nrange = [Nl,Nh]

            print('%d <= N200 <= %d' % tuple(Nrange))

            print("    reducing and jackknifing by lens")
            comb,w = reduce_from_ranges(data,'ngals_r200',Nrange, range_type=self.range_type,
                                        getind=True)
        
            print("    found",w.size,"in bin")
            # first copy all common tags
            for n in comb.dtype.names:
                bs[n][i] = comb[n][0]

            # now the things we are averaging by lens weight
            mn,err,sdev = lens_wmom(data,'photoz_cts',ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = data['photoz_cts'][w].min(), data['photoz_cts'][w].max()

            mn,err,sdev = lens_wmom(data,'ngals_r200',ind=w, sdev=True)
            bs['ngals200_mean'][i] = mn
            bs['ngals200_err'][i] = err
            bs['ngals200_sdev'][i] = sdev
            bs['ngals200_range'][i] = Nrange

            i+=1

        return bs

    def select_bin(self, data, binnum):
        """

        Although not used by bin(), this is useful for other programs such as
        the random hist matching and correction code

        """
        if binnum > self['nbin']:
            raise ValueError("bin number must be in [0,%d]" % (self['nbin']-1,))

        nlow, nhigh = self.bin_ranges()
        logic = get_range_logic(data, 'ngals_r200', [nlow[binnum], nhigh[binnum]], self.range_type)
        return where1(logic)


    def set_bin_ranges(self):
        if self['nbin'] == 12:
            self.lowlim =  numpy.array([3, 4, 5, 6, 7, 8,  9, 12, 18,  26, 41,  71], dtype='i8')
            self.highlim = numpy.array([3, 4, 5, 6, 7, 8, 11, 17, 25,  40, 70, 220], dtype='i8')
        else:
            raise ValueError("Unsupported nbin: %d\n", self['nbin'])

    def bin_ranges(self, binnum=None):
        if binnum is not None:
            if (binnum < 0) or (binnum > (self['nbin']-1)):
                raise ValueError("binnum out of bounds: [%d,%d]" % (0,self['nbin']-1))
            return self.lowlim[binnum], self.highlim[binnum]

        return self.lowlim, self.highlim


    def plot_dsig_byrun(self, run, show=False):
        """

        Run is a run of the lensing code.
        We must generalize this
        The different z bins are overplotted

        """

        name = self.get_name()
        data=lensing.files.sample_read(type='binned',sample=run,name=name)

        biggles.configure('screen','width', 1140)
        biggles.configure('screen','height', 1140)
        biggles.configure('fontsize_min',1.0)

        if self['nbin'] == 12:
            nrow = 3
            ncol = 4
        else:
            raise ValueError("Unsupported nbin: %s" % self['nmass'])

        pa = FramedArray(nrow, ncol)
        pa.aspect_ratio = 1.0/1.5

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'

        xrnge = [0.01,60.0]
        yrnge = [1.e-2,8000]
        pa.xrange = xrnge
        pa.yrange = yrnge
        pa.xlog = True
        pa.ylog = True


        row = -1

        Nlow, Nhigh = self.bin_ranges()

        i = 0
        for i in xrange(self['nbin']):
            col = i % ncol
            if col == 0:
                row += 1

            Nrange = Nlow[i], Nhigh[i]

            eu.plotting.bscatter(data['r'][i],
                                 data['dsig'][i],
                                 yerr=data['dsigerr'][i],
                                 xrange=xrnge,
                                 yrange=yrnge,
                                 xlog=True,ylog=True,
                                 show=False, 
                                 plt=pa[row,col])
            label = self.bin_label(i)
            pl = PlotLabel(.85, .85, label, halign='right')

            pa[row,col].add(pl)

        if show:
            pa.show()

        d = lensing.files.sample_dir(type='binned',sample=run,name=name)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        epsfile=lensing.files.sample_file(type='binned-plots',sample=run,name=name,extra='-allplot',ext='eps')
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)


    def bin_label(self, binnum):
        Nrange = self.bin_ranges(binnum)
        if Nrange[0] == Nrange[1]:
            label = r'$N_{200}$ = %d' % Nrange[0]
        else:
            label = r'%d $\le N_{200} \le$ %d' % tuple(Nrange)
        return label


class MZBinner(dict):
    """
    Mass and z binning
    """
    def __init__(self, nmass, nz):
        self['nmass'] = nmass
        self['nz']    = nz
    def get_name(self):
        return 'm%02iz%01i' % (self['nmass'],self['nz'])

    def bin_byrun(self, run):
        """
        Do the binning and write out a file
        """

        name=self.get_name()
        d = lensing.files.sample_read(type='collated',sample=run)
        res = self.bin(d)

        lensing.files.sample_write(data=res,
                                   type='binned',
                                   sample=run,
                                   name=name)

    def bin(self, data):

        mlow,mhigh = self.m200bins()
        zlow,zhigh = self.zbins()

        nrbin = data['rsum'][0].size
        dt = lensbin_dtype(nrbin, ['m200','r200','z'])

        bs = numpy.zeros(self['nmass']*self['nz'], dtype=dt)

        i=0
        for ml,mh in zip(mlow,mhigh):
            mrange = [ml,mh]

            ll = (log10(mrange[0]),log10(mrange[1]))
            print('%0.2f < log M200 < %0.2f' % ll)
            for zl,zh in zip(zlow,zhigh):
                zrange = (zl,zh)
                print('    %0.2f < z < %0.2f' % zrange)
                comb,w = reduce_from_ranges(data,'m200',mrange,tag2='z', range2=zrange,
                                            getind=True)
        
                bs['r'][i] = comb['r']
                bs['dsig'][i] = comb['dsig']
                bs['dsigerr'][i] = comb['dsigerr']
                bs['osig'][i] = comb['osig']
                bs['npair'][i] = comb['npair']

                mn,err,sdev = lens_wmom(data,'z',ind=w, sdev=True)
                bs['z_mean'][i] = mn
                bs['z_err'][i] = err
                bs['z_sdev'][i] = sdev
                bs['z_range'][i] = zrange

                mn,err,sdev = lens_wmom(data,'m200',ind=w, sdev=True)
                bs['m200_mean'][i] = mn
                bs['m200_err'][i] = err
                bs['m200_sdev'][i] = sdev
                bs['m200_range'][i] = mrange

                mn,err,sdev = lens_wmom(data,'r200',ind=w, sdev=True)
                bs['r200_mean'][i] = mn
                bs['r200_err'][i] = err
                bs['r200_sdev'][i] = sdev
                bs['r200_range'][i] = mrange



                i+=1

        return bs


    def plot_dsig_byrun(self, run, show=False):
        """

        Run is a run of the lensing code.
        We must generalize this
        The different z bins are overplotted

        """

        name = self.get_name()
        data=lensing.files.lensbin_read(run,name)

        biggles.configure('screen','width', 1140)
        biggles.configure('screen','height', 1140)
        biggles.configure('fontsize_min',1.0)

        if self['nmass'] == 12:
            nrow = 3
            ncol = 4
        else:
            raise ValueError("Unsupported nmass: %s" % self['nmass'])

        pa = FramedArray(nrow, ncol)
        #pa.aspect_ratio = 1.0/1.61803399
        pa.aspect_ratio = 1.0/1.5

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'
        pa.xrange = [0.01,60.0]
        pa.yrange = [1.e-2,8000]
        pa.xlog = True
        pa.ylog = True


        fmt = "%-2i %-2i %-2i %8i %11.5f %11.5f\n"
        row = -1
        colors = ['blue','magenta','red']

        i = 0
        for mi in xrange(self['nmass']):
            col = i % ncol
            if col == 0:
                row += 1

            mrange = data['m200_range'][i]

            pplist = []
            for zi in xrange(self['nz']):
                zr = data['z_range'][i]
                pp = lensing.plotting.add_to_log_plot(pa[row,col],
                                                      data['r'][i],
                                                      data['dsig'][i],
                                                      data['dsigerr'][i],
                                                      color=colors[zi])
                #pp.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
                pp['p'].label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
                pplist.append(pp['p'])
                i+=1

            if mi == 0:
                key = PlotKey(0.1,0.2,pplist, fontsize=0.1)
                pa[row,col].add(key)


            ll = (log10(mrange[0]),log10(mrange[1]))
            lab = r'$%0.2f < logM_{200} < %0.2f$' % ll
            pl = PlotLabel(.5, .85, lab)
            pa[row,col].add(pl)


        if show:
            pa.show()

        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        epsfile = path_join(d, 'lensbin-%s-%s-allplot.eps' % (run,name))
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)


    def plot_nfwmass_byrun(self,run, 
                           withlin=True,
                           residual=False,
                           noerr=False,
                           fudge=None,tag='m200',
                           yrange=None,
                           show=False):
        """

        M200 fit vs M200 true, in z bins

        You have to create these files first

        """
        if withlin:
            ex='lin'
            nex='-lin'
        else:
            nex=''
            ex=None

        name=self.get_name()
        d = lensing.files.sample_read(type='fit', sample=run, name=name, extra=ex)

        plt = FramedPlot()
        colors = ['blue','magenta','red']
        i=0

        if fudge is not None:
            d[tag+'_fit'] *= fudge
            d[tag+'_fit_err'] *= fudge
        for mi in xrange(self['nmass']):
            plist=[]
            for zi in xrange(self['nz']):
                zr = d['z_range'][i]

                mt = d['m200_mean'][i]
                mterr = d['m200_err'][i]
                mfit = d[tag+'_fit'][i]
                mfiterr = d[tag+'_fit_err'][i]

                if residual:
                    if True:
                        yp = (mt - mfit)/mt
                        yperr = sqrt( (mterr/mt)**2 + (mfiterr/mfit)**2 )
                    else:
                        yp = mt - mfit
                        yperr = sqrt(mterr**2 + mfiterr**2)
                else:
                    yp = mfit 
                    yperr = mfiterr

                p=Points([mt],[yp],type='filled circle',color=colors[zi])
                p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
                plist.append(p)

                plt.add(p)
                if not noerr:
                    plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
                    plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

                i+=1

            if mi == 0:
                key = PlotKey(0.1,0.9,plist)
                plt.add(key)

        plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
        if not residual:
            plt.yrange = [0.5*d[tag+'_fit'].min(), 1.5*d[tag+'_fit'].max()]
            plt.ylabel = r'$M_{200}^{fit} [h^{-1} M_{\odot}]$'
            plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
        else:
            #plt.ylabel = r'$M_{true}-M_{fit} [h^{-1} M_{\odot}]$'
            #plt.ylabel = r'$(M_{200}^{true}-M_{200}^{fit})/M_{200}^{true}$'
            if yrange is not None:
                plt.yrange=yrange
            plt.ylabel = '(true-fit)/true'
            plt.add(Curve([0.1,1.e15],[0.0,0.0]))
        plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{\odot}]$'
        plt.xlog = True
        if not residual:
            plt.ylog = True

        plt.aspect_ratio=1

        if not withlin:
            plt.add(PlotLabel(0.7,0.1,"no linear bias"))

        if show:
            plt.show()

        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'fit-m200-residual-vs-true-%s-%s%s.eps' % (run,name,nex)
        else:
            f = 'fit-m200-vs-true-%s-%s%s.eps' % (run,name,nex)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)


    def plot_invmass_byrun(self, run, type='',
                           residual=False,
                           yrange=None,
                           show=False):

        """

        Inverted M200 vs M200 true, in z bins. You must create the inverted
        files first of course

        Inputs
        ------
        run: string
            run identifier
        type: string, optional
            The type of mass to use, either 
                ''    for mean of in and out
                'in'  for mass determined from points inside radius
                'out' for mass determined from points outside radius
        residual: boolean, optional
            Plot the residual from true value.
        """
        name=self.get_name()
        d = lensing.files.sample_read(type='invert',sample=run,name=name)

        plt = FramedPlot()
        colors = ['blue','magenta','red']
        i=0

        tag = 'm200'+type+'_inv'
        for mi in xrange(self['nmass']):
            plist=[]
            for zi in xrange(self['nz']):
                zr = d['z_range'][i]

                mt = d['m200_mean'][i]
                mterr = d['m200_err'][i]
                minv = d[tag][i]
                minverr = d[tag+'_err'][i]

                if residual:
                    if True:
                        yp = (mt - minv)/mt
                        yperr = sqrt( (mterr/mt)**2 + (minverr/minv)**2 )
                    else:
                        yp = mt - minv
                        yperr = sqrt(mterr**2 + minverr**2)
                else:
                    yp = minv 
                    yperr = minverr

                p=Points([mt],[yp],type='filled circle',color=colors[zi])
                p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
                plist.append(p)

                plt.add(p)
                plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
                plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

                i+=1

            if mi == 0:
                key = PlotKey(0.1,0.9,plist)
                plt.add(key)

        plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
        if not residual:
            plt.yrange = [0.5*d[tag].min(), 1.5*d[tag].max()]
            plt.ylabel = r'$M_{200}^{inv} [h^{-1} M_{\dot}]$'
            plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
        else:
            if yrange is not None:
                plt.yrange=yrange
            plt.ylabel = '(true-inv)/true'
            plt.add(Curve([0.1,1.e15],[0.0,0.0]))
        plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{\odot}]$'
        plt.xlog = True
        if not residual:
            plt.ylog = True

        plt.aspect_ratio=1

        if show:
            plt.show()

        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'invert-m200'+type+'-residual-vs-true-%s-%s.eps' % (run,name)
        else:
            f = 'invert-m200'+type+'-vs-true-%s-%s.eps' % (run,name)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)




    def zbins(self):
        if self['nz'] == 3:
            low  = [0.0,  0.34, 0.53]
            high = [0.34, 0.53, 1.4]
        elif self['nz'] == 2:
            low  = [0.0, 0.44]
            high = [0.44, 1.4]
        else:
            raise ValueError("Unsupported nz: %s\n" % self['nz'])
        
        return low,high


    def m200bins(self):
        if self['nmass'] == 12:
            # I ran logbin_linear_edges with self['nmass']=14 and combined last three
            low = [5.03e+12,
                   7.28568896657e+12,
                   1.05369194771e+13,
                   1.52390079477e+13,
                   2.20393981122e+13,
                   3.18744547422e+13,
                   4.60983943365e+13,
                   6.66697509837e+13,
                   9.64210524077e+13,
                   1.39448838645e+14,
                   2.01677726116e+14,
                   2.91676184662e+14]

            high = [7.28568896657e+12,
                    1.05369194771e+13,
                    1.52390079477e+13,
                    2.20393981122e+13,
                    3.18744547422e+13,
                    4.60983943365e+13,
                    6.66697509837e+13,
                    9.64210524077e+13,
                    1.39448838645e+14,
                    2.01677726116e+14,
                    2.91676184662e+14,
                    8.83e+14]
        else:
            raise ValueError("Unsupported nmass: %s\n" % self['nmass'])

        low=numpy.array(low,dtype='f8')
        high=numpy.array(high,dtype='f8')

        return low,high


#
# Old procedural Mass and z binning
#

def mzbin_name(nmass,nz):
    name='m%02iz%01i' % (nmass,nz)
    return name

def zbins(nbin):
    if nbin == 3:
        low  = [0.0,  0.34, 0.53]
        high = [0.34, 0.53, 1.4]
    elif nbin == 2:
        low  = [0.0, 0.44]
        high = [0.44, 1.4]
    else:
        raise ValueError("Unsupported nbin: %s\n" % nbin)
    
    return low,high


def m200bins(nbin):
    if nbin == 12:
        # I ran logbin_linear_edges with nbin=14 and combined last three
        low = [5.03e+12,
               7.28568896657e+12,
               1.05369194771e+13,
               1.52390079477e+13,
               2.20393981122e+13,
               3.18744547422e+13,
               4.60983943365e+13,
               6.66697509837e+13,
               9.64210524077e+13,
               1.39448838645e+14,
               2.01677726116e+14,
               2.91676184662e+14]

        high = [7.28568896657e+12,
                1.05369194771e+13,
                1.52390079477e+13,
                2.20393981122e+13,
                3.18744547422e+13,
                4.60983943365e+13,
                6.66697509837e+13,
                9.64210524077e+13,
                1.39448838645e+14,
                2.01677726116e+14,
                2.91676184662e+14,
                8.83e+14]
    else:
        raise ValueError("Unsupported nbin: %s\n" % nbin)

    low=numpy.array(low,dtype='f8')
    high=numpy.array(high,dtype='f8')

    return low,high

def mzbin(data, nmass, nz):

    mlow,mhigh = m200bins(nmass)
    zlow,zhigh = zbins(nz)

    nrbin = data['rsum'][0].size
    dt = lensbin_dtype(nrbin, ['m200','r200','z'])

    bs = numpy.zeros(nmass*nz, dtype=dt)

    i=0
    for ml,mh in zip(mlow,mhigh):
        mrange = [ml,mh]

        ll = (log10(mrange[0]),log10(mrange[1]))
        print('%0.2f < log M200 < %0.2f' % ll)
        for zl,zh in zip(zlow,zhigh):
            zrange = (zl,zh)
            print('    %0.2f < z < %0.2f' % zrange)
            comb,w = combine_mzbin_lensum_from_ranges(data,'m200',mrange,zrange=zrange,
                                                      getind=True)
    
            bs['r'][i] = comb['r']
            bs['dsig'][i] = comb['dsig']
            bs['dsigerr'][i] = comb['dsigerr']
            bs['osig'][i] = comb['osig']
            bs['npair'][i] = comb['npair']

            mn,err,sdev = lens_wmom(data,'z',ind=w, sdev=True)
            bs['z_mean'][i] = mn
            bs['z_err'][i] = err
            bs['z_sdev'][i] = sdev
            bs['z_range'][i] = zrange

            mn,err,sdev = lens_wmom(data,'m200',ind=w, sdev=True)
            bs['m200_mean'][i] = mn
            bs['m200_err'][i] = err
            bs['m200_sdev'][i] = sdev
            bs['m200_range'][i] = mrange

            mn,err,sdev = lens_wmom(data,'r200',ind=w, sdev=True)
            bs['r200_mean'][i] = mn
            bs['r200_err'][i] = err
            bs['r200_sdev'][i] = sdev
            bs['r200_range'][i] = mrange



            i+=1

    return bs

def plot_mzbin_byrun(run, nmass, nz, show=False):
    """

    Run is a run of the lensing code.
    We must generalize this
    The different z bins are overplotted

    """

    name = 'm%02iz%01i' % (nmass,nz)
    data=lensing.files.lensbin_read(run,name)

    biggles.configure('screen','width', 1140)
    biggles.configure('screen','height', 1140)
    biggles.configure('fontsize_min',1.0)

    if nmass == 12:
        nrow = 3
        ncol = 4
    else:
        raise ValueError("Unsupported nmass: %s" % nmass)

    pa = FramedArray(nrow, ncol)
    #pa.aspect_ratio = 1.0/1.61803399
    pa.aspect_ratio = 1.0/1.5

    pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
    pa.ylabel = r'$\Delta\Sigma ~ [M_{\odot} pc^{-2}]$'
    pa.xrange = [0.01,60.0]
    pa.yrange = [1.e-2,8000]
    pa.xlog = True
    pa.ylog = True


    fmt = "%-2i %-2i %-2i %8i %11.5f %11.5f\n"
    row = -1
    colors = ['blue','magenta','red']

    i = 0
    for mi in xrange(nmass):
        col = i % ncol
        if col == 0:
            row += 1

        mrange = data['m200_range'][i]

        pplist = []
        for zi in xrange(nz):
            zr = data['z_range'][i]
            pp = lensing.plotting.add_to_log_plot(pa[row,col],
                                                   data['r'][i],
                                                   data['dsig'][i],
                                                   data['dsigerr'][i],
                                                   color=colors[zi])
            pp.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            pplist.append(pp)
            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.2,pplist, fontsize=0.1)
            pa[row,col].add(key)


        ll = (log10(mrange[0]),log10(mrange[1]))
        lab = r'$%0.2f < logM_{200} < %0.2f$' % ll
        pl = PlotLabel(.5, .85, lab)
        pa[row,col].add(pl)


    if show:
        pa.show()

    d = lensing.files.lensbin_plot_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    epsfile = path_join(d, 'lensbin-%s-%s-allplot.eps' % (run,name))
    stdout.write("Plotting to file: %s\n" % epsfile)
    pa.write_eps(epsfile)

def plot_mzbin_invmass_byrun(run, nmass, nz, type='',
                             residual=False,
                             yrange=None,
                             show=False):

    """
    Inverted M200 vs M200 true, in z bins
    """
    name=mzbin_name(nmass,nz)
    conf = lensing.files.read_config(run)
    d = lensing.files.sample_read(type='invert',sample=run,name=name)

    plt = FramedPlot()
    colors = ['blue','magenta','red']
    i=0

    tag = 'm200'+type+'_inv'
    for mi in xrange(nmass):
        plist=[]
        for zi in xrange(nz):
            zr = d['z_range'][i]

            mt = d['m200_mean'][i]
            mterr = d['m200_err'][i]
            minv = d[tag][i]
            minverr = d[tag+'_err'][i]

            if residual:
                if True:
                    yp = (mt - minv)/mt
                    yperr = sqrt( (mterr/mt)**2 + (minverr/minv)**2 )
                else:
                    yp = mt - minv
                    yperr = sqrt(mterr**2 + minverr**2)
            else:
                yp = minv 
                yperr = minverr

            p=Points([mt],[yp],type='filled circle',color=colors[zi])
            p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            plist.append(p)

            plt.add(p)
            plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
            plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.9,plist)
            plt.add(key)

    plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
    if not residual:
        plt.yrange = [0.5*d[tag].min(), 1.5*d[tag].max()]
        plt.ylabel = r'$M_{200}^{inv} [h^{-1} M_{\odot}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-inv)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{\odot}]$'
    plt.xlog = True
    if not residual:
        plt.ylog = True

    plt.aspect_ratio=1

    if show:
        plt.show()

    d = lensing.files.lensbin_plot_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    if residual:
        f = 'invert-m200'+type+'-residual-vs-true-%s-%s.eps' % (run,name)
    else:
        f = 'invert-m200'+type+'-vs-true-%s-%s.eps' % (run,name)
    epsfile = path_join(d, f)
    stdout.write("Plotting to file: %s\n" % epsfile)
    plt.write_eps(epsfile)



def plot_mzbin_mass_byrun(run, nmass, nz, 
                          withlin=True,
                          residual=False,
                          noerr=False,
                          fudge=None,tag='m200',
                          yrange=None,
                          show=False):
    """

    M200 fit vs M200 true, in z bins

    """
    if withlin:
        ex='lin'
        nex='-lin'
    else:
        nex=''
        ex=None

    name=mzbin_name(nmass,nz)
    conf = lensing.files.read_config(run)
    d = lensing.files.sample_read(type='fit', sample=run, name=name, extra=ex)

    plt = FramedPlot()
    colors = ['blue','magenta','red']
    i=0

    if fudge is not None:
        d[tag+'_fit'] *= fudge
        d[tag+'_fit_err'] *= fudge
    for mi in xrange(nmass):
        plist=[]
        for zi in xrange(nz):
            zr = d['z_range'][i]

            mt = d['m200_mean'][i]
            mterr = d['m200_err'][i]
            mfit = d[tag+'_fit'][i]
            mfiterr = d[tag+'_fit_err'][i]

            if residual:
                if True:
                    yp = (mt - mfit)/mt
                    yperr = sqrt( (mterr/mt)**2 + (mfiterr/mfit)**2 )
                else:
                    yp = mt - mfit
                    yperr = sqrt(mterr**2 + mfiterr**2)
            else:
                yp = mfit 
                yperr = mfiterr

            p=Points([mt],[yp],type='filled circle',color=colors[zi])
            p.label = '%0.2f < z < %0.2f' % (zr[0],zr[1])
            plist.append(p)

            plt.add(p)
            if not noerr:
                plt.add(SymmetricErrorBarsX([mt],[yp],[mterr],color=colors[zi]))
                plt.add(SymmetricErrorBarsY([mt],[yp],[yperr],color=colors[zi]))

            i+=1

        if mi == 0:
            key = PlotKey(0.1,0.9,plist)
            plt.add(key)

    plt.xrange = [0.5*d['m200_mean'].min(), 1.5*d['m200_mean'].max()]
    if not residual:
        plt.yrange = [0.5*d[tag+'_fit'].min(), 1.5*d[tag+'_fit'].max()]
        plt.ylabel = r'$M_{200}^{fit} [h^{-1} M_{\odot}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        #plt.ylabel = r'$M_{true}-M_{fit} [h^{-1} M_{\odot}]$'
        #plt.ylabel = r'$(M_{200}^{true}-M_{200}^{fit})/M_{200}^{true}$'
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-fit)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{\odot}]$'
    plt.xlog = True
    if not residual:
        plt.ylog = True

    plt.aspect_ratio=1

    if not withlin:
        plt.add(PlotLabel(0.7,0.1,"no linear bias"))

    if show:
        plt.show()

    d = lensing.files.lensbin_plot_dir(run,name)
    if not os.path.exists(d):
        os.makedirs(d)
    if residual:
        f = 'fit-m200-residual-vs-true-%s-%s%s.eps' % (run,name,nex)
    else:
        f = 'fit-m200-vs-true-%s-%s%s.eps' % (run,name,nex)
    epsfile = path_join(d, f)
    stdout.write("Plotting to file: %s\n" % epsfile)
    plt.write_eps(epsfile)

def combine_mzbin_lensum_from_ranges(data, tag, trange, zrange=None,getind=False):
    logic = (data[tag] >= trange[0]) & (data[tag] < trange[1])
    if zrange is not None:
        logic = logic & \
            (data['z'] >= zrange[0]) & (data['z'] < zrange[1])
    w=where1(logic)

    comb = average_lensums(data[w])

    if getind:
        return comb, w
    else:
        return comb
 



def reduce_from_ranges(data, 
                       tag, range1,
                       range_type = '[)', 
                       tag2=None,
                       range2=None,
                       range2_type = '[)',
                       getind=False):
    """

    Range 1 and range 2 are anded.

    Can add more ranges if needed

    """
    #logic = (data[tag] >= range1[0]) & (data[tag] < range1[1])

    logic = get_range_logic(data, tag, range1, range_type)

    if range2 is not None and tag2 is not None:
        logic = logic & get_range_logic(data, tag2, range2, range2_type)

    w=where1(logic)

    comb = lensing.outputs.average_lensums(data[w])

    if getind:
        return comb, w
    else:
        return comb
 
def reduce_from_ranges_many(data, tags_and_ranges, getind=False):
    """

    parameters
    ----------
    data: ndarray with fields
        The lensum data
    tags_and_ranges: list
        each element in the list is of the form
            (tagname, range, range_type)
    getind: bool
        If True, get the indices from data that are in range
    """

    logic = numpy.ones(data.size, dtype='bool')
    for tag_range in tags_and_ranges:
        print("    ",tag_range)
        tag, range, range_type = tag_range
        logic = logic & get_range_logic(data, tag, range, range_type)

    w=where1(logic)

    comb = lensing.outputs.average_lensums(data[w])

    if getind:
        return comb, w
    else:
        return comb
 
def get_range_logic(data, tag, brange, type):
    minval,maxval = brange
    if minval is None:
        minval = data[tag].min()
    if maxval is None:
        maxval = data[tag].max()

    print(minval,maxval)
    if type == '[]':
        logic = (data[tag] >= minval) & (data[tag] <= maxval)
    elif type == '[)':
        logic = (data[tag] >= minval) & (data[tag] < maxval)
    elif type == '(]':
        logic = (data[tag] > minval) & (data[tag] <= maxval)
    elif type == '()':
        logic = (data[tag] > minval) & (data[tag] < maxval)
    else:
        raise ValueError("Bad range type: '%s'" % type)

    return logic

def logbin_linear_edges_struct(data, nbin, tag):
    """

    Do a log binning and print out the corresponding bin edges in
    the linear variable.  
    
    For binning by mass estimators, probably want to combine the last 2 or 3
    bins because of the exponential cutoff
    
    """

    low, high, rev=logbin_linear_edges(data[tag], nbin, getrev=True) 

    mean_data = numpy.zeros(nbin,dtype='f8')
    h=numpy.zeros(nbin,dtype='i8')
    print("Getting mean of data")
    for i in xrange(nbin):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]
            mean_data[i],err = lens_wmom(data, tag, ind=w)
            h[i] = w.size
        else:
            mean_data[i] = -9999

    low=10.0**hdict['low']
    high=10.0**hdict['high']

    fmt = "%-2i %8i %11.5e %11.5e %11.5e\n"
    for i in xrange(nbin):
        stdout.write(fmt % (i+1,h[i],low[i],high[i],mean_data[i]))

    print(eu.numpy_util.arr2str(low))
    print(eu.numpy_util.arr2str(high))
 
def logbin_linear_edges(data, nbin, getrev=False):
    """

    Do a log binning and print out the corresponding bin edges in
    the linear variable.  
    
    For binning by mass estimators, probably want to combine the last 2 or 3
    bins because of the exponential cutoff
    
    """

    #log_data = numpy.log10(data[tag])
    log_data = numpy.log10(data)

    # we can't just take the returned means because we want
    # the mean in the linear value of the parameter
    hdict = histogram(log_data, nbin=nbin, more=True, rev=True)

    low=10**hdict['low']
    high=10**hdict['high']

    if getrev:
        return low, high, hdict['rev']
    else:
        return low, high


    
def lensbin_struct(nrbin, bintags=None, n=1):
    dt = lensbin_dtype(nrbin, bintags=bintags)
    return numpy.zeros(n, dtype=dt)

def lensbin_dtype(nrbin, bintags=None):
    """
    This is the same as lensing.outputs.averaged_dtype
    but with the averages added
    """
    dt=[('nlenses','i8')]
    if bintags is not None:
        if not isinstance(bintags,list):
            bintags = [bintags]

        for bt in bintags:
            tn = bt+'_range'
            dt.append( (tn,'f8',2) )

            tn = bt+'_minmax'
            dt.append( (tn,'f8',2) )

            tn = bt+'_mean'
            dt.append( (tn,'f8') )

            tn = bt+'_err'
            dt.append( (tn,'f8') )

            tn = bt+'_sdev'
            dt.append( (tn,'f8') )

    dt += lensing.outputs.averaged_dtype(nrbin)

    return numpy.dtype(dt)

