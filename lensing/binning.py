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
from sys import stdout
import numpy
from numpy import log10,sqrt,linspace,where
import esutil as eu
from esutil.ostools import path_join
from esutil.stat import histogram
from esutil.numpy_util import where1

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

def instantiate_binner(type, nbin):
    if type == 'lambda':
        b = LambdaBinner(nbin)
    elif type == 'n200':
        b = N200Binner(nbin)
    elif type == 'mz':
        b = MZBinner(nbin)
    else:
        raise ValueError("unsupported binner type: '%s'" % type)
    return b


def bin_lenses_byrun(run, type, nbin):
    b=instantiate_binner(type, nbin)
    b.bin_byrun(run)

class BinnerBase(dict):
    def __init__(self, nbin):
        self['nbin'] = nbin
        self.set_bin_ranges()

    def name(self):
        raise RuntimeError("override this method")

    def bin_byrun(self, run):
        """
        Do the binning and write out a file
        """

        name=self.name()
        d = lensing.files.sample_read('collated',run, fs='hdfs')
        res = self.bin(d)

        outdir = lensing.files.sample_dir('binned',run,name=name)
        if not os.path.exists(outdir):
            print("Making output dir:",outdir)
            os.makedirs(outdir)
        lensing.files.sample_write(res,'binned',run,name=name,clobber=True)

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

    def bin_label(binnum):
        raise RuntimeError("override this method")

    def plot_dsig_osig_byrun_bin(self, run, type, binnum, **keys):
        """

        Plot a single bin, with dsig and osig.

        See lensing.plotting.plot2dsig
        """
        name = self.name()
        data=lensing.files.sample_read(type,run,name=name)

        max_binnum = self['nbin']-1
        if (binnum < 0) or (binnum > max_binnum):
            raise ValueError("binnum is out of range [0,%d]" % max_binnum)
        data = data[binnum]

        keys['plot_label'] = self.bin_label(binnum)
        lensing.plotting.plot2dsig(data['r'], 
                                   data['dsig'], data['dsigerr'],
                                   data['osig'], data['dsigerr'],
                                   **keys)

    def plot_dsig_osig_byrun(self, run, type, **keys):
        """

        Plot all bins dsig+osig but one at a time, not in
        an array

        """
        for binnum in xrange(self['nbin']):
            self.plot_dsig_osig_byrun_bin(run, type, binnum, **keys)

    def plot_dsig_byrun_1var(self, run, type, show=False):
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

        name = self.name()
        data=lensing.files.sample_read(type,run,name=name)

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
        else:
            raise ValueError("Unsupported nbin: %s" % self['nbin'])

        pa = FramedArray(nrow, ncol)
        pa.aspect_ratio = aspect_ratio

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

        xrnge = [0.01,60.0]
        yrnge = [1.e-2,8000]
        pa.xrange = xrnge
        pa.yrange = yrnge
        pa.xlog = True
        pa.ylog = True


        row = -1

        low, high = self.bin_ranges()

        i = 0
        for i in xrange(self['nbin']):
            col = i % ncol
            if col == 0:
                row += 1

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

        epsfile=lensing.files.sample_file(type+'-plots',run,name=name,extra='allplot',ext='eps')
        d = os.path.dirname(epsfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
        converter.convert(epsfile, dpi=120, verbose=True)


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

        name = self.name()
        data=lensing.files.sample_read(type,run,name=name)

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
        else:
            raise ValueError("Unsupported nbin: %s" % self['nbin'])

        pa = FramedArray(nrow, ncol)
        pa.aspect_ratio = aspect_ratio

        pa.xlabel = r'$r$ [$h^{-1}$ Mpc]'
        pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

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
                                 plt=pa[row,col])
            label = self.bin_label(i)
            pl = PlotLabel(.85, .85, label, halign='right')

            pa[row,col].add(pl)

            line=biggles.Curve([xrnge[0],xrnge[1]],[0,0])
            pa[row,col].add(line)


        if show:
            pa.show()

        epsfile=lensing.files.sample_file(type+'-plots',run,name=name,extra='osig-allplot',ext='eps')
        d = os.path.dirname(epsfile)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
        converter.convert(epsfile, dpi=120, verbose=True)




def define_lambda_bins(sample, lastmin=58.):
    l=lensing.files.read_original_catalog('lens',sample)
    w=where1(l['z_lambda'] < 0.4)
    l=l[w]

    # 12 is unimportant here
    lb = instantiate_binner('lambda', 12)
    lb.define_bins(l['lambda_chisq'], lastmin=lastmin)

class LambdaBinner(BinnerBase):
    """
    RedMapper binner on lambda
    """
    range_type='()'

    lambda_field = 'lambda_chisq'
    z_field = 'z_lambda'
    def name(self):
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
            comb,w = reduce_from_ranges(data,
                                        lambda_field,
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

            # ran logbin_linear_edges with nbin=14 and combined last three
            # also made upper exactly 189
            """
            lowlim = [  10.00009882,   12.52914267,   15.69778647,   19.66778627,
                      24.64180651,   30.87376584,   38.68179945,   48.46449947,
                      60.72126278,   76.07778464,   95.31799984,  119.42410174]

            highlim = [  12.52914267,   15.69778647,   19.66778627,   24.64180651,
                       30.87376584,   38.68179945,   48.46449947,   60.72126278,
                       76.07778464,   95.31799984,  119.42410174,  234.87844989]
            """

        elif self['nbin'] == 16:
            # ran logbin_linear_edges with nbin=18 and combined last three
            # also made upper exactly 189

            lowlim  = [10.00, 11.77, 13.86, 16.31, 19.20, 22.59, 26.60, 31.30, 36.85, 43.37, 51.05, 60.09, 70.73, 83.26, 98.00, 115.35]
            highlim = [11.77, 13.86, 16.31, 19.20, 22.59, 26.60, 31.30, 36.85, 43.37, 51.05, 60.09, 70.73, 83.26, 98.00, 115.35, 189.00]

        else:
            raise ValueError("Unsupported nbin: %d\n", self['nbin'])

        #self.lowlim = numpy.array(lowlim,dtype='f8')
        #self.highlim = numpy.array(highlim,dtype='f8')
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
    def name(self):
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

        name = self.name()
        data=lensing.files.sample_read('binned',run,name=name)

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
        pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'

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

        d = lensing.files.sample_dir('binned',run,name=name)
        if not os.path.exists(d):
            print("making dir:",d)
            os.makedirs(d)
        epsfile=lensing.files.sample_file('binned-plots',run,name=name,extra='-allplot',ext='eps')
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
    def name(self):
        return 'm%02iz%01i' % (self['nmass'],self['nz'])

    def bin_byrun(self, run):
        """
        Do the binning and write out a file
        """

        name=self.name()
        d = lensing.files.sample_read('collated',run)
        res = self.bin(d)

        lensing.files.sample_write(res,'binned',run,name=name,clobber=True)

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

        name = self.name()
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
        pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'
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

        name=self.name()
        d = lensing.files.sample_read('fit', run, name=name, extra=ex)

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
            plt.ylabel = r'$M_{200}^{fit} [h^{-1} M_{sun}]$'
            plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
        else:
            #plt.ylabel = r'$M_{true}-M_{fit} [h^{-1} M_{sun}]$'
            #plt.ylabel = r'$(M_{200}^{true}-M_{200}^{fit})/M_{200}^{true}$'
            if yrange is not None:
                plt.yrange=yrange
            plt.ylabel = '(true-fit)/true'
            plt.add(Curve([0.1,1.e15],[0.0,0.0]))
        plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
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
        name=self.name()
        d = lensing.files.sample_read('invert',run,name=name)

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
            plt.ylabel = r'$M_{200}^{inv} [h^{-1} M_{sun}]$'
            plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
        else:
            if yrange is not None:
                plt.yrange=yrange
            plt.ylabel = '(true-inv)/true'
            plt.add(Curve([0.1,1.e15],[0.0,0.0]))
        plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
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
    pa.ylabel = r'$\Delta\Sigma ~ [M_{sun} pc^{-2}]$'
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
    d = lensing.files.sample_read('invert',run,name=name)

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
        plt.ylabel = r'$M_{200}^{inv} [h^{-1} M_{sun}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-inv)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
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
    d = lensing.files.sample_read('fit', run, name=name, extra=ex)

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
        plt.ylabel = r'$M_{200}^{fit} [h^{-1} M_{sun}]$'
        plt.add(Curve([0.1,1.e15],[0.1,1.e15]))
    else:
        #plt.ylabel = r'$M_{true}-M_{fit} [h^{-1} M_{sun}]$'
        #plt.ylabel = r'$(M_{200}^{true}-M_{200}^{fit})/M_{200}^{true}$'
        if yrange is not None:
            plt.yrange=yrange
        plt.ylabel = '(true-fit)/true'
        plt.add(Curve([0.1,1.e15],[0.0,0.0]))
    plt.xlabel = r'$M_{200}^{true} [h^{-1} M_{sun}]$'
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
                       tag, range1, range_type = '[)', 
                       tag2=None, range2=None, range2_type = '[)',
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
 

def get_range_logic(data, tag, brange, type):
    minval,maxval = brange
    if minval is None:
        minval = data[tag].min()
    if maxval is None:
        maxval = data[tag].max()

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
    dt=[]
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

