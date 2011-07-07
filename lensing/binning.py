"""
Binner classes
    MZBinner: bin by true mass and redshift

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

def reduce_from_ranges(data, tag, range1, tag2=None, range2=None, getind=False):
    """
    Can add more ranges if needed
    """
    logic = (data[tag] >= range1[0]) & (data[tag] < range1[1])
    if range2 is not None and tag2 is not None:
        logic = logic & \
            (data[tag2] >= range2[0]) & (data[tag2] < range2[1])

    w=where1(logic)

    comb = lensing.outputs.reduce_lensums(data[w])

    if getind:
        return comb, w
    else:
        return comb
 



def logbin_linear_edges(data, nbin, tag):
    """

    Do a log binning and print out the corresponding bin edges in
    the linear variable.  
    
    For binning by mass estimators, probably want to combine the last 2 or 3
    bins because of the exponential cutoff
    
    """

    log_data = numpy.log10(data[tag])

    # we can't just take the returned means because we want
    # the mean in the linear value of the parameter
    hdict = histogram(log_data, nbin=nbin, more=True, rev=True)

    rev = hdict['rev']
    mean_data = numpy.zeros(nbin,dtype='f8')
    stdout.write("Getting mean of '%s'\n" % tag)
    for i in xrange(nbin):
        if rev[i] != rev[i+1]:
            w=rev[ rev[i]:rev[i+1] ]
            mean_data[i],err = lens_wmom(data, tag, ind=w)

    h = hdict['hist']
    low=10.0**hdict['low']
    high=10.0**hdict['high']
    fmt = "%-2i %8i %11.5e %11.5e %11.5e\n"
    for i in xrange(nbin):
        stdout.write(fmt % (i+1,h[i],low[i],high[i],mean_data[i]))

    print(eu.numpy_util.arr2str(low))
    print(eu.numpy_util.arr2str(high))
    
def lensbin_dtype(nrbin, bintags):
    if not isinstance(bintags,list):
        bintags = [bintags]

    dt = []
    for bt in bintags:
        tn = bt+'_range'
        dt.append( (tn,'f8',2) )
        tn = bt+'_mean'
        dt.append( (tn,'f8') )
        tn = bt+'_err'
        dt.append( (tn,'f8') )
        tn = bt+'_sdev'
        dt.append( (tn,'f8') )

    nrbin = int(nrbin)
    dt += [('r','f8',nrbin),
           ('dsig','f8',nrbin),
           ('dsigerr','f8',nrbin),
           ('osig','f8',nrbin),
           ('npair','i8',nrbin)]
    return numpy.dtype(dt)



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
        d = lensing.files.lensred_collated_read(run)
        res = self.bin(d)
        lensing.files.lensbin_write(res,run,name) 

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


    def plot_dsig_byrun(self, run, dops=False):
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


        if dops:
            d = lensing.files.lensbin_plot_dir(run,name)
            if not os.path.exists(d):
                os.makedirs(d)
            epsfile = path_join(d, 'lensbin-%s-%s-allplot.eps' % (run,name))
            stdout.write("Plotting to file: %s\n" % epsfile)
            pa.write_eps(epsfile)
        else:
            pa.show()

    def plot_nfwmass_byrun(self,run, 
                           withlin=True,
                           residual=False,
                           noerr=False,
                           fudge=None,tag='m200',
                           yrange=None,
                           dops=False):
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
        d = lensing.files.lensfit_read(run,name,extra=ex)

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

        if dops:
            d = lensing.files.lensbin_plot_dir(run,name)
            if not os.path.exists(d):
                os.makedirs(d)
            if residual:
                f = 'lensfit-m200-residual-vs-true-%s-%s%s.eps' % (run,name,nex)
            else:
                f = 'lensfit-m200-vs-true-%s-%s%s.eps' % (run,name,nex)
            epsfile = path_join(d, f)
            stdout.write("Plotting to file: %s\n" % epsfile)
            plt.write_eps(epsfile)

        else:
            plt.show()


    def plot_invmass_byrun(self, run, type='',
                           residual=False,
                           yrange=None,
                           dops=False):

        """
        Inverted M200 vs M200 true, in z bins

        You must create the inverted files first of course
        """
        name=self.name()
        d = lensing.files.lensinv_read(run,name)

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

        if dops:
            d = lensing.files.lensbin_plot_dir(run,name)
            if not os.path.exists(d):
                os.makedirs(d)
            if residual:
                f = 'lensinv-m200'+type+'-residual-vs-true-%s-%s.eps' % (run,name)
            else:
                f = 'lensinv-m200'+type+'-vs-true-%s-%s.eps' % (run,name)
            epsfile = path_join(d, f)
            stdout.write("Plotting to file: %s\n" % epsfile)
            plt.write_eps(epsfile)

        else:
            plt.show()




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

def mzbin_byrun(run, nmass, nz):
    """
    Do the binning and write out a file
    """
    name='m%02iz%01i' % (nmass,nz)
    d = lensing.files.lensout_collate(run)
    res = mzbin(d, nmass, nz)
    lensing.files.lensbin_write(res,run,name) 

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

def plot_mzbin_byrun(run, nmass, nz, dops=False):
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


    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        epsfile = path_join(d, 'lensbin-%s-%s-allplot.eps' % (run,name))
        stdout.write("Plotting to file: %s\n" % epsfile)
        pa.write_eps(epsfile)
    else:
        pa.show()

def plot_mzbin_invmass_byrun(run, nmass, nz, type='',
                             residual=False,
                             yrange=None,
                             dops=False):

    """
    Inverted M200 vs M200 true, in z bins
    """
    name=mzbin_name(nmass,nz)
    conf = lensing.files.read_config(run)
    d = lensing.files.lensinv_read(run,name)

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

    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'lensinv-m200'+type+'-residual-vs-true-%s-%s.eps' % (run,name)
        else:
            f = 'lensinv-m200'+type+'-vs-true-%s-%s.eps' % (run,name)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)

    else:
        plt.show()



def plot_mzbin_mass_byrun(run, nmass, nz, 
                          withlin=True,
                          residual=False,
                          noerr=False,
                          fudge=None,tag='m200',
                          yrange=None,
                          dops=False):
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
    d = lensing.files.lensfit_read(run,name,extra=ex)

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

    if dops:
        d = lensing.files.lensbin_plot_dir(run,name)
        if not os.path.exists(d):
            os.makedirs(d)
        if residual:
            f = 'lensfit-m200-residual-vs-true-%s-%s%s.eps' % (run,name,nex)
        else:
            f = 'lensfit-m200-vs-true-%s-%s%s.eps' % (run,name,nex)
        epsfile = path_join(d, f)
        stdout.write("Plotting to file: %s\n" % epsfile)
        plt.write_eps(epsfile)

    else:
        plt.show()

def combine_mzbin_lensum_from_ranges(data, tag, trange, zrange=None,getind=False):
    logic = (data[tag] >= trange[0]) & (data[tag] < trange[1])
    if zrange is not None:
        logic = logic & \
            (data['z'] >= zrange[0]) & (data['z'] < zrange[1])
    w=where1(logic)

    comb = reduce_lensout(data[w])

    if getind:
        return comb, w
    else:
        return comb
 


