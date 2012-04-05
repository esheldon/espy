"""
Run the tester to see trends of <e> etc. with various parameters.

The most important for gal ellip is <e> vs R.  A polynomial
will be fit for this case and written to a file.

You can then run the Tester.detrend(rmag_max) function to 
create a new column that is detrended.

"""
from __future__ import print_function
import os
import columns
import numpy
from numpy import sqrt
import esutil as eu
from esutil.numpy_util import where1
from esutil.ostools import path_join, expand_path

import biggles
from biggles import FramedPlot, PlotKey, Table, PlotLabel, Points, \
        SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX, \
        Curve,FramedArray

import pprint

import fimage
from fimage.conversions import mom2sigma,mom2fwhm

from . import regauss

import converter

class Tester(dict):
    """
    testing a single run
    """
    def __init__(self, procrun, sweeptype, band, 
                 run='any', camcol='any', coldir=None):
        self.procrun = procrun
        self.sweeptype=sweeptype
        self.run = run
        self.camcol=camcol
        self.band = band
        self.coldir=coldir

        if self.procrun in ['01','02','03']:
            self.old=True
        else:
            self.old=False

    def open_columns(self):
        if self.coldir is not None:
            print("using coldir:",self.coldir)
            c = columns.Columns(self.coldir)
        else:
            c = regauss.open_columns(self.procrun,self.sweeptype)
            print("using coldir:",c.dir)

        return c

    def load_data(self, field):

        # loading basic columns
        if 'e1' not in self or field not in self:
            print("Opening columns")
            c = self.open_columns()

            meta = c['meta'].read()
            if meta['sweeptype'] != self.sweeptype:
                raise ValueError("Requested sweeptype %s but "
                    "columns dir has %s " % (self.sweeptype,meta['sweeptype']))

            if self.run != 'any' and self.camcol != 'any':
                print('Getting indices of run:',self.run,'camcol:',self.camcol)
                w=c['run'].match(self.run) & (c['camcol'] == self.camcol)

                if w.size == 0:
                    raise ValueError("Run",self.run,"not found")
            elif self.run != 'any':
                print('Getting indices of run:',self.run)
                w=(c['run'].match(self.run))

                if w.size == 0:
                    raise ValueError("Run",self.run,"not found")
            elif self.camcol != 'any':

                print('Getting indices of camcol:',self.camcol)
                w=(c['camcol'] == self.camcol)
                if w.size == 0:
                    raise ValueError("camcol",self.camcol,"not found")
            else:
                w=None


            if 'e1' not in self:
                fdict={}
                if self.sweeptype == 'gal':
                    fdict = {'R_rg':'R',
                             'corrflags_rg':'amflags',
                             'e1_rg':'e1',
                             'e2_rg':'e2',
                             'uncer_rg':'uncer'}
                else:
                    fdict = {'amflags':'amflags',
                             'e1':'e1','e2':'e2',
                             'uncer':'uncer'}
                
                for f in ['Irr','Icc','Irr_psf','Icc_psf','e1_psf','e2_psf']:
                    fdict[f] = f

                for name in fdict:
                    uname = fdict[name]

                    cname = name+'_'+self.band
                    print("    loading col",cname,"as",uname)
                    if w is None:
                        self[uname] = c[cname][:]
                    else:
                        self[uname] = c[cname][w]

                Tpsf = self['Irr_psf'] + self['Icc_psf']
                self['sigma_psf'] = mom2sigma(Tpsf)
                self['fwhm_psf'] = mom2fwhm(Tpsf, pixscale=0.396)

                Tobj = self['Irr'] + self['Icc']
                self['sigma'] = mom2sigma(Tobj)

                fdict={'field':'field'}
                if self.sweeptype == 'gal':
                    fdict['cmodelmag_dered_r'] = 'rmag'
                else:
                    fdict['psfmag_dered_r'] = 'rmag'

                for name in fdict:
                    uname = fdict[name]
                    print("    loading col",name,"as",uname)
                    if w is None:
                        self[uname] = c[name][:]
                    else:
                        self[uname] = c[name][w]

            if field not in self:
                print("loading plot field:",field)
                if w is None:
                    self[field] = c[field][:]
                else:
                    self[field] = c[field][w]


    def plot_vs_field(self, field, plot_type, 
                      rmag_min=None, 
                      rmag_max=None, 
                      fmin=None, fmax=None, 
                      nbin=20, nperbin=50000,
                      xrng=None, yrng=None, show=True):

        
        allowed=['meane','residual']
        if plot_type not in allowed:
            raise ValueError("plot_type should be in [%s]" % ','.join(allowed))
        if plot_type == 'residual' and self.sweeptype != 'star':
            raise ValueError("residuals only supported for stars")

        if rmag_min is None:
            if self.sweeptype == 'gal':
                rmag_min=18.0
            else:
                rmag_min=15.0

        if rmag_max is None:
            if self.sweeptype == 'gal':
                rmag_max=21.8
            else:
                rmag_max=19.0

        # this will only load the main data once.
        self.load_data(field)

        print("Using rmag range: [%0.2f,%0.2f]" % (rmag_min,rmag_max))
        # notes 
        #  - amflags here is really corrflags_rg for gals
        #  - e1 is really e1_rg for gals
        logic = ((self['amflags'] == 0)
                 & (self['e1'] < 4)
                 & (self['e1'] > -4)
                 & (self['rmag'] > rmag_min)
                 & (self['rmag'] < rmag_max) )
        
        if self.sweeptype == 'gal':
            logic = logic & (self['R'] > 1.0/3.0) & (self['R'] < 1.0)

        w=where1(logic)
        print("Number passing cuts:",w.size)
        minnum=31000
        if w.size < minnum:
            print("want %d good objects, found %d" % (minnum,w.size))
            return

        weights = 1.0/(0.32**2 + self['uncer'][w]**2)

        # we can try to get nice labels for some fields
        if field == 'fwhm_psf':
            field_data = self['fwhm_psf'][w]
            fstr = 'PSF FWHM (arcsec)'
        elif field == 'sigma_psf':
            field_data = self['sigma_psf'][w]
            fstr = r'\sigma_{PSF}'
        elif field == 'sigma':
            field_data = self['sigma'][w]
            fstr = r'\sigma_{obj+PSF}'
        else:
            field_data = self[field][w]
            fstr=field
            fstr = fstr.replace('_','\_')

        print("Plotting for field:",field)

        if plot_type == 'residual':
            print('  doing: residual')
            be1 = eu.stat.Binner(field_data, self['e1'][w]-self['e1_psf'][w], 
                                 weights=weights)
            be2 = eu.stat.Binner(field_data, self['e2'][w]-self['e2_psf'][w], 
                                 weights=weights)
            ylabel = r'$<e_{star}-e_{PSF}>$'
        else:
            print('  doing meane')
            be1 = eu.stat.Binner(field_data, self['e1'][w], weights=weights)
            be2 = eu.stat.Binner(field_data, self['e2'][w], weights=weights)
            ylabel = r'$<e>$'


        # regular hist for display
        print("  regular fixed binsize hist")
        xm,xe,xstd=eu.stat.wmom(field_data, weights, sdev=True)
        #hxmin = xm-4.0*xstd
        #hxmax = xm+4.0*xstd
        bsize = xstd/5.
        hist = eu.stat.histogram(field_data, binsize=bsize, 
                                 weights=weights, more=True)


        print("  hist  e1, nperbin: ",nperbin)
        be1.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be1.dohist(nbin=nbin, min=fmin, max=fmax)
        print("  stats e1")
        be1.calc_stats()
        print("  hist  e2, nperbin: ",nperbin)
        be2.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be2.dohist(nbin=nbin, min=fmin, max=fmax)
        print("  stats e2")
        be2.calc_stats()



        plt = FramedPlot()

        if xrng is not None:
            plt.xrange=xrng
        else:
            if field == 'R':
                plt.xrange=[0.29,1.01]

        if yrng is not None:
            plt.yrange=yrng
            ymin = yrng[0]
            ymax = 0.8*yrng[1]
        else:
            ymin = min( be1['wymean'].min(),be2['wymean'].min() )
            ymax = 0.8*max( be1['wymean'].max(),be2['wymean'].max() )

        # this is a histogram-like object
        ph = eu.plotting.make_hist_curve(hist['low'], hist['high'], hist['whist'], 
                                         ymin=ymin, ymax=ymax, color='grey50')
        plt.add(ph)

        p1 = Points( be1['wxmean'], be1['wymean'], 
                    type='filled circle', color='blue')
        p1err = SymErrY( be1['wxmean'], be1['wymean'], be1['wyerr2'], color='blue')
        p1.label = r'$e_1$'

        p2 = Points( be2['wxmean'], be2['wymean'], 
                    type='filled circle', color='red')
        p2.label = r'$e_2$'
        p2err = SymErrY( be2['wxmean'], be2['wymean'], be2['wyerr2'], color='red')

        key = PlotKey(0.1,0.9, [p1,p2])
        plt.add(p1, p1err, p2, p2err, key)

        if self.camcol != 'any' and field == 'R' and plot_type=='meane':
            order=3
            print("  getting poly order",order)
            coeff1 = numpy.polyfit(be1['wxmean'], be1['wymean'], order)
            poly1=numpy.poly1d(coeff1)
            coeff2 = numpy.polyfit(be2['wxmean'], be2['wymean'], order)
            poly2=numpy.poly1d(coeff2)

            ps1 = Curve( be1['wxmean'], poly1(be1['wxmean']), color='blue')
            ps2 = Curve( be2['wxmean'], poly2(be2['wxmean']), color='red')
            plt.add(ps1,ps2)

            polyf = self.R_polyfile(self.camcol, rmag_max)
            out={'coeff_e1':list([float(c) for c in coeff1]),
                 'coeff_e2':list([float(c) for c in coeff2])}
            print("    -- Writing poly coeffs to:",polyf)
            eu.io.write(polyf,out)
        if field != 'rmag':
            rmag_lab = \
                PlotLabel(0.1,0.05,'%0.2f < rmag < %0.2f' % (rmag_min,rmag_max), 
                          halign='left')
            plt.add(rmag_lab)

        procrun_lab = PlotLabel(0.1,0.1,
                            'procrun: %s filter: %s' % (self.procrun, self.band), 
                            halign='left')
        plt.add(procrun_lab)
        cy=0.9
        if self.run != 'any':
            run_lab = PlotLabel(0.9,0.9, 'run: %06i' % self.run, halign='right')
            plt.add(run_lab)
            cy=0.8
        if self.camcol != 'any':
            run_lab = PlotLabel(0.9,cy, 'camcol: %i' % self.camcol, halign='right')
            plt.add(run_lab)



        plt.xlabel = r'$'+fstr+'$'
        plt.ylabel = ylabel


        if show:
            plt.show()
        epsfile = self.plotfile(field, rmag_max, plot_type=plot_type)
        eu.ostools.makedirs_fromfile(epsfile, verbose=True)
        print("  Writing eps file:",epsfile)
        plt.write_eps(epsfile)

        converter.convert(epsfile, verbose=True)

    def detrend(self, rmag_min, rmag_max, show=False):
        """
        Use the polynomial fit to <e1,2> vs R to detrend
        the data
        """

        biggles.configure('fontsize_min', 1)
        band=self.band
        c = self.open_columns()

        cols=['R_rg','e1_rg','e2_rg','corrflags_rg', 'e1_psf','e2_psf','uncer_rg']
        cols = [col+'_'+band for col in cols]
        cols += ['camcol','cmodelmag_dered_r']

        print("Reading columns:")
        pprint.pprint(cols)
        data = c.read_columns(cols)
        
        dtflag=numpy.ones(data.size,dtype='u1')
        e1dt=numpy.zeros(data.size,dtype='f8') - 9999.
        e2dt=numpy.zeros(data.size,dtype='f8') - 9999.

        for camcol in [1,2,3,4,5,6]:
            print("camcol:",camcol)
            pfile=self.R_polyfile(camcol, rmag_max)
            print("    ",pfile)
            coeffs=eu.io.read(pfile)
            poly_e1=numpy.poly1d(coeffs['coeff_e1'])
            poly_e2=numpy.poly1d(coeffs['coeff_e2'])

            # this is the same logic as in
            # Tester.plot_vs_field
            print("    getting logic")
            w = where1((data['camcol'] == camcol)
                       & (data['R_rg_'+band] > 1.0/3.0)
                       & (data['R_rg_'+band] < 1.0)
                       & (data['corrflags_rg_'+band] == 0)
                       & (data['e1_rg_'+band] < 4)
                       & (data['e1_rg_'+band] > -4)
                       & (data['cmodelmag_dered_r'] > rmag_min)
                       & (data['cmodelmag_dered_r'] < rmag_max) )
            print("    de-trending")
            Rsub = data['R_rg_'+band][w]
            trend1 = poly_e1(Rsub)
            e1dt[w] = data['e1_rg_'+band][w] - trend1
            del trend1
            trend2 = poly_e2(Rsub)
            e2dt[w] = data['e2_rg_'+band][w] - trend2

            dtflag[w] = 0

            self.plot_detrend(data,e1dt,e2dt,w,camcol,rmag_min,rmag_max,
                              coeffs, show=show)

        
        ls='%0.1f' % rmag_max
        ls = ls.replace('.','')
        e1name='e1_rg_dt'+ls+'_'+band
        e2name='e2_rg_dt'+ls+'_'+band
        dtflagname='dt'+ls+'_flag_'+band
        print("writing:",e1name)
        c.write_column(e1name, e1dt, create=True)
        print("writing:",e2name)
        c.write_column(e2name, e2dt, create=True)
        print("writing:",dtflagname)
        c.write_column(dtflagname, dtflag, create=True)

    def plot_detrend(self, data, e1dt, e2dt, w, camcol, rmag_min, rmag_max, 
                     coeffs, show=False):
        """
        Make a table of plots.  
        """
        band=self.band
        d=self.plotdir()
        f = 'rg-%(run)s%(band)s-meane-vs-R-%(rmag_max)0.2f-%(camcol)s-detrend.eps'
        f = f % {'run':self.procrun,'band':self.band,
                 'rmag_max':rmag_max,'camcol':camcol}
        epsfile=path_join(d,'bycamcol',f)

        # one column for each type (R,e1psf,e2psf)
        tab=Table(1,3)
        gratio=1.61803399
        #tab.aspect_ratio = 1/gratio
        tab.aspect_ratio = 1./2.

        weights = 1.0/(0.32**2 + data['uncer_rg_'+band][w]**2)
        nperbin = 200000

        # trends vs R
        print("\n* R")
        xmin=0.29
        xmax=1.01
        if band=='r':
            ymin = -0.03
            ymax =  0.03
        else:
            ymin = -0.03
            ymax =  0.06
        pe1,pe1err,pe2,pe2err,ph,t1,t2 = \
            self.make_plot_components(data['R_rg_'+band][w],
                                      data['e1_rg_'+band][w], 
                                      data['e2_rg_'+band][w], 
                                      weights, nperbin, ymin, ymax,
                                      fmin=xmin,fmax=xmax,
                                      coeffs=coeffs)

        pe1dt,pe1dterr,pe2dt,pe2dterr,phdt = \
            self.make_plot_components(data['R_rg_'+band][w],
                                      e1dt[w], 
                                      e2dt[w], 
                                      weights, nperbin, ymin, ymax,
                                      fmin=xmin,fmax=xmax)



        # one row for before and after detrend
        Ra=FramedArray(2,1)
        Ra.xrange = [xmin,xmax]
        Ra.yrange = [ymin,ymax]
        Ra.xlabel = 'R'
        Ra.ylabel = '<e>'

        rmag_lab = \
            PlotLabel(0.1,0.07,'%0.2f < rmag < %0.2f' % (rmag_min,rmag_max), 
                      halign='left')
        cflab = PlotLabel(0.1,0.15,
                          'camcol: %s filter: %s' % (camcol, self.band), 
                          halign='left')

        key = PlotKey(0.1,0.9, [pe1,pe2])
        keydt = PlotKey(0.1,0.9, [pe1dt,pe2dt])
        Rz=Curve([xmin,xmax],[0,0])
        Ra[0,0].add( Rz,pe1,pe1err,pe2,pe2err,ph,t1,t2,key)
        Ra[1,0].add( Rz,pe1dt,pe1dterr,pe2dt,pe2dterr,phdt,keydt,rmag_lab,cflab)

        if show:
            Ra.show()
        tab[0,0] = Ra

        # trends vs e1 psf
        print("\n* e1_psf")
        xmin = -0.275
        xmax =  0.275
        #ymin = -0.015
        #ymax =  0.015

        pe1psf_e1,pe1psf_e1err,pe1psf_e2,pe1psf_e2err,pe1psf_ph = \
                self.make_plot_components(data['e1_psf_'+band][w],
                                          data['e1_rg_'+band][w], 
                                          data['e2_rg_'+band][w], 
                                          weights, nperbin, ymin, ymax,
                                          fmin=xmin,fmax=xmax)

        pe1psf_e1_dt,pe1psf_e1_dterr,pe1psf_e2_dt,pe1psf_e2_dterr,pe1psf_ph_dt = \
                self.make_plot_components(data['e1_psf_'+band][w],
                                          e1dt[w], 
                                          e2dt[w], 
                                          weights, nperbin, ymin, ymax,
                                          fmin=xmin,fmax=xmax)



        # one row for before and after detrend
        e1psf_a=FramedArray(2,1)
        e1psf_a.xrange = [xmin,xmax]
        e1psf_a.yrange = [ymin,ymax]
        e1psf_a.xlabel = r'$e1_{PSF}$'
        #e1psf_a.ylabel = '<e>'

        e1psf_key = PlotKey(0.1,0.9, [pe1psf_e1,pe1psf_e2])
        e1psf_keydt = PlotKey(0.1,0.9, [pe1psf_e1_dt,pe1psf_e2_dt])
        e1psf_z=Curve([xmin,xmax],[0,0])
        e1psf_a[0,0].add(e1psf_z, 
                         pe1psf_e1,pe1psf_e1err, 
                         pe1psf_e2, pe1psf_e2err,
                         pe1psf_ph, e1psf_key)
        e1psf_a[1,0].add(e1psf_z, 
                         pe1psf_e1_dt, pe1psf_e1_dterr,
                         pe1psf_e2_dt, pe1psf_e2_dterr,
                         pe1psf_ph_dt, 
                         e1psf_keydt)

        if show:
            e1psf_a.show()
        tab[0,1] = e1psf_a


        # trends vs e2 psf
        print("\n* e2_psf")
        pe2psf_e1,pe2psf_e1err,pe2psf_e2,pe2psf_e2err,pe2psf_ph = \
                self.make_plot_components(data['e2_psf_'+band][w],
                                          data['e1_rg_'+band][w], 
                                          data['e2_rg_'+band][w], 
                                          weights, nperbin, ymin, ymax,
                                          fmin=xmin,fmax=xmax)

        pe2psf_e1_dt,pe2psf_e1_dterr,pe2psf_e2_dt,pe2psf_e2_dterr,pe2psf_ph_dt = \
                self.make_plot_components(data['e2_psf_'+band][w],
                                          e1dt[w], 
                                          e2dt[w], 
                                          weights, nperbin, ymin, ymax,
                                          fmin=xmin,fmax=xmax)



        # one row for before and after detrend
        e2psf_a=FramedArray(2,1)
        e2psf_a.xrange = [xmin,xmax]
        e2psf_a.yrange = [ymin,ymax]
        e2psf_a.xlabel = r'$e2_{PSF}$'
        #e2psf_a.ylabel = '<e>'

        e2psf_key = PlotKey(0.1,0.9, [pe2psf_e1,pe2psf_e2])
        e2psf_keydt = PlotKey(0.1,0.9, [pe2psf_e1_dt,pe2psf_e2_dt])
        e2psf_z=Curve([xmin,xmax],[0,0])
        e2psf_a[0,0].add(e2psf_z, 
                         pe2psf_e1, pe2psf_e1err,
                         pe2psf_e2, pe2psf_e2err,
                         pe2psf_ph, e2psf_key)
        e2psf_a[1,0].add(e2psf_z, 
                         pe2psf_e1_dt, pe2psf_e1_dterr,
                         pe2psf_e2_dt, pe2psf_e2_dterr,
                         pe2psf_ph_dt, 
                         e2psf_keydt)

        if show:
            e2psf_a.show()
        tab[0,2] = e2psf_a

        if show:
            tab.show()

        tab.write_eps(epsfile)
        converter.convert(epsfile, dpi=150, verbose=True)

    def R_polyfile(self, camcol, rmag_max):

        f = 'rg-%(run)s%(band)s-meane-vs-R-%(rmag_max)0.2f-%(camcol)s-polyfit.yaml'
        f = f % {'run':self.procrun,'band':self.band,
                 'rmag_max':rmag_max,'camcol':camcol}
        d=self.plotdir()
        f=path_join(d,'bycamcol',f)
        return f

    def plotfile(self, field, rmag_max, plot_type='meane'):
        camcol=self.camcol
        d=self.plotdir()

        f = 'rg-%(run)s%(band)s-%(type)s-vs-%(field)s-%(rmag_max)0.2f'
        f = f % {'run':self.procrun,'band':self.band,'type':plot_type,
                 'field':field,'rmag_max':rmag_max}
        if self.run != 'any':
            f += '-%06i' % self.run
        if camcol != 'any':
            d = path_join(d,'bycamcol')
            f += '-%i' % camcol
        f += '.eps'
        f=path_join(d,f)
        return f

    def plotdir(self):
        d=os.environ['LENSDIR']
        d=path_join(d,'regauss-tests',self.procrun)
        if self.run != 'any':
            d = path_join(d,'%06i' % self.run)
        #if self.camcol != 'any':
        #    d = path_join(d,'bycamcol')
        return d


    def make_plot_components(self, fdata, e1, e2, weights, nperbin, 
                             ymin,ymax,
                             fmin=None, fmax=None,
                             coeffs=None):
        be1 = eu.stat.Binner(fdata, e1, weights=weights)
        be2 = eu.stat.Binner(fdata, e2, weights=weights)

        print("  hist  e1, nperbin: ",nperbin)
        be1.dohist(nperbin=nperbin, min=fmin, max=fmax)
        print("  stats e1")
        be1.calc_stats()
        print("  hist  e2, nperbin: ",nperbin)
        be2.dohist(nperbin=nperbin, min=fmin, max=fmax)
        print("  stats e2")
        be2.calc_stats()

        # regular hist for display
        print("  regular fixed binsize hist")
        xm,xe,xstd=eu.stat.wmom(fdata, weights, sdev=True)
        bsize = xstd/5.
        hist = eu.stat.histogram(fdata, binsize=bsize, 
                                 weights=weights, more=True)

        ph = eu.plotting.make_hist_curve(hist['low'], hist['high'], hist['whist'], 
                                         ymin=ymin, ymax=ymax, color='grey50')

        p1 = Points( be1['wxmean'], be1['wymean'], 
                    type='filled circle', color='blue')
        p1err = SymErrY( be1['wxmean'], be1['wymean'], be1['wyerr2'], color='blue')
        p1.label = r'$e_1$'

        p2 = Points( be2['wxmean'], be2['wymean'], 
                    type='filled circle', color='red')
        p2.label = r'$e_2$'
        p2err = SymErrY( be2['wxmean'], be2['wymean'], be2['wyerr2'], color='red')

        if coeffs is not None:
            poly1=numpy.poly1d(coeffs['coeff_e1'])
            poly2=numpy.poly1d(coeffs['coeff_e2'])
            t1 = Curve( be1['wxmean'], poly1(be1['wxmean']), color='blue')
            t2 = Curve( be2['wxmean'], poly2(be2['wxmean']), color='red')
            return p1,p1err,p2,p2err,ph,t1,t2
        else:
            return p1,p1err,p2,p2err,ph


def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = numpy.linspace(-2,2,0.1)
    x = numpy.sin(t)+numpy.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = numpy.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = numpy.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]
