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
        SymmetricErrorBarsY as SymErrY, SymmetricErrorBarsX as SymErrX

import fimage
from fimage.conversions import mom2sigma,mom2fwhm

from . import regauss

import converter

class Tester(dict):
    """
    testing a single run
    """
    def __init__(self,procrun,sweeptype,band, run='any', camcol='any', coldir=None):
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
                             'amflags_rg':'amflags',
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
        #  - amflags here is really rg for gals
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
            be1 = eu.stat.Binner(field_data, self['e1'][w]-self['e1_psf'][w], weights=weights)
            be2 = eu.stat.Binner(field_data, self['e2'][w]-self['e2_psf'][w], weights=weights)
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
        hist = eu.stat.histogram(field_data, binsize=bsize, weights=weights, more=True)


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

        p1 = Points( be1['wxmean'], be1['wymean'], type='filled circle', color='blue')
        p1err = SymErrY( be1['wxmean'], be1['wymean'], be1['wyerr2'], color='blue')
        p1.label = r'$e_1$'

        p2 = Points( be2['wxmean'], be2['wymean'], type='filled circle', color='red')
        p2.label = r'$e_2$'
        p2err = SymErrY( be2['wxmean'], be2['wymean'], be2['wyerr2'], color='red')

        key = PlotKey(0.1,0.9, [p1,p2])
        plt.add(p1, p1err, p2, p2err, key)

        if field != 'rmag':
            rmag_lab = PlotLabel(0.1,0.05,'%0.2f < rmag < %0.2f' % (rmag_min,rmag_max), halign='left')
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

       

    def plotfile(self, field, rmag_max, plot_type='meane'):
        f = 'rg-%(run)s%(band)s-%(type)s-vs-%(field)s-%(rmag_max)0.2f'
        f = f % {'run':self.procrun,'band':self.band,'type':plot_type,
                 'field':field,'rmag_max':rmag_max}
        if self.run != 'any':
            f += '-%06i' % self.run
        if self.camcol != 'any':
            f += '-%i' % self.camcol
        f += '.eps'
        d=self.plotdir()
        f=path_join(d,f)
        return f

    def plotdir(self):
        d=os.environ['LENSDIR']
        d=path_join(d,'regauss-tests',self.procrun)
        if self.run != 'any':
            d = path_join(d,'%06i' % self.run)
        if self.camcol != 'any':
            #d = path_join(d,'%i' % self.camcol)
            d = path_join(d,'bycamcol')
        return d
