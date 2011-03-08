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

class Tester(dict):
    """
    testing a single run
    """
    def __init__(self,procrun,band, run='any'):
        self.procrun = procrun
        self.run = run
        self.band = band

        if self.procrun in ['01','02','03']:
            self.old=True
        else:
            self.old=False


    def load_data(self):
        if 'field' not in self:
            print("Opening columns")
            c = regauss.open_columns(self.procrun)

            if self.run != 'any':
                print('Getting indices of run:',self.run)
                w=(c['run'].match(self.run))

                if w.size == 0:
                    raise ValueError("Run",self.run,"not found")
            else:
                w=None

            # basic tags to load
            if self.old:
                cn_byband={'R_rg':'r_rg',
                           'amflags_rg':'whyflag_rg',
                           'e1_rg':'e1_rg','e2_rg':'e2_rg',
                           'uncer_rg':'momerr_rg',
                           'Irr_psf':'iyy_rg','Icc_psf':'ixx_rg',
                           'e1_psf':'e1_psf','e2_psf':'e2_psf'}
                for n in cn_byband:
                    print("    loading col: ",n)
                    nn = cn_byband[n] + '_'+self.band
                    if w is None:
                        self[n] = c[nn][:]
                    else:
                        self[n] = c[nn][w]

            else:
                cn_byband=['R_rg',
                           'amflags_rg',
                           'e1_rg','e2_rg',
                           'uncer_rg',
                           'Irr_psf','Icc_psf',
                           'e1_psf','e2_psf']

                for name in cn_byband:
                    cname = name+'_'+self.band
                    print("    loading col: ",cname)
                    if w is None:
                        self[name] = c[cname][:]
                    else:
                        self[name] = c[cname][w]

            Tpsf = self['Irr_psf'] + self['Icc_psf']
            self['psf_sigma'] = mom2sigma(Tpsf)
            self['psf_fwhm'] = mom2fwhm(Tpsf, pixscale=0.396)

            cn = ['field','cmodelmag_dered_r']
            for name in cn:
                print("    loading col: ",name)
                if w is None:
                    self[name] = c[name][:]
                else:
                    self[name] = c[name][w]



    def plot_ellip_vs_field(self, field, rmag_max=21.8, fmin=None, fmax=None, nbin=20, nperbin=50000,
                            yrange=None, show=True):
        self.load_data()

        w=where1(  (self['R_rg'] > 1.0/3.0) 
                 & (self['R_rg'] < 1.0)
                 & (self['amflags_rg'] == 0)
                 & (self['e1_rg'] < 4)
                 & (self['e1_rg'] > -4)
                 & (self['cmodelmag_dered_r'] > 18.0)
                 & (self['cmodelmag_dered_r'] < rmag_max) )

        if w.size == 0:
            print("no good objects")
            return

        weights = 1.0/(0.32**2 + self['uncer_rg'][w]**2)

        if field == 'psf_fwhm':
            field_data = self['psf_fwhm'][w]
            fstr = 'PSF FWHM (arcsec)'
        elif field == 'psf_sigma':
            field_data = self['psf_sigma'][w]
            fstr = r'$\sigma_{PSF}$'
        else:
            field_data = self[field][w]
            fstr=field

        print("Plotting mean e for field:",field)

        fstr = fstr.replace('_','\_')

        be1 = eu.stat.Binner(field_data, self['e1_rg'][w], weights=weights)
        be2 = eu.stat.Binner(field_data, self['e2_rg'][w], weights=weights)

        print("  hist  e1")
        be1.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be1.dohist(nbin=nbin, min=fmin, max=fmax)
        print("  stats e1")
        be1.calc_stats()
        print("  hist  e2")
        be2.dohist(nperbin=nperbin, min=fmin, max=fmax)
        #be2.dohist(nbin=nbin, min=fmin, max=fmax)
        print("  stats e2")
        be2.calc_stats()

        plt = FramedPlot()
        p1 = Points( be1['wxmean'], be1['wymean'], type='filled circle', color='blue')
        p1err = SymErrY( be1['wxmean'], be1['wymean'], be1['wyerr2'], color='blue')
        p1.label = r'$e_1$'

        p2 = Points( be2['wxmean'], be2['wymean'], type='filled circle', color='red')
        p2.label = r'$e_2$'
        p2err = SymErrY( be2['wxmean'], be2['wymean'], be2['wyerr2'], color='red')

        key = PlotKey(0.1,0.9, [p1,p2])
        plt.add(p1, p1err, p2, p2err, key)

        if field != 'cmodelmag_dered_r':
            rmag_lab = PlotLabel(0.1,0.05,'rmag < %0.2f' % rmag_max, halign='left')
            plt.add(rmag_lab)

        procrun_lab = PlotLabel(0.1,0.1,
                                'procrun: %s filter: %s' % (self.procrun, self.band), 
                                halign='left')
        plt.add(procrun_lab)

        plt.xlabel = r'$'+fstr+'$'
        plt.ylabel = r'$<e>$'

        if yrange is not None:
            plt.yrange=yrange

        if show:
            plt.show()
        epsfile = self.plotfile(field, rmag_max)
        print("  Writing eps file:",epsfile)
        plt.write_eps(epsfile)

       

    def plotfile(self, field, rmag_max):
        f = 'rg-meane-%s-vs-%s-rmag%0.2f' % (self.band, field, rmag_max)
        if self.run != 'any':
            f += '-%06i' % self.run
        f += '.eps'
        d=self.plotdir()
        f=path_join(d,f)
        return f
    def plotdir(self):
        d=os.environ['LENSDIR']
        d=path_join(d,'regauss-tests')
        return d
