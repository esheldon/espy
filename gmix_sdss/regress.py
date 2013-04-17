"""
Regressions vs s2n or epsf.

For epsf, the mean shape vs s2n is subtracted

All regressions are done on the pre-sensitivity corrected values.
"""
import os
import time
import numpy
import sdsspy
import esutil as eu

from . import collate
from . import files

from . import cuts
from .cuts import SRATIO_MIN,PSF_EMIN,PSF_EMAX,S2N_MIN,S2N_MAX

# value used for epsf regressions
SRATIO_MIN_EPSF=1.0
S2N_MIN_EPSF=20.0

def get_plot_dir(**keys):
    if 'gmix_run' not in keys:
        raise ValueError("send gmix_run=")

    d=files.get_basedir()
    d=os.path.join(d,keys['gmix_run'],'regress')
    return d

def get_plot_file(**keys):
    d=get_plot_dir(**keys)
    type=keys['type']

    f='regress-%s' % keys['gmix_run']
    f += '-%d' % keys['camcol']

    f += '-%s.eps' % type
    f=os.path.join(d,f)
    return f

def get_fit_file(**keys):
    """
    same args as for get_plot_file with sratio= also,
    type set to s2n
    """
    keys['type'] = 's2n'
    epsfile=get_plot_file(**keys)

    sstr='%.2f' % keys['sratio']
    fitfile=epsfile.replace('.eps','-fits-%s.pickle' % sstr)

    return fitfile

def write_fit_file(data, **keys):
    """
    same args as for get_plot_file with sratio= also,
    type set to s2n
    """
    import pickle
    fitfile=get_fit_file(**keys)

    try:
        eu.ostools.makedirs_fromfile(fitfile)
    except:
        pass

    print 'writing poly fit file:',fitfile
    with open(fitfile,'w') as fobj:
        pickle.dump(data,fobj)


def read_fit_file(**keys):
    """
    same args as for get_plot_file with sratio= also,
    type set to s2n
    """
    import pickle
    fitfile=get_fit_file(**keys)
    print 'reading poly fit file:',fitfile
    with open(fitfile) as fobj:
        data=pickle.load(fobj)
    return data


def detrend_g_vs_s2n(**keys):
    """
    gmix_run=, camcol=, sratio=
    s2n=, g1=, g2=
    sratio=
    """
    fdata=read_fit_file(**keys)

    s2n=keys['s2n']
    g1in = keys['g1']
    g2in = keys['g2']

    lx = numpy.log10(s2n)
    sub1 = fdata['ply1'](lx)
    sub2 = fdata['ply2'](lx)
    g1 = g1in - sub1
    g2 = g2in - sub2

    return g1,g2


class S2NRegressor(dict):
    def __init__(self, **keys):
        self.update(keys)


        if ('gmix_run' not in self 
                or 'camcol' not in self
                or 'nbin' not in self):
            raise ValueError("send gmix_run=,nbin=")

        self._conf=files.read_config(self['gmix_run'])

        self['s2n_label'] = 'S/N'

    def _dofit(self, x, y, yerr):
        import mcmc
        nwalkers=100
        burnin=40
        nstep=40
        degree=2
        fitter=mcmc.PolyFitter(degree, x, y, nwalkers, burnin, nstep,
                               yerr=yerr)
        return fitter

    def _plot_fit(self, plt, x, y, yerr, color):
        import biggles
        fitter=self._dofit(x,y,yerr)
        print fitter
        line=biggles.Curve(x, fitter(x), color=color)
        plt.add(line)

        return fitter.get_poly()

    def doplot(self):
        import biggles
        import pcolors

        data=self.get_data()
    
        arr=biggles.FramedArray(2,1)
        arr.xlabel=r'$log_{10}(%s)$' % self['s2n_label']
        arr.ylabel=r'$\gamma$'
        arr.yrange=[-0.015,0.015]
        arr.uniform_limits = 1
        arr.title=self._get_title()

        colors=['tan3','green4','blue','red']
        types=['filled triangle','filled circle','filled square',
               'filled diamond']
        pts_arr1=[]

        xmin=9999.e9
        xmax=-9999e9
        thresholds=self.get_sratio_thresh()
        for i,srthresh in enumerate(thresholds):

            print 'sratio thresh:',srthresh
            w,=numpy.where(data['sratio'] > srthresh)

            
            x1,y1,y1err,x2,y2,y2err=self._get_binned_data(w)

            pts1=biggles.Points(x1, y1, color=colors[i], type=types[i])
            err1=biggles.SymmetricErrorBarsY(x1, y1, y1err,
                   color=colors[i])


            pts2=biggles.Points(x2, y2, color=colors[i], type=types[i])
            err2=biggles.SymmetricErrorBarsY(x2,y2,y2err,
                   color=colors[i])

            labstr=r'$\sigma_{gal}/\sigma_{PSF} > %0.2f' % srthresh
            pts1.label=labstr

            pts_arr1.append(pts1)

            arr[0,0].add(pts1,err1)
            arr[1,0].add(pts2,err2)

            xmin0=min(x1.min(), x2.min())
            xmax0=min(x1.max(), x2.max())
            if xmin0 < xmin:
                xmin=xmin0
            if xmax0 > xmax:
                xmax=xmax0

            if srthresh == 1.0:
                ply1=self._plot_fit(arr[0,0], x1, y1, y1err, colors[i])
                ply2=self._plot_fit(arr[1,0], x2, y2, y2err, colors[i])
                self._write_fit_file(ply1,ply2,srthresh)


        key=biggles.PlotKey(0.9,0.9,pts_arr1, halign='right')
        arr[0,0].add(key)

        lab1=biggles.PlotLabel(0.1,0.9,r'$\gamma_1$')
        lab2=biggles.PlotLabel(0.1,0.9,r'$\gamma_2$')
        arr[0,0].add(lab1)
        arr[1,0].add(lab2)

        z=biggles.Curve([xmin,xmax],[0,0])
        arr[0,0].add(z)
        arr[1,0].add(z)

        arr.xrange=[0.9*xmin, 1.1*xmax]

        self._write_plots(arr)

    def _get_binned_data(self, w):

        data=self.get_data()

        lx = numpy.log10(data['s2n'][w])
        lxmin = numpy.log10(S2N_MIN)
        lxmax = numpy.log10(S2N_MAX)

        bs1=eu.stat.Binner(lx, data['g'][w,0])
        bs1.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
        bs1.calc_stats()

        bs2=eu.stat.Binner(lx, data['g'][w,1])
        bs2.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
        bs2.calc_stats()

        x1=bs1['xmean']
        y1=bs1['ymean']
        y1err=bs1['yerr']

        x2=bs2['xmean']
        y2=bs2['ymean']
        y2err=bs2['yerr']

        return x1,y1,y1err,x2,y2,y2err
            


        

    def _write_plots(self, plt):
        import converter
        epsfile=self._get_epsfile()
        print epsfile
        plt.write_eps(epsfile)
        print epsfile.replace('.eps','.png')
        converter.convert(epsfile,dpi=100)

    def get_data(self):
        if not hasattr(self,'_data'):
            self._load_data()
        return self._data


    def get_sratio_thresh(self):
        #return [0.5,1.0,1.5,2.0]
        return [2.0,1.5,1.0,0.5]

    def _get_epsfile(self):
        type='s2n'
        epsfile=get_plot_file(type=type, **self)
        try:
            eu.ostools.makedirs_fromfile(epsfile)
        except:
            pass
        return epsfile

    def _write_fit_file(self, ply1, ply2, sratio):
        data={'ply1':ply1, 'ply2':ply2}
        write_fit_file(data, sratio=sratio, **self)

    def _get_title(self):
        title='%s' % self['gmix_run']
        if 'run' in self:
            title += '-%06d' % self['run']
        title += '-%d' % self['camcol']

        return title

    def _load_data(self):
        cols=collate.open_columns(self['gmix_run'])

        selector=cuts.Selector(cols, do_sratio_cut=False)
        # using defaults with broad cuts
        selector.do_select(camcol=self['camcol'])

        ind=selector.indices
        self._data=selector.data

        print 'loading g'
        self._data['g'] = cols['g'][ind]


class EPSFRegressor(dict):
    def __init__(self, **keys):
        self.update(keys)

        if ('gmix_run' not in self 
                or 'camcol' not in self):
            raise ValueError("send gmix_run=,camcol=")

        if self['camcol'] is None:
            raise ValueError("send camcol for PSF ellip regress")

        self['ebinsize'] = 0.01

        self._conf=files.read_config(self['gmix_run'])

        self._load_data()
        self._correct_vs_s2n()

    def doplot(self):
        import biggles
        x1,y1,y1err,x2,y2,y2err=self._get_binned_data()

        arr=biggles.FramedArray(2,1)
        arr.xlabel=r'$e^{PSF}$'
        arr.ylabel=r'$\gamma$'

        arr.uniform_limits = 1
        arr.title=self._get_title()

        #data=self._data
        #arr.xrange=self._get_sym_range(data['psf_e1'],data['psf_e2'])
        arr.xrange=[PSF_EMIN,PSF_EMAX]
        arr.yrange=[-0.015,0.015]

        pts1=biggles.Points(x1, y1, type='filled circle')
        err1=biggles.SymmetricErrorBarsY(x1, y1, y1err)

        pts2=biggles.Points(x2, y2, type='filled circle')
        err2=biggles.SymmetricErrorBarsY(x2, y2, y2err)

        arr[0,0].add(pts1,err1)
        arr[1,0].add(pts2,err2)

        lab1=biggles.PlotLabel(0.1,0.9,r'$\gamma_1, e_1^{PSF}$')
        lab2=biggles.PlotLabel(0.1,0.9,r'$\gamma_2, e_2^{PSF}$')

        dlab=biggles.PlotLabel(0.95,0.2,'detrended',halign='right')
        cutstr=self._get_cut_string()
        cutlab=biggles.PlotLabel(0.95,0.1,cutstr,
                                halign='right')

        arr[0,0].add(lab1)
        arr[1,0].add(lab2)
        arr[1,0].add(dlab)
        arr[1,0].add(cutlab)

        z=biggles.Curve([-1,1],[0,0])
        arr[0,0].add(z)
        arr[1,0].add(z)

        self._write_plots(arr)

    def _get_cut_string(self):
        s=r'$S/N > %d \sigma_{gal}/\sigma_{PSF} > %.2f'
        s = s % (S2N_MIN_EPSF,SRATIO_MIN_EPSF)
        return s

    def _get_sym_range(self, x1, x2):
        data=self._data
        xmin=min(x1.min(), x2.min())
        xmax=max(x1.max(), x2.max())

        if abs(xmin) > abs(xmax):
            return [-abs(xmin), abs(xmin)]
        else:
            return [-abs(xmax), abs(xmax)]

    def _get_binned_data(self):
        data=self._data

        bs11=eu.stat.Binner(data['psf_e1'], data['g'][:,0])
        bs11.dohist(binsize=self['ebinsize'], 
                    min=PSF_EMIN, max=PSF_EMAX)
        bs11.calc_stats()

        bs22=eu.stat.Binner(data['psf_e2'], data['g'][:,1])
        bs22.dohist(binsize=self['ebinsize'], 
                    min=PSF_EMIN, max=PSF_EMAX)
        bs22.calc_stats()

        w,=numpy.where(bs11['hist'] > 0)
        x1=bs11['xmean'][w]
        y1=bs11['ymean'][w]
        y1err=bs11['yerr'][w]

        w,=numpy.where(bs22['hist'] > 0)
        x2=bs22['xmean'][w]
        y2=bs22['ymean'][w]
        y2err=bs22['yerr'][w]

        return x1,y1,y1err,x2,y2,y2err


    def _dofit(self, x, y, yerr):
        import mcmc
        nwalkers=100
        burnin=40
        nstep=40
        degree=2
        fitter=mcmc.PolyFitter(degree, x, y, nwalkers, burnin, nstep,
                               yerr=yerr)
        return fitter

    def _correct_vs_s2n(self):
        data=self._data

        print 'de-trending g vs s2n'
        g1,g2=detrend_g_vs_s2n(s2n=data['s2n'], 
                               g1=data['g'][:,0],
                               g2=data['g'][:,1],
                               sratio=SRATIO_MIN_EPSF,
                               **self)

        data['g'][:,0] = g1
        data['g'][:,1] = g2

    def _load_data(self):
        cols=collate.open_columns(self['gmix_run'])

        selector=cuts.Selector(cols)
        # using default broad cuts
        selector.do_select(camcol=self['camcol'],
                           sratio_min=SRATIO_MIN_EPSF,
                           s2n_min=S2N_MIN_EPSF)



        ind=selector.indices
        self._data=selector.data

        print 'loading g'
        self._data['g'] = cols['g'][ind]

    def _get_title(self):
        title='%s-%d' % (self['gmix_run'], self['camcol'])
        return title


    def _write_plots(self, plt):
        import converter
        epsfile=self._get_epsfile()
        print epsfile
        plt.write_eps(epsfile)
        print epsfile.replace('.eps','.png')
        converter.convert(epsfile,dpi=100)

    def _get_epsfile(self):
        type='epsf'
        epsfile=get_plot_file(type=type.replace('_','-'), 
                              **self)
        try:
            eu.ostools.makedirs_fromfile(epsfile)
        except:
            pass
        return epsfile

