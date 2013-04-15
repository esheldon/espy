import os
import time
import numpy
import sdsspy
import esutil as eu

from . import collate
from . import files

def get_plot_dir(**keys):
    if 'gmix_run' not in keys:
        raise ValueError("send gmix_run=")

    d=files.get_basedir()
    d=os.path.join(d,keys['gmix_run'],'plots')
    return d

def get_plot_file(**keys):
    d=get_plot_dir(**keys)
    type=keys['type']

    f='regress-%s' % keys['gmix_run']
    if 'run' in keys and keys['run'] is not None:
        f += '-%06d' % keys['run']
    if 'camcol' in keys and keys['camcol'] is not None:
        f += '-%d' % keys['camcol']

    f += '-%s.eps' % type
    f=os.path.join(d,f)
    return f

class S2NRegressor(dict):
    def __init__(self, **keys):
        self.update(keys)


        if ('gmix_run' not in self 
                or 'nbin' not in self):
            raise ValueError("send gmix_run=,nbin=")

        self._conf=files.read_config(self['gmix_run'])

        self['s2n_field'] = self.get('s2n_field','s2n')
        self['objtype']=self.get('objtype',None)

        # range allowed for PSF shape
        self['emin'] = -0.2
        self['emax'] =  0.2


        if self['s2n_field'] == 'Ts2n':
            self['s2n_label'] = r'$(S/N)_T$'
        else:
            self['s2n_label'] = 'S/N'
        self._set_fields()
        self._set_sratio_thresh()

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

        #colors=pcolors.rainbow(len(self['sratio_thresh']))
        colors=['red','tan3','green4','blue']
        types=['filled triangle','filled circle','filled square',
               'filled diamond']
        pts_arr1=[]

        xmin=9999.e9
        xmax=-9999e9
        for i,srthresh in enumerate(self['sratio_thresh']):

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
                self._plot_fit(arr[0,0], x1, y1, y1err, colors[i])
                self._plot_fit(arr[1,0], x2, y2, y2err, colors[i])


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

        s2nf=self['s2n_field']

        if s2nf == 'Ts2n':
            xmin=1.0
            xmax=100.0
        else:
            xmin=10.0
            xmax=300.0

        lx = numpy.log10(data[s2nf][w])
        lxmin = numpy.log10(xmin)
        lxmax = numpy.log10(xmax)

        bs1=eu.stat.Binner(lx, data['g'][w,0])
        bs1.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
        bs1.calc_stats()

        bs2=eu.stat.Binner(lx, data['g'][w,1])
        bs2.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
        bs2.calc_stats()

        x1=bs1['xmean']
        #x1=10**(bs1['xmean'])
        y1=bs1['ymean']
        y1err=bs1['yerr']

        x2=bs2['xmean']
        #x2=10**(bs2['xmean'])
        y2=bs2['ymean']
        y2err=bs2['yerr']

        if self._conf['obj_fitter'] == 'mcmc':
            sbs1=eu.stat.Binner(lx, data['gsens'][w,0])
            sbs1.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
            sbs1.calc_stats()

            sbs2=eu.stat.Binner(lx, data['gsens'][w,1])
            sbs2.dohist(nbin=self['nbin'], min=lxmin, max=lxmax)
            sbs2.calc_stats()

            sbs1['ymean'] = sbs1['ymean'].clip(0.001, 1.0)
            sbs2['ymean'] = sbs2['ymean'].clip(0.001, 1.0)
            #print 'gsens1',sbs1['ymean']
            #print 'gsens2',sbs2['ymean']

            y1 /= sbs1['ymean']
            y1err /= sbs1['ymean']
            y2 /= sbs2['ymean']
            y2err /= sbs2['ymean']

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

    def _set_fields(self):
        self['fields']=[self['s2n_field'],'g','sratio']

        if self._conf['obj_fitter']=='mcmc':
            self['fields'] += ['gsens']

        if self['objtype'] is not None:
            self['fields'] += ['model']
        print "self['fields']:",self['fields']

    def _set_sratio_thresh(self):
        self['sratio_thresh']=[0.5,1.0,1.5,2.0]

    def _get_epsfile(self):
        type=self['s2n_field']
        if self['objtype'] is not None:
            type += '-%s' % self['objtype']
        epsfile=get_plot_file(type=type.replace('_','-'), 
                              **self)
        try:
            eu.ostools.makedirs_fromfile(epsfile)
        except:
            pass
        return epsfile

    def _get_title(self):
        title='%s' % self['gmix_run']
        if 'run' in self:
            title += '-%06d' % self['run']
        if 'camcol' in self:
            title += '-%d' % self['camcol']

        return title

    def _load_data(self):
        cols=collate.open_columns(self['gmix_run'])

        print 'epsf read'
        ppars=cols['psf_pars'][:]
        print 'epsf logic'
        logic = (  (ppars[:,2] > self['emin'])
                 & (ppars[:,2] < self['emax'])
                 & (ppars[:,3] > self['emin'])
                 & (ppars[:,3] < self['emax']) )
        wkeep,=numpy.where(logic)
        print '%d/%d passed epsf' % (wkeep.size,ppars.size)

        if 'camcol' in self:
            print 'camcol logic:',self['camcol']
            camcol=cols['camcol'][:]
            logic = logic & (camcol==self['camcol'])

        print 'calling where'
        w,=numpy.where(logic)
        print '%d/%d passed epsf+' % (w.size,ppars.size)

        data=cols.read_columns(self['fields'], rows=w, verbose=True)

        if self['objtype'] is not None:
            from esutil.numpy_util import strmatch
            print 'selecting model:',self['objtype']
            reg='.*'+self['objtype']+'.*'
            w,=numpy.where( strmatch(data['model'], reg) )
            if w.size==0:
                raise ValueError("no objects of type '%s'" % self['objtype'])
            print '    ',w.size
            data=data[w]
        self._data=data



class RunS2NRegressor(S2NRegressor):
    def __init__(self, **keys):
        super(RunS2NRegressor,self).__init__(**keys)

        if 'run' not in self:
            raise ValueError("send run=[, camcol=]")

    def _load_data(self):
        cols=collate.open_columns(self['gmix_run'])


        w=(cols['run']==self['run'])
        if 'camcol' in self and self['camcol'] is not None:
            camcol=cols['camcol'][w]
            w2,=numpy.where(camcol==self['camcol'])
            w=w[w2]

        data=cols.read_columns(self['fields'],
                               rows=w,
                               verbose=True)

        self._data=data


class EPSFRegressor(dict):
    def __init__(self, **keys):
        self.update(keys)

        if ('gmix_run' not in self 
                or 'camcol' not in self):
            raise ValueError("send gmix_run=,nbin=,camcol=")

        if self['camcol'] is None:
            raise ValueError("send camcol for PSF ellip regress")

        self['ebinsize'] = 0.02
        self['emin'] = -0.2
        self['emax'] =  0.2

        self._conf=files.read_config(self['gmix_run'])

        self['objtype']=self.get('objtype',None)

        self['s2n_min'] = 20.
        self['sratio_min'] = 1.0

        self._load_data()

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
        arr.xrange=[self['emin'], self['emax']]
        arr.yrange=[-0.1,0.1]

        pts1=biggles.Points(x1, y1, type='filled circle')
        err1=biggles.SymmetricErrorBarsY(x1, y1, y1err)

        pts2=biggles.Points(x2, y2, type='filled circle')
        err2=biggles.SymmetricErrorBarsY(x2, y2, y2err)

        arr[0,0].add(pts1,err1)
        arr[1,0].add(pts2,err2)

        lab1=biggles.PlotLabel(0.1,0.9,r'$\gamma_1, e_1^{PSF}$')
        lab2=biggles.PlotLabel(0.1,0.9,r'$\gamma_2, e_2^{PSF}$')
        cutstr=self._get_cut_string()
        cutlab=biggles.PlotLabel(0.9,0.1,cutstr,
                                halign='right')

        arr[0,0].add(lab1)
        arr[1,0].add(lab2)
        arr[1,0].add(cutlab)

        z=biggles.Curve([-1,1],[0,0])
        arr[0,0].add(z)
        arr[1,0].add(z)

        self._write_plots(arr)

    def _get_cut_string(self):
        s=r'$S/N > %d \sigma_{gal}/\sigma_{PSF} > %.2f'
        s = s % (self['s2n_min'],self['sratio_min'])
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
        bs11.dohist(binsize=self['ebinsize'], min=self['emin'], max=self['emax'])
        bs11.calc_stats()

        bs22=eu.stat.Binner(data['psf_e2'], data['g'][:,1])
        bs22.dohist(binsize=self['ebinsize'], min=self['emin'], max=self['emax'])
        bs22.calc_stats()

        w,=numpy.where(bs11['hist'] > 0)
        x1=bs11['xmean'][w]
        y1=bs11['ymean'][w]
        y1err=bs11['yerr'][w]

        w,=numpy.where(bs22['hist'] > 0)
        x2=bs22['xmean'][w]
        y2=bs22['ymean'][w]
        y2err=bs22['yerr'][w]

        if self._conf['obj_fitter'] == 'mcmc':
            sbs1=eu.stat.Binner(data['psf_e1'], data['gsens'][:,0])
            sbs1.dohist(nbin=self['nbin'])
            sbs1.calc_stats()

            sbs2=eu.stat.Binner(data['psf_e2'], data['gsens'][:,1])
            sbs2.dohist(nbin=self['nbin'])
            sbs2.calc_stats()

            sbs1['ymean'] = sbs1['ymean'].clip(0.001, 1.0)
            sbs2['ymean'] = sbs2['ymean'].clip(0.001, 1.0)

            y1 /= sbs1['ymean']
            y1err /= sbs1['ymean']
            y2 /= sbs2['ymean']
            y2err /= sbs2['ymean']

        return x1,y1,y1err,x2,y2,y2err


    def _load_data(self):
        cols=collate.open_columns(self['gmix_run'])

        print 'loading select columns'
        print '    camcol'
        camcol = cols['camcol'][:]
        print '    sratio'
        sratio = cols['sratio'][:]
        print '    s2n'
        s2n    = cols['s2n'][:]

        print 'getting logic'
        logic=(  (camcol==self['camcol']) 
               & (s2n > self['s2n_min'])
               & (sratio > self['sratio_min']) )

        if self['objtype'] is not None:
            print '    (model)'
            model=cols['model'][:]
            logic = logic & (model == self['objtype'])

        print 'where'
        w,=numpy.where(logic)

        print 'psf pars'
        ppars_rec=cols['psf_pars'][w]
        psf_e1=ppars_rec[:,2]
        psf_e2=ppars_rec[:,3]

        w2,=numpy.where(  (psf_e1 > self['emin'])
                        & (psf_e1 < self['emax'])
                        & (psf_e2 > self['emin'])
                        & (psf_e2 < self['emax']) )

        psf_e1 = psf_e1[w2]
        psf_e2 = psf_e2[w2]
        w=w[w2]

        data={}

        print 'loading data'
        print '    g'
        data['g'] = cols['g'][w]
        if self._conf['obj_fitter'] == 'mcmc':
            print '    gsens'
            data['gsens'] = cols['gsens'][w]


        data['psf_e1'] = psf_e1
        data['psf_e2'] = psf_e2

        self._data=data

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
        if self['objtype'] is not None:
            type += '-%s' % self['objtype']
        epsfile=get_plot_file(type=type.replace('_','-'), 
                              **self)
        try:
            eu.ostools.makedirs_fromfile(epsfile)
        except:
            pass
        return epsfile



