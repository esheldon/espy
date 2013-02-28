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
                or 'nperbin' not in self):
            raise ValueError("send gmix_run=,nperbin=")

        self._conf=files.read_config(self['gmix_run'])

        self['s2n_field'] = self.get('s2n_field','s2n')
        self['objtype']=self.get('objtype',None)

        self['s2n_label'] = 'S/N'
        self._set_fields()
        self._set_sratio_thresh()

    def doplot(self):
        import biggles
        import pcolors

        data=self.get_data()
    
        arr=biggles.FramedArray(2,1)
        arr.xlabel=self['s2n_label']
        arr.ylabel=r'$\gamma$'
        arr.xlog=True
        arr.uniform_limits = 1
        arr.title=self._get_title()
        #arr.yrange=[-

        colors=pcolors.rainbow(len(self['sratio_thresh']))
        pts_arr1=[]

        xmin=9999.e9
        xmax=-9999e9
        for i,srthresh in enumerate(self['sratio_thresh']):

            print 'sratio thresh:',srthresh
            w,=numpy.where(data['sratio'] > srthresh)

            
            x1,y1,y1err,x2,y2,y2err=self._get_binned_data(w)

            pts1=biggles.Points(x1, y1, color=colors[i], type='filled circle')
            err1=biggles.SymmetricErrorBarsY(x1, y1, y1err,
                                             color=colors[i], type='filled circle')


            pts2=biggles.Points(x2, y2, color=colors[i], type='filled circle')
            err2=biggles.SymmetricErrorBarsY(x2,y2,y2err,
                                             color=colors[i], type='filled circle')

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

        s2nf=self['s2n_field']
        data=self.get_data()

        bs1=eu.stat.Binner(data[s2nf][w], data['g'][w,0])
        bs1.dohist(nperbin=self['nperbin'])
        bs1.calc_stats()

        bs2=eu.stat.Binner(data[s2nf][w], data['g'][w,1])
        bs2.dohist(nperbin=self['nperbin'])
        bs2.calc_stats()

        x1=bs1['xmean']
        y1=bs1['ymean']
        y1err=bs1['yerr']

        x2=bs2['xmean']
        y2=bs2['ymean']
        y2err=bs2['yerr']

        if self._conf['obj_fitter'] == 'mcmc':
            sbs1=eu.stat.Binner(data[s2nf][w], data['gsens'][w,0])
            sbs1.dohist(nperbin=self['nperbin'])
            sbs1.calc_stats()

            sbs2=eu.stat.Binner(data[s2nf][w], data['gsens'][w,1])
            sbs2.dohist(nperbin=self['nperbin'])
            sbs2.calc_stats()

            print 'gsens1',sbs1['ymean']
            print 'gsens2',sbs2['ymean']

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
        print "self['fields']:",self['fields']

        if self._conf['obj_fitter']=='mcmc':
            self['fields'] += ['gsens']

        if self['objtype'] is not None:
            self['fields'] += ['model']

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

        if 'camcol' in self:
            print 'getting camcol:',self['camcol']
            w=cols['camcol']==self['camcol']
            data=cols.read_columns(self['fields'], rows=w, verbose=True)
        else:
            data=cols.read_columns(self['fields'], verbose=True)

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


