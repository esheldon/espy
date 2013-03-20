import os
from sys import stderr
import numpy
import esutil as eu
from esutil.numpy_util import where1
import deswl
import converter
try:
    import lensing
except:
    pass

from . import collate


def get_match_dir(run, mock_catalog):
    d=deswl.files.collated_dir(run)
    d=os.path.join(d,'match-%s' % mock_catalog)
    if not os.path.exists(d):
        os.makedirs(d)
    return d

def get_plot_dir(run, mock_catalog):
    d=get_match_dir(run, mock_catalog)
    d=os.path.join(d,'plots')
    return d

def get_shear_compare_html_url(run, mock_catalog, type):
    f=get_shear_compare_plot_url(run, mock_catalog, type)
    f=f.replace('.eps','.html')
    return f
def get_shear_compare_plot_url(run, mock_catalog, type, num=None):
    d=get_plot_dir(run, mock_catalog)
    fname='shearcmp-%s-%s-%s' % (run,mock_catalog,type)
    if num is not None:
        fname = fname+'-%02d' % num
    fname += '.eps'
    fname = os.path.join(d,fname)
    eu.ostools.makedirs_fromfile(fname)
    return fname

def get_match_files(run,mock_catalog):
    cdir=deswl.files.collated_dir(run)
    runf=os.path.join(cdir,'%s-radec.dat' % run)

    d=get_match_dir(run, mock_catalog)
    mockf=os.path.join(d,'%s-radec.dat' % mock_catalog)
    matchf='matches-%s-%s.dat' % (run,mock_catalog)
    matchf=os.path.join(d,matchf)
    scriptf=os.path.join(d,'domatch.sh')

    return {'dir':d,
            'runf':runf,
            'mockf':mockf,
            'matchf':matchf,
            'scriptf':scriptf}

def match_colname(mock_catalog):
    colname = 'match_%s' % mock_catalog
    colname = colname.replace('.','_')
    colname = colname.replace('-','_')
    return colname

class Comparator(dict):
    def __init__(self, run, mock_catalog):
        self['run'] = run
        self['mock_catalog'] = mock_catalog

    def plot_meanshear_vs_trueshear(self, nperbin, indices=None, label=None,
                       show=False, num=None):
        import biggles
        from biggles import PlotLabel, Curve

        type='vs-shear'
        # will use cache on successive calls
        data=self.get_data()
        epsfile=get_shear_compare_plot_url(self['run'],self['mock_catalog'],
                                           type,
                                           num=num)
        print >>stderr,'will write to file:',epsfile

        if indices is None:
            b1=eu.stat.Binner(data['shear1_true'],
                              data['shear1'],
                              weights=data['weight'])
            b2=eu.stat.Binner(data['shear2_true'],
                              data['shear2'],
                              weights=data['weight'])
        else:
            w=indices # alias
            b1=eu.stat.Binner(data['shear1_true'][w],
                              data['shear1'][w],
                              weights=data['weight'][w])
            b2=eu.stat.Binner(data['shear2_true'][w],
                              data['shear2'][w],
                              weights=data['weight'][w])

        tab=biggles.Table(1,2)

        print >>stderr,'shear1'
        print >>stderr,'    hist'
        b1.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b1.calc_stats()
        xlabel=r'$\gamma_1^{true}$'
        ylabel=r'$\gamma_1^{meas}$'
        plt1=eu.plotting.bscatter(b1['wxmean'],b1['wymean'],yerr=b1['wyerr'],
                                  show=False, xlabel=xlabel, ylabel=ylabel)
        
        if label is not None:
            lab1=PlotLabel(0.9,0.9,label,halign='right')
            plt1.add(lab1)

        coeff1 = numpy.polyfit(b1['wxmean'], b1['wymean'], 1)
        poly1=numpy.poly1d(coeff1)
        print coeff1
        ps1 = Curve( b1['wxmean'], poly1(b1['wxmean']), color='blue')
        plt1.add(ps1)

        flabt1='m: %0.2f b: %0.3f' % (coeff1[0],coeff1[1])
        flab1=PlotLabel(0.1,0.2,flabt1,halign='left')
        plt1.add(flab1)

        plt1.aspect_ratio=1

        print >>stderr,'shear2'
        print >>stderr,'    hist'
        b2.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b2.calc_stats()
        xlabel=r'$\gamma_2^{true}$'
        ylabel=r'$\gamma_2^{meas}$'
        plt2=eu.plotting.bscatter(b2['wxmean'],b2['wymean'],yerr=b2['wyerr'],
                                  show=False, xlabel=xlabel, ylabel=ylabel)

        if label is not None:
            lab2=PlotLabel(0.9,0.9,label,halign='right')
            plt2.add(lab2)

        coeff2 = numpy.polyfit(b2['wxmean'], b2['wymean'], 1)
        poly2=numpy.poly1d(coeff2)
        print coeff2
        ps2 = Curve( b2['wxmean'], poly2(b2['wxmean']), color='blue')
        plt2.add(ps2)

        flabt2='m: %0.2f b: %0.3f' % (coeff2[0],coeff2[1])
        flab2=PlotLabel(0.1,0.2,flabt2,halign='left')
        plt2.add(flab2)

        plt2.aspect_ratio=1

        tab[0,0] = plt1
        tab[0,1] = plt2
        if show:
            tab.show()
        print >>stderr,'writing plot:',epsfile
        tab.write_eps(epsfile)
        converter.convert(epsfile,dpi=120,verbose=True)

        res={'coeff1':coeff1,'coeff2':coeff2,'plt1':plt1,'plt2':plt2}
        return res,epsfile.replace('.eps','.png')

    def plot_sheardiff_vs_field(self, nperbin, field, show=False):
        import biggles
        biggles.configure('default','fontsize_min',1.5)
        biggles.configure('_HalfAxis','ticks_size',2.5)
        biggles.configure('_HalfAxis','subticks_size',1.25)
        data=self.get_data()
        type='vs_'+field

        epsfile=get_shear_compare_plot_url(self['run'],self['mock_catalog'],
                                           type)

        tab=biggles.Table(1,2)
        names={'shear1_diff':r'$\gamma_1-\gamma_1^{true}',
               'shear2_diff':r'$\gamma_2-\gamma_2^{true}'}
        if field == 'shear_s2n':
            xlog=True
        else:
            xlog=False
        plots=eu.plotting.bhist_vs(data, field, 'shear1_diff','shear2_diff', 
                                   nperbin=nperbin,names=names, 
                                   size=2.5,
                                   xlog=xlog,
                                   show=False)
        tab[0,0] = plots[0]
        tab[0,1] = plots[1]
        tab[0,0].aspect_ratio=1
        tab[0,1].aspect_ratio=1
        tab.write_eps(epsfile)
        converter.convert(epsfile,dpi=120,verbose=True)
        pngfile=epsfile.replace('.eps','.png')
        if show:
            tab.show()
        return pngfile

    def plot_sheardiff_vs(self, nperbin, show=False):
        type='vs'
        html_file = get_shear_compare_html_url(self['run'], 
                                               self['mock_catalog'], 
                                               type)
        pngfiles=[]
        fields=['shear_s2n']
        data=self.get_data()
        if 'mag_model' in data:
            fields+=['mag_model']

        for field in fields:
            png=self.plot_sheardiff_vs_field(nperbin, field, show=show)
            pngfiles.append(png)

        self.write_html(pngfiles, html_file)

    def plot_meanshear_vs_trueshear_bys2n(self, nperbin, nperbin_sub, 
                                          show=False):
        from biggles import FramedPlot, Curve, Points, Table, PlotKey

        type='vs_shear'

        html_file = get_shear_compare_html_url(self['run'], 
                                               self['mock_catalog'], 
                                               type)
        print 'html_file:',html_file

        data=self.get_data()
        
        epsfile=get_shear_compare_plot_url(self['run'],self['mock_catalog'],
                                           type)
        print >>stderr,'Will write summary plot:',epsfile

        print >>stderr,'histogramming shear_s2n',nperbin
        hdict=eu.stat.histogram(data['shear_s2n'],nperbin=nperbin,
                                rev=True,more=True)
        rev=hdict['rev']
        cdlist=[]
        pngfiles=[]
        for i in xrange(hdict['hist'].size):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]
                label=r'$%0.3g < S/N < %0.3g$' \
                    % (hdict['low'][i],hdict['high'][i])
                print >>stderr,label,'mean:',hdict['mean'][i],'num:',w.size
                cd,png=self.plot_meanshear_vs_trueshear(nperbin_sub,
                                                    indices=w,label=label,
                                                    num=i,
                                                    show=show)
                cdlist.append(cd)
                pngfiles.append(png)
        
        slopes1=[cd['coeff1'][0] for cd in cdlist]
        offsets1=[cd['coeff1'][1] for cd in cdlist]
        slopes2=[cd['coeff2'][0] for cd in cdlist]
        offsets2=[cd['coeff2'][1] for cd in cdlist]

        tab=Table(2,1)

        plt_slope=FramedPlot()
        cslope1 = Curve(hdict['mean'],slopes1,color='red')
        cslope2 = Curve(hdict['mean'],slopes2,color='blue')
        pslope1 = Points(hdict['mean'],slopes1,color='red',
                         type='filled circle')
        pslope2 = Points(hdict['mean'],slopes2,color='blue',
                         type='filled circle')

        cslope1.label = r'$\gamma_1$'
        cslope2.label = r'$\gamma_2$'
        key=PlotKey(0.1,0.2,[cslope1,cslope2],halign='left')

        plt_slope.add(cslope1,cslope2,pslope1,pslope2,key)
        plt_slope.xlabel = 'Shear S/N'
        plt_slope.ylabel = 'Slope'
        plt_slope.xlog=True
        plt_slope.xrange = eu.plotting.get_log_plot_range(hdict['mean'])

        plt_offset=FramedPlot()
        coffset1 = Curve(hdict['mean'],offsets1,color='red')
        coffset2 = Curve(hdict['mean'],offsets2,color='blue')
        poffset1 = Points(hdict['mean'],offsets1,color='red',
                          type='filled circle')
        poffset2 = Points(hdict['mean'],offsets2,color='blue',
                          type='filled circle')
        plt_offset.add(coffset1,coffset2,poffset1,poffset2)
        plt_offset.xlabel = 'Shear S/N'
        plt_offset.ylabel = 'Offset'
        plt_offset.xlog=True
        plt_offset.xrange = eu.plotting.get_log_plot_range(hdict['mean'])


        tab[0,0] = plt_slope
        tab[1,0] = plt_offset
        if show:
            tab.show()

        print >>stderr,'Writing summary plot:',epsfile
        tab.write_eps(epsfile)
        converter.convert(epsfile,dpi=90,verbose=True)
        png=epsfile.replace('.eps','.png')
        pngfiles = [png] + pngfiles
        
        self.write_html(pngfiles, html_file)

    def get_data(self):
        if not hasattr(self,'_data'):
            dmc=lensing.scat.DESMockCatalog(self['mock_catalog'])
            mock=dmc.read()

            c=collate.open_columns(self['run'])
            colname = match_colname(self['mock_catalog'])

            print >>stderr,'Reading match column'
            matches=c[colname][:]
            print >>stderr,'Reading flags'
            flags=c['shear_flags'][:]
            print >>stderr,'Getting good matches and flags'
            w=where1( (matches >= 0) & (flags==0) )

            d={}
            # deal with old names in multi epoch
            if 'shear_cov22' in c:
                covnames = ['shear_cov11','shear_cov22']
                fixnames=False
            else:
                covnames = ['shear_cov00','shear_cov11']
                fixnames=True

            names2read = \
                ['shear1','shear2','shear_s2n'] + covnames

            if self['run'][0:2] == 'me':
                names2read+= ['mag_model']

            for k in names2read:
                print >>stderr,'Reading',k
                tmp = c[k][:]
                tmp = tmp[w]
                if fixnames:
                    if k == 'shear_cov00':
                        k = 'shear_cov11'
                    elif k == 'shear_cov11':
                        k = 'shear_cov22'
                d[k] = tmp
                
            matches=matches[w]

            d['shear1_true'] = mock['gamma1'][matches]
            d['shear2_true'] = mock['gamma2'][matches]
            d['shear1_diff'] = d['shear1']-d['shear1_true']
            d['shear2_diff'] = d['shear2']-d['shear2_true']
            d['matches'] = matches
            d['mock'] = mock
            d['weight'] = 1.0/(0.32 + d['shear_cov11']+d['shear_cov22'])
            self._data = d
        return self._data

    def write_html(self, pngfiles, html_file):
        print 'writing html:',html_file
        with open(html_file,'w') as fobj:
            fobj.write("<html>\t<body>\n")
            for png in pngfiles:
                png = os.path.basename(png)
                fobj.write('\t\t<img src="%s"><b>\n' % png)
                            
            fobj.write("\t</body>\n</html>\n")

