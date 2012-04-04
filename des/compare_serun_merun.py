import os
import sys
from sys import stderr
import deswl
import des

import lensing
import esutil as eu
from esutil.numpy_util import where1
import numpy

import converter
import recfile

def get_match_dir(serun,merun):
    sedir=deswl.files.collated_dir(serun)
    match_dir=os.path.join(sedir, 'match-%s' % merun)
    return match_dir
def get_match_files(serun, merun):
    """
    Matches go in the se collated dir
    """
    sedir=deswl.files.collated_dir(serun)
    medir=deswl.files.collated_dir(merun)
    se_radec=os.path.join(sedir,'%s-radec.dat' % serun)
    me_radec=os.path.join(medir,'%s-radec.dat' % merun)

    match_dir=get_match_dir(serun,merun)
    matchf = 'matches-%s-%s.dat' % (serun,merun)
    matchf = os.path.join(match_dir, matchf)

    scriptf=os.path.join(match_dir,'domatch.sh')

    return {'dir':match_dir,
            'sef':se_radec,
            'mef':me_radec,
            'matchf':matchf,
            'scriptf':scriptf}

def get_plot_dir(serun, merun):
    d=get_match_dir(serun, merun)
    d=os.path.join(d,'plots')
    return d


def get_shear_compare_plot_url(serun, merun, num=None):
    d=get_plot_dir(serun,merun)
    fname='shearcmp-%s-%s' % (serun,merun)
    if num is not None:
        fname = fname+'-%02d' % num
    fname += '.eps'
    fname = os.path.join(d,fname)
    eu.ostools.makedirs_fromfile(fname)
    return fname


def get_match_colname(merun):
    colname = 'match_%s' % merun
    return colname

class Comparator(dict):
    def __init__(self, serun, merun):
        self['serun'] = serun
        self['merun'] = merun

    def plot_sheardiff(self, nperbin, indices=None, label=None,
                       show=False, num=None):
        import biggles
        from biggles import PlotLabel, Curve

        # will use cache on successive calls
        data=self.get_sheardiff_data()
        epsfile=get_shear_compare_plot_url(self['serun'],self['merun'],
                                           num=num)
        print >>stderr,'will write to file:',epsfile

        if indices is None:
            b1=eu.stat.Binner(data['shear1me'],
                              data['shear1'],
                              weights=data['weight'])
            b2=eu.stat.Binner(data['shear2me'],
                              data['shear2'],
                              weights=data['weight'])
        else:
            w=indices # alias
            b1=eu.stat.Binner(data['shear1me'][w],
                              data['shear1'][w],
                              weights=data['weight'][w])
            b2=eu.stat.Binner(data['shear2me'][w],
                              data['shear2'][w],
                              weights=data['weight'][w])

        tab=biggles.Table(1,2)

        print >>stderr,'shear1'
        print >>stderr,'    hist'
        b1.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b1.calc_stats()
        xlabel=r'$\gamma_1^{ME}$'
        ylabel=r'$\gamma_1^{SE}$'
        plt1=eu.plotting.bscatter(b1['wxmean'],b1['wymean'],yerr=b1['wyerr'],
                                  show=False, xlabel=xlabel, ylabel=ylabel)
        
        if label is not None:
            lab1=PlotLabel(0.1,0.9,label,halign='left')
            plt1.add(lab1)

        coeff1 = numpy.polyfit(b1['wxmean'], b1['wymean'], 1)
        poly1=numpy.poly1d(coeff1)
        print coeff1
        ps1 = Curve( b1['wxmean'], poly1(b1['wxmean']), color='blue')
        plt1.add(ps1)

        flabt1='m: %0.2f b: %0.3f' % (coeff1[0],coeff1[1])
        flab1=PlotLabel(0.9,0.2,flabt1,halign='right')
        plt1.add(flab1)

        plt1.aspect_ratio=1

        print >>stderr,'shear2'
        print >>stderr,'    hist'
        b2.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b2.calc_stats()
        xlabel=r'$\gamma_2^{ME}$'
        ylabel=r'$\gamma_2^{SE}$'
        plt2=eu.plotting.bscatter(b2['wxmean'],b2['wymean'],yerr=b2['wyerr'],
                                  show=False, xlabel=xlabel, ylabel=ylabel)

        if label is not None:
            lab2=PlotLabel(0.1,0.9,label,halign='left')
            plt2.add(lab2)

        coeff2 = numpy.polyfit(b2['wxmean'], b2['wymean'], 1)
        poly2=numpy.poly1d(coeff2)
        print coeff2
        ps2 = Curve( b2['wxmean'], poly2(b2['wxmean']), color='blue')
        plt2.add(ps2)

        flabt2='m: %0.2f b: %0.3f' % (coeff2[0],coeff2[1])
        flab2=PlotLabel(0.9,0.2,flabt2,halign='right')
        plt2.add(flab2)

        plt2.aspect_ratio=1

        tab[0,0] = plt1
        tab[0,1] = plt2
        if show:
            tab.show()
        print >>stderr,'writing plot:',epsfile
        tab.write_eps(epsfile)
        converter.convert(epsfile,dpi=120,verbose=True)

        return {'coeff1':coeff1,'coeff2':coeff2,'plt1':plt1,'plt2':plt2}


    def plot_sheardiff_bys2n(self, nperbin, nperbin_sub, show=False):
        import biggles
        from biggles import FramedPlot,PlotLabel, Curve, Points, Table, PlotKey

        data=self.get_sheardiff_data()
        
        epsfile=get_shear_compare_plot_url(self['serun'],self['merun'])
        print >>stderr,'Will write summary plot:',epsfile

        print >>stderr,'histogramming shear_s2n',nperbin
        hdict=eu.stat.histogram(data['shear_s2n'],nperbin=nperbin,
                                rev=True,more=True)
        rev=hdict['rev']
        cdlist=[]
        for i in xrange(hdict['hist'].size):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]
                label=r'$%0.3g < S/N < %0.3g$' \
                    % (hdict['low'][i],hdict['high'][i])
                print >>stderr,label,'mean:',hdict['mean'][i],'num:',w.size
                cd=self.plot_sheardiff(nperbin_sub,indices=w,label=label,
                                       num=i,
                                       show=show)
                cdlist.append(cd)
        
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


    def get_sheardiff_data(self):
        if not hasattr(self,'_data'):
            sec=des.collate.open_columns(self['serun'])
            mec=des.collate.open_columns(self['merun'])
            colname = get_match_colname(self['merun'])

            print >>stderr,'Reading match column'
            matches=sec[colname][:]
            print >>stderr,'Reading se flags'
            flags=sec['shear_flags'][:]
            print >>stderr,'Getting good matches and se flags'
            w=where1( (matches >= 0) & (flags==0) )
            matches=matches[w]

            print >>stderr,'Getting good me flags'
            meflags=mec['shear_flags'][:]
            meflags = meflags[matches]
            
            wm=where1(meflags == 0)
            matches=matches[wm]
            w=w[wm]

            d={}
            for k in ['shear1','shear2','shear_cov11','shear_cov22','shear_s2n']:
                print >>stderr,'Reading',k
                tmp = sec[k][:]
                tmp = tmp[w]
                d[k] = tmp
                
            # could be dups, need to extract all first
            shear1me = mec['shear1'][:]
            shear2me = mec['shear2'][:]
            d['shear1me'] = shear1me[matches]
            d['shear2me'] = shear2me[matches]
            d['matches'] = matches
            d['weight'] = 1.0/(0.32 + d['shear_cov11']+d['shear_cov22'])
            self._data = d
        return self._data

