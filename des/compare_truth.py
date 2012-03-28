import os
from sys import stderr
import numpy
import esutil as eu
from esutil.numpy_util import where1
import lensing
import deswl
from . import collate


def get_match_dir(run, mock_catalog):
    d=deswl.files.collated_dir(run)
    d=os.path.join(d,'match-%s' % mock_catalog)
    if not os.path.exists(d):
        os.makedirs(d)
    return d

def get_match_files(run,mock_catalog):
    d=get_match_dir(run, mock_catalog)
    runf=os.path.join(d,'%s-radec.dat' % run)
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

    def plot_sheardiff(self):
        import biggles
        from biggles import PlotLabel, Curve

        nperbin=100000
        data=self.get_sheardiff_data()

        b1=eu.stat.Binner(data['shear1_true'],
                          data['shear1'],
                          weights=data['weight'])
        b2=eu.stat.Binner(data['shear2_true'],
                          data['shear2'],
                          weights=data['weight'])

        print >>stderr,'shear1'
        print >>stderr,'    hist'
        b1.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b1.calc_stats()
        plt1=eu.plotting.bscatter(b1['wxmean'],b1['wymean'],yerr=b1['yerr'],
                                  show=False)
        lab1=PlotLabel(0.9,0.9,r'$\gamma_1',halign='right')
        plt1.add(lab1)

        coeff1 = numpy.polyfit(b1['wxmean'], b1['wymean'], 1)
        poly1=numpy.poly1d(coeff1)
        print coeff1
        ps1 = Curve( b1['wxmean'], poly1(b1['wxmean']), color='blue')
        plt1.add(ps1)

        plt1.show()

        print >>stderr,'shear2'
        print >>stderr,'    hist'
        b2.dohist(nperbin=nperbin)
        print >>stderr,'    calc_stats'
        b2.calc_stats()
        plt2=eu.plotting.bscatter(b2['wxmean'],b2['wymean'],yerr=b2['yerr'],
                                  show=False)

        lab2=PlotLabel(0.9,0.9,r'$\gamma_2',halign='right')
        plt2.add(lab2)

        coeff2 = numpy.polyfit(b2['wxmean'], b2['wymean'], 1)
        poly2=numpy.poly1d(coeff2)
        print coeff2
        ps2 = Curve( b2['wxmean'], poly2(b2['wxmean']), color='blue')
        plt2.add(ps2)

        plt2.show()
        
    def get_sheardiff_data(self):
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
        for k in ['shear1','shear2','shear_cov11','shear_cov22']:
            print >>stderr,'Reading',k
            tmp = c[k][:]
            tmp = tmp[w]
            d[k] = tmp
            
        matches=matches[w]

        d['shear1_true'] = mock['gamma1'][matches]
        d['shear2_true'] = mock['gamma2'][matches]
        d['matches'] = matches
        d['mock'] = mock
        d['weight'] = 1.0/(0.32 + d['shear_cov11']+d['shear_cov22'])
        return d
